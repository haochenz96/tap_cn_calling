import numpy as np
import pandas as pd
import math

import statsmodels
import statsmodels.api as sm
from statsmodels.base.model import GenericLikelihoodModel

from scipy.stats import nbinom
from scipy.special import logsumexp

def get_responsibilities_and_marginal(y, mixing_props, cn, total_reads, phi):
    nobservations = len(y)
    nclones = len(cn)
    # coeffs = np.zeros((nobservations, nclones))
    logcoeffs = np.zeros((nobservations, nclones))
    
    for clone_idx in range(nclones):
        mu = cn[clone_idx] * total_reads
        prob  = phi / (phi + mu)
        # coeffs[:, clone_idx] = mixing_props[clone_idx] * nbinom.pmf(y, phi, prob)
        logcoeffs[:, clone_idx] =  np.log(mixing_props[clone_idx]) + nbinom.logpmf(y, phi, prob)

    # marginal = np.sum(np.log(np.sum(coeffs, axis = 1)))
    # responsibilities = coeffs / np.sum(coeffs, axis = 1)[:, np.newaxis]
    
    marginal = np.sum(logsumexp(logcoeffs, axis=1))    
    responsibilities = np.exp(logcoeffs - logsumexp(logcoeffs, axis=1)[:, np.newaxis])            
    
    return responsibilities, marginal

# fixed dispersion
class NBin_mixture_Q_fixed_phi(GenericLikelihoodModel):
    def __init__(self, endog, coeffs, total_reads, phi, seed=0, **kwds):
        super(NBin_mixture_Q_fixed_phi, self).__init__(endog, exog = np.zeros(endog.shape), **kwds)
        self.coeffs = coeffs
        self.total_reads = total_reads
        self.nclones = coeffs.shape[1]
        self.seed = seed
        self.phi = phi

    def nloglikeobs(self, params):
        # phi = params[-1]
        cn = params
        ll = _ll_mixture_nb(self.endog, self.coeffs, cn, self.total_reads, self.phi)
        return -ll

    def fit(self, start_params=None, maxiter=10000, maxfun=5000, **kwds):
        for clone_idx in range(self.nclones):
            self.exog_names.append(f'cn_clone_{clone_idx}')

        if start_params is None:
            # Reasonable starting values
            cn_guess = np.arange(self.nclones) + 1
            start_params = cn_guess
        
        return super(NBin_mixture_Q_fixed_phi, self).fit(start_params=start_params,
                                               maxiter=maxiter, maxfun=maxfun,
                                               **kwds)
    
def mixed_NB_EM_fixed_dispersion(y, total_reads, nclones=1, phi=100, cn_max=8, maxiter=10, seed=0, tol = 1e-6):
    
    np.random.seed(seed)
    nobservations = len(y)
    
    # initial guess
    mixing_props = [1/nclones]*nclones
    cn = np.arange(nclones) + 1
    # cn = np.random.randint(cn_max - 1, size=nclones) + 1

    relative_marginal_gain = np.inf
    old_marginal = np.inf
    iter_count = 0
    em_data = []
    
    while (iter_count < maxiter) & (relative_marginal_gain >= tol):
    # while (iter_count < maxiter):
        
        # E-step
        coeffs, new_marginal = get_responsibilities_and_marginal(y, mixing_props, cn, total_reads, phi)
        
        print(new_marginal, mixing_props)
        
        # M-step
        new_mixing_props = coeffs.sum(axis=0) / nobservations

        curr_model = NBin_mixture_Q_fixed_phi(y, coeffs, total_reads, phi, seed)
        # res = curr_model.fit(disp=0)
        
        res = curr_model.fit(start_params=cn, disp=0, xtol=1e-9, ftol=1e-9)
        
        # optimization solution
        params = res.params
        new_cn = params[:nclones]
        
        # summary
        inner_iter = 0
        # inner_iter = res.mle_retvals['iterations']
        inner_converged = res.mle_retvals['converged']
        fopt = res.mle_retvals['fopt']
        fcalls = res.mle_retvals['fcalls']
        
        if iter_count > 0:
            relative_marginal_gain = (new_marginal - old_marginal) / np.abs(old_marginal)
            
        em_data.append([iter_count, old_marginal, new_marginal, relative_marginal_gain, fopt, fcalls, inner_iter, inner_converged])
        
        if (relative_marginal_gain > 0) | (np.abs(new_marginal) == np.inf) | (np.abs(old_marginal) == np.inf):
            old_marginal = new_marginal
            cn = new_cn
            mixing_props = new_mixing_props
            
            iter_count += 1
            if np.isnan(relative_marginal_gain):
                relative_marginal_gain = np.inf
        else:
            break
        
    df_EM = pd.DataFrame(em_data, columns = ['iterations', 'old_marginal', 'marginal', 'relative_marginal_gain', 'fopt', 'fcalls', 'inner_iter', 'inner_converged'])
    
    return mixing_props, cn, phi, df_EM


# unknown dispersion
def _ll_mixture_nb(y, coeffs, cn, total_reads, phi):
    ll = 0
    nclones = len(cn)
    for clone_idx in range(nclones):
        mu = cn[clone_idx] * total_reads
        prob = phi / (phi + mu)
        ll += coeffs[:, clone_idx] * nbinom.logpmf(y, phi, prob)
    return ll

class NBin_mixture_Q(GenericLikelihoodModel):
    def __init__(self, endog, coeffs, total_reads, seed=0, **kwds):
        super(NBin_mixture_Q, self).__init__(endog, exog = np.zeros(endog.shape), **kwds)
        self.coeffs = coeffs
        self.total_reads = total_reads
        self.nclones = coeffs.shape[1]
        self.seed = seed

    def nloglikeobs(self, params):
        phi = params[-1]
        cn = params[:-1]
        ll = _ll_mixture_nb(self.endog, self.coeffs, cn, self.total_reads, phi)
        return -ll

    def fit(self, start_params=None, maxiter=10000, maxfun=5000, **kwds):
        for clone_idx in range(self.nclones):
            self.exog_names.append(f'cn_clone_{clone_idx}')
        self.exog_names.append('phi')

        if start_params is None:
            # Reasonable starting values
            cn_guess = np.arange(self.nclones) + 1
            phi_guess = 100
            
            start_params = np.append(cn_guess, phi_guess)
        
        return super(NBin_mixture_Q, self).fit(start_params=start_params,
                                               maxiter=maxiter, maxfun=maxfun,
                                               **kwds)

def mixed_NB_EM(y, total_reads, nclones=1, maxiter=10, seed=0, cn_max=8, tol = 1e-6):
    
    np.random.seed(seed)
    nobservations = len(y)
    
    # initial guess
    mixing_props = [1/nclones]*nclones
    cn = np.arange(nclones) + 1
    # cn = np.random.randint(cn_max - 1, size=nclones) + 1
    phi = 10
    relative_marginal_gain = np.inf
    old_marginal = -np.inf
    iter_count = 0
    em_data = []
    
    while (iter_count < maxiter) & (relative_marginal_gain >= tol):
    # while (iter_count < maxiter):
        
        # E-step
        coeffs, old_marginal = get_responsibilities_and_marginal(y, mixing_props, cn, total_reads, phi)
        
        # print(new_marginal, mixing_props)
        
        # M-step
        new_mixing_props = coeffs.sum(axis=0) / nobservations
        curr_model = NBin_mixture_Q(y, coeffs, total_reads, seed)
        # res = curr_model.fit(disp=0)
        
        res = curr_model.fit(start_params=np.append(cn, phi), disp=1)
        
        # optimization solution
        params = res.params
        new_cn = params[:nclones]
        new_phi = params[-1]
        
        # summary
        inner_iter = 0
        # inner_iter = res.mle_retvals['iterations']
        inner_converged = res.mle_retvals['converged']
        fopt = res.mle_retvals['fopt']
        fcalls = res.mle_retvals['fcalls']
        
        _, new_marginal = get_responsibilities_and_marginal(y, new_mixing_props, new_cn, total_reads, new_phi)        
        
        if iter_count > 0:
            # print(new_marginal, old_marginal, sep='\t')
            relative_marginal_gain = (new_marginal - old_marginal) / np.abs(old_marginal)
            
        em_data.append([iter_count, old_marginal, new_marginal, relative_marginal_gain, fopt, fcalls, inner_iter, inner_converged])
        
        if (relative_marginal_gain > 0) | (np.abs(new_marginal) == np.inf) | (np.abs(old_marginal) == np.inf):
            # old_marginal = new_marginal
            cn = new_cn
            phi = new_phi
            mixing_props = new_mixing_props
            
            iter_count += 1
            if np.isnan(relative_marginal_gain):
                relative_marginal_gain = np.inf
        else:
            break
        
    df_EM = pd.DataFrame(em_data, columns = ['iterations', 'old_marginal', 'marginal', 'relative_marginal_gain', 'fopt', 'fcalls', 'inner_iter', 'inner_converged'])
    
    return mixing_props, cn, phi, df_EM


# def main(args):
    

# def str2bool(v):
#     if isinstance(v, bool):
#         return v
#     if v.lower() in ('yes', 'true', 't', 'y', '1'):
#         return True
#     elif v.lower() in ('no', 'false', 'f', 'n', '0'):
#         return False
#     else:
#         raise argparse.ArgumentTypeError('Boolean value expected.')

# if __name__ == "__main__":
#     parser = argparse.ArgumentParser()
#     parser.add_argument('--sample', type=str, help='loom file name')
#     parser.add_argument('--dir', type=int, help='read depth threshold [10]', default=10)
#     parser.add_argument('--vaf', type=float, help='variant allele frequency threshold [0.33]', default=0.33)
#     parser.add_argument('--presence', type=float, help='prevalence threshold [0.2]', default=0.2)
#     parser.add_argument('--amplicon', type=int, help='amplicon depth threshold [0]', default=0) 
#     parser.add_argument('--hetero', type=float, help='vaf threshold to call a mutation heterozygous [0.2]', default=0.2)
#     parser.add_argument('--homo', type=float, help='vaf threshold to call a mutation homozygous [0.8]', default=0.8)
    
#     parser.add_argument('--prefix', type=str, help='prefix for output files of plotting')

#     args = parser.parse_args(None if sys.argv[1:] else ['-h'])

#     main(args)
