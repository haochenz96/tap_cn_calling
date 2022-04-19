import numpy as np
import pandas as pd
import math
import argparse
import sys

import statsmodels
import statsmodels.api as sm
from statsmodels.base.model import GenericLikelihoodModel

from scipy.stats import nbinom

def get_responsibilities_and_marginal_gene_level(df_observed_read_counts, df_gene_amplicons, cell_total_reads, mixing_props, cn):    
    ncells = len(df_observed_read_counts)
    nclones = len(cn)

    coeffs = np.zeros((ncells, nclones))    
    
    for clone_idx in range(nclones):
        mu = cn[clone_idx] * cell_total_reads.values[:, np.newaxis] * df_gene_amplicons['amplicon_factor'].values[np.newaxis, :]
        phi_matrix = np.array([1]*ncells)[:, np.newaxis] * df_gene_amplicons['phi'].values[np.newaxis, :]

        prob = phi_matrix / (phi_matrix + mu)
        coeff_ij = nbinom.pmf(df_observed_read_counts.values, phi_matrix, prob)

        coeffs[:, clone_idx] =  mixing_props[clone_idx] * np.prod(coeff_ij, axis=1)    
        
    marginal = np.sum(np.log(np.sum(coeffs, axis = 1)))
    responsibilities = coeffs / np.sum(coeffs, axis = 1)[:, np.newaxis]        
    
    return responsibilities, marginal

def _ll_mixture_nb_gene_level(responsibilities, df_observed_read_counts, df_gene_amplicons, cell_total_reads, cn):
    ll = 0
    ncells = len(df_observed_read_counts)
    nclones = len(cn)
    for clone_idx in range(nclones):
        mu = cn[clone_idx] * cell_total_reads.values[:, np.newaxis] * df_gene_amplicons['amplicon_factor'].values[np.newaxis, :]
        phi_matrix = np.array([1]*ncells)[:, np.newaxis] * df_gene_amplicons['phi'].values[np.newaxis, :]

        prob = phi_matrix / (phi_matrix + mu)

        ll += responsibilities[:, clone_idx] * np.sum(nbinom.logpmf(df_observed_read_counts.values, phi_matrix, prob), axis=1)
    
    return ll

class NBin_mixture_Q_fixed_phi_gene_level(GenericLikelihoodModel):
    def __init__(self, endog, responsibilities, df_observed_read_counts, df_gene_amplicons, cell_total_reads, seed=0, **kwds):
        super(NBin_mixture_Q_fixed_phi_gene_level, self).__init__(endog, exog = [0], **kwds)
        self.responsibilities = responsibilities
        self.cell_total_reads = cell_total_reads
        self.nclones = responsibilities.shape[1]
        self.df_observed_read_counts = df_observed_read_counts
        self.df_gene_amplicons = df_gene_amplicons
        self.seed = seed

    def nloglikeobs(self, params):
        # phi = params[-1]
        cn = params
        ll = _ll_mixture_nb_gene_level(self.responsibilities, self.df_observed_read_counts, self.df_gene_amplicons, self.cell_total_reads, cn)
        return -ll

    def fit(self, start_params=None, maxiter=10000, maxfun=5000, **kwds):
        for clone_idx in range(self.nclones):
            self.exog_names.append(f'cn_clone_{clone_idx}')

        if start_params is None:
            # Reasonable starting values
            cn_guess = np.arange(self.nclones) + 1
            start_params = cn_guess
        
        return super(NBin_mixture_Q_fixed_phi_gene_level, self).fit(start_params=start_params, maxiter=maxiter, maxfun=maxfun, **kwds)

def mixed_NB_EM_fixed_dispersion_gene_level(df_observed_read_counts, df_gene_amplicons, cell_total_reads, nclones=1, cn_max=8, maxiter=20, seed=0, tol = 1e-6):
    
    np.random.seed(seed)
    ncells = len(df_observed_read_counts)
    
    # initial guess
    mixing_props = [1/nclones]*nclones
    # cn = np.arange(nclones) + 1
    cn = np.random.randint(cn_max - 1, size=nclones) + 1

    relative_marginal_gain = np.inf
    old_marginal = np.inf
    iter_count = 0
    em_data = []
    
    while (iter_count < maxiter) & (relative_marginal_gain >= tol):
    # while (iter_count < maxiter):
        
        # E-step
        responsibilities, new_marginal = get_responsibilities_and_marginal_gene_level(df_observed_read_counts, df_gene_amplicons, cell_total_reads,
                                                                            mixing_props, cn)
        
        print(new_marginal, mixing_props)
        
        # M-step
        new_mixing_props = responsibilities.sum(axis=0) / ncells

        curr_model = NBin_mixture_Q_fixed_phi_gene_level([1], responsibilities, df_observed_read_counts, df_gene_amplicons, cell_total_reads, seed)
        # res = curr_model.fit(disp=0)
        
        res = curr_model.fit(start_params=cn, disp=0)
        
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
    
    return mixing_props, cn, df_EM

def main(args):
    
    np.random.seed(args.seed)
    
    # prepare the input data
    df_tsv = pd.read_csv(args.readcounts, sep='\t', index_col = 0)
    df_nb_residual = pd.read_csv(args.amplicon, index_col = 0)
    df_nb_residual['phi'] = 1 / df_nb_residual['alpha']
    curr_gene = args.gene
    nclones = args.nclones
    
    curr_selected_amplicons = list(df_nb_residual[df_nb_residual['gene'] == curr_gene].index)
    df_observed_read_counts = df_tsv[curr_selected_amplicons]
    cell_total_read_counts = df_tsv.sum(axis = 1)
    df_selected_amplicons = df_nb_residual.loc[curr_selected_amplicons][['amplicon_factor', 'phi']]    
    
    # multiple restarts of EM
    nrestarts = args.nrestarts
    max_marginal = -np.inf
    
    seed_list = np.random.permutation(np.arange(100))[:nrestarts]
    for restart_idx in range(nrestarts):
        inferred_mixing_props, inferred_cn, df_EM = mixed_NB_EM_fixed_dispersion_gene_level(df_observed_read_counts, df_selected_amplicons, cell_total_read_counts,
                                                                                            nclones=nclones, seed=seed_list[restart_idx], cn_max = args.maxcn)
        
        curr_max_marginal = df_EM.iloc[-1]['marginal']
        
        if curr_max_marginal > max_marginal:
            final_df_EM = df_EM
            final_mixing_props = inferred_mixing_props
            final_cn = inferred_cn
            max_marginal = curr_max_marginal
        
    # compute bic and aic
    nparams = 2 * nclones
    nsamples = np.prod(df_observed_read_counts.values.shape)
    bic = nparams * np.log(nsamples) - 2 * max_marginal
    aic = 2 * nparams - 2 * max_marginal
    
    # write the output
    prefix = args.prefix
    sample = args.sample
    
    df_result = pd.DataFrame([[sample, curr_gene, nclones, args.seed, max_marginal, bic, aic]],
                             columns = ['sample', 'gene', 'nclones', 'seed', 'marginal', 'BIC', 'AIC'])
    df_result.to_csv(f'{prefix}_result.csv', index=False)
    
    df_clone_info = pd.DataFrame({'clone_id': list(np.arange(nclones)),
                                  'cn': list(final_cn),
                                  'props': list(final_mixing_props)})
    df_clone_info.to_csv(f'{prefix}_clone_info.csv', index=False)
    
    final_df_EM.to_csv(f'{prefix}_EM_info.csv', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample', type=str, help='sample name')
    parser.add_argument('--readcounts', type=str, help='read count file')
    parser.add_argument('--amplicon', type=str, help='amplicon dataframe')
    parser.add_argument('--gene', type=str, help='gene of interest')
    parser.add_argument('--nclones', type=int, help='number of clones', default=1)
    parser.add_argument('--nrestarts', type=int, help='number of restarts', default=1)
    parser.add_argument('--seed', type=int, help='seed', default=0)
    parser.add_argument('--maxcn', type=int, help='maximum possible copy number', default=8)
        
    parser.add_argument('--prefix', type=str, help='prefix for output files')

    args = parser.parse_args(None if sys.argv[1:] else ['-h'])

    main(args)
