from cgi import print_arguments
from pprint import pprint
import numpy as np
import pandas as pd
import math
import argparse
import sys

import statsmodels
import statsmodels.api as sm
from statsmodels.base.model import GenericLikelihoodModel

from scipy.stats import nbinom
from scipy.special import logsumexp

def get_responsibilities_and_marginal_panel_level(df_observed_read_counts, df_amplicon_params, cell_total_reads, genelist, mixing_props, cn_profiles_df, pi, homdel_params):
    '''
    E-step: get responsibilities of each cluster for each cell, assuming all other parameters are known. Also get marginal which is the total mixing props of each cluster in the population.

    Parameters
    ----------
    df_observed_read_counts : pandas.DataFrame
        ncells x namplicons
        Observed read counts for each amplicon in each cell.

    df_amplicon_params : pandas.DataFrame
        namplicons x nparams
        Pretrained amplicon parameters.
    
    cell_total_reads : pandas.Series
        Total read counts for each cell.

    mixing_props : numpy.array
        nclones x 1
        Known mixing proportions for each cluster.

    cn_profiles_df : pandas.DataFrame
        nclones x ngenes
        Known copy number profiles for each cluster.
    
    pi : pandas.DataFrame
        nclones x namplicons
        Known amplicon-specific dropout/homdel probabilities for each cluster.
    
    homdel_nb_params : Tuple
        (n0, p0) for NB distribution of read counts for a homdel/dropout amplicon.

    Returns:
    --------
    responsibilities: ncells x nclones
        Responsibilities of each cluster for each cell.

    marginal: nclones x 1
        Marginal proportion of each cluster in the population.
    '''
    ncells = len(df_observed_read_counts) # N
    nclones = cn_profiles_df.shape[0] # K
    namplicons = df_amplicon_params.shape[0] # M

    # convert:
    # gene level read counts (K x ngenes) 
    # to
    # amplicon level read counts (K x M)
    cn_profiles_df = pd.DataFrame(cn_profiles_df, index = range(nclones), columns = genelist)
    expanded_cn_profile = pd.DataFrame(
        index = range(nclones), 
        columns = df_observed_read_counts.columns
        )
    # embed()
    expanded_cn_profile = expanded_cn_profile.apply(
        lambda x: cn_profiles_df[df_amplicon_params.loc[x.name, 'gene']],
    )
    
    n0, p0 = homdel_params

    # logcoeffs = np.zeros((ncells, nclones)) # N x K

    # --- calculate likelihood through matrix operation --- 
    mu = (expanded_cn_profile.values.flatten() * cell_total_reads.values[:, np.newaxis]).reshape(ncells, nclones, namplicons).transpose(1, 0, 2) * df_amplicon_params['amplicon_factor'].values[np.newaxis, :] # [K x N x M]

    phi_matrix = np.array([1]*ncells)[:, np.newaxis] * df_amplicon_params['phi'].values[np.newaxis, :] # [N x M]
    prob = phi_matrix / (phi_matrix + mu) # [K x N x M]
    
    # coeff = pi * nb_likelihood + (1 - pi) * beta_likelihood
    coeff_ijk = nbinom.pmf(df_observed_read_counts.values, phi_matrix, prob) * pi.values[:, None, :] + nbinom.pmf(df_observed_read_counts.values, n0, p0) * (1 - pi.values)[:, None, :] # this should guarantee that coeff_ijk is non-zero, but still replace just in case
    
    coeff_ijk = np.where(coeff_ijk == 0, 1e-300, coeff_ijk) # @HZ: replace 0 with a small number

    logcoeffs = (np.log(mixing_props)[:, None] + np.sum(np.log(coeff_ijk), axis=2))

    # embed()
    marginal = np.sum(logsumexp(logcoeffs, axis=0))
    # np.exp(logcoeffs - logsumexp(logcoeffs, axis=1)[:, np.newaxis])  
    responsibilities = np.exp(logcoeffs - logsumexp(logcoeffs, axis=0))
    responsibilities = responsibilities.T
    print(np.sum(responsibilities, axis=0))
    
    return responsibilities, marginal

def get_optimal_cn_profile(df_observed_read_counts, df_amplicons, cell_total_reads, genelist, responsibilities, homdel_params, maxcn = 8):
    '''
    M-step: given known responsibilities, calculate the optimal copy number profile for each cluster.
    '''
    ncells = len(df_observed_read_counts)
    nclones = responsibilities.shape[1]
    ngenes = len(genelist)
    namplicons = len(df_amplicons)

    cn_profile = np.zeros((nclones, ngenes))
    pi = pd.DataFrame(np.ones((nclones, namplicons)), index = range(nclones), columns = df_observed_read_counts.columns)
    n0, p0 = homdel_params

    for gene_idx, gene in enumerate(genelist):
        # iterate over each gene
        curr_amplicons = list(df_amplicons[df_amplicons['gene'] == gene].index)
        df_gene_observed_read_counts = df_observed_read_counts[curr_amplicons]
        df_gene_amplicons = df_amplicons[df_amplicons['gene'] == gene]
        
        for clone_idx in range(nclones - 1):
            cn_nb_lls = pd.Series(range(1, maxcn + 1)).apply(
                lambda cn_i: _evaluate_coeff(
                    responsibilities[:, clone_idx],
                    df_gene_observed_read_counts,
                    df_gene_amplicons,
                    cell_total_reads, 
                    cn_i,
                    )
                ) # 1 x maxcn
            best_cn = np.argmax(cn_nb_lls) + 1
            cn_profile[clone_idx, gene_idx] = best_cn
            if best_cn == 1:
                mu = 1 * cell_total_reads.values[:, np.newaxis] * df_gene_amplicons['amplicon_factor'].values[np.newaxis, :]

                phi_matrix = np.array([1]*ncells)[:, np.newaxis] * df_gene_amplicons['phi'].values[np.newaxis, :]

                prob = phi_matrix / (phi_matrix + mu)
                psi = nbinom.logpmf(df_gene_observed_read_counts.values, phi_matrix, prob)
                amplicon_nb_lls = np.sum(psi, axis = 0)
                # print(f"amplicon_nb_lls: {amplicon_nb_lls}")
                # print(f"amplicon_nb_lls shape: {amplicon_nb_lls.shape}")
                # 1 x M (namplicons selected): likelihood of the NB distribution with cn = 1
                # print(df_gene_observed_read_counts)
                psi0 = nbinom.logpmf(df_gene_observed_read_counts.values, n0, p0)
                amplicon_homdel_lls = np.sum(psi0, axis = 0)
                # 1 x M (namplicons selected): likelihood of amplicon being homdel
                homdel_amps = np.array(curr_amplicons)[
                        np.where(amplicon_homdel_lls > amplicon_nb_lls)[0]
                        ]
                pi.loc[clone_idx, homdel_amps] = 0

    cn_profile[-1,:] = 2
    
    return cn_profile, pi

def _evaluate_coeff(responsibilities, df_gene_observed_read_counts, df_gene_amplicons, cell_total_reads, cn):
    '''
    M-step helper function: given known responsibilities, and a particular integer cn value, calculate the parameters for each negative binomial distribution.

    Parameters

    Returns:
    --------
    coeff: float
        The coefficient for the given cn value.
    '''
    ll = 0
    ncells = len(df_gene_observed_read_counts)
    
    mu = cn * cell_total_reads.values[:, np.newaxis] * df_gene_amplicons['amplicon_factor'].values[np.newaxis, :]
    phi_matrix = np.array([1]*ncells)[:, np.newaxis] * df_gene_amplicons['phi'].values[np.newaxis, :]
    prob = phi_matrix / (phi_matrix + mu)
    
    psi = nbinom.logpmf(df_gene_observed_read_counts.values, phi_matrix, prob)
    
    # print(psi.shape, responsibilities.shape, np.sum(psi, axis=1).shape, sum(responsibilities * np.sum(psi, axis=1)).shape)
    
    # print(np.sum(np.isnan(psi)), np.sum(np.isnan(responsibilities)))
    return sum(responsibilities * np.sum(psi, axis=1))

def mixed_NB_EM_fixed_dispersion_panel_level(df_observed_read_counts, df_amplicons, cell_total_reads, genelist, nclones=1, cn_max=8, maxiter=20, seed=0, tol = 1e-6):
    
    df_amplicons = df_amplicons.loc[df_observed_read_counts.columns]
    
    np.random.seed(seed)
    ncells = len(df_observed_read_counts)
    ngenes = len(genelist)
    namplicons = len(df_amplicons)
    
    # initial guess
    # initialize the CN profiles with a diploid clone and N-1 random ploidy clones
    # initialize the mixing proportions with uniform probability
    mixing_props = [1/nclones]*nclones
    cn_profiles = np.random.randint(cn_max - 1, size=(nclones - 1, ngenes)) + 1
    cn_profiles = np.vstack([cn_profiles, np.ones((1, ngenes))*2]) 
    pi = pd.DataFrame(
            np.random.choice([0, 1], size=(nclones, namplicons), p=[0.01, 0.99]), 
            # np.ones((nclones, namplicons)),
            index = range(nclones), 
            columns = df_observed_read_counts.columns
            )
    homdel_params = (1e-10, 0.01)
    print('='*20)
    print('----- initial guess -----')
    print(cn_profiles)
    print(' '*20)
    print(' '*20)
    relative_marginal_gain = np.inf
    old_marginal = np.inf
    iter_count = 0
    em_data = []
    
    while (iter_count < maxiter) & (relative_marginal_gain >= tol):
    # while (iter_count < maxiter):
        
        # E-step
        responsibilities, new_marginal = get_responsibilities_and_marginal_panel_level(
            df_observed_read_counts, df_amplicons, cell_total_reads, genelist, mixing_props, cn_profiles, pi, homdel_params
            )
        print(f'number of unique clusters: {len(np.unique(responsibilities.argmax(axis=1)))}')
        
        # M-step
        new_mixing_props = responsibilities.sum(axis=0) / ncells

        new_cn_profiles, new_pi = get_optimal_cn_profile(df_observed_read_counts, df_amplicons, cell_total_reads, genelist, responsibilities, homdel_params, cn_max)
        
        print('='*20)
        print(f'----- iter_count = {iter_count} updated cn_profiles -----')
        print(new_cn_profiles)
        print('-'*5)
        print(f'New marginal = {new_marginal}')
        print(f'Mixing props = {new_mixing_props}')
        print(' '*20)
        print(' '*20)

        # print(np.sum(np.isnan(responsibilities)))
        # print(responsibilities.shape)
        
        if iter_count > 0:
            relative_marginal_gain = (new_marginal - old_marginal) / np.abs(old_marginal)
            
        em_data.append([iter_count, old_marginal, new_marginal, relative_marginal_gain])
        
        if (relative_marginal_gain > 0) | (np.abs(new_marginal) == np.inf) | (np.abs(old_marginal) == np.inf):
            old_marginal = new_marginal
            cn_profiles = new_cn_profiles
            mixing_props = new_mixing_props
            pi = new_pi
            
            iter_count += 1
            if np.isnan(relative_marginal_gain):
                relative_marginal_gain = np.inf
        else:
            break
        
    df_EM = pd.DataFrame(em_data, columns = ['iterations', 'old_marginal', 'marginal', 'relative_marginal_gain'])
    
    return mixing_props, cn_profiles, pi, df_EM


def main(args):
    '''
    Main execution function
    '''
    np.random.seed(args.seed)
    
    # prepare the input data
    cn_calling_mode = args.cn_calling_mode
    if cn_calling_mode not in ['single-sample', 'cohort']:
        print('Error: cn_calling_mode must be either single-sample or cohort')
        exit(1)
    elif cn_calling_mode == 'single-sample':
        sample_name = args.sample_name
        if type(sample_name) != str:
            raise ValueError('for single-sample mode, sample_name must be a string')
            exit(1)

        df_tsv = pd.read_csv(args.readcounts, sep='\t', index_col = 0)
    else: # cn_calling_mode == ‘cohort’
        sample_names = args.sample_name
        cohort_name = args.cohort_name
        print(sample_names)
        print(len(args.readcounts))
        if type(sample_names) != list:
            raise ValueError('for cohort mode, sample_names must be a list')
            exit(1)
        elif len(sample_names) != len(args.readcounts):
            raise ValueError('for cohort mode, sample_names must be a list mapping one-to-one to readcount matricies')
            exit(1)

        df_tsv = pd.concat(
            [pd.read_csv(f, sep='\t', index_col = 0) for f in args.readcounts], 
            keys = sample_names,
            )
        
    df_selected_amplicons = pd.read_csv(args.amplicon_parameters_f, index_col = 0)
    num_total_amplicons = df_selected_amplicons.shape[0]
    df_selected_amplicons = df_selected_amplicons.loc[df_selected_amplicons['converged'] == True] # filter out unconverged amplicons
    num_converge_amplicons = df_selected_amplicons.shape[0]
    print(f'[INFO] using {num_converge_amplicons}/{num_total_amplicons} converged amplicons'
    )
    df_selected_amplicons['phi'] = 1 / df_selected_amplicons['alpha']
    nclones = args.nclones
    
    genelist = list(df_selected_amplicons['gene'].unique())
    ngenes = len(genelist)
    
    curr_selected_amplicons = list(df_selected_amplicons.index)
    df_observed_read_counts = df_tsv[curr_selected_amplicons]
    cell_total_read_counts = df_tsv.sum(axis = 1)
    
    # multiple restarts of EM
    nrestarts = args.nrestarts
    max_marginal = -np.inf
    
    seed_list = np.random.permutation(np.arange(100))[:nrestarts]
    for restart_idx in range(nrestarts):
        inferred_mixing_props, inferred_cn_profiles, inferred_pi, df_EM = mixed_NB_EM_fixed_dispersion_panel_level(df_observed_read_counts, df_selected_amplicons, cell_total_read_counts,
                                                                                                      genelist, nclones=nclones, cn_max = args.maxcn,
                                                                                                      seed=seed_list[restart_idx])
        
        curr_max_marginal = df_EM.iloc[-1]['marginal']
        
        if curr_max_marginal > max_marginal:
            final_df_EM = df_EM
            final_mixing_props = inferred_mixing_props
            final_cn_profiles = inferred_cn_profiles
            final_pi = inferred_pi
            max_marginal = curr_max_marginal
        
    # compute bic and aic
    # nparams = nclones + (nclones - 1) * ngenes 
    nparams = nclones + (nclones - 1) * ngenes * 2 # @HZ: this is doubled because of pi
    nsamples = np.prod(df_observed_read_counts.values.shape)
    bic = nparams * np.log(nsamples) - 2 * max_marginal
    aic = 2 * nparams - 2 * max_marginal
    
    # write the output
    prefix = args.prefix
    # sample = args.sample
    if cn_calling_mode == 'single-sample':
        df_result = pd.DataFrame([[sample_name, nclones, args.seed, max_marginal, bic, aic]],
                             columns = ['sample', 'nclones', 'seed', 'marginal', 'BIC', 'AIC'])
    else:
        df_result = pd.DataFrame([[cohort_name, nclones, args.seed, max_marginal, bic, aic]],
                             columns = ['cohort', 'nclones', 'seed', 'marginal', 'BIC', 'AIC'])
    df_result.to_csv(f'{prefix}_result.csv', index=False)
    
    with open(f'{prefix}_clone_info.csv', 'w') as out:
        out.write('clone_idx,gene,cn,prop\n')
        for clone_idx in range(nclones):
            for gene_idx, gene in enumerate(genelist):
                out.write(f'{clone_idx},{gene},{final_cn_profiles[clone_idx][gene_idx]},{final_mixing_props[clone_idx]}\n')
    inferred_pi.to_csv(
        f"{prefix}_pi.csv",
    )
    
    # df_clone_info = pd.DataFrame({'clone_id': list(np.arange(nclones)),
    #                               'cn': list(final_cn),
    #                               'props': list(final_mixing_props)})
    # df_clone_info.to_csv(f'{prefix}_clone_info.csv', index=False)
    
    final_df_EM.to_csv(f'{prefix}_EM_info.csv', index=False)
    
    with open(f'{prefix}_genelist.txt', 'w') as out:
        for gene in genelist:
            out.write(f'{gene}\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--cn_calling_mode', type=str, choices= ['single-sample', 'cohort'], help='whether to do NB_EM on each sample, or a cohort')
    parser.add_argument('--sample_name', nargs="*", help='sample name; if cohort-level run, this needs to be a list of sample names, matching the order of tsvs.')
    parser.add_argument('--cohort_name', type=str, help='cohort name', default='')
    parser.add_argument('--readcounts', nargs="*", help='read count file(s); if cohort-level run, this should be a list of rc files for all samples in cohort.')
    # parser.add_argument('--amplicon', type=str, help='amplicon dataframe containing pre-trained amplicon factors')
    parser.add_argument('--nclones', type=int, help='number of clones', default=1)
    parser.add_argument('--amplicon_parameters_f', type=str, help='''
    amplicon parameters dataframe containing the following necessary columns: 
        - amplicon_ID (AMPL41099') in the first column to be used as index;
        - 'gene', the corresponding gene's Hugo_Symbol for each amplicon;
        - 'amplicon_factor' & 'alpha' & 'beta_zero' & 'beta_one' & 'method' & 'mean' & 'variance' trained NB parameters specific for each amplicon.
    ''')
    parser.add_argument('--nrestarts', type=int, help='number of restarts', default=1)
    parser.add_argument('--seed', type=int, help='seed', default=0)
    parser.add_argument('--maxcn', type=int, help='maximum possible copy number', default=8)     
    parser.add_argument('--prefix', type=str, help='prefix for output files')

    args = parser.parse_args(None if sys.argv[1:] else ['-h'])

    main(args)
