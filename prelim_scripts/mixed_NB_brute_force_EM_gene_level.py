import numpy as np
import pandas as pd
import math
import argparse
import sys
import itertools

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

def get_cn_string(cn_list, max_cn):
    cn_string = ['0']*(max_cn+1)
    for cn in cn_list:
        cn_string[cn] = '1'    
    return ''.join(cn_string)

def main(args):
    
    seed = args.seed
    np.random.seed(seed)
    
    # prepare the input data
    sample = args.sample    
    df_tsv = pd.read_csv(args.readcounts, sep='\t', index_col = 0)
    df_nb_residual = pd.read_csv(args.amplicon, index_col = 0)
    df_nb_residual['phi'] = 1 / df_nb_residual['alpha']
    curr_gene = args.gene
    
    curr_selected_amplicons = list(df_nb_residual[df_nb_residual['gene'] == curr_gene].index)
    df_observed_read_counts = df_tsv[curr_selected_amplicons]
    cell_total_read_counts = df_tsv.sum(axis = 1)
    df_selected_amplicons = df_nb_residual.loc[curr_selected_amplicons][['amplicon_factor', 'phi']]
    
    
    # enumerate all possible combination of copy number clones
    max_cn = args.maxcn
    max_nclones = args.maxclones
    ncells = len(df_observed_read_counts)

    # multiple restarts of EM
    nsamples = np.prod(df_observed_read_counts.values.shape)
    nrestarts = args.nrestarts
    seed_list = np.random.permutation(np.arange(100))[:nrestarts]
    tol = 1e-6
    maxiter = 20
    cn_list_idx = 0
    result_data = []
    
    for nclones in range(1, max_nclones+1):
        for cn_list in itertools.combinations(np.arange(max_cn+1), nclones):
            prop_string = ''
            cn_string = get_cn_string(cn_list, max_cn)
            final_mixing_props = [None]
            for restart_idx in range(nrestarts):
                np.random.seed(seed_list[restart_idx])

                # initial guess
                mixing_props = np.random.dirichlet([1]*nclones)

                relative_marginal_gain = np.inf
                old_marginal = -np.inf
                iter_count = 0
                em_data = []
                max_marginal = -np.inf

                while (iter_count < maxiter) & (relative_marginal_gain >= tol):
                    responsibilities, new_marginal = get_responsibilities_and_marginal_gene_level(df_observed_read_counts, df_selected_amplicons, cell_total_read_counts,
                                                                                                  mixing_props, cn_list)
                    new_mixing_props = responsibilities.sum(axis=0) / ncells

                    if iter_count > 0:
                        relative_marginal_gain = (new_marginal - old_marginal) / np.abs(old_marginal)

                    em_data.append([sample, curr_gene, nclones,
                                    iter_count, old_marginal, new_marginal, relative_marginal_gain])

                    if (relative_marginal_gain > 0) | (np.abs(new_marginal) == np.inf) | (np.abs(old_marginal) == np.inf):
                        old_marginal = new_marginal
                        mixing_props = new_mixing_props

                        iter_count += 1
                        if np.isnan(relative_marginal_gain):
                            relative_marginal_gain = np.inf
                    else:
                        break

                df_EM = pd.DataFrame(em_data, columns = ['sample', 'gene', 'nclones',
                                                         'iterations', 'old_marginal', 'marginal', 'relative_marginal_gain'])

                curr_max_marginal = df_EM.iloc[-1]['marginal']

                if curr_max_marginal > max_marginal:
                    final_df_EM = df_EM
                    final_mixing_props = mixing_props
                    final_cn = cn_list
                    max_marginal = curr_max_marginal

            prop_string = '|'.join(map(str, final_mixing_props))

            bic = 2 * nclones * np.log(nsamples) - 2 * max_marginal
            aic = 2 * 2 * nclones - 2 * max_marginal
            result_data.append([sample, curr_gene, nclones, seed, cn_string, prop_string, max_marginal, bic, aic])
                
    # write the output
    prefix = args.prefix
    
    df_result = pd.DataFrame(result_data,
                             columns = ['sample', 'gene', 'nclones', 'seed',
                                        'cn_string', 'prop_string', 'marginal', 'BIC', 'AIC'])
    df_result.to_csv(f'{prefix}_result.csv', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample', type=str, help='sample name')
    parser.add_argument('--readcounts', type=str, help='read count file')
    parser.add_argument('--amplicon', type=str, help='amplicon dataframe')
    parser.add_argument('--gene', type=str, help='gene of interest')
    parser.add_argument('--maxcn', type=int, help='maximum value of copy number [5]', default=5)
    parser.add_argument('--maxclones', type=int, help='maximum number of clones [6]', default=6)
    parser.add_argument('--nrestarts', type=int, help='number of restarts', default=10)
    parser.add_argument('--seed', type=int, help='seed', default=0)
        
    parser.add_argument('--prefix', type=str, help='prefix for output files')

    args = parser.parse_args(None if sys.argv[1:] else ['-h'])

    main(args)
