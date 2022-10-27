# gather NB EM results

import argparse
import sys
from pathlib import Path
import shutil
import glob

import pandas as pd
import numpy as np
from scipy.stats import nbinom
from scipy.special import logsumexp

import seaborn as sns
import matplotlib.pyplot as plt
# from IPython import embed

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
    # embed()
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

def main(args):

    # ----- prepare inputs ----- 
    # <<< inputs >>> 
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
    else:
        sample_name = args.cohort_name
        sample_names = args.sample_name
        # print(sample_names)
        # print(len(args.readcounts))
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
    num_unique_amplicons = len(df_tsv.columns)
    df_cell_total_read_counts = df_tsv.sum(axis = 1)

    df_selected_amplicons = pd.read_csv(args.amplicon_parameters_f, index_col = 0)
    df_selected_amplicons = df_selected_amplicons.loc[df_selected_amplicons['converged'] == True] # filter out unconverged amplicons
    df_selected_amplicons['phi'] = 1 / df_selected_amplicons['alpha']

    nclones = args.nclones
    genelist = list(df_selected_amplicons['gene'].unique())
    ngenes = len(genelist)
    
    # subset to amplicons of interest
    curr_selected_amplicons = list(df_selected_amplicons.index)
    print(f'[INFO] number of selected amplicons: {len(curr_selected_amplicons)} / {num_unique_amplicons}')
    df_observed_read_counts = df_tsv[curr_selected_amplicons]

    # <<< intermediate outputs >>>
    inputs_dir = Path(args.inputs_dir)

    # <<< outputs >>>
    output_prefix = args.prefix
    outputs_dir = Path(output_prefix).parents[0]
    if not outputs_dir.exists():
        # automatically done with snakemake so should be redundant
        outputs_dir.mkdir(parents=True, exist_ok=True)

    # ----- summarize EM results -----
    # print(inputs_dir)
    # print(sample_name)
    # print(nclones)
    print( str(inputs_dir / f'{sample_name}_nclones={nclones}_seed=*_result.csv') )
    EM_results_fs = glob.glob(str(inputs_dir / f'{sample_name}_nclones={nclones}_seed=*_result.csv')) # note we are restricting to a particular nclones
    # print(EM_results_fs)
    EM_results_dfs = [pd.read_csv(f) for f in EM_results_fs]
    EM_summary_df = pd.concat(EM_results_dfs, ignore_index=True)

    # plot ranked BIC's
    bic_sorted = EM_summary_df['BIC'].sort_values(ascending = True).reset_index(drop = True)
    plt.figure(figsize=(7, 7))
    ax = sns.lineplot(x = bic_sorted.index, y = bic_sorted.values)
    ax.set_xlabel('Rank')
    ax.set_ylabel('BIC')
    plt.grid()
    plt.savefig(
        f'{output_prefix}_BIC_ranked.png', dpi=300
        )

    # identify the best solution from all seeds
    best_idx = EM_summary_df['BIC'].idxmin()
    seed = int(EM_summary_df.iloc[best_idx]['seed'])
    nclones = int(EM_summary_df.iloc[best_idx]['nclones'])
    print(f'[INFO] best solution: seed = {seed}, nclones = {nclones}')

    # copy over the best solution
    shutil.copyfile(
        EM_results_fs[best_idx], 
        f'{output_prefix}_solution-EM_info.csv'
        )
    solution_clone_info_f = (inputs_dir / f'{sample_name}_nclones={nclones}_seed={seed}_clone_info.csv')
    solution_pi = (inputs_dir / f'{sample_name}_nclones={nclones}_seed={seed}_pi.csv')
    if not solution_clone_info_f.is_file() or not solution_pi.is_file():
        print('[ERROR} clone info file or homdel file not found')
    else:
        shutil.copyfile(
            f'{inputs_dir}/{sample_name}_nclones={nclones}_seed={seed}_clone_info.csv', 
            f'{output_prefix}_solution-clone_info.csv'
            )
        shutil.copyfile(
            f'{inputs_dir}/{sample_name}_nclones={nclones}_seed={seed}_pi.csv', 
            f'{output_prefix}_solution-pi.csv'
            )
    
    # ----- get sc-amplicon ploidy with best solution
    # << prepare inputs >>
    df_solution_clone_info = pd.read_csv(solution_clone_info_f)
    df_solution_pi = pd.read_csv(solution_pi, index_col = 0)
    # -- the clone_profiles numpy array need to be clone_idx by gene_idx; and the genes need to be in the same order as the genelist
    df_wide_solution_clone_info = df_solution_clone_info.pivot(index=['clone_idx','prop'], columns='gene', values='cn').reset_index()
    array_solution_clone_cn_profiles = df_wide_solution_clone_info.loc[:, genelist].values
    # -- the clone_props needs to be a dict mapping from clone_idx to clone_prop
    dict_clone_props = df_wide_solution_clone_info.loc[:, ['clone_idx', 'prop']].set_index('clone_idx')['prop'].to_dict()
    
    # get responsibilities
    # def get_responsibilities_and_marginal_panel_level(df_observed_read_counts, df_amplicon_params, cell_total_reads, genelist, mixing_props, cn_profiles_df, pi, homdel_params):
    responsibilities, _ = get_responsibilities_and_marginal_panel_level(
        df_observed_read_counts, 
        df_selected_amplicons,
        df_cell_total_read_counts, 
        genelist, 
        # dict_clone_props, 
        list(dict_clone_props.values()),
        array_solution_clone_cn_profiles,
        df_solution_pi,
        homdel_params = (1e-10, 0.1)
        )
    cell_assignments = np.argmax(responsibilities, axis=1)
    df_cell_assignments = pd.DataFrame(cell_assignments, index = df_cell_total_read_counts.index, columns = ['clone_id'])
    df_cell_assignments.to_csv(f'{output_prefix}_solution-cell_assignments.csv', index=True, header=True)

    # final output: single-cell x amplicon 
    df_sc_amplicon_ploidy = pd.DataFrame(
        index = df_observed_read_counts.index, 
        columns = df_observed_read_counts.columns
        )
    # embed()
    for amplicon_idx, amplicon in enumerate(df_observed_read_counts):
        curr_gene = df_selected_amplicons.loc[amplicon]['gene']
        gene_idx = genelist.index(curr_gene)
        df_sc_amplicon_ploidy.iloc[:, amplicon_idx] = array_solution_clone_cn_profiles[
            df_cell_assignments['clone_id'].values, 
            gene_idx
            ] * df_solution_pi.loc[df_cell_assignments['clone_id'].values, amplicon].values # homdel
    df_sc_amplicon_ploidy.to_csv(f'{output_prefix}_solution-sc_amplicon_ploidy.csv', index=True, header=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--cn_calling_mode', type=str, choices= ['single-sample', 'cohort'], help='whether to do NB_EM on each sample, or a cohort')
    parser.add_argument('--sample_name', nargs="*", help='sample name; if cohort-level run, this needs to be a list of sample names, matching the order of tsvs.')
    parser.add_argument('--cohort_name', type=str, help='cohort name', default='')
    parser.add_argument('--readcounts', nargs="*", help='read count file(s); if cohort-level run, this should be a list of rc files for all samples in cohort.')
    parser.add_argument('--nclones', type=int, help='number of clones', default=1)
    parser.add_argument('--amplicon_parameters_f', type=str, help='''
    amplicon parameters dataframe containing the following necessary columns: 
        - amplicon_ID (AMPL41099') in the first column to be used as index;
        - 'gene', the corresponding gene's Hugo_Symbol for each amplicon;
        - 'amplicon_factor' & 'alpha' & 'beta_zero' & 'beta_one' & 'method' & 'mean' & 'variance' trained NB parameters specific for each amplicon.
    ''')
    parser.add_argument('--inputs_dir', type=str, help='directory containing the EM results')
    # parser.add_argument('--nclones', type=int, help='number of clones', default=1)
    # parser.add_argument('--nrestarts', type=int, help='number of restarts', default=1)
    # parser.add_argument('--seed', type=int, help='seed', default=0)
    # parser.add_argument('--maxcn', type=int, help='maximum possible copy number', default=8)
    parser.add_argument('--prefix', type=str, help='prefix for output files')

    args = parser.parse_args(None if sys.argv[1:] else ['-h'])

    main(args)