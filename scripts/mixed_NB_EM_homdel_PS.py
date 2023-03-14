import numpy as np
import pandas as pd
import argparse
import sys

from scipy.stats import nbinom
from scipy.special import logsumexp

# from IPython import embed

# @HZ: 01-23-2023:
# PS homdel method
zero_cn_mean = 0.1

def get_responsibilities_and_marginal_panel_level(df_observed_read_counts, df_amplicons, cell_total_reads, genelist, mixing_props, cn_profiles, homdel_profiles):
    '''
    E-step: get responsibilities of each cluster for each cell, assuming all other parameters are known. Also get marginal which is the total mixing props of each cluster in the population.

    Parameters
    ----------
    df_observed_read_counts : pandas.DataFrame
        Observed read counts for each amplicon in each cell.
    df_amplicons : pandas.DataFrame
        Amplicon information.
    cell_total_reads : pandas.Series
        Total read counts for each cell.
    genelist : list
        List of genes.
    mixing_props : numpy.array
        Known mixing proportions for each cluster.
    cn_profiles : numpy.array
        Known copy number profiles for each cluster.
    homdel_profiles : numpy.array
        Known homdel profiles for each cluster.

    Returns:
    --------
    responsibilities: ncells x nclones
        Responsibilities of each cluster for each cell.

    marginal: nclones x 1
        Marginal proportion of each cluster in the population.
    '''
    ncells = len(df_observed_read_counts)
    nclones = cn_profiles.shape[0]
    ngenes = cn_profiles.shape[1]
    namplicons = df_amplicons.shape[0] # M
    amplicon_list = list(df_amplicons.index)
    print(f'DEBUG: length of amplicon_list {len(amplicon_list)}')
    
    expanded_cn_profile = np.zeros((nclones, namplicons))
    for amplicon_idx, amplicon in enumerate(amplicon_list):
        curr_gene = df_amplicons.loc[amplicon]['gene']
        gene_idx = genelist.index(curr_gene)
        expanded_cn_profile[:, amplicon_idx] = cn_profiles[:, gene_idx]    
        for clone_idx in range(nclones):
            if cn_profiles[clone_idx, gene_idx] == 0:
                expanded_cn_profile[clone_idx, amplicon_idx] = zero_cn_mean
            elif homdel_profiles[clone_idx, amplicon_idx] == 1:
                expanded_cn_profile[clone_idx, amplicon_idx] = zero_cn_mean
            else:
                expanded_cn_profile[clone_idx, amplicon_idx] = cn_profiles[clone_idx, gene_idx]    
    
    logcoeffs = np.zeros((ncells, nclones))
    
    for clone_idx in range(nclones):
        mu = expanded_cn_profile[clone_idx, :] * cell_total_reads.values[:, np.newaxis] * df_amplicons['amplicon_factor'].values[np.newaxis, :]
        phi_matrix = np.array([1]*ncells)[:, np.newaxis] * df_amplicons['phi'].values[np.newaxis, :]

        prob = phi_matrix / (phi_matrix + mu)
        coeff_ij = nbinom.pmf(df_observed_read_counts.values, phi_matrix, prob)
        
        logcoeffs[:, clone_idx] =  np.log(mixing_props[clone_idx]) + np.sum(np.log(coeff_ij), axis=1)

    marginal = np.sum(logsumexp(logcoeffs, axis=1))
    responsibilities = np.exp(logcoeffs - logsumexp(logcoeffs, axis=1)[:, np.newaxis])
    
    return responsibilities, marginal

def get_optimal_cn_profile(df_observed_read_counts, df_amplicons, cell_total_reads, genelist, responsibilities, maxcn = 8):
    '''
    M-step: given known responsibilities, calculate the optimal copy number profile for each cluster.
    '''
    ncells = len(df_observed_read_counts)
    nclones = responsibilities.shape[1]
    ngenes = len(genelist)
    namplicons = len(df_amplicons)

    # amplicon_to_index_dict = {amplicon_list[idx]: idx for idx in range(namplicons)}

    cn_profile = np.zeros((nclones, ngenes))
    homdel_profiles = pd.DataFrame(np.ones((nclones, namplicons)), index = range(nclones), columns = df_observed_read_counts.columns)
    
    for gene_idx, gene in enumerate(genelist):
        
        # iterate over each gene
        curr_amplicons = list(df_amplicons[df_amplicons['gene'] == gene].index)
        df_gene_observed_read_counts = df_observed_read_counts[curr_amplicons]
        df_gene_amplicons = df_amplicons[df_amplicons['gene'] == gene]

        for clone_idx in range(nclones - 1):
            clone_responsibilities = responsibilities[:, clone_idx]
            
            # ----- #HOMDEL get zero cn values -----
            mu_zero = zero_cn_mean * cell_total_reads.values[:, np.newaxis] * df_gene_amplicons['amplicon_factor'].values[np.newaxis, :]
            phi_matrix = np.array([1]*ncells)[:, np.newaxis] * df_gene_amplicons['phi'].values[np.newaxis, :]
            prob_zero = phi_matrix / (phi_matrix + mu_zero)
            psi_zero = nbinom.logpmf(df_gene_observed_read_counts.values, phi_matrix, prob_zero)
            homdel_coeffs = np.sum(clone_responsibilities[:, np.newaxis] * psi_zero, axis=0)

            cn_nb_lls_plus_per_amplicon_homdel_mask = pd.Series(range(1, maxcn + 1)).apply(
                lambda cn_i: _evaluate_coeff(
                    responsibilities[:, clone_idx],
                    df_gene_observed_read_counts,
                    df_gene_amplicons,
                    cell_total_reads, 
                    cn_i,
                    homdel_coeffs,
                    )
                ) # 1 x maxcn
            # embed()
            # sys.exit(1)
            best_cn = np.argmax(cn_nb_lls_plus_per_amplicon_homdel_mask.apply(lambda row_i: row_i[0])) + 1
            cn_profile[clone_idx, gene_idx] = best_cn
            homdel_profiles.loc[clone_idx, df_gene_amplicons.index] = cn_nb_lls_plus_per_amplicon_homdel_mask.iloc[best_cn-1][1]

            # max_coeff = -np.inf
            # for cn in range(maxcn + 1):
            #     # try out all the integer copy number values
            #     curr_coeff = _evaluate_coeff(
            #         responsibilities[:, clone_idx],
            #         df_gene_observed_read_counts,
            #         df_gene_amplicons,
            #         cell_total_reads, 
            #         cn
            #         )
                
            #     # print(curr_coeff, gene_idx, clone_idx, cn)
            #     if curr_coeff > max_coeff:
            #         max_coeff = curr_coeff
            #         cn_profile[clone_idx, gene_idx] = cn
            
            # print('updated', gene_idx, clone_idx, cn_profile[clone_idx, gene_idx])

    cn_profile[-1,:] = 2
    homdel_profiles.iloc[-1, :] = 0
    
    return cn_profile, homdel_profiles

def _evaluate_coeff(responsibilities, df_gene_observed_read_counts, df_gene_amplicons, cell_total_reads, cn, homdel_coeffs):
    '''
    M-step helper function: given known responsibilities, and a particular integer cn value, calculate the parameters for each negative binomial distribution.

    Parameters

    Returns:
    --------
    total_coeff: float
        The coefficient for the given cn value.
    per_amplicon_homdel_mask: np.array
        A mask of amplicons that are homozygously deleted.
    '''
    ll = 0
    ncells = len(df_gene_observed_read_counts)
    namplicons = len(df_gene_amplicons)
    
    mu = cn * cell_total_reads.values[:, np.newaxis] * df_gene_amplicons['amplicon_factor'].values[np.newaxis, :]
    phi_matrix = np.array([1]*ncells)[:, np.newaxis] * df_gene_amplicons['phi'].values[np.newaxis, :]
    prob = phi_matrix / (phi_matrix + mu)
    
    psi = nbinom.logpmf(df_gene_observed_read_counts.values, phi_matrix, prob)
    cn_coeffs = np.sum(responsibilities[:, np.newaxis] * psi, axis = 0)

    total_coeff = 0
    per_amplicon_homdel_mask = np.zeros(namplicons)
    for amplicon_idx in range(namplicons):
        if cn_coeffs[amplicon_idx] < homdel_coeffs[amplicon_idx]:
            per_amplicon_homdel_mask[amplicon_idx] = 1
            total_coeff += homdel_coeffs[amplicon_idx]
        else:
            total_coeff += cn_coeffs[amplicon_idx]

    return total_coeff, per_amplicon_homdel_mask    

def mixed_NB_EM_fixed_dispersion_panel_level(df_observed_read_counts, df_amplicons, cell_total_reads, genelist, nclones=1, cn_max=8, maxiter=20, seed=0, tol = 1e-6, predefined_cn_clone_info = None, predefined_cn_clone_props = None, init_guess_maxcn = None):
    '''
    EM algorithm for mixed negative binomial model with fixed dispersion.

    Parameters
    ----------
    df_observed_read_counts: pd.DataFrame

    ...
    predefined_cn_clone_info: np.array
        A predefined cn profile for each clone. If not None, the algorithm will use this as the initial guess.
    predefined_cn_clone_props: list
        A predefined mixing proportion for each clone. If not None, the algorithm will use this as the initial guess.
    init_guess_maxcn: int
        The maximum copy number to use for the initial guess. If None, it will use cn_max.

    '''
    df_amplicons = df_amplicons.loc[df_observed_read_counts.columns]
    
    np.random.seed(seed)
    ncells = len(df_observed_read_counts)
    ngenes = len(genelist)
    namplicons = len(df_amplicons)

    if init_guess_maxcn is None:
        init_guess_maxcn = cn_max
    elif init_guess_maxcn > cn_max:
        raise('[WARNING] init_guess_maxcn cannot be larger than cn_max, setting it to cn_max')
        init_guess_maxcn = cn_max
    else:
        pass

    # #FOR HOMDEL pre-defined initial guess
    if predefined_cn_clone_info is not None:
        cn_profiles = predefined_cn_clone_info
        # nclones = cn_profiles.shape[0]
        if predefined_cn_clone_props is not None:
            mixing_props = predefined_cn_clone_props
        else:
            mixing_props = [1/nclones]*nclones
        
        if nclones > cn_profiles.shape[0]:
            # if nclones greater than previous solution's nclones, add random clones, with equal mixing proportions
            cn_profiles = np.vstack([cn_profiles, np.random.randint(init_guess_maxcn - 1, size=(nclones - cn_profiles.shape[0], ngenes)) + 1])
            mixing_props = [ prop_i * cn_profiles.shape[0] / nclones for prop_i in mixing_props] + [(1-cn_profiles.shape[0]/nclones)/(nclones-cn_profiles.shape[0])]* (nclones-cn_profiles.shape[0])
        elif nclones < cn_profiles.shape[0]:
            nclones = cn_profiles.shape[0]
        
    else:
        # randomly initialize the CN profiles and mixing proportions
        # initialize the CN profiles with a diploid clone and N-1 random ploidy clones
        # initialize the mixing proportions with uniform probability
        cn_profiles = np.random.randint(cn_max - 1, size=(nclones - 1, ngenes)) + 1 # note to avoid initialize with a CN=0
        cn_profiles = np.vstack([cn_profiles, np.ones((1, ngenes))*2])
        mixing_props = [1/nclones]*nclones
    
    # #HOMDEL
    homdel_profiles = np.zeros((nclones, namplicons))
        

    print('='*20)
    print('----- initialized CN clone profiles: -----')
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
        # df_observed_read_counts, df_amplicons, cell_total_reads, genelist, amplicon_list, mixing_props, cn_profiles, homdel_profiles
        responsibilities, new_marginal = get_responsibilities_and_marginal_panel_level(
            df_observed_read_counts, df_amplicons, cell_total_reads, genelist, mixing_props, cn_profiles, homdel_profiles
            )
        print(f'number of unique clusters: {len(np.unique(responsibilities.argmax(axis=1)))}')
        
        # M-step
        new_mixing_props = responsibilities.sum(axis=0) / ncells
        
        new_cn_profiles, homdel_profiles = get_optimal_cn_profile(df_observed_read_counts, df_amplicons, cell_total_reads, genelist, responsibilities, cn_max)
        
        print('='*20)
        print(f'----- iter_count = {iter_count} updated cn_profiles -----')
        print(new_cn_profiles)
        print(f'----- homdel_profiles -----')
        print(homdel_profiles)
        print(f"[INFO] number of homdel amplicons: {(homdel_profiles==1).any(axis=0).sum()}")
        print(f"[INFO] genes involved: {df_amplicons.loc[homdel_profiles.columns[(homdel_profiles==1).any(axis=0)]]['gene'].values}")
        homdel_profiles = homdel_profiles.values # prepare for next iteration
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
            
            iter_count += 1
            if np.isnan(relative_marginal_gain):
                relative_marginal_gain = np.inf
        else:
            break
        
    df_EM = pd.DataFrame(em_data, columns = ['iterations', 'old_marginal', 'marginal', 'relative_marginal_gain'])
    
    return mixing_props, cn_profiles, homdel_profiles, df_EM


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
    df_observed_read_counts = df_tsv[curr_selected_amplicons] # select to converged amplicons
    cell_total_read_counts = df_tsv.sum(axis = 1)

    # parse predefined cn profiles
    if args.predefined_cn_clone_info is not None:
        df_predefined_cn_profiles = pd.read_csv(args.predefined_cn_clone_info, header=0, index_col=None)
        # -- the clone_profiles numpy array need to be clone_idx by gene_idx; and the genes need to be in the same order as the genelist
        df_wide_solution_clone_info = df_predefined_cn_profiles.pivot(index=['clone_idx','prop'], columns='gene', values='cn').reset_index()
        array_predefined_cn_profiles = df_wide_solution_clone_info.loc[:, genelist].values
        # -- the clone_props needs to be an ordered list of clone_prop of each CN clone
        array_clone_props = df_wide_solution_clone_info.loc[:, ['prop']].values
    else:
        array_predefined_cn_profiles = None
        array_clone_props = None

    init_guess_maxcn = args.maxcn

    # multiple restarts of EM @HZ- DEPRECATED?
    nrestarts = args.nrestarts
    max_marginal = -np.inf
    
    seed_list = np.random.permutation(np.arange(100))[:nrestarts]
    for restart_idx in range(nrestarts):
        inferred_mixing_props, inferred_cn_profiles, homdel_profiles,df_EM = mixed_NB_EM_fixed_dispersion_panel_level(
            df_observed_read_counts, 
            df_selected_amplicons, 
            cell_total_read_counts, 
            genelist, 
            nclones=nclones, 
            cn_max = args.maxcn, 
            seed=seed_list[restart_idx], predefined_cn_clone_info=array_predefined_cn_profiles, predefined_cn_clone_props=array_clone_props, init_guess_maxcn=init_guess_maxcn
            )
        
        curr_max_marginal = df_EM.iloc[-1]['marginal']
        
        if curr_max_marginal > max_marginal:
            final_df_EM = df_EM
            final_mixing_props = inferred_mixing_props
            final_cn_profiles = inferred_cn_profiles
            max_marginal = curr_max_marginal
        
    # compute bic and aic
    nparams = nclones + (nclones - 1) * ngenes
    nsamples = np.prod(df_observed_read_counts.values.shape)
    bic = nparams * np.log(nsamples) - 2 * max_marginal
    aic = 2 * nparams - 2 * max_marginal
    
    # ----- write the output -----
    prefix = args.prefix
    # sample = args.sample
    if cn_calling_mode == 'single-sample':
        df_result = pd.DataFrame([[sample_name, nclones, args.seed, max_marginal, bic, aic]],
                             columns = ['sample', 'nclones', 'seed', 'marginal', 'BIC', 'AIC'])
    else:
        df_result = pd.DataFrame([[cohort_name, nclones, args.seed, max_marginal, bic, aic]],
                             columns = ['cohort', 'nclones', 'seed', 'marginal', 'BIC', 'AIC'])
    df_result.to_csv(f'{prefix}.result.csv', index=False)
    
    with open(f'{prefix}.clone_info.csv', 'w') as out:
        out.write('clone_idx,gene,cn,prop\n')
        for clone_idx in range(nclones):
            for gene_idx, gene in enumerate(genelist):
                out.write(f'{clone_idx},{gene},{final_cn_profiles[clone_idx][gene_idx]},{final_mixing_props[clone_idx]}\n')
    
    pd.DataFrame(homdel_profiles, index = range(nclones), columns=df_observed_read_counts.columns).to_csv(f'{prefix}.homdel_profiles.csv', index=True)
    # df_clone_info.to_csv(f'{prefix}_clone_info.csv', index=False)
    
    final_df_EM.to_csv(f'{prefix}.EM_info.csv', index=False)
    
    # with open(f'{prefix}_genelist.txt', 'w') as out:
    #     for gene in genelist:
    #         out.write(f'{gene}\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--cn_calling_mode', type=str, choices= ['single-sample', 'cohort'], help='whether to do NB_EM on each sample, or a cohort')
    parser.add_argument('--sample_name', nargs="*", help='sample name; if cohort-level run, this needs to be a list of sample names, matching the order of tsvs.')
    parser.add_argument('--cohort_name', type=str, help='cohort name', default='')
    parser.add_argument('--readcounts', nargs="*", help='read count file(s); if cohort-level run, this should be a list of rc files for all samples in cohort.')
    # parser.add_argument('--amplicon', type=str, help='amplicon dataframe containing pre-trained amplicon factors')
    parser.add_argument('--nclones', type=int, help='number of clones to run the clustering', default=1)
    parser.add_argument('--amplicon_parameters_f', type=str, help='''
    amplicon parameters dataframe containing the following necessary columns: 
        - amplicon_ID (AMPL41099') in the first column to be used as index;
        - 'gene', the corresponding gene's Hugo_Symbol for each amplicon;
        - 'amplicon_factor' & 'alpha' & 'beta_zero' & 'beta_one' & 'method' & 'mean' & 'variance' trained NB parameters specific for each amplicon.
    ''')
    parser.add_argument('--nrestarts', type=int, help='number of restarts', default=1)
    parser.add_argument('--seed', type=int, help='seed', default=0)
    parser.add_argument('--maxcn', type=int, help='maximum possible copy number', default=8)   
    parser.add_argument('--predefined_cn_clone_info', type=str, help='for homdel training, can start from previous optimal solution. This is the `{prefix}-clone_info.csv` file output by the NB-EM_CN_caller without homdel. Necessary columns are `clone_idx`,`prop`,`gene`, `cn`. Note that if the nclones from the previous solution is smaller than the currently defined nclones, random clones will be appended with uniform distribution.')
    parser.add_argument('--init_guess_maxcn', type=int, help='limit the maximum possible copy number for the initial guess', default=3)
    parser.add_argument('--prefix', type=str, help='prefix for output files')

    args = parser.parse_args(None if sys.argv[1:] else ['-h'])

    main(args)
