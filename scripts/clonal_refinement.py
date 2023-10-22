# 2023-08 Akhil
# 2023-09 Haochen

"""
This script refines selected FALCON solutions
by the following steps:
"""

import os, sys, argparse
import pandas as pd
import numpy as np
from copy import deepcopy
from collections import defaultdict
import mosaic.io as mio
from pathlib import Path

import warnings 
warnings.filterwarnings("ignore", category=RuntimeWarning) # disable division by zero warnings
warnings.filterwarnings("ignore", message=".*The 'nopython' keyword.*") # disable numba warnings

def __binarize_NGT_matrix(NGT_df):
    """
    takes as input a numbered genotype matrix, where the last two columns are sample ID and cell barcode, and binarize the matrix by setting HET to 0 and HOM to 1
    """
    NGT_df[NGT_df.columns[:-2]] = np.where(NGT_df[NGT_df.columns[:-2]] == 1, 0, NGT_df[NGT_df.columns[:-2]])
    NGT_df[NGT_df.columns[:-2]] = np.where(NGT_df[NGT_df.columns[:-2]] == 2, 1, NGT_df[NGT_df.columns[:-2]])

    return NGT_df

def profile_distance(medians, cell_profile):
    distances = {}
    for median_clone, median_profile in medians.items():
        mask = median_profile != 3
        mask2 = cell_profile != 3
        mask = np.logical_and(mask, mask2)
        series1 = median_profile[mask]
        series2 = cell_profile[mask]
        absolute_sum = np.sqrt(((series1 - series2)**2).sum())
        distances[median_clone] = absolute_sum/series2.shape[0]

    return distances


def median_distance(medians, clone1, clone2):
    """
    get L2 distance between two SNV profiles provided in `medians`. Omit SNVs with missing values (3) in either profile.
    """
    median_profile1 = medians[clone1]
    median_profile2 = medians[clone2]
    mask = (median_profile1 != 3 & median_profile2 != 3)
    absolute_sum = np.sqrt(((median_profile1[mask] - median_profile2[mask])**2).sum())
    distance = absolute_sum/median_profile2[mask].shape[0]

    return distance

def generate_median_cluster_profiles(NGT_df):
    median_clusters = {}
    median_clusters_all = {}

    for clone in NGT_df['clone_id'].unique():
        NGT_curr_clone = NGT_df[NGT_df['clone_id'] == clone][NGT_df.columns[:-2]]
        NGT_curr_clone_mode = NGT_curr_clone.mode(axis=0)
        if NGT_curr_clone.shape[0] > 0:
            if NGT_curr_clone_mode.shape[0] > 1: 
            # if multiple modes exist, choose the first one
                median_clusters_all[clone] = NGT_curr_clone_mode.iloc[0].squeeze(axis=0)
            else:
                median_clusters_all[clone] = NGT_curr_clone_mode.squeeze(axis=0)

            if NGT_curr_clone.shape[0] > 0.05 * NGT_df.shape[0]:
                median_clusters[clone] = median_clusters_all[clone]

    return median_clusters, median_clusters_all

def copy_number_distance(df_tapestri_clones, clone1, clone2):
    c1_profile = np.array(df_tapestri_clones.loc[clone1].tolist())
    c2_profile = np.array(df_tapestri_clones.loc[clone2].tolist())
    absolute_sum = np.abs(c1_profile - c2_profile).sum()
    absolute_sum = absolute_sum/len(c2_profile)
    return absolute_sum

def combine_similar_final_clusters(NGT_df, cn_assignment_df, median_clusters, df_tapestri_clones):
    normal_idx = -1
    if len(df_tapestri_clones.index[df_tapestri_clones.eq(2.0).all(axis=1)].to_list()) > 0:
        normal_idx = df_tapestri_clones.index[df_tapestri_clones.eq(2.0).all(axis=1)].to_list()[0]
    for c1 in NGT_df['clone_id'].unique():
        for c2 in NGT_df['clone_id'].unique():
            if c1 != c2:
                if median_distance(median_clusters, c1, c2) <= 0.10:
                    if copy_number_distance(df_tapestri_clones, c1, c2) <= 0.40:
                        c1_size = NGT_curr_clone = NGT_df[NGT_df['clone_id'] == c1][NGT_df.columns[:-1]].shape[0]
                        c2_size = NGT_curr_clone = NGT_df[NGT_df['clone_id'] == c2][NGT_df.columns[:-1]].shape[0]
                            
                        if c1 != normal_idx and c2 != normal_idx:
                            if c1_size >= c2_size:
                                cn_assignment_df['clone_id'] = cn_assignment_df['clone_id'].replace(c2, c1)
                                NGT_df['clone_id'] = NGT_df['clone_id'].replace(c2, c1)
                            else:
                                cn_assignment_df['clone_id'] = cn_assignment_df['clone_id'].replace(c1, c2)
                                NGT_df['clone_id'] = NGT_df['clone_id'].replace(c1, c2)

                        else:
                            if c1 == normal_idx:
                                cn_assignment_df['clone_id'] = cn_assignment_df['clone_id'].replace(c2, c1)
                                NGT_df['clone_id'] = NGT_df['clone_id'].replace(c2, c1)
                            else:
                                cn_assignment_df['clone_id'] = cn_assignment_df['clone_id'].replace(c1, c2)
                                NGT_df['clone_id'] = NGT_df['clone_id'].replace(c1, c2)


    return NGT_df, cn_assignment_df, df_tapestri_clones

def merge_large_clusters(NGT_df, NGT_df_homvaf, cn_assignment_df, median_clusters,  count):
    for clone in NGT_df['clone_id'].unique():
        NGT_curr_clone = NGT_df[NGT_df['clone_id'] == clone][NGT_df.columns[:-1]]
        if NGT_curr_clone.shape[0] > 0.05 * NGT_df.shape[0]:
            for cell, cell_profile in NGT_curr_clone.iterrows():
                cell_barcode = cell_profile['cell_barcode']
                cell_profile = cell_profile.drop('cell_barcode')
                distances = profile_distance(median_clusters, cell_profile)
                distances = {k:v/distances[clone] for k,v in distances.items()}
                k,v = min(distances.items(), key=lambda x: x[1])
                if v < 0.85:
                    count += 1
                    cn_assignment_df.loc[(cn_assignment_df.index == cell) & (cn_assignment_df['cell_barcode'] == cell_barcode), 'clone_id'] = k
                    NGT_df_homvaf.loc[(NGT_df_homvaf.index == cell) & (NGT_df_homvaf['cell_barcode'] == cell_barcode), 'clone_id'] = k
                    NGT_df.loc[(NGT_df.index == cell) & (NGT_df['cell_barcode'] == cell_barcode), 'clone_id'] = k

    return NGT_df, cn_assignment_df, count

def dissolve_small_clusters(NGT_df, NGT_df_homvaf, cn_assignment_df, median_clusters,  count):
    for clone in NGT_df['clone_id'].unique():
        NGT_curr_clone = NGT_df[NGT_df['clone_id'] == clone][NGT_df.columns[:-1]]
        if NGT_curr_clone.shape[0] > 0 and NGT_curr_clone.shape[0] <= 0.05 * NGT_df.shape[0]:
            for cell, cell_profile in NGT_curr_clone.iterrows():
                cell_barcode = cell_profile['cell_barcode']
                cell_profile = cell_profile.drop('cell_barcode')
                distances = profile_distance(median_clusters, cell_profile)
                distances = {k:v for k,v in distances.items()}
                distances.pop(clone, None)
                k,v = min(distances.items(), key=lambda x: x[1])
                count += 1
                cn_assignment_df.loc[(cn_assignment_df.index == cell) & (cn_assignment_df['cell_barcode'] == cell_barcode), 'clone_id'] = k
                NGT_df_homvaf.loc[(NGT_df_homvaf.index == cell) & (NGT_df_homvaf['cell_barcode'] == cell_barcode), 'clone_id'] = k
                NGT_df.loc[(NGT_df.index == cell) & (NGT_df['cell_barcode'] == cell_barcode), 'clone_id'] = k

    return NGT_df, cn_assignment_df, count

def split_sample_clones(sample_name, NGT_df, NGT_df_homvaf, df_tapestri_clones, cn_assignment_df, added_clone_mut_dict):
    count = 0
    added_clone = False
    
    median_clusters, median_clusters_all = generate_median_cluster_profiles(NGT_df)
    NGT_df, cn_assignment_df, count = merge_large_clusters(NGT_df, NGT_df_homvaf, cn_assignment_df, median_clusters,  count)
    NGT_df, cn_assignment_df, count = dissolve_small_clusters(NGT_df, NGT_df_homvaf,  cn_assignment_df, median_clusters,  count)
    #split

    max_cn = cn_assignment_df['clone_id'].max()
    for clone in NGT_df['clone_id'].unique():
        candidates_not_found = False
        while(candidates_not_found == False):
            candidate_muts = {}
            NGT_curr_clone = NGT_df_homvaf[NGT_df_homvaf['clone_id'] == clone][NGT_df_homvaf.columns[:-1]]
            non_missing_NGT_curr_clone = NGT_curr_clone.replace(3, np.NaN)

            for mut in NGT_curr_clone.columns:
                if mut != 'cell_barcode':
                    mean = non_missing_NGT_curr_clone[mut].mean()
                    if mean < 0.75:
                        if 2 in NGT_curr_clone[mut].unique():
                            if NGT_curr_clone[mut].value_counts()[2.0] > 0.10 * NGT_df_homvaf.shape[0] and NGT_curr_clone[mut].value_counts()[3.0] < 0.30 * NGT_curr_clone.shape[0]:
                                candidate_muts[mut] = NGT_curr_clone[mut].value_counts()[2.0]/NGT_curr_clone.shape[0]

            candidate_muts = {k: v for k, v in sorted(candidate_muts.items(), key=lambda item: item[1], reverse=True)}
            if len(candidate_muts) < 2:
                candidates_not_found = 1
            else:
                added_clone = True
                max_cn = max_cn + 1
                candidate = list(candidate_muts.keys())[0]
                print(candidate, clone, sample_name)
                added_clone_mut_dict[max_cn] = candidate
                new_clone_profile = df_tapestri_clones.loc[clone].copy()
                df_tapestri_clones = pd.concat([df_tapestri_clones, pd.DataFrame([new_clone_profile])])
                #df_tapestri_clones = df_tapestri_clones.append(new_clone_profile)
                df_tapestri_clones.index.values[-1] = max_cn
                for cell, cell_profile in NGT_curr_clone.iterrows():
                    cell_barcode = cell_profile['cell_barcode']
                    cell_profile = cell_profile.drop('cell_barcode')
                    if cell_profile.loc[candidate] == 2.0:
                        count += 1
                        NGT_df.loc[(NGT_df.index == cell) & (NGT_df['cell_barcode'] == cell_barcode), 'clone_id'] = max_cn
                        NGT_df_homvaf.loc[(NGT_df_homvaf.index == cell) & (NGT_df_homvaf['cell_barcode'] == cell_barcode), 'clone_id'] = max_cn
                        cn_assignment_df.loc[(cn_assignment_df.index == cell) & (cn_assignment_df['cell_barcode'] == cell_barcode), 'clone_id'] = max_cn

    if added_clone == True:
        median_clusters, median_clusters_all = generate_median_cluster_profiles(NGT_df)
        NGT_df, cn_assignment_df, count = merge_large_clusters(NGT_df, NGT_df_homvaf, cn_assignment_df, median_clusters,  count)
        NGT_df, cn_assignment_df, count = dissolve_small_clusters(NGT_df, NGT_df_homvaf, cn_assignment_df, median_clusters,  count)
    
    print(sample_name, 'number cells changed assignment', count)
    return cn_assignment_df, NGT_df, df_tapestri_clones, added_clone_mut_dict

def main(args):
    ####################################             
    # ##### 0. prepare inputs ##########
    ####################################
    amplicon_parameters = pd.read_csv(args.w) # for homdel. get gene connectivity.
    dataset = args.d
    snv_dir = Path(args.l) # 
    print('Current dataset:', dataset)
    added_clone_mut_dict = defaultdict(str)
    

    cn_assignment_df = pd.read_csv(args.i, index_col=0)
    cn_assignment_df.index.name = 'sample_name'
    cn_assignment_df = cn_assignment_df.sort_values(by=['sample_name', 'cell_barcode'])   
    
    df_tapestri_clones = pd.read_csv(args.n, index_col=0)
    
    raw_reads = []
    raw_reads_dir = Path(args.r)
    for f in raw_reads_dir.glob('*read_counts.tsv'):
        raw_read = pd.read_csv(f, sep='\t')
        raw_read.index = [f.stem] * len(raw_read)
        raw_reads.append(raw_read)

    raw_reads_df = pd.concat(raw_reads, axis=0)
    raw_reads_df.index.name = 'sample_name'

    # read in annotated SNV file
    if args.s is not None and os.stat(args.s).st_size != 0:
        snvs = pd.read_csv(args.s, sep='\t', comment='#')
        snvs.replace(np.nan, '', inplace=True)
        snvs = snvs[snvs['annotation'].str.startswith('germline')]['condensed_format'].to_list()
        germline_mutations = list(set(snvs))

    else:
        raise ValueError('No annotated SNV file provided. Please provide a file with filtered, annotated SNVs.')
    
    # ----- 1. dissolve clones that take up >5% cells in the normal sample -----
    # @HZ: this might be too strong a filter to put here, but ok for now
    normal_idx = None
    if len(df_tapestri_clones.index[df_tapestri_clones.eq(2.0).all(axis=1)].to_list()) > 0:
        normal_idx =  df_tapestri_clones.index[df_tapestri_clones.eq(2.0).all(axis=1)].to_list()[0]
    
    if normal_idx is not None:
        cn_normal_df = cn_assignment_df[cn_assignment_df.index == 'RA18_18-11_1']

        for c in cn_assignment_df[cn_assignment_df.index == 'RA18_18-11_1']['clone_id'].unique():

            if len(cn_normal_df[cn_normal_df['clone_id'] == c])/cn_normal_df.shape[0] > 0.05:
                cn_assignment_df['clone_id'] = cn_assignment_df['clone_id'].replace(c, normal_idx)
    # ------------------------------------------

    samples = cn_assignment_df.index.to_series()
    samples = samples.unique()

    patient_fillout_h5 = snv_dir.glob(f'{dataset}.patient*.h5')[0]
    sample_obj = mio.load(patient_fillout_h5)
    sample_obj_homvaf = deepcopy(sample_obj)
    sample_obj.dna.genotype_variants(min_dp = 8, min_alt_read = 3, assign_low_conf_genotype=False)
    sample_obj_homvaf.dna.genotype_variants(min_dp = 8, het_vaf=5, hom_vaf = 99)# assign_low_conf_genotype=False)

    # Initialize NGT dataframe
    NGT_df = sample_obj.dna.get_attribute('NGT', constraint='row')
    NGT_df['sample_name'] = NGT_df.index.str.split('-').str[1] + '-' + NGT_df.index.str.split('-').str[2]
    NGT_df['cell_barcode'] = NGT_df.index.str.split('-').str[0] 
    NGT_df = NGT_df.reset_index(drop=True)
    NGT_df.set_index('sample_name', inplace=True)
    NGT_df = pd.merge(NGT_df, cn_assignment_df,on=['sample_name', 'cell_barcode'], how='inner')
    NGT_df = NGT_df[NGT_df['clone_id'].notna()]

    NGT_df_homvaf = sample_obj_homvaf.dna.get_attribute('NGT', constraint='row')
    NGT_df_homvaf['sample_name'] = NGT_df_homvaf.index.str.split('-').str[1] + '-' + NGT_df_homvaf.index.str.split('-').str[2]
    NGT_df_homvaf['cell_barcode'] = NGT_df_homvaf.index.str.split('-').str[0] 
    NGT_df_homvaf = NGT_df_homvaf.reset_index(drop=True)
    NGT_df_homvaf.set_index('sample_name', inplace=True)
    NGT_df_homvaf = pd.merge(NGT_df_homvaf, cn_assignment_df,on=['sample_name', 'cell_barcode'], how='inner')
    NGT_df_homvaf = NGT_df_homvaf[NGT_df_homvaf['clone_id'].notna()]
    bin_NGT_df = __binarize_NGT_matrix(NGT_df)
    ####################################             
    # ##### A. clone splitting #########
    ####################################
    
    name = ""
    cn_assignment_df, NGT_df, df_tapestri_clones, added_clone_mut_dict = split_sample_clones(
        name, 
        bin_NGT_df[germline_mutations + ['cell_barcode', 'clone_id']], 
        NGT_df_homvaf[germline_mutations + ['cell_barcode', 'clone_id']],
        df_tapestri_clones, 
        cn_assignment_df, 
        added_clone_mut_dict
        )

    ####################################             
    # ##### B. clone merging ###########
    ####################################

    _, median_clusters_all = generate_median_cluster_profiles(NGT_df) #bin_NGT_df?

    combined_NGT_df, cn_assignment_df, df_tapestri_clones = combine_similar_final_clusters(
        NGT_df, 
        cn_assignment_df, 
        median_clusters_all, 
        df_tapestri_clones
        )
    _, median_clusters_all = generate_median_cluster_profiles(combined_NGT_df)

    #################################################
    # ##### C. organize clones to write outputs #####
    #################################################
    n_idx = 0
    new_idx_dict = {}
    for c_id in df_tapestri_clones.index:
        if c_id not in list(cn_assignment_df['clone_id'].unique()):
            df_tapestri_clones = df_tapestri_clones.drop(c_id)
        else:
            new_idx_dict[c_id] = n_idx
            n_idx += 1

    df_tapestri_clones = df_tapestri_clones.reset_index(drop=True)
    cn_assignment_df['clone_id'] = cn_assignment_df['clone_id'].map(new_idx_dict)

    ampl_cols = [c for c in df_tapestri_clones.columns if c.startswith('AMPL')]
    matching_rows = df_tapestri_clones.loc[(df_tapestri_clones[ampl_cols] == len(ampl_cols) * [2.0]).all(axis=1)]

    if len(list(matching_rows.index)) > 0:
        row_index1 = 0
        row_index2 = list(matching_rows.index)[0]
        cn_assignment_df['clone_id'].replace({row_index1: row_index2, row_index2: row_index1}, inplace=True)

        df_tapestri_clones.iloc[row_index1], df_tapestri_clones.iloc[row_index2] = df_tapestri_clones.iloc[row_index2].copy(), df_tapestri_clones.iloc[row_index1].copy()
    
    cn_assignment_df.to_csv(args.a)
    df_tapestri_clones.to_csv(args.p)

    return 0



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', type=str, help='dataset name')
    parser.add_argument('-i', type=str, help='input path to optimal falcon clone assignment output')
    parser.add_argument('-n', type=str, help='input path to optimal falcon clone profile output')
    parser.add_argument('-r', type=str, help="input path to the patient's raw read count directory. The files should be {sample_name}.*_read_counts.tsv}")
    parser.add_argument('-l', type=str, help='input path to hdf5 files directory')
    parser.add_argument('-w', type=str, help='input path to amplicon parameters')
    parser.add_argument('-s', type=str, default=None, help='input path to manually annotated SNVs')
    parser.add_argument('-a', type=str, help='output refined clone assignment csv file')
    parser.add_argument('-p', type=str, help='output refined copy number profiles csv file')
    args = parser.parse_args(None if sys.argv[1:] else ['-h'])
    main(args)
