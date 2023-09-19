import os
import sys
import pandas as pd
import numpy as np
import math
import shutil
import argparse
import itertools
import yaml

import h5py
import mosaic.io as mio
import plotly.express as px
from collections import Counter
from tea.cravat import get_technical_artifact_mask
from tea.format import isNaN
from tea.plots import plot_snv_clone
from tea.cravat import NONFUNC_SO

def get_filtered_mutations(cravat_df, sample_obj, sample_name, analysis_config, manually_annotated_snvs):
    num_cells = sample_obj.dna.shape[0]
    
    snv_selection_params = analysis_config['snv_selection_params']
        
    # 1. mutational prevalence
    mut_prev_threshold = snv_selection_params['mut_prev_threshold']
    if not type(mut_prev_threshold) is list:
        mut_prev_threshold = [mut_prev_threshold]

    # 2. technical artifact filters
    bq_prev_threshold = snv_selection_params['bq_prev_threshold']
    if bq_prev_threshold is not None and type(bq_prev_threshold) is not float: # single-value
        raise ValueError(f"bq_prev_threshold must be float, not {type(bq_prev_threshold)}")

    normals_occurences = snv_selection_params['normals_occurences'] if 'normals_occurences' in snv_selection_params else 3 # <--- default to 3

    if 'ado_threshold' in snv_selection_params and snv_selection_params['ado_threshold'] is not None:
        print(f"[INFO] filtering out SNVs with ADO > {snv_selection_params['ado_threshold']} in ANY sample...")
        ado_threshold = snv_selection_params['ado_threshold']
    else:
        ado_threshold = None

    # 3. functional SNVs
    topic = snv_selection_params['topic']
    if not type(topic) is str: # single-value
        raise ValueError(f"topic must be str, not {type(topic)}")
    try: 
        func_only = snv_selection_params['func_only']
        func_only = bool(func_only)
    except KeyError:
        func_only = False
    except TypeError:
        func_only = False

    # 4. germline snvs
    germline_attrs = {}
    for germline_attr_i in ['remove_hom_germline_snps_af', 'rescue_1000genome_af']:
        if not germline_attr_i in snv_selection_params:
            germline_attrs[germline_attr_i] = False
        else:
            germline_attrs[germline_attr_i] = snv_selection_params[germline_attr_i]

    
    if 'filter_TtoC_artifact' in snv_selection_params:
        try:
            filter_TtoC_artifact = snv_selection_params['filter_TtoC_artifact']['filter'] 
            filter_TtoC_artifact_lower_thres = snv_selection_params['filter_TtoC_artifact']['lower_thres']
            filter_TtoC_artifact_upper_thres = snv_selection_params['filter_TtoC_artifact']['upper_thres']
        except KeyError:
            raise ValueError(f"[ERROR] filter_TtoC_artifact must have keys ['filter', 'lower_thres', 'upper_thres']")
    else:
        if filter_TtoC_artifact_lower_thres >= filter_TtoC_artifact_upper_thres:
            raise ValueError(f"[ERROR] filter_TtoC_artifact_lower_thres must be strictly smaller than filter_TtoC_artifact_upper_thres")
        else:
            if mut_prev_i < 0.01:
                print(f"[INFO] mut_prev_i is lower than default upper_thres (0.01) for T>C filter. The filter will be applied.")
                filter_TtoC_artifact = True
                filter_TtoC_artifact_lower_thres = mut_prev_threshold[0]
                filter_TtoC_artifact_upper_thres = 0.01 
            else:
                print(f"[WARNING] mut_prev_i is higher than default upper_thres (0.01) for T>C filter. The filter will not be applied.")
                filter_TtoC_artifact = False

    voi_union = set()
    voi_count_union = {}
    ann_map_union = {}
    bulk_germline_vars = set()
    bulk_somatic_vars = set()
    TtoC_artifact_blacklist = set()
    ado_blacklist = set()
    
    mask = get_technical_artifact_mask(
            cravat_df,
            num_cells = num_cells, 
            bq_prev_threshold = bq_prev_threshold,
            normals_pon_occurence=normals_occurences,
            rescue_1000genome_af = germline_attrs['rescue_1000genome_af'],
            filter_broad_wes_pon = False)
    
    if filter_TtoC_artifact:
        tc_mask = (cravat_df.index.str.endswith('T/C')) & (cravat_df[('Tapestri_result', 'sc_mut_prev')] >= filter_TtoC_artifact_lower_thres * num_cells) & (cravat_df[('Tapestri_result', 'sc_mut_prev')] <= filter_TtoC_artifact_upper_thres * num_cells)
        mask = mask & ~tc_mask
        TtoC_artifact_blacklist = TtoC_artifact_blacklist.union(
        set(cravat_df.index[tc_mask].tolist())
        )

    # filters on functional SNVs
    if func_only:
        mask = mask & ~cravat_df[('Variant Annotation', 'Sequence Ontology')].isin(NONFUNC_SO)

    voi = cravat_df.index[mask].tolist()

    # filters on mut_prev_threshold
    prev_filtered_vars = sample_obj.dna.ids()[
        sample_obj.dna.get_attribute("mut_filtered", constraint="row").sum(axis=0) >= (mut_prev_threshold[0] * num_cells)
    ]
    # take intersection
    voi = [ v for v in voi if v in prev_filtered_vars ]
    

    if ado_threshold is not None:
        ado_high_vars = sample_obj.dna.ids()[
            (sample_obj.dna.get_attribute('NGT',constraint='row') == 3).sum(axis=0) > (ado_threshold*num_cells)
        ]
        voi = [ v for v in voi if v not in ado_high_vars ]
        ado_blacklist = ado_blacklist.union(set(ado_high_vars))

    voi_mut_prev = Counter(cravat_df.loc[voi, ('Tapestri_result', 'sc_mut_prev')].to_dict())

    ann = cravat_df.loc[voi, :].index.map(
        lambda x: 
        cravat_df.loc[x, ('Variant Annotation', 'Gene')] + ' ' + cravat_df.loc[x, ('Variant Annotation', 'Protein Change')] if not isNaN(cravat_df.loc[x, ('Variant Annotation','Protein Change')])
        else cravat_df.loc[x, ('Variant Annotation','Gene')] + ' ' + cravat_df.loc[x, ('Variant Annotation','Sequence Ontology')]
    )
    ann_map = dict(zip(voi, ann))


    print(len(voi))

    try:
        bulk_normal_vars = cravat_df.index[(cravat_df[('bulk_comparison', 'bulk-matched_bulk_normal-AF')] > 0)]
    except KeyError:
        print(f'bulk normal annotation not found in CRAVAT DF')
    else:
        if len(bulk_normal_vars) == 0:
            print(f'[WARNING] No bulk normal SNVs detected')
        bulk_germline_vars.update(bulk_normal_vars)

    # select SNVs detected in matched bulk cohort
    try:
        bulk_cohort_vars = cravat_df.index[(cravat_df[('bulk_comparison', 'bulk-matched_bulk_cohort-AF')] > 0)]
    except KeyError:
        print(f'bulk tumor annotation not found in CRAVAT DF')
    else:
        if len(bulk_cohort_vars) == 0:
            print(f'[WARNING] No bulk cohort SNVs detected')
        bulk_somatic_vars.update(bulk_cohort_vars)
    

    voi_union = voi_union.union(set(voi))
    voi_count_union.update(voi_mut_prev)
    ann_map_union.update(ann_map)

    # remove SNVs that are blacklisted
    print(f"[DEBUG] {len(voi_union)} SNVs before blacklist filtering")
    voi_union = voi_union.difference(TtoC_artifact_blacklist)
    print(f"[DEBUG] {len(voi_union)} SNVs after TtoC blacklist filtering")
    voi_union = voi_union.difference(ado_blacklist)
    print(f"[DEBUG] {len(voi_union)} SNVs after ADO blacklist filtering")
    
    total_mutations = []
    germline_mutations = []
    for var_i in ann_map:
        total_mutations.append(var_i)
        # germline
        if var_i in bulk_germline_vars:
            germline_mutations.append(var_i)
        else:
            pass

    germline_var_col = '#00cc66' # dark green
    somatic_var_col = '#ff0000' # red

    for var_i in ann_map_union:
        # germline
        overall_af = sample_obj.dna.get_attribute('alt_read_count', constraint='row').sum(axis=0) / sample_obj.dna.get_attribute('DP', constraint='row').sum(axis=0)
        if overall_af[var_i] >= snv_selection_params['remove_hom_germline_snps_af']:
            germline_mutations.append(var_i)
        if var_i in bulk_germline_vars:
            ann_map_union[var_i] = f'<span style="color:{germline_var_col};">' + ann_map_union[var_i] + '</span>'
        elif var_i in bulk_somatic_vars:
            ann_map_union[var_i] = f'<span style="color:{somatic_var_col};">' + ann_map_union[var_i] + '</span>'
        else:
            pass

    voi_sorted = sorted(voi_union, key=voi_count_union.get, reverse=True)
    
    if manually_annotated_snvs is not None and os.stat(manually_annotated_snvs).st_size != 0:
        snvs = pd.read_csv(manually_annotated_snvs, sep='\t', comment='#')
        snvs.replace(np.nan, '', inplace=True)
        germline_mutations = snvs[(snvs['annotation'] == 'germline_HET')]['condensed_format'].to_list()
        total_mutations = snvs[(snvs['annotation'] != 'likely_artifact') & (snvs['annotation'] != 'germline_HOM')]['condensed_format'].to_list()
        somatic_mutations = snvs[snvs['annotation'] == 'bulk_somatic']['condensed_format'].to_list()
 

    germline_mutations = list(set(germline_mutations))
    somatic_mutations = list(set(somatic_mutations))
    total_mutations = list(set(total_mutations))

    print(len(germline_mutations), len(somatic_mutations), len(total_mutations))

    print(total_mutations)
    return total_mutations, germline_mutations, somatic_mutations

    
def generate_condor_input(sample_name, cn_assignment_df, args, analysis_config, bin_thres=0.5):
    merged_cn_assignment_df = cn_assignment_df.copy()
    print(merged_cn_assignment_df['clone_id'].value_counts())

    merged_cn_assignment_df['cell_barcode'] =  merged_cn_assignment_df['idx'].astype(str) + ':' + merged_cn_assignment_df['cell_barcode'].astype(str)
    merged_cn_assignment_df.drop(['idx'], axis=1, inplace=True)

    df_total_samples = []
    df_alt_samples = []
    common_mutations = []
    germline_mutations = []
    somatic_mutations = []
    cravat_f = args.l + sample_name + '/' + sample_name + '_CRAVAT_output_cleaned.txt'
    cravat_df = pd.read_csv(cravat_f, sep='\t', index_col=0, header=[0,1])
    for file in os.listdir(args.l + sample_name):
        if file.startswith(sample_name):
            if file.endswith('.h5'):
                print(file)
                hdf5_f = args.l + sample_name + '/' + file
                sample = file.split('.')[0]
                sample_obj = mio.load(hdf5_f)
                sample_obj.dna.genotype_variants(
                    min_dp = 8,
                    min_alt_read = 3,
                    assign_low_conf_genotype = True,
                    )            
                df_alt_snv = sample_obj.dna.get_attribute('alt_read_count', constraint='row')
                df_total_snv = sample_obj.dna.get_attribute('DP', constraint='row')
                
                print(df_alt_snv.shape)
                cmuts, gmuts, smuts = get_filtered_mutations(cravat_df, sample_obj, sample_name, analysis_config, args.snvs)
                common_mutations = cmuts
                germline_mutations = gmuts
                somatic_mutations = smuts
                
                df_alt_snv.reset_index(inplace=True)
                df_alt_snv = df_alt_snv.rename(columns = {'index':'cell_barcode'})

                df_alt_samples.append(df_alt_snv)
                    
                df_total_snv.reset_index(inplace=True)
                df_total_snv = df_total_snv.rename(columns = {'index':'cell_barcode'})
                
                df_total_samples.append(df_total_snv)
    
    def mut_replace(x):
        x = x.replace(":", "_").replace("/" , "_").split('_')
        x[2], x[3] = x[3], x[2]
        return "_".join(x)

    common_mutations = list(map(mut_replace, common_mutations))
    germline_mutations_list = list(map(mut_replace, germline_mutations))
    somatic_mutations_list = list(map(mut_replace, somatic_mutations))
    print(len(common_mutations))
    
    df_total = pd.concat(df_total_samples, ignore_index=True)
    df_total = df_total.rename(columns = {'DP':'cell_barcode'})
    df_total = df_total.rename(columns={c: mut_replace(c) for c in df_total.columns if c not in ['cell_barcode']})

    df_total = df_total[['cell_barcode'] + common_mutations]
    df_total = df_total.set_index('cell_barcode')
    df_total = df_total.fillna(0)
    
    df_alt = pd.concat(df_alt_samples)
    df_alt = df_alt.rename(columns = {'alt_read_count':'cell_barcode'})
    df_alt = df_alt.rename(columns={c: mut_replace(c) for c in df_alt.columns if c not in ['cell_barcode']})

    df_alt = df_alt[['cell_barcode'] + common_mutations]
    df_alt = df_alt.set_index('cell_barcode')
    df_alt = df_alt.fillna(0)
    
    print(df_total.shape)
    print(df_alt.shape)

    df_character_mat = df_total.copy()
    df_character_mat =  df_alt.divide(df_total)


    print(df_character_mat.index)
    print(merged_cn_assignment_df['cell_barcode'])
    
    def rename_barcode(s):
        return s.split(':')[1] + '-' + s.split(':')[0]
    
    merged_cn_assignment_df['cell_barcode'] = merged_cn_assignment_df['cell_barcode'].apply(rename_barcode)
    print(merged_cn_assignment_df['clone_id'].value_counts())
    print(set([c.split('-')[2] for c in df_character_mat.index.tolist()]), merged_cn_assignment_df.shape)
    
    df_character_mat = pd.merge(df_character_mat, merged_cn_assignment_df,left_on=df_character_mat.index, right_on='cell_barcode', how='left')
    df_character_mat = df_character_mat.set_index('cell_barcode')
    print(df_character_mat['clone_id'].value_counts())
    
    df_character_mat.rename(columns={'clone_id': 'cluster_id'}, inplace=True)
        
    l_ids = list(df_character_mat['cluster_id'].unique())
    print(l_ids)
    l_ids.sort()
    l_ids_dict = {}
    index = 0
    for i in l_ids:
        l_ids_dict[i] = index
        index += 1

    df_character_mat['cluster_id'] =  df_character_mat['cluster_id'].replace(l_ids_dict)

    df_character_vaf = df_character_mat.copy()
    df_character_mat[df_character_mat.columns[:-1]] = df_character_mat[df_character_mat.columns[:-1]].applymap(lambda x: -1 if pd.isna(x) else 0 if x < bin_thres else 1)

    return df_total, df_alt, df_character_mat, df_character_vaf, germline_mutations_list, somatic_mutations_list

def main(args):
    analysis_config_yaml = args.c
    with open(analysis_config_yaml, 'r') as f:
        analysis_config = yaml.safe_load(f)

    dataset = args.d
    cn_assignment_df = pd.read_csv(args.i)
    df_total, df_alt, df_character_mat, df_character_vaf, germline_mutations, somatic_mutations = generate_condor_input(dataset, cn_assignment_df, args, analysis_config)
    with open(args.g, 'w') as fp:
        for item in germline_mutations:
            fp.write("%s\n" % item)
    
    with open(args.s, 'w') as fp:
        for item in somatic_mutations:
            fp.write("%s\n" % item)


    df_character_vaf.to_csv(args.v, index_label='')
    df_character_mat.to_csv(args.m, index_label='')
    df_alt.to_csv(args.a,index_label='')
    df_total.to_csv(args.t,index_label='')

    return 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', type=str, help='dataset name')
    parser.add_argument('-l', type=str, help='input path for hdf5 and CRAVAT files')
    parser.add_argument('-i', type=str, help='input path to refined cell assignments csv file')
    parser.add_argument('-c', type=str, help='input path to snv selection config parameter file')
    parser.add_argument('-snvs', type=str, default=None, help='input path to manually annotated SNVs')
    parser.add_argument('-v', type=str, help='output path for ConDoR VAF matrix')
    parser.add_argument('-m', type=str, help='output path for ConDoR binary matrix')
    parser.add_argument('-a', type=str, help='output path for ConDoR alternate readcount matrix')
    parser.add_argument('-t', type=str, help='output path for ConDoR total readcount matrix ')
    parser.add_argument('-g', type=str, help='output path for ConDoR germline mutation list')
    parser.add_argument('-s', type=str, help='output path for somatic mutation list')

    args = parser.parse_args(None if sys.argv[1:] else ['-h'])
    main(args)
