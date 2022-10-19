import allel
import loompy
import numpy as np
import pandas as pd
import pysam
import seaborn as sns
import math
import argparse
import sys

def get_condensed_vaf_matrix(ds, df_pos, amplicon_list,
                             read_depth_threshold = 10, vaf_threshold = 0.33, presence_threshold = 0.2, amplicon_threshold=0):
    
    namplicon = len(amplicon_list)
    vaf_mat = None
    alt_mat = None
    read_depth_mat = None
    read_depth = None
    cum_pos_indices = None
    cum_cell_indices = None
    for idx, amplicon in enumerate(amplicon_list):
        selected_indices = df_pos[df_pos['amplicon'] == amplicon]['index'].values
        
        b = ds.layers['AD'][selected_indices, :]
        a = ds.layers['DP'][selected_indices, :]
        
        # build read depth mat
        if read_depth is not None:
            read_depth = np.vstack((read_depth, np.median(a, axis = 0)))
        else:
            read_depth = np.median(a, axis = 0)        
        
        # remove positions with median read depth less than read_depth_threshold
        selected_indices = selected_indices[np.median(a, axis = 1) >= read_depth_threshold][:,None]
        b = b[np.median(a, axis = 1) >= read_depth_threshold, :]
        a = a[np.median(a, axis = 1) >= read_depth_threshold, :]        

        if vaf_mat is not None:
            vaf_mat = np.vstack((vaf_mat, b / a))
            read_depth_mat = np.vstack((read_depth_mat, a))
            alt_mat = np.vstack((alt_mat, b))
        else:
            vaf_mat = b / a
            read_depth_mat = a
            alt_mat = b
        
        if cum_pos_indices is not None:
            cum_pos_indices = np.vstack((cum_pos_indices, selected_indices))
        else:
            cum_pos_indices = selected_indices.copy()

    ncells = vaf_mat.shape[1]
    cum_cell_indices = np.arange(ncells)
    
    # selecting the cells in which the read depth is more than read depth threshold for at least amplicon threshold percent of the amplicons
    if amplicon_threshold > 0:
        curr_selected_indices = np.sum(read_depth >= read_depth_threshold, axis = 0) >= amplicon_threshold*namplicon
        vaf_mat = vaf_mat[:, curr_selected_indices]
        read_depth_mat = read_depth_mat[:, curr_selected_indices]
        alt_mat = alt_mat[:, curr_selected_indices]
        cum_cell_indices = cum_cell_indices[curr_selected_indices]

    # remove positions with max vaf less than vaf_threshold
    curr_selected_indices = np.nanmax(vaf_mat, axis = 1) >= vaf_threshold
    cum_pos_indices = cum_pos_indices[curr_selected_indices]
    vaf_mat = vaf_mat[curr_selected_indices, :]
    read_depth_mat = read_depth_mat[curr_selected_indices, :]
    alt_mat = alt_mat[curr_selected_indices, :]
    
    # remove positions with vaf less than vaf_threshold in less than presense_threshold fraction of the cells
    ncells = vaf_mat.shape[1]
    vaf_threshold_cell_count = (np.nan_to_num(vaf_mat) >= vaf_threshold).sum(axis = 1)
    curr_selected_indices = vaf_threshold_cell_count >= presence_threshold*ncells
    cum_pos_indices = cum_pos_indices[curr_selected_indices]
    vaf_mat = vaf_mat[curr_selected_indices, :]
    read_depth_mat = read_depth_mat[curr_selected_indices, :]
    alt_mat = alt_mat[curr_selected_indices, :]
    
    return vaf_mat, read_depth_mat, alt_mat, cum_pos_indices.flatten(), cum_cell_indices

def write_output_files(dest_prefix, vaf_mat, read_depth_mat, pos_indices, cell_indices, df_selected_pos,
                       mutation_presence_threshold=0.2, homozygous_mutation_threshold=0.8):
    snv_mat = np.zeros(vaf_mat.shape)
    snv_mat[vaf_mat < mutation_presence_threshold] = 0
    snv_mat[vaf_mat >= homozygous_mutation_threshold] = 2
    snv_mat[(vaf_mat >= mutation_presence_threshold) & (vaf_mat < homozygous_mutation_threshold)] = 1
    snv_mat[np.isnan(vaf_mat)] = 3
    snv_mat = snv_mat.astype(int)
    
    np.savetxt(f"{dest_prefix}_scite_snv_mat.txt", snv_mat, delimiter=" ", fmt='%d')
    np.savetxt(f"{dest_prefix}_sifit_snv_mat.txt", np.hstack((np.arange(snv_mat.shape[0])[:,None], snv_mat)), delimiter=" ", fmt='%d')
#     df_pos.set_index('index').loc[pos_indices].reset_index().drop(['ref_len', 'alt_len', 'normal'], axis=1).to_csv(f'{dest_prefix}_pos_indices.csv', sep=',', index=False)
    
    snv_mat_str = snv_mat.astype(str)
    snv_mat_str = np.char.replace(snv_mat_str, '3', '?').T
    snv_mat_str = np.char.replace(snv_mat_str, '2', '1')
    
#     snv_mat_phiscs = np.vstack((np.hstack((np.array([['cell_idx/mut_idx']]), pos_indices[None, :].astype(str))), np.hstack((cell_indices[:,None], snv_mat_str))))
    snv_mat_phiscs = np.vstack((np.hstack((np.array([['cell_idx/mut_idx']]), df_selected_pos['gene'].values[None, :].astype('str'))), np.hstack((cell_indices[:,None], snv_mat_str))))
    np.savetxt(f"{dest_prefix}_phiscs_snv_mat.txt", snv_mat_phiscs, delimiter="\t", fmt='%s')
    
#     return snv_mat

def main(args):
    
    # connect to the loom file
    ds = loompy.connect(args.loom)
    
    # construct position dataframe
    df_pos = pd.DataFrame({'pos': list(ds.ra['POS']),
                           'chrom': list(ds.ra['CHROM']),
                           'amplicon': list(ds.ra['amplicon']),
                           'ref': list(ds.ra['REF']),
                           'alt': list(ds.ra['ALT'])})

    df_pos['ref_len'] = df_pos['ref'].apply(len)
    df_pos['alt_len'] = df_pos['alt'].apply(len)
    df_pos['normal'] = (df_pos['ref_len'] == 1) & (df_pos['alt_len'] == 1)

    df_pos['index'] = df_pos.index
    
    amplicon_list = list(df_pos.groupby(['pos', 'amplicon']).last().reset_index()['amplicon'].unique())
    
    df_gene = pd.read_csv('/n/fs/ragr-research/projects/starch/STARCH/hgTables_hg19.txt', sep='\t', index_col=0)

    df_gene = df_gene[df_gene['chrom'].apply(len) <= 5]
    
    # construct the VAF and read depth matrices
    vaf_mat, read_depth_mat, alt_mat, pos_indices, cell_indices = get_condensed_vaf_matrix(ds, df_pos, amplicon_list, 
                                                                                           read_depth_threshold = args.depth,
                                                                                           vaf_threshold = args.vaf,
                                                                                           presence_threshold = args.presence,
                                                                                           amplicon_threshold= args.amplicon)
    
    # construct and write dataframes
    df_selected_pos = df_pos.set_index('index').loc[pos_indices].reset_index().drop(['ref_len', 'alt_len', 'normal'], axis=1)
    
    gene_list = []
    for idx, row in df_selected_pos.iterrows():
        pos = row['pos']
        chrom = row['chrom']

        df_select = df_gene[(df_gene['chrom'] == ''.join(['chr', str(chrom)])) & (df_gene['cdsStart'] <= pos) & (df_gene['cdsEnd'] >= pos)]
        if len(df_select) == 0:
            gene_list.append(row['index'])
        else:
            gene = list(df_select['name2'])[0]
            if gene in gene_list:
                gene_list.append('_'.join([gene, str(len([x for x in gene_list if str(x).startswith(gene)]))]))
            else:
                gene_list.append(gene)

    df_selected_pos['gene'] = gene_list    
    
    df_selected_pos.to_csv(f'{args.prefix}_pos_indices.csv', index=False)
    pd.concat([df_selected_pos, pd.DataFrame(vaf_mat, columns=cell_indices)], axis=1).to_csv(f'{args.prefix}_vaf.csv', index=False)
    pd.concat([df_selected_pos, pd.DataFrame(alt_mat, columns=cell_indices)], axis=1).to_csv(f'{args.prefix}_alt.csv', index=False)
    pd.concat([df_selected_pos, pd.DataFrame(read_depth_mat, columns=cell_indices)], axis=1).to_csv(f'{args.prefix}_total.csv', index=False)    
    
    # write the output files
    write_output_files(args.prefix, vaf_mat, read_depth_mat, pos_indices, cell_indices, df_selected_pos,
                       mutation_presence_threshold = args.hetero, homozygous_mutation_threshold = args.homo)
    
    # disconnect
    ds.close()
    
def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--loom', type=str, help='loom file name')
    parser.add_argument('--depth', type=int, help='read depth threshold [10]', default=10)
    parser.add_argument('--vaf', type=float, help='variant allele frequency threshold [0.33]', default=0.33)
    parser.add_argument('--presence', type=float, help='prevalence threshold [0.2]', default=0.2)
    parser.add_argument('--amplicon', type=int, help='amplicon depth threshold [0]', default=0) 
    parser.add_argument('--hetero', type=float, help='vaf threshold to call a mutation heterozygous [0.2]', default=0.2)
    parser.add_argument('--homo', type=float, help='vaf threshold to call a mutation homozygous [0.8]', default=0.8)
    
    parser.add_argument('--prefix', type=str, help='prefix for output files of plotting')

    args = parser.parse_args(None if sys.argv[1:] else ['-h'])

    main(args)
