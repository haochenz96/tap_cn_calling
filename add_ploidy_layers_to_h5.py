# run in mosaic-custom

import mosaic.io as mio
import pandas as pd
import numpy as np
from pathlib import Path
import re
import os, sys
import argparse
import plotly.express as px
import json
from tea.format import TAPESTRI_BARCODE_FORMAT
from IPython import embed

def main(args):
    cn_calling_mode = args.cn_calling_mode
    # cohort_name = args.cohort_name
    sample_name = args.sample_name
    sc_amplicon_ploidy_csvs = args.sc_amplicon_ploidy_csvs
    cell_assignments_csvs = args.cell_assignments_csvs
    EM_info_csvs = args.EM_info_csvs
    output_EM_summary = args.output_EM_summary
    input_h5 = args.input_h5
    add_ploidy_method = args.add_ploidy_for_all_or_best
    output_h5 = args.output_h5

    # infer output_dir
    output_dir = Path(output_EM_summary).parent
    output_dir.mkdir(parents=True, exist_ok=True)

    # ----- summarize EM results -----
    EM_info_dfs = [pd.read_csv(f) for f in EM_info_csvs]
    EM_summary_df = pd.concat(EM_info_dfs, ignore_index=True).sort_values(by='nclones', ignore_index=True)
    EM_summary_df.to_csv(output_EM_summary, index=False)
    # embed()
    best_idx = EM_summary_df['BIC'].idxmin()
    best_nclones = int(EM_summary_df.iloc[best_idx]['nclones'])
    print(f'[INFO] best solution: nclones = {best_nclones}')
    fig = px.line(
        EM_summary_df,
        x = 'nclones',
        y = 'BIC',
        title = f'sample {sample_name}: BIC_vs_nclones',
        markers = True,
        )
    fig.write_image(file = str(output_dir / f'{sample_name}_BIC_vs_nclones.png'),
            format="png", width=500, height=500, scale=2)

    # ----- add ploidy layers to H5 -----
    if not Path(input_h5).exists():
        print(input_h5)
        print(f'[ERROR] some of input H5 file(s)  does not exist')
        exit(1)

    # sanity check inputs:
    # print(sample_name)
    # print(input_h5)
    assert type(sample_name) == str
    assert type(input_h5) == str

    sample = mio.load(input_h5, name = sample_name, raw = False)
    sample.cnv.metadata['NB_EM_info'] = json.dumps(EM_summary_df.to_dict())

    idx = 0
    for sc_amplicon_ploidy_csv_iclone in sc_amplicon_ploidy_csvs:
        
        # infer iclone from filename
        if re.findall("nclones=\d+", sc_amplicon_ploidy_csv_iclone) is not None:
            iclone_string = re.findall("nclones=\d+", sc_amplicon_ploidy_csv_iclone)[0]
            print(f'[INFO] processing solution for {iclone_string}')
        else:
            print(f'[ERROR] could not infer iclone from filename {sc_amplicon_ploidy_csv_iclone}; skipping.')
            continue

        if add_ploidy_method == 'best' and iclone_string != f'nclones={best_nclones}':
            # skip if not the best solution
            idx += 1
            continue
        else:
            # add ploidy layer to sample
            if cn_calling_mode == 'single-sample':
                sc_amplicon_ploidy_df = pd.read_csv(
                    sc_amplicon_ploidy_csv_iclone, 
                    index_col = 0, 
                    header = 0
                    )
                cell_assignment_df = pd.read_csv(
                    cell_assignments_csvs[idx],
                    index_col = 0,
                    header = 0
                )
            else:
            # cn_calling_mode == 'cohort':
                sc_amplicon_ploidy_df = pd.read_csv(
                    sc_amplicon_ploidy_csv_iclone, 
                    index_col = [0,1], 
                    header = 0
                    )
                cell_assignment_df = pd.read_csv(
                    cell_assignments_csvs[idx],
                    index_col = [0,1],
                    header = 0
                )
                
                sc_amplicon_ploidy_df.index.names = ['sample', 'cell_barcode']
                cell_assignment_df.index.names = ['sample', 'cell_barcode']
                sc_amplicon_ploidy_df = sc_amplicon_ploidy_df.loc[sample_name]
                cell_assignment_df = cell_assignment_df.loc[sample_name] 
                sc_amplicon_ploidy_df.index = sc_amplicon_ploidy_df.index.drop_duplicates() # @HZ: I forgot what this does
                cell_assignment_df.index = cell_assignment_df.index.drop_duplicates()
            # embed()

            # sanity check to make sure the number of clones in the sc_amplicon_ploidy_df matches the number of clones in the cell_assignment_df
            nclones_sc_ploidy_df = sc_amplicon_ploidy_df.drop_duplicates().shape[0]
            nclones_cell_assignment_df = cell_assignment_df.drop_duplicates().shape[0]
            assert nclones_sc_ploidy_df == nclones_cell_assignment_df, f'[ERROR] number of clones in sc_amplicon_ploidy_df does not match number of clones in cell_assignment_df'
            
            ploidy_layer = sc_amplicon_ploidy_df.reindex(
                        index = sample.cnv.barcodes(),
                        columns = sample.cnv.ids()
                    ).fillna(0) # make sure the barcodes and amplicons are in the same order as in the input H5 
                    # reindexing fills in non-existing amplicons with NaN
                    # we then fill them with 0's to avoid plotting error later
            sample.cnv.add_layer(f'ploidy-NB_EM_{iclone_string}', ploidy_layer.values)
            print('-'*50)
            print(f'[INFO] added ploidy layer ploidy-NB_EM_{iclone_string} to {sample_name} input H5')

            cn_clone_label = cell_assignment_df.reindex(
                        index = sample.cnv.barcodes(),
                    )['clone_id']
            assert cn_clone_label.isna().any() == False, f'[ERROR] some barcodes are not assigned to any clone'
            sample.cnv.add_row_attr(f'clone_id-NB_EM_{iclone_string}', cn_clone_label.values)
            print(f'[INFO] added single-cell CN clone-id assignment clone_id-NB_EM_{iclone_string} to {sample_name} input H5')
            idx += 1

    if Path(output_h5).exists():
        print(f'[WARNING] output H5 file {output_h5} already exists, removing & overwriting')
        os.remove(output_h5)
    mio.save(sample, output_h5)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--cn_calling_mode', type=str, choices= ['single-sample', 'cohort'], help='whether to do NB_EM on each sample, or a cohort')
    # parser.add_argument('--cohort_name', type=str, help='cohort name', default='')
    parser.add_argument('--sample_name', type=str, help='sample name')
    parser.add_argument('--sc_amplicon_ploidy_csvs', nargs='+', type=str, help='dataframe containing sc-per-amplicon ploidy', required=True)
    parser.add_argument('--cell_assignments_csvs', nargs='+', type=str, help='dataframe containing cell assignments', required=True)
    parser.add_argument('--EM_info_csvs', nargs='+', type=str, help='EM info for all nclones', required=True)
    parser.add_argument('--input_h5', type=str, help='input h5 file')
    parser.add_argument('--add_ploidy_for_all_or_best', type=str, help='whether to add all nclones ploidy or only the best one', choices=['all', 'best'], default='all')
    parser.add_argument('--output_h5', type=str, help='output h5 file')
    parser.add_argument('--output_EM_summary', type=str, help='output EM summary file')

    args = parser.parse_args(None if sys.argv[1:] else ['-h'])

    main(args)