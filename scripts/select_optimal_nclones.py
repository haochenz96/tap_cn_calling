import os
import sys
import pandas as pd
import numpy as np
import math
import shutil
import argparse
import itertools
import plotly.express as px
# import matplotlib.pyplot as plt
from pathlib import Path
from IPython import embed

def average_distance(df_tapestri_clones, cn_assignment_df, amps_of_interest=None):
    if amps_of_interest is None:
        amps_of_interest = df_tapestri_clones.columns
    sample_prop = cn_assignment_df.groupby(by='clone_id').size()/len(cn_assignment_df)
    found_clones = list(cn_assignment_df.groupby(by='clone_id').size().index)
    # print(found_clones)
    # print(df_tapestri_clones.index)
    distances = []
    for c_id in range(len(found_clones)):
        clone_id = found_clones[c_id]
        s1 = len(cn_assignment_df[cn_assignment_df['clone_id'] == clone_id])
        for c_id2 in range(c_id + 1, len(found_clones)):
            clone_id2 = found_clones[c_id2]
            s2 = len(cn_assignment_df[cn_assignment_df['clone_id'] == clone_id2])
 #           if s1 > 10 and s2 > 10:
            l1 = df_tapestri_clones.loc[clone_id].values
            l2 = df_tapestri_clones.loc[clone_id2].values
            dist = np.linalg.norm(l1 - l2, ord=2)
            dist = dist * sample_prop.loc[clone_id] * sample_prop.loc[clone_id2]
            distances.append(dist)
    return np.sum(distances)

def main(args):
    if args.c is not None:
        clones = [int(''.join(i)) for i in args.c]
    else:
        # need to parse out the nclone string from the file names
        clones = [int(i.split('=')[1].split('.')[0]) for i in os.listdir(args.i) if i.endswith('.csv') and "unique_cn_clone_profiles.csv" in i]

    if args.amplicon_gene_map is not None:
        amplicon_gene_map_df = pd.read_csv(args.amplicon_gene_map, index_col=0,)
        if not amplicon_gene_map_df.index[0].startswith('AMPL'):
            # embed()
            raise ValueError("Amplicon gene map file is not formatted correctly")
        # use either gene or gene_name column
        found = False
        for col in ['gene', 'gene_name']:
            if col in amplicon_gene_map_df.columns:
                found = True
                break
        if found is False:
            raise ValueError("Amplicon gene map file is not formatted correctly")
        genes_of_interest = amplicon_gene_map_df[col].value_counts()[amplicon_gene_map_df[col].value_counts() >= args.min_amps_per_gene].index
        amplicons_of_interest = amplicon_gene_map_df[amplicon_gene_map_df[col].isin(genes_of_interest)].index
        print(f"[INFO] Found {len(amplicons_of_interest)} amplicons with at least {args.min_amps_per_gene} amplicons targeting each gene")
    else:
        amplicons_of_interest = None
        

    avg_distance = []
    proper_clones = []
    for nclone in clones:
        # print(args.i, args.d)
        df_tapestri_clones = pd.read_csv(Path(args.i) / f"{args.d}_homdel_nclones={nclone}.unique_cn_clone_profiles.csv", index_col = 0)
        cn_assignment_df = pd.read_csv(Path(args.i) / f"{args.d}_homdel_nclones={nclone}.sample_sc_clone_assignment.updated.csv", index_col = 0)
        df_tapestri_clones = df_tapestri_clones.drop(df_tapestri_clones.columns[df_tapestri_clones.isin([-1]).any()], axis=1)
        threshold = 0.10 * df_tapestri_clones.shape[1]
        malignant_rows = ((df_tapestri_clones == 0.0).sum(axis=1) > threshold)
        
        proper_result = True
        if args.no_solution_filter:
            print("Not filtering solutions!")
        else:
            for mal_row in malignant_rows.index[malignant_rows]:
                if len(cn_assignment_df[cn_assignment_df['clone_id'] == mal_row]) > 10:
                    proper_result = False

        if proper_result:
            proper_clones.append(nclone)
            avg_val = average_distance(df_tapestri_clones, cn_assignment_df, amplicons_of_interest)
            avg_distance.append(avg_val)
        else:
            print(f'Removed solution with {nclone} clones due to presence of at least one nontrivial clone with high number of amplicons called as homdel')
            # print(nclone)
            # print(malignant_rows)

    distances_df = pd.DataFrame({'nclone': proper_clones, 'avg_distance': avg_distance})
    distances_df['distance_from_prev_clone'] = distances_df['avg_distance'].diff()

    Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    # plot the distance from previous clone using plotly.express
    # sort the distances_df by nclone
    distances_df = distances_df.sort_values(by='nclone')
    fig = px.line(distances_df, x='nclone', y='avg_distance')
    fig.update_layout(
        title=f"Average inter-clone distances {args.d}",
        xaxis_title="Nclone",
        yaxis_title="Average inter-clone distance",
        font=dict(
            family="Courier New, monospace",
            size=12,
            color="#7f7f7f"
        )
    )
    # save a png
    fig.write_image(Path(args.output_dir) / f"{args.d}_solution.avg_distance.png")

    # plt.figure(figsize=(10, 10))
    # sns.lineplot(x='nclone', y='avg_distance', data=distances_df)
    # plt.savefig(Path(args.output_dir) / f"{args.d}_solution.avg_distance.png")

    # # @HZ: optimal nclone is the one with the largest distance from previous clone
    # opt_nclone = distances_df['nclone'].iloc[distances_df['distance_from_prev_clone'].idxmax()]
    # embed()
    # @PS:
    # optimal nclone is the one with the largest between-clone average distance
    opt_nclone = distances_df.loc[distances_df['avg_distance'].idxmax()]['nclone'].astype(int)


    print(f"[INFO] sample {args.d} -- opt_nclone: {opt_nclone}")
    # copy over
    for fi in os.listdir(args.i):
        if f"nclones={opt_nclone}" in fi:
            shutil.copyfile(Path(args.i) / fi, Path(args.output_dir) / fi)
    # cn_profile_fp = Path(args.i) / f"{args.d}_homdel_nclones={opt_nclone}.unique_cn_clone_profiles.csv"
    # clone_assignment_fp = Path(args.i) / f"{args.d}_homdel_nclones={opt_nclone}.sample_sc_clone_assignment.updated.csv"

    # # shutil.copyfile(clone_assignment_fp, args.a) 
    # # shutil.copyfile(cn_profile_fp, args.p) 

    # shutil.copyfile(cn_profile_fp, Path(args.output_dir) / f"{args.d}_homdel_nclones={opt_nclone}.unique_cn_clone_profiles.csv")
    # shutil.copyfile(clone_assignment_fp, Path(args.output_dir) / f"{args.d}_homdel_nclones={opt_nclone}.sample_sc_clone_assignment.updated.csv")

    return 0

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', type=str, help='dataset name')
    parser.add_argument('-i', type=str, help='input path to falcon solution')
    parser.add_argument('-c', type=list, nargs='+', help='nclone numbers to consider. If None, all nclones in solutions directory will be considered', default=None)
    parser.add_argument('--amplicon_gene_map', type=str, help='input path to amplicon gene map', default=None)
    parser.add_argument('--min_amps_per_gene', type=int, help='minimum number of amplicons per gene', default=3)
    parser.add_argument('--output_dir', type=str, help='output directory')
    # parser.add_argument('-a', type=str, help='output optimal clone assignment csv file')
    # parser.add_argument('-p', type=str, help='output optimal copy number profiles csv file')
    # add option to turn off filtering of solutions with too many clones with high number of amplicons called as homdel
    parser.add_argument('--no_solution_filter', action='store_true', help='do not filter out solutions with too many clones with high number of amplicons called as homdel')
    args = parser.parse_args(None if sys.argv[1:] else ['-h'])
    main(args)
