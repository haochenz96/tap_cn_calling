import pandas as pd
import numpy as np
from pathlib import Path
import re
import os, sys, glob
import argparse
import plotly.express as px

def main(args):

    # ----- io -----
    sample_name = args.sample_name
    EM_info_dir = Path(args.EM_info_dir)
    output_dir = Path(args.output_dir)

    # ----- summarize EM results -----
    EM_info_csvs = glob.glob(str(EM_info_dir / f'{sample_name}_nclones=*solution-EM_info.csv'))
    EM_info_dfs = [pd.read_csv(f) for f in EM_info_csvs]
    EM_summary_df = pd.concat(EM_info_dfs, ignore_index=True).sort_values(by='nclones', ignore_index=True)
    EM_summary_df.to_csv(output_dir / f"{sample_name}.NB_EM_summary.csv", index=False)
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

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample_name', type=str, help='sample name', required=True)
    parser.add_argument('--EM_info_dir', type=str, help='path to gathered EM info for each model (nclones)', required=True)
    parser.add_argument('--output_dir', type=str, help='output_dir ')

    args = parser.parse_args(None if sys.argv[1:] else ['-h'])

    main(args)