import mosaic.io as mio
import pandas as pd
import numpy as np
from pathlib import Path
import re
import os, sys
import argparse


def main(args):
    sample_name = args.sample_name
    sc_amplicon_ploidy_csvs = args.sc_amplicon_ploidy_csvs
    input_H5 = args.input_h5
    output_H5 = args.output_h5

    if not Path(input_H5).exists():
        print(f'[ERROR] output H5 file {output_H5} does not exist')
        exit(1)
    else:
        sample = mio.load(input_H5, name = sample_name, raw = False)

    for sc_amplicon_ploidy_csv_iclone in sc_amplicon_ploidy_csvs:
        # infer iclone from filename
        if re.findall("nclones=\d+", sc_amplicon_ploidy_csv_iclone) is not None:
            iclone_string = re.findall("nclones=\d+", sc_amplicon_ploidy_csv_iclone)[0]
        else:
            print(f'[ERROR] could not infer iclone from filename {sc_amplicon_ploidy_csv_iclone}; skipping.')
            continue
        # add ploidy layer to sample
        sc_amplicon_ploidy_df = pd.read_csv(
            sc_amplicon_ploidy_csv_iclone, 
            index_col = 0, 
            header = 0
            )
        ploidy_layer = sc_amplicon_ploidy_df.reindex(
                    index = sample.cnv.barcodes(),
                    columns = sample.cnv.ids()
                ) # make sure the barcodes and amplicons are in the same order as in the input H5 
                # reindexing fills in non-existing amplicons with NaN
        sample.cnv.add_layer(f'NB_EM_{iclone_string}', ploidy_layer.values)
        print('-'*50)
        print(f'[INFO] added ploidy layer {iclone_string} to {sample_name} input H5')

    if Path(output_H5).exists():
        print(f'[WARNING] output H5 file {output_H5} already exists, removing & overwriting')
        os.remove(output_H5)
    mio.save(sample, output_H5)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample_name', type=str, help='sample name')
    parser.add_argument('--sc_amplicon_ploidy_csvs', nargs='+', type=str, help='number of clones', required=True)
    parser.add_argument('--input_h5', type=str, help='input h5 file')
    parser.add_argument('--output_h5', type=str, help='output h5 file')

    args = parser.parse_args(None if sys.argv[1:] else ['-h'])

    main(args)