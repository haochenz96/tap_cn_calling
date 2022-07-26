import glob
import argparse
import plotly.express as px



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample_name', type=str, help='sample name')
    parser.add_argument('--nclones', type=int, help='number of clones')
    parser.add_argument('--readcounts', type=str, help='read count file')
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