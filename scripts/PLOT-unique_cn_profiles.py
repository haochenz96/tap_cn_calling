# run in mosaic-custom

import argparse
import sys
from pathlib import Path
import shutil
import glob
import logging

import pandas as pd
import numpy as np
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import plotly.express as px
import yaml

from IPython import embed

def main(args):

    ######### --LOG-- ###########
    LOG_FILE = args.log_file
    if LOG_FILE is None:
        logging.basicConfig(
            stream=sys.stdout,
            level=logging.INFO, 
            filemode='a',
            force=True,
            format='%(asctime)s %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S',
            )
    else:
        logging.basicConfig(
            filename = LOG_FILE,
            level = logging.INFO,
            filemode='a',
            force=True,
            format='%(asctime)s %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S',
            )
    # ----- IO ----- #
    cohort_name = args.cohort_name
    try:
        amp_gene_map_df = pd.read_csv(args.amp_gene_map_f, index_col = 'amplicon_number', sep=None, engine='python')
    except ValueError:
        try:
            amp_gene_map_df = pd.read_csv(args.amp_gene_map_f, index_col = 'amplicon', sep=None, engine='python')
        except:
            raise ValueError(f'amp_gene_map_f {args.amp_gene_map_f} is not in the correct format!')
    cn_clone_profiles_df = pd.read_csv(args.cn_clone_profiles_csv, index_col=0, header=0)
    sample_sc_clone_assignment_df = pd.read_csv(args.sample_sc_clone_assignment_csv, index_col=False, header=0)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_f_prefix = args.output_f_prefix

    # _____check input_____
    assert amp_gene_map_df.index.astype(str)[0].startswith('AMPL'), '[WARNING] amp_gene_map_df.index should start with "AMPL"'
    if ('chr' not in amp_gene_map_df.columns) and ('chromosome' not in amp_gene_map_df.columns):
        raise ValueError(f'amp_gene_map_f {args.amp_gene_map_f} is not in the correct format!')
    else:
        CHROM_COL_NAME = 'chr' if 'chr' in amp_gene_map_df.columns else 'chromosome'

    if ('gene' not in amp_gene_map_df.columns) and ('gene_name' not in amp_gene_map_df.columns):
        raise ValueError('amp_gene_map_df must have "gene" or "gene_name" column(s)')
    else:
        GENE_COL_NAME = 'gene' if 'gene' in amp_gene_map_df.columns else 'gene_name'
    
    assert 'clone_id' in sample_sc_clone_assignment_df.columns, '[WARNING] sample_sc_clone_assignment_df must have "cluster_id" column'
    sample_sc_clone_assignment_df.columns = ['sample', 'sc_idx', 'cluster_id']

    # by design, clusters with identical CN profiles should not be in sample_sc_clone_assignment_df, but if they do, we remove them:
    # @HZ TODO: this is still faulty
    duplicated_clusters = []
    if not cn_clone_profiles_df.shape[0] == sample_sc_clone_assignment_df['cluster_id'].unique().shape[0]:
        logging.warning('sample_sc_clone_assignment_df has duplicated clusters, removing them now...')
        
        for cluster_i in cn_clone_profiles_df.index.unique():
            if cluster_i in duplicated_clusters:
                continue
            for cluster_j in cn_clone_profiles_df.index.unique():
                if cluster_i == cluster_j:
                    continue
                if (cn_clone_profiles_df.loc[cluster_i] == cn_clone_profiles_df.loc[cluster_j]).all():
                    # in the sample_sc_clone_assignment_df, replace cluster_j with cluster_i
                    logging.info(f'{cluster_i} and {cluster_j} are identical')
                    print(f'[WARNING] {cluster_i} and {cluster_j} are identical')
                    duplicated_clusters.append(cluster_j)
                    sample_sc_clone_assignment_df.loc[sample_sc_clone_assignment_df['cluster_id'] == cluster_j, 'cluster_id'] = cluster_i
        # remove duplicated clusters
    unique_cn_clone_profiles_df = cn_clone_profiles_df.drop(duplicated_clusters, axis=0) # delete cluster_j in cn_clone_profiles_df

    # cluser_name_rename_map = dict(zip(unique_cn_clone_profiles_df.index, np.arange(unique_cn_clone_profiles_df.shape[0])))
    # unique_cn_clone_profiles_df.index = np.arange(unique_cn_clone_profiles_df.shape[0])
    # # remember to rename the sample_sc_clone_assignment_df
    # sample_sc_clone_assignment_df['cluster_id'] = sample_sc_clone_assignment_df['cluster_id'].map(cluser_name_rename_map)

    # non-converged amplicons
    unique_cn_clone_profiles_df = unique_cn_clone_profiles_df.reindex(columns = amp_gene_map_df.index, fill_value=-1)
    
    unique_cn_clone_profiles_df.to_csv(
        output_dir / f'{cohort_name}{output_f_prefix}.unique_cn_clone_profiles.csv',
        index=True, header=True
    )

    # _____________________

    # ----- plotting params -----
    # Draw cluster labels
    # cluster_labels = np.arange(unique_cn_clone_profiles_df.shape[0])
    cluster_labels = unique_cn_clone_profiles_df.index.values
    logging.info(f"identified {len(cluster_labels)} unique clusters")
    # embed()
    cn_clone_palette = dict(zip(cluster_labels, np.array(px.colors.qualitative.Set3)[cluster_labels]))
    cluster_colors = [cn_clone_palette[i] for i in cluster_labels]

    # embed()

    # gene names to label
    amplicon_order = amp_gene_map_df.index.values
    gene_names = amp_gene_map_df[GENE_COL_NAME].values
    gene_names_to_plot = pd.Series(gene_names).value_counts()[pd.Series(gene_names).value_counts() >= 3].index # plot gene names covered by at least 3 amplicons

    ###########  ----- CN cluster profiles -----  ################
    # draw subplots
    fig = make_subplots(
            rows=2, cols=2,
            shared_yaxes=True, shared_xaxes=True,
            horizontal_spacing=0.01,
            vertical_spacing=0.01,
            column_widths=[1 / 25, 24 / 25],
            row_heights=[1 / 25, 24 / 25],
            )

    # get labels
    labs = go.Heatmap(z=cluster_labels,
                        y=np.arange(unique_cn_clone_profiles_df.shape[0]),
                        x=[0] * unique_cn_clone_profiles_df.shape[0],
                        customdata=cluster_labels,
                        colorscale=cluster_colors,
                        hovertemplate='label: %{customdata}<extra></extra>',
                        showlegend=False,
                        showscale=False)
    fig.add_trace(labs, row=2, col=1)

    # @HZ 12/20/2022: both ticktext and tickvals are needed to manually draw the ticklabels, otherwise it won't show up
    fig.layout.yaxis3.ticktext = cluster_labels
    fig.layout.yaxis3.tickvals = np.arange(unique_cn_clone_profiles_df.shape[0])

    # Draw gene names
    un_genes, idx_genes, inv_genes, cnts_genes = np.unique(gene_names, return_index=True, return_inverse=True, return_counts=True)
    gene_col_binary = [0]
    for i in np.arange(1,len(inv_genes)):
        # iterate through inv_genes, get connectivity of amplicons
        if inv_genes[i] == inv_genes[i-1]:
            gene_col_binary.append(gene_col_binary[i-1])
        else:
            gene_col_binary.append(abs(1-gene_col_binary[i-1]))

    ticks = (idx_genes + cnts_genes / 2).astype(int)
    gene_names_subplot = go.Heatmap(
        z = gene_col_binary,
        x = unique_cn_clone_profiles_df.columns,
        y = [0] * unique_cn_clone_profiles_df.shape[1],
        colorscale = [[0, 'rgb(0,0,0)'], [1, 'rgb(144,144,144)']],
        showlegend=False,
        showscale=False,
    )
    fig.add_trace(gene_names_subplot, row=1, col=2)

    # Draw main heatmap
    # labels = np.tile(labels[:, None], (1, unique_cn_clone_profiles_df.shape[1]))
    vals = go.Heatmap(
        z=unique_cn_clone_profiles_df,
        y=np.arange(unique_cn_clone_profiles_df.shape[0]),
        x=unique_cn_clone_profiles_df.columns,
        # customdata=labels,
        coloraxis='coloraxis',
        hovertemplate='%{z:.2f}<br>%{x}<extra>%{customdata}</extra>',
        showlegend=False,
        showscale=False
    )
    fig.add_trace(vals, row=2, col=2)

    # draw gene names
    genes_ticks = (idx_genes + cnts_genes / 2).astype(int)
        
    # embed()

    fig.layout.xaxis2.ticktext = [ i if i in gene_names_to_plot else "" for i in gene_names[genes_ticks] ]
    fig.layout.xaxis2.tickvals = genes_ticks
    fig.layout.xaxis2.tickfont = {
        'size': 8,
    }
    fig.update_layout({'xaxis2': {'ticklen': 4, 'side': 'top', 'tickangle': -60, 'showticklabels': True}})

    # draw chromosome numbers
    chromosome_ordered = amp_gene_map_df.loc[amplicon_order, CHROM_COL_NAME].values
    un, ind, cnts = np.unique(chromosome_ordered, return_index=True, return_counts=True)
    ticks = (ind + cnts / 2).astype(int)

    fig.layout.xaxis4.ticktext = chromosome_ordered[ticks]
    fig.layout.xaxis4.tickvals = ticks
    fig.update_layout({'xaxis4': {'tickangle': -45, 'showticklabels': True}})

    for i in ind:
        fig.add_vline(i - 0.5, line_color='lightcyan', line_width=1, row=2, col=2)
        
    ######################################################
    # update color schemes
    num_vals = 7
    colorscale = [
        (0, 'rgb(163, 163, 163)'), (1/num_vals, 'rgb(163, 163, 163)'),
        (1/num_vals, 'rgb(0,0,0)'), (2/num_vals, 'rgb(0,0,0)'),
        (2/num_vals, 'rgb(7, 95, 237)'), (3/num_vals, 'rgb(7, 95, 237)'),
        (3/num_vals, 'rgb(146, 170, 209)'), (4/num_vals, 'rgb(146, 170, 209)'), 
        (4/num_vals, 'rgb(237, 156, 57)'), (5/num_vals, 'rgb(237, 156, 57)'),
        (5/num_vals, 'rgb(242, 29, 22)'), (6/num_vals, 'rgb(242, 29, 22)'),
        (6/num_vals, 'rgb(163, 8, 60)'), (1, 'rgb(163, 8, 60)')
        ]

    colorbar_ticktext=[str(i) for i in list(range(num_vals-1))]
    colorbar_ticktext.insert(0, 'NA')
    colorbar_ticktext[-1] += '+'
    colorbar_tickvals = [(num_vals-1)/(num_vals*2) * (2*i + 1) - 1 for i in range(num_vals)] 
    fig.update_layout(
            coloraxis=dict(
                colorscale=colorscale,
                colorbar_tickvals = colorbar_tickvals,
                colorbar_ticktext = colorbar_ticktext,
                colorbar_title=dict(
                    font_size = 10,
                    text = 'total_copy_number',
                ),
                cmax=num_vals-2,
                cmin=-1
            ),
            font_family = 'Arial'
            )

    fig.write_image(
        output_dir / f'{cohort_name}{output_f_prefix}.cn_clone_profiles.png',
        format='png',
        width=1600,
        height=400,
        scale = 3,
    )
    logging.info(f"[SUCCESS] finished plotting unique CN clone profiles for {cohort_name}")
    ######################################################

    ############ ---- sample cluster composition ---- ################
    fig = px.histogram(
        sample_sc_clone_assignment_df,
        x = 'sample',
        color = 'cluster_id',
        color_discrete_map = cn_clone_palette,
        barnorm = 'percent',
        )
    fig.update_layout(
        legend = {
            'title': 'cluster ID',
            'traceorder': 'normal',
        },
        xaxis = {
            'title': 'sample name',
        },
        yaxis = {
            'title': 'percent of cells in each sample',
        },
        width = 600,
        height = 600,
    )

    fig.write_image(
        output_dir / f'{cohort_name}{output_f_prefix}.sample_CN-cluster_composition.png',
        format='png',
        scale = 3,
    )
    logging.info(f"[SUCCESS] finished plotting cluster compositions of samples in cohort {cohort_name}")

    ############ ---- homdel amps' raw rc distribution ---- ################
    plot_homdel = args.plot_homdel
    cn_call_yaml = args.cn_call_yaml
    
    if plot_homdel:
        if cn_call_yaml is None:
            logging.error("Please provide the rc_tsvs argument and sample names if you want to plot raw read count distributions for homdel amplicons.")
        else:
            (output_dir / 'gene_rc_distributions').mkdir(exist_ok=True, parents=True)

            # _____ read from YAML _____
            with open(cn_call_yaml, 'r') as f:
                cn_call_config = yaml.safe_load(f)
            sample_names = list(cn_call_config['tsv_file_dict'].keys())
            rc_tsvs = list(cn_call_config['tsv_file_dict'].values())
            amp_params_df = pd.read_csv(
                cn_call_config['panel_amplicon_parameters'],
                index_col= 0,
            )
            # _________________________

            homdel_amps = unique_cn_clone_profiles_df.columns[(unique_cn_clone_profiles_df == 0).any(axis=0)].tolist()
            homdel_genes = amp_params_df.loc[homdel_amps, 'gene'].unique()
            control_genes = ["KRAS", 'TP53', 'SMAD4', 'CDKN2A', 'TGFBR2'] # <------- @TODO: make more customizable
            if args.homdel_genes_oi is not None:
                control_genes += args.homdel_genes_oi
            goi = list(set(homdel_genes).union(set(control_genes)))

            # read in raw read count data
            concat_rc_df = pd.concat(
                [pd.read_csv(f, sep='\t', index_col = 0) for f in rc_tsvs], 
                keys = sample_names,
            )
            sample_sc_clone_assignment_df.set_index(['sample', 'sc_idx'], inplace=True)
            for sample_i in sample_names:
                assert concat_rc_df.loc[sample_i, :].shape[0] == sample_sc_clone_assignment_df.loc[sample_i, :].shape[0], f"Number of rows in raw read count data for sample {sample_i} does not match the number of cells in the sample."
            # embed()
            rcs_melted = pd.melt(concat_rc_df[homdel_amps], ignore_index=False, var_name = 'amplicon', value_name = 'read_counts')
            rcs_melted['cn_cluster'] = rcs_melted.index.map(lambda x: sample_sc_clone_assignment_df.loc[x, 'cluster_id'])

            rc_dist_plots = {}
            rc_box_plots = {}
            
            
            for gene_i in goi:
                amps_i = amp_params_df.index[amp_params_df['gene'] == gene_i]
                logging.info(f"for gene {gene_i}, the amplicons are: {amps_i.tolist()}")

                # normalize by amplicon factor and cell total read counts
                concat_rc_df_normed = concat_rc_df[amps_i] / amp_params_df.loc[amps_i]['amplicon_factor'] / concat_rc_df.sum(axis=1).values[:, None]
                rcs_i = pd.melt(concat_rc_df_normed, ignore_index=False, var_name = 'amplicon', value_name = 'read_counts')
                rcs_i['cn_cluster'] = rcs_i.index.map(lambda x: sample_sc_clone_assignment_df.loc[x, 'cluster_id'])

                # rc_dist_plots[gene_i] = px.histogram(
                #     data_frame = rcs_i,
                #     x = 'read_counts',
                #     facet_col = 'amplicon',
                #     facet_row = 'cn_cluster',
                # )
                # rc_dist_plots[gene_i].update_layout(
                #     title = f"{gene_i}",
                #     title_x = 0.5,
                #     font_family = 'Arial',
                #     font_size = 15,
                # ).write_image(
                #     output_dir / 'gene_rc_distributions' / f'{sample_i}_{gene_i}_rc_histogram.png',
                # )
                rc_box_plots[gene_i] = px.box(
                    data_frame = rcs_i,
                    x = 'amplicon',
                    y = 'read_counts',
                    color = 'cn_cluster',
                    color_discrete_map = cn_clone_palette,
                    points = False
                )
                rc_box_plots[gene_i].update_yaxes(
                    title_text = 'normalized read counts',
                    range = [0,10]
                )

                rc_box_plots[gene_i].update_layout(
                    title = f"{gene_i}",
                    title_x = 0.5,
                    font_family = 'Arial',
                    font_size = 15,
                ).write_image(
                    output_dir / 'gene_rc_distributions' / f'{cohort_name}{output_f_prefix}_{gene_i}.rc_boxplot.png',
                )
            logging.info(f"[SUCCESS] finished plotting raw read count distributions for homdel amplicons in cohort {cohort_name}")


    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--cohort_name', type=str, required=True)
    parser.add_argument('--amp_gene_map_f', type=str, required=True, help = 'CSV/TSV file with amplicon IDs in the first column. Required columns: `chr` and `gene`. Please make sure the amplicon IDs are ordered by genomic coordinates for best plotting results.')
    parser.add_argument('--cn_clone_profiles_csv', type=str, help='amplicon-level CN clone profile', required=True)
    parser.add_argument('--sample_sc_clone_assignment_csv', type=str, help='df assigning each sample, each single cell to each CN cluster. Therefore the 3 columns, in order, must be `sample_name`, `single-cell ID`, `cluster_id`. Only the `cluster_id` column needs to be named.', required=True)
    parser.add_argument('--output_dir', type=str, help='output directory', required=True)
    parser.add_argument('--output_f_prefix', type=str, help='output file prefix', default='')
    parser.add_argument('--plot_homdel', type=bool, help='for genes likely affected by homdel, plot distribution of all its amplicons', default=False)
    parser.add_argument('--homdel_genes_oi', type = str, nargs='+', help = 'genes of interest for plotting homdel amplicon distributions', default = None)
    parser.add_argument('--cn_call_yaml', type=str, help='Required when plotting raw rc distribution. YAML file for this cn-call run', default=None)
    parser.add_argument('--log_file', type=str, default=None)


    args = parser.parse_args(None if sys.argv[1:] else ['-h'])

    main(args)