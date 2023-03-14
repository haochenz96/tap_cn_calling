from pathlib import Path
from glob import glob
import sys, re
import pandas as pd
import plotly.express as px


# ----- fetch run params -----
top_dir = Path(config['top_dir'])
cohort_name = config['cohort_name']
sample_names = config['samples']
run_suffix = config['run_suffix']
cn_calling_mode = config['cn_calling_mode'] 

# @HZ 03/13/2023: only allows for cohort mode now
if not cn_calling_mode == 'cohort':
    print(f'[ERROR] only cohort mode is supported now')
    sys.exit(1)

output_prefix = f'{cohort_name}-cohort-cn_calling-{run_suffix}'
working_dir = top_dir / cohort_name / output_prefix
print(f'[INFO] --- working directory ---- {working_dir}')
workdir: working_dir

panel_seeds=[i for i in range(config['panel_nseeds'])]

# ----- check run scripts exist ----- 
for job_i in config['scripts']:
    f_path = config['scripts'][job_i]
    if not Path(f_path).is_file():
        print(f'[ERROR] script for job: {job_i} not found')
        sys.exit(1)
    else:
        print('-'*50)
        print(f'[INFO] script for job: {job_i} found')
        print(f'---> {f_path}')

rule pass1_outputs:
    input:
        # ===== pass 1 =====
        # # call results
        # expand('intermediate_results/{cohort_name}_nclones={nclones}_seed={seed}_result.csv', sample=samples, nclones=config['nclones'], seed=panel_seeds, output_dir = output_dir),
        # # # gather NB_EM results
        expand(
            'cn_call_no_homdel/solutions/{cohort_name}_nclones={nclones}_solution-EM_info.csv', 
            cohort_name = cohort_name,
            nclones = config['nclones']
        ),
        expand(
            'cn_call_no_homdel/solutions/{cohort_name}_nclones={nclones}_solution-cell_assignments.csv', 
            cohort_name = cohort_name, 
            nclones=config['nclones'], 
        ),
        expand(
            'cn_call_no_homdel/outputs/{cohort_name}_BIC_vs_nclones.png',
            cohort_name = cohort_name,
        ),
        

rule cn_calling_panel_level_no_homdel:
# @HZ 03/13/2023: first run without homdel
# scatter by nclones, seed
    input:
        input_rc_tsv_file = [config['tsv_file_dict'][sample_i] for sample_i in sample_names],
    output:
        result_file = 'cn_call_no_homdel/intermediate_results/{cohort_name}_nclones={nclones}_seed={seed}_result.csv',
    params:
        cn_calling_mode = config['cn_calling_mode'],
        input_cohort_name = cohort_name,
        input_sample_name = sample_names,
        prefix = 'cn_call_no_homdel/intermediate_results/{cohort_name}_nclones={nclones}_seed={seed}',
        python_script = config['scripts']['cn_calling_panel_level_no_homdel'],
        amplicon_df = config['panel_amplicon_parameters'],
        maxcn = config['panel_maxcn'],
    log:
        std = 'cn_call_no_homdel/intermediate_results/std/{cohort_name}_nclones={nclones}_seed={seed}.log',
        err = 'cn_call_no_homdel/intermediate_results/std/{cohort_name}_nclones={nclones}_seed={seed}.err.log',
    conda: 
        "envs/sc_cn_calling.yaml",
    threads: lambda wildcards, attempt: attempt * 4
    resources:
        # mem_mb = lambda wildcards, attempt: attempt * 2000,
        mem_mb = 2000,
        time_min = lambda wildcards, attempt: attempt * 179,
    retries: 3
    shell:
        '''
        python {params.python_script} \
            --cn_calling_mode {params.cn_calling_mode} \
            --sample_name {params.input_sample_name} \
            --cohort_name {params.input_cohort_name} \
            --readcounts {input.input_rc_tsv_file} \
            --amplicon_parameters_f {params.amplicon_df} \
            --nclones {wildcards.nclones} \
            --seed {wildcards.seed} \
            --prefix {params.prefix} \
            --maxcn {params.maxcn} \
            1> {log.std} 2> {log.err}
        '''

rule gather_NB_EM_results_no_homdel:
# scatter by nclones
    input:
        EM_result_files = expand('cn_call_no_homdel/intermediate_results/{{cohort_name}}_nclones={{nclones}}_seed={seed}_result.csv', seed=panel_seeds),
    output:
        solution_cell_assignments = 'cn_call_no_homdel/solutions/{cohort_name}_nclones={nclones}_solution-cell_assignments.csv',
        solution_EM_info = 'cn_call_no_homdel/solutions/{cohort_name}_nclones={nclones}_solution-EM_info.csv',
        solution_clone_info = 'cn_call_no_homdel/solutions/{cohort_name}_nclones={nclones}_solution-clone_info.csv',
        sc_amplicon_ploidy = 'cn_call_no_homdel/solutions/{cohort_name}_nclones={nclones}_solution-sc_amplicon_ploidy.csv',
    params:
        cn_calling_mode = cn_calling_mode,
        input_cohort_name = cohort_name,
        input_sample_name = sample_names,
        input_rc_tsv_file = [config['tsv_file_dict'][sample_i] for sample_i in sample_names],
        python_script = config['scripts']['gather_NB_EM_results_no_homdel'],
        intermediate_results_dir = 'cn_call_no_homdel/intermediate_results',
        prefix = 'cn_call_no_homdel/solutions/{cohort_name}_nclones={nclones}',
        amplicon_parameters_f = config['panel_amplicon_parameters'],
    log:
        std = 'cn_call_no_homdel/std/gather_NB_EM_results-{cohort_name}_nclones={nclones}.log',
        err = 'cn_call_no_homdel/std/gather_NB_EM_results-{cohort_name}_nclones={nclones}.err.log',
    conda: 
        "envs/sc_cn_calling.yaml",
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2000,
        time_min = lambda wildcards, attempt: attempt * 29,
    shell:
        '''
        python {params.python_script} \
            --cn_calling_mode {params.cn_calling_mode} \
            --cohort_name {params.input_cohort_name} \
            --sample_name {params.input_sample_name} \
            --readcounts {params.input_rc_tsv_file} \
            --nclones {wildcards.nclones} \
            --amplicon_parameters_f {params.amplicon_parameters_f} \
            --inputs_dir {params.intermediate_results_dir} \
            --prefix {params.prefix} \
            1> {log.std} 2> {log.err}
        '''
        
checkpoint get_best_nclones:
    input:
        # gather all nclones
        solution_EM_infos = expand('cn_call_no_homdel/solutions/{{cohort_name}}_nclones={nclones}_solution-EM_info.csv', nclones = config['nclones']),
        solution_cell_assignments = expand('cn_call_no_homdel/solutions/{{cohort_name}}_nclones={nclones}_solution-cell_assignments.csv', nclones = config['nclones']),
    output:
        # # output directory
        # cn_call_no_homdel_output_dir = directory('cn_call_no_homdel/outputs'),
        # plot: final BIC vs number of clones
        BIC_vs_nclones_lineplot = 'cn_call_no_homdel/outputs/{cohort_name}_BIC_vs_nclones.png',
    params:
        cohort_name = cohort_name,
        output_dir = 'cn_call_no_homdel/outputs'
    log:
        std = 'cn_call_no_homdel/std/EM_results_analysis-{cohort_name}.log',
        err = 'cn_call_no_homdel/std/EM_results_analysis-{cohort_name}.err.log',
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2000,
        time_min = lambda wildcards, attempt: attempt * 29,
    run:
        # ----- IO -----
        EM_info_csvs = input.solution_EM_infos
        EM_info_dfs = [pd.read_csv(f) for f in EM_info_csvs]
        output_dir = Path(params.output_dir)
        # ----- summarize EM results -----
        EM_summary_df = pd.concat(EM_info_dfs, ignore_index=True).sort_values(by='nclones', ignore_index=True)
        EM_summary_df.to_csv(output_dir / f"{params.cohort_name}_NB_EM_summary.csv", index=False)
        best_idx = EM_summary_df['BIC'].idxmin()
        best_nclones = int(EM_summary_df.iloc[best_idx]['nclones']) # <-- global variable defined here
        # print(f'[INFO] best solution: nclones = {best_nclones}')
        # ----- symlink best nclones solution
        fnames = ['cell_assignments', 'amp_clone_profiles', 'clone_info']
        best_solution = {}
        for fi in fnames:
            best_solution[fi] = f'cn_call_no_homdel/solutions/{params.cohort_name}_nclones={best_nclones}_solution-{fi}.csv'
            (output_dir / f'{params.cohort_name}_nclones={best_nclones}_solution-{fi}.csv').symlink_to(best_solution[fi])
        # ----- plot -----
        fig = px.line(
            EM_summary_df,
            x = 'nclones',
            y = 'BIC',
            title = f'sample {params.cohort_name}: BIC_vs_nclones',
            markers = True,
            )
        fig.write_image(file = str(output_dir / f'{params.cohort_name}_BIC_vs_nclones.png'),
                format="png", width=500, height=500, scale=2)

# rule add_ploidy_layers_to_h5:
# # for both ss and cohort-mode, scatter by sample; merge all nclones
# # for ss mode, `sample` and `sample_i` are the same thing- it is each individual sample
# # for cohort mode, `sample` is cohort_name; `sample_i` is each individual sample
#     input:
#         # gather for each nclones
#         sc_amplicon_ploidy_all_nclones = expand('{{cohort_name}}/solutions/{{cohort_name}}_nclones={nclones}_solution-sc_amplicon_ploidy.csv', nclones=config['nclones']),
#         EM_cell_assignments_all_nclones = expand('{{cohort_name}}/solutions/{{cohort_name}}_nclones={nclones}_solution-cell_assignments.csv', nclones=config['nclones']),
#         EM_info_csv_all_nclones = expand('{{cohort_name}}/solutions/{{cohort_name}}_nclones={nclones}_solution-EM_info.csv', nclones=config['nclones']),

#         # each individua sample's H5
#         input_H5 = lambda wildcards: config['input_H5'][wildcards.sample_i],
#     output:
#         H5_ploidy_added = 'outputs/{sample_i}_m2_f.ploidy_added.h5',
#         EM_summary = 'outputs/{sample_i}.NB_EM_summary.csv',
#     params:
#         cn_calling_mode = config['cn_calling_mode'],
#         input_sample_name = lambda wildcards: wildcards.sample_i,
#         python_script = config['scripts']['add_ploidy_layers_to_h5'],
#         add_all_or_best = config['add_all_or_best']
#     conda: 
#         "envs/mosaic-custom.yaml",
#     log:
#         std = 'std/add_ploidy_layers_to_h5-{sample_i}.log',
#         err = 'std/add_ploidy_layers_to_h5-{sample_i}.err.log',
#     threads: 4
#     resources:
#         mem_mb = lambda wildcards, attempt: attempt * 2000,
#         time_min = lambda wildcards, attempt: attempt * 29,
#     shell:
#         '''
#         python {params.python_script} \
#             --cn_calling_mode {params.cn_calling_mode} \
#             --sample_name {params.input_sample_name} \
#             --sc_amplicon_ploidy_csvs {input.sc_amplicon_ploidy_all_nclones} \
#             --cell_assignments_csvs {input.EM_cell_assignments_all_nclones} \
#             --EM_info_csvs {input.EM_info_csv_all_nclones} \
#             --input_h5 {input.input_H5} \
#             --add_ploidy_for_all_or_best {params.add_all_or_best} \
#             --output_h5 {output.H5_ploidy_added} \
#             --output_EM_summary {output.EM_summary} \
#             1> {log.std} 2> {log.err}
#         '''