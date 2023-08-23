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

################# homdel snake-specific params #################
if config['start_from_best_sol'] == 'yes':
    # @HZ 03/13/2023: next, run with homdel, starting from the best solution without homdel
    # scatter by nclones, seed
    try:
        nclones = config['best_nclones'] # <--- this is the best nclones from pass1. Can be a list.
        random_seeds = [0]
        no_homdel_nclones_clone_info = 'cn_call_no_homdel/solutions/{cohort_name}_nclones={nclones}_solution-clone_info.csv'
    except KeyError:
        print(f'[ERROR] best_nclones is not specified in the config file. Please specify it.')
        sys.exit(1)
elif config['start_from_best_sol'] == 'no':
    # start from scratch
    nclones = config['nclones']
    random_seeds=[i for i in range(config['panel_nseeds'])]
    no_homdel_nclones_clone_info = None
else:
    print(f'[ERROR] start_from_best_sol can only be `yes` or `no`')
    sys.exit(1)

# @HZ: enforce `init_maxcn`
if not 'init_maxcn' in config:
    print(f'[ERROR] init_maxcn is not specified in the config file. Please specify it.')
    sys.exit(1)
#################


# @HZ 03/13/2023: only allows for cohort mode now
if not cn_calling_mode == 'cohort':
    print(f'[ERROR] only cohort mode is supported now')
    sys.exit(1)

output_prefix = f'{cohort_name}-cohort-cn_calling-{run_suffix}'
working_dir = top_dir / cohort_name / output_prefix
print(f'[INFO] --- working directory ---- {working_dir}')
workdir: working_dir

rule pass2_outputs:
    input:
        expand(
            'cn_call_with_homdel/outputs/{cohort_name}_homdel_nclones={nclones}.unique_cn_clone_profiles.csv',
            cohort_name = cohort_name,
            nclones = nclones,
        ),
        expand(
            'cn_call_with_homdel/outputs/{cohort_name}_BIC_vs_nclones.png', 
            cohort_name = cohort_name,
        )

rule cn_calling_panel_level_with_homdel:
# @HZ 03/13/2023: next, run with homdel, starting from the best solution without homdel
# scatter by nclones, seed
    input:
        input_rc_tsv_files = [config['tsv_file_dict'][sample_i] for sample_i in sample_names],
        # no_homdel_nclones_clone_info = 'cn_call_no_homdel/solutions/{cohort_name}_nclones={nclones}_solution-clone_info.csv',
    output:
        result_homdel_profiles = 'cn_call_with_homdel/intermediate_results/{cohort_name}-homdel-nclones={nclones}_seed={seed}.homdel_profiles.csv',
    params:
        # ----- homdel params -----
        cn_call_with_homdel_script = config['scripts']['cn_calling_panel_level_with_homdel'],
        start_from_best_sol = config['start_from_best_sol'],
        no_homdel_nclones_clone_info = no_homdel_nclones_clone_info, # None if start from scratch
        output_prefix = 'cn_call_with_homdel/intermediate_results/{cohort_name}-homdel-nclones={nclones}_seed={seed}',
        seed = lambda wildcards: wildcards.seed,
        # ----- common params -----
        cn_calling_mode = config['cn_calling_mode'],
        cohort_name = cohort_name,
        sample_names = sample_names,
        amplicon_df = config['panel_amplicon_parameters'],
        nclones = lambda wildcards: wildcards.nclones,
        maxcn = config['panel_maxcn'] if 'panel_maxcn' in config else 8, # default to 8
        init_maxcn = config['init_maxcn'] if 'init_maxcn' in config else 3, # default to 3
        min_num_amps_per_gene = config['min_num_amps_per_gene'] if 'min_num_amps_per_gene' in config else 1,
    log:
        std = 'cn_call_with_homdel/std/call/{cohort_name}-homdel-nclones={nclones}_seed={seed}.call.log',
        err = 'cn_call_with_homdel/std/call/{cohort_name}-homdel-nclones={nclones}_seed={seed}.call.err.log',
    conda: 
        "envs/sc_cn_calling.yaml",
    threads: lambda wildcards, attempt: attempt * 4
    resources:
        mem_mb = 8000,
        time_min = lambda wildcards, attempt: attempt * 179,
    retries: 2
    shell:
        '''
        python {params.cn_call_with_homdel_script} \
            --cn_calling_mode {params.cn_calling_mode} \
            --cohort_name {params.cohort_name} \
            --sample_name {params.sample_names} \
            --readcounts {input.input_rc_tsv_files} \
            --amplicon_parameters_f {params.amplicon_df} \
            --nclones {params.nclones} \
            --start_from_best_sol {params.start_from_best_sol} \
            --predefined_cn_clone_info {params.no_homdel_nclones_clone_info} \
            --maxcn {params.maxcn} \
            --init_maxcn {params.init_maxcn} \
            --seed {params.seed} \
            --prefix {params.output_prefix} \
            --min_num_amps_per_gene {params.min_num_amps_per_gene} \
            1> {log.std} 2> {log.err}        
        '''
        
rule gather_NB_EM_results_with_homdel:
    input:
        input_rc_tsv_files = [config['tsv_file_dict'][sample_i] for sample_i in sample_names],
        result_homdel_profiles = expand('cn_call_with_homdel/intermediate_results/{{cohort_name}}-homdel-nclones={{nclones}}_seed={seed}.homdel_profiles.csv', seed=random_seeds),
    output:
        solution_EM_info = 'cn_call_with_homdel/solutions/{cohort_name}-homdel-nclones={nclones}_solution.EM_info.csv',
        solution_clone_profiles = 'cn_call_with_homdel/solutions/{cohort_name}-homdel-nclones={nclones}_solution.amp_clone_profiles.csv',
        solution_sc_assignments = 'cn_call_with_homdel/solutions/{cohort_name}-homdel-nclones={nclones}_solution.cell_assignments.csv',  
    params:
        gather_with_homdel_script = config['scripts']['gather_NB_EM_results_with_homdel'],
        intermediate_results_dir = 'cn_call_with_homdel/intermediate_results',
        output_prefix = 'cn_call_with_homdel/solutions/{cohort_name}-homdel-nclones={nclones}',
        # common params
        cn_calling_mode = config['cn_calling_mode'],
        cohort_name = cohort_name,
        sample_names = sample_names,
        amplicon_df = config['panel_amplicon_parameters'],
        nclones = lambda wildcards: wildcards.nclones,
    log:
        std = 'cn_call_with_homdel/std/gather/{cohort_name}-homdel-nclones={nclones}.gather.log',
        err = 'cn_call_with_homdel/std/gather/{cohort_name}-homdel-nclones={nclones}.gather.err.log',
    conda: 
        "envs/sc_cn_calling.yaml",
    threads: lambda wildcards, attempt: attempt * 4
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 8000,
        time_min = 179,
    retries: 2
    shell:
        '''
        python {params.gather_with_homdel_script} \
        --cn_calling_mode cohort \
        --cohort_name {params.cohort_name} \
        --sample_name {params.sample_names} \
        --readcounts {input.input_rc_tsv_files} \
        --amplicon_parameters_f {params.amplicon_df} \
        --nclones {params.nclones} \
        --inputs_dir {params.intermediate_results_dir} \
        --prefix {params.output_prefix} \
        1> {log.std} 2> {log.err}
        '''

rule plot_cn_clone_profiles_compositions:
    input:
        solution_clone_profiles = 'cn_call_with_homdel/solutions/{cohort_name}-homdel-nclones={nclones}_solution.amp_clone_profiles.csv',   
        solution_sc_assignments = 'cn_call_with_homdel/solutions/{cohort_name}-homdel-nclones={nclones}_solution.cell_assignments.csv',     
    output:
        # solution_clone_profiles_plot = 'cn_call_with_homdel/outputs/{cohort_name}_homdel_nclones={nclones}.cn_clone_profiles.png',
        # result_clone_compo_plot = 'cn_call_with_homdel/outputs/{cohort_name}_homdel_nclones={nclones}.sample_CN-cluster_composition.png',
        unique_clone_profiles_csv = 'cn_call_with_homdel/outputs/{cohort_name}_homdel_nclones={nclones}.unique_cn_clone_profiles.csv',
    params:
        plot_cn_clone_profiles_script = config['scripts']['plot_cn_clone_profiles'],
        amp_gene_map_f = config['amplicon_gene_map_f'],
        output_dir = f'cn_call_with_homdel/outputs',
        output_f_prefix = lambda wildcards: f'_homdel_nclones={wildcards.nclones}',
        # common params
        cohort_name = cohort_name,
        sample_names = sample_names,
    log:
        std = 'cn_call_with_homdel/std/plot/{cohort_name}-homdel-nclones={nclones}.plot.log',
        err = 'cn_call_with_homdel/std/plot/{cohort_name}-homdel-nclones={nclones}.plot.err.log',
    conda: 
        "envs/mosaic-custom.yaml",
    shell:
        '''
        python {params.plot_cn_clone_profiles_script} \
            --cohort_name {params.cohort_name} \
            --amp_gene_map_f {params.amp_gene_map_f} \
            --cn_clone_profiles_csv {input.solution_clone_profiles} \
            --sample_sc_clone_assignment_csv {input.solution_sc_assignments} \
            --output_dir {params.output_dir} \
            --output_f_prefix {params.output_f_prefix} \
            1> {log.std} 2> {log.err}
        '''

checkpoint get_best_nclones:
    input:
        # gather all nclones
        solution_EM_infos = expand('cn_call_with_homdel/solutions/{{cohort_name}}-homdel-nclones={nclones}_solution.EM_info.csv', nclones = config['nclones']),
        solution_cell_assignments = expand('cn_call_with_homdel/solutions/{{cohort_name}}-homdel-nclones={nclones}_solution.cell_assignments.csv', nclones = config['nclones']),
    output:
        # # output directory
        # cn_call_no_homdel_output_dir = directory('cn_call_no_homdel/outputs'),
        # plot: final BIC vs number of clones
        BIC_vs_nclones_lineplot = 'cn_call_with_homdel/outputs/{cohort_name}_BIC_vs_nclones.png',
    params:
        cohort_name = cohort_name,
        output_dir = 'cn_call_with_homdel/outputs'
    log:
        std = 'cn_call_with_homdel/std/EM_results_analysis-{cohort_name}.log',
        err = 'cn_call_with_homdel/std/EM_results_analysis-{cohort_name}.err.log',
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
            best_solution[fi] = f'cn_call_with_homdel/solutions/{params.cohort_name}_nclones={best_nclones}_solution-{fi}.csv'
            # (output_dir / f'{params.cohort_name}_nclones={best_nclones}_solution-{fi}.csv').symlink_to(best_solution[fi])
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
# # ----- output-dependent final output generation -----
# # need to determine the best nclones by run result
# def get_final_output_amp_clone_profiles(wildcards):
#     checkpoint_output = checkpoints.get_best_nclones.get(**wildcards).output['BIC_vs_nclones_lineplot']
#     best_nclones_clone_info = glob(str(Path(checkpoint_output).parent / '*solution-clone_info.csv'))[0]
#     best_nclones = re.findall("nclones=\d+", best_nclones_clone_info)[0]
#     return f'cn_call_with_homdel/outputs/{wildcards.cohort_name}-homdel-nclones={best_nclones}.sample_CN-cluster_composition.png',