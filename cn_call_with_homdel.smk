# ----- fetch run params -----
top_dir = Path(config['top_dir'])
cohort_name = config['cohort_name']
sample_names = config['samples']
run_suffix = config['run_suffix']
cn_calling_mode = config['cn_calling_mode'] 

#################
best_nclones = config['best_nclones'] # <--- this is the best nclones from pass1. Can be a list.
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
            'cn_call_with_homdel/outputs/{cohort_name}_homdel_nclones={best_nclones}.sample_CN-cluster_composition.png',
            cohort_name = cohort_name,
            best_nclones = best_nclones,
        )

rule cn_calling_panel_level_with_homdel:
# @HZ 03/13/2023: next, run with homdel, starting from the best solution without homdel
# scatter by nclones, seed
    input:
        input_rc_tsv_files = [config['tsv_file_dict'][sample_i] for sample_i in sample_names],
        no_homdel_best_nclones_clone_info = 'cn_call_no_homdel/solutions/{cohort_name}_nclones={best_nclones}_solution-clone_info.csv',
    output:
        result_homdel_profiles = 'cn_call_with_homdel/intermediate_results/{cohort_name}-homdel-nclones={best_nclones}_start_from_best_sol.homdel_profiles.csv',
    params:
        cn_call_with_homdel_script = config['scripts']['cn_calling_panel_level_with_homdel'],
        output_prefix = 'cn_call_with_homdel/intermediate_results/{cohort_name}-homdel-nclones={best_nclones}_start_from_best_sol',
        # common params
        cn_calling_mode = config['cn_calling_mode'],
        cohort_name = cohort_name,
        sample_names = sample_names,
        amplicon_df = config['panel_amplicon_parameters'],
        best_nclones = lambda wildcards: wildcards.best_nclones,
    log:
        std = 'cn_call_with_homdel/std/{cohort_name}-homdel-nclones={best_nclones}_start_from_best_sol.call.log',
        err = 'cn_call_with_homdel/std/{cohort_name}-homdel-nclones={best_nclones}_start_from_best_sol.call.err.log',
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
            --nclones {params.best_nclones} \
            --predefined_cn_clone_info {input.no_homdel_best_nclones_clone_info} \
            --prefix {params.output_prefix} \
            1> {log.std} 2> {log.err}        
        '''
        
rule gather_NB_EM_results_with_homdel:
    input:
        input_rc_tsv_files = [config['tsv_file_dict'][sample_i] for sample_i in sample_names],
        result_homdel_profiles = 'cn_call_with_homdel/intermediate_results/{cohort_name}-homdel-nclones={best_nclones}_start_from_best_sol.homdel_profiles.csv',
    output:
        result_clone_profiles = 'cn_call_with_homdel/solutions/{cohort_name}-homdel-nclones={best_nclones}_start_from_best_sol_solution.amp_clone_profiles.csv',
        result_sc_assignments = 'cn_call_with_homdel/solutions/{cohort_name}-homdel-nclones={best_nclones}_start_from_best_sol_solution.cell_assignments.csv',  
    params:
        gather_with_homdel_script = config['scripts']['gather_NB_EM_results_with_homdel'],
        intermediate_results_dir = 'cn_call_with_homdel/intermediate_results',
        output_prefix = 'cn_call_with_homdel/solutions/{cohort_name}-homdel-nclones={best_nclones}_start_from_best_sol',
        # common params
        cn_calling_mode = config['cn_calling_mode'],
        cohort_name = cohort_name,
        sample_names = sample_names,
        amplicon_df = config['panel_amplicon_parameters'],
        best_nclones = lambda wildcards: wildcards.best_nclones,
    log:
        std = 'cn_call_with_homdel/std/{cohort_name}-homdel-nclones={best_nclones}_start_from_best_sol.gather.log',
        err = 'cn_call_with_homdel/std/{cohort_name}-homdel-nclones={best_nclones}_start_from_best_sol.gather.err.log',
    conda: 
        "envs/sc_cn_calling.yaml",
    shell:
        '''
        python {params.gather_with_homdel_script} \
        --cn_calling_mode cohort \
        --cohort_name {params.cohort_name} \
        --sample_name {params.sample_names} \
        --readcounts {input.input_rc_tsv_files} \
        --amplicon_parameters_f {params.amplicon_df} \
        --nclones {params.best_nclones} \
        --inputs_dir {params.intermediate_results_dir} \
        --prefix {params.output_prefix} \
        1> {log.std} 2> {log.err}
        '''

rule plot_cn_clone_profiles_compositions:
    input:
        result_clone_profiles = 'cn_call_with_homdel/solutions/{cohort_name}-homdel-nclones={best_nclones}_start_from_best_sol_solution.amp_clone_profiles.csv',   
        result_sc_assignments = 'cn_call_with_homdel/solutions/{cohort_name}-homdel-nclones={best_nclones}_start_from_best_sol_solution.cell_assignments.csv',     
    output:
        result_clone_profiles_plot = 'cn_call_with_homdel/outputs/{cohort_name}_homdel_nclones={best_nclones}.cn_clone_profiles.png',
        result_clone_compo_plot = 'cn_call_with_homdel/outputs/{cohort_name}_homdel_nclones={best_nclones}.sample_CN-cluster_composition.png',
    params:
        plot_cn_clone_profiles_script = config['scripts']['plot_cn_clone_profiles'],
        amp_gene_map_f = config['amplicon_gene_map_f'],
        output_dir = f'cn_call_with_homdel/outputs',
        output_f_prefix = lambda wildcards: f'_homdel_nclones={wildcards.best_nclones}',
        # common params
        cohort_name = cohort_name,
        sample_names = sample_names,
    conda: 
        "envs/mosaic-custom.yaml",
    shell:
        '''
        python {params.plot_cn_clone_profiles_script} \
            --cohort_name {params.cohort_name} \
            --amp_gene_map_f {params.amp_gene_map_f} \
            --cn_clone_profiles_csv {input.result_clone_profiles} \
            --sample_sc_clone_assignment_csv {input.result_sc_assignments} \
            --output_dir {params.output_dir} \
            --output_f_prefix {params.output_f_prefix} 
        '''

# # ----- output-dependent final output generation -----
# # need to determine the best nclones by run result
# def get_final_output_amp_clone_profiles(wildcards):
#     checkpoint_output = checkpoints.get_best_nclones.get(**wildcards).output['BIC_vs_nclones_lineplot']
#     best_nclones_clone_info = glob(str(Path(checkpoint_output).parent / '*solution-clone_info.csv'))[0]
#     best_nclones = re.findall("nclones=\d+", best_nclones_clone_info)[0]
#     return f'cn_call_with_homdel/outputs/{wildcards.cohort_name}-homdel-nclones={best_nclones}.sample_CN-cluster_composition.png',