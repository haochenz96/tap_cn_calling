from pathlib import Path
import sys

# to be overriden by command line
configfile: "sc_cn_calling_config.yaml"

# ----- specify directories -----
top_dir = Path(config['top_dir'])
cohort_name = config['cohort_name']
working_dir = top_dir / cohort_name
print(f'[INFO] --- working directory ---- {working_dir}')
workdir: working_dir

# ----- fetch run params ----- 
# ******* pay special attention to the wildcard naming here *******
# the difference between `sample` and `samples` should be noted. Although confusing, this ensures the flexibility to do both cohort-level run and single-sample run
# *****************************************************************
try:
    if config['cn_calling_mode'] not in ['single-sample', 'cohort']:
        raise ValueError(f'cn_calling_mode must be either single-sample or cohort , not {config["cn_calling_mode"]}')
    if config['cn_calling_mode'] == 'single-sample':
        cn_calling_mode = 'ss'
        samples = config['samples']
    else:
        cn_calling_mode = 'cohort'
        samples = cohort_name
    output_dir = f'{cn_calling_mode}-cn_calling'
except KeyError:
    print('[ERROR] cn_calling_mode not specified in config file')
    sys.exit(1)

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

rule all:
    input:
        # call results
        #expand('sc_cn_calling/{sample}_{gene}_{nclones}_{seed}_result.csv', sample=config['samples'], gene=config['genes'], nclones=config['nclones'], seed=config['seeds']),
        expand('{output_dir}/{sample}/intermediate_results/{sample}_nclones={nclones}_seed={seed}_result.csv', sample=samples, nclones=config['nclones'], seed=panel_seeds, output_dir = output_dir),
        # # gather NB_EM results
        expand('{output_dir}/{sample}/solutions/{sample}_nclones={nclones}_solution-cell_assignments.csv', sample=samples, nclones=config['nclones'], output_dir = output_dir),
        expand('{output_dir}/{sample}/solutions/{sample}_nclones={nclones}_solution-sc_amplicon_ploidy.csv', sample = samples, nclones = config['nclones'], output_dir = output_dir),
        # add ploidy layer to H5
        expand('{output_dir}/{sample}/outputs/{sample_i}_m2_f.ploidy_added.h5', sample=samples, sample_i=config['samples'], output_dir = output_dir),

rule cn_calling_panel_level:
# scatter by sample (if in ss mode), nclones, seed
    input:
        input_rc_tsv_file = lambda wildcards: config['tsv_file_dict'][wildcards.sample] if cn_calling_mode == 'ss' else 
        [config['tsv_file_dict'][sample_i] for sample_i in config['samples']],
    output:
        result_file = '{output_dir}/{sample}/intermediate_results/{sample}_nclones={nclones}_seed={seed}_result.csv',
    benchmark: 
        '{output_dir}/{sample}/benchmark/{sample}_nclones={nclones}_seed={seed}_benchmark.log'
    params:
        cn_calling_mode = cn_calling_mode,
        input_sample_name = lambda wildcards: wildcards.sample if cn_calling_mode == 'ss' else list(config['tsv_file_dict'].keys()),
        prefix = '{output_dir}/{sample}/intermediate_results/{sample}_nclones={nclones}_seed={seed}',
        python_script = config['scripts']['cn_calling_panel_level'],
        amplicon_df = config['panel_amplicon_parameters'],
        maxcn = config['panel_maxcn'],
    log:
        std = '{output_dir}/{sample}/intermediate_results/std/{sample}_nclones={nclones}_seed={seed}.log',
        err = '{output_dir}/{sample}/intermediate_results/std/{sample}_nclones={nclones}_seed={seed}.err.log',
    conda: 
        "envs/sc_cn_calling.conda.yaml",
    threads: 4
    resources:
        # mem_mb = lambda wildcards, attempt: attempt * 2000,
        mem_mb = 2000,
        time_min = lambda wildcards, attempt: attempt * 179,
    shell:
        '''
        python {params.python_script} \
            --cn_calling_mode {params.cn_calling_mode} \
            --sample_name {params.input_sample_name} \
            --readcounts {input.input_rc_tsv_file} \
            --amplicon_parameters_f {params.amplicon_df} \
            --nclones {wildcards.nclones} \
            --seed {wildcards.seed} \
            --prefix {params.prefix} \
            --maxcn {params.maxcn} \
            1> {log.std} 2> {log.err}
        '''

rule gather_NB_EM_results:
# scatter by sample and nclones
    input:
        EM_result_files = expand('{{output_dir}}/{{sample}}/intermediate_results/{{sample}}_nclones={{nclones}}_seed={seed}_result.csv', seed=panel_seeds),
    output:
        solution_cell_assignments = '{output_dir}/{sample}/solutions/{sample}_nclones={nclones}_solution-cell_assignments.csv',
        solution_EM_info = '{output_dir}/{sample}/solutions/{sample}_nclones={nclones}_solution-EM_info.csv',
        solution_clone_info = '{output_dir}/{sample}/solutions/{sample}_nclones={nclones}_solution-clone_info.csv',
        sc_amplicon_ploidy = '{output_dir}/{sample}/solutions/{sample}_nclones={nclones}_solution-sc_amplicon_ploidy.csv',
    params:
        cn_calling_mode = cn_calling_mode,
        input_cohort_name = lambda wildcards: wildcards.sample,
        input_sample_name = lambda wildcards: wildcards.sample if cn_calling_mode == 'ss' else list(config['tsv_file_dict'].keys()),
        input_rc_tsv_file = lambda wildcards: config['tsv_file_dict'][wildcards.sample] if cn_calling_mode == 'ss' else 
        [config['tsv_file_dict'][sample_i] for sample_i in config['samples']],
        python_script = config['scripts']['gather_NB_EM_results'],
        intermediate_results_dir = '{output_dir}/{sample}/intermediate_results',
        prefix = '{output_dir}/{sample}/solutions/{sample}_nclones={nclones}',
        amplicon_parameters_f = config['panel_amplicon_parameters'],
    log:
        std = '{output_dir}/{sample}/std/gather_NB_EM_results-{sample}_nclones={nclones}.log',
        err = '{output_dir}/{sample}/std/gather_NB_EM_results-{sample}_nclones={nclones}.err.log',
    conda: 
        "envs/sc_cn_calling.conda.yaml",
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

rule add_ploidy_layers_to_h5:
# for both ss and cohort-mode, scatter by sample; merge all nclones
# for ss mode, `sample` and `sample_i` are the same thing- it is each individual sample
# for cohort mode, `sample` is cohort_name; `sample_i` is each individual sample
    input:
        # gather for each nclones
        sc_amplicon_ploidy_all_nclones = expand('{{output_dir}}/{{sample}}/solutions/{{sample}}_nclones={nclones}_solution-sc_amplicon_ploidy.csv', nclones=config['nclones']),
        # EM_cell_assignments_all_nclones = expand('{{output_dir}}/{{sample}}/solutions/{{sample}}_nclones={nclones}_solution-cell_assignments.csv', nclones=config['nclones']),
        EM_info_csv_all_nclones = expand('{{output_dir}}/{{sample}}/solutions/{{sample}}_nclones={nclones}_solution-EM_info.csv', nclones=config['nclones']),

        # each individua sample's H5
        input_H5 = lambda wildcards: config['input_H5'][wildcards.sample_i],
    output:
        H5_ploidy_added = '{output_dir}/{sample}/outputs/{sample_i}_m2_f.ploidy_added.h5',
        EM_summary = '{output_dir}/{sample}/outputs/{sample_i}.NB_EM_summary.csv',
    params:
        cn_calling_mode = cn_calling_mode,
        input_sample_name = lambda wildcards: wildcards.sample_i,
        python_script = config['scripts']['add_ploidy_layers_to_h5'],
    conda: 
        "envs/mosaic-custom.yaml",
    log:
        std = '{output_dir}/{sample}/std/add_ploidy_layers_to_h5-{sample_i}.log',
        err = '{output_dir}/{sample}/std/add_ploidy_layers_to_h5-{sample_i}.err.log',
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2000,
        time_min = lambda wildcards, attempt: attempt * 29,
    shell:
        '''
        python {params.python_script} \
            --cn_calling_mode {params.cn_calling_mode} \
            --sample_name {params.input_sample_name} \
            --sc_amplicon_ploidy_csvs {input.sc_amplicon_ploidy_all_nclones} \
            --EM_info_csvs {input.EM_info_csv_all_nclones} \
            --input_h5 {input.input_H5} \
            --add_ploidy_for_all_or_best best \
            --output_h5 {output.H5_ploidy_added} \
            --output_EM_summary {output.EM_summary} \
            1> {log.std} 2> {log.err}
        '''

rule analyze_EM_results:
    # scatter by sample
    input:
        # gather all nclones
        optimal_solution = expand('ss-cn_calling/{{sample}}/solution/{{sample}}_nclones={nclones}_solution-cell_assignments.csv', nclones = config['nclones']),
    output:
        # plot1: number of iterations vs number of clones
        num_itr_vs_num_clones = 'ss-cn_calling/{sample}/analysis/{sample}_NB_EM-num_itr_vs_num_clones.png',
        # plot2: final BIC vs number of clones
        final_BIC_vs_num_clones = 'ss-cn_calling/{sample}/analysis/{sample}_NB_EM-final_BIC_vs_num_clones.png',
    params:
        solution_dir = 'ss-cn_calling/{sample}/solution',
        output_dir = 'ss-cn_calling/{sample}/analysis',
    conda:
        "envs/mosaic-custom.yaml",
    log:
        std = 'ss-cn_calling/{sample}/std/EM_results_analysis-{sample}.log',
        err = 'ss-cn_calling/{sample}/std/EM_results_analysis-{sample}.err.log',
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2000,
        time_min = lambda wildcards, attempt: attempt * 29,
    shell:
        '''
        python {params.python_script} \
            --sample_name {wildcards.sample} \
            --solution_dir {params.solution_dir} \
            --output_dir {output.dir} \
            1> {log.std} 2> {log.err}
        '''
