from pathlib import Path

# to be overriden by command line
configfile: "sc_cn_calling_config.yaml"

# ----- specify directories -----
top_dir = Path(config['top_dir'])
cohort_name = config['cohort_name']
working_dir = top_dir / cohort_name
print(f'[INFO] --- working directory ---- {working_dir}')
workdir: working_dir

# ----- fetch panel seeds ----- 
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
        expand('sc_cn_calling_panel_level/{sample}/intermediate_results/{sample}_nclones={nclones}_seed={seed}_result.csv', sample=config['samples'], nclones=config['nclones'], seed=panel_seeds),
        # gather NB_EM results
        expand('sc_cn_calling_panel_level/{sample}/solution/{sample}_nclones={nclones}_solution-cell_assignments.csv', sample=config['samples'], nclones=config['nclones']),
        expand('sc_cn_calling_panel_level/{sample}/solution/{sample}_nclones={nclones}_solution-sc_amplicon_ploidy.csv', sample = config['samples'], nclones = config['nclones']),
        # add ploidy layer to H5
        expand('sc_cn_calling_panel_level/{sample}/ploidy_added_H5/{sample}_m2_f.ploidy_added.h5', sample=config['samples']),


rule cn_calling_gene_level:
    input:
        tsv_file=lambda wildcards: config['tsv_file_dict'][wildcards.sample],
    output:
        result_file = 'sc_cn_calling/{sample}_{gene}_{nclones}_{seed}_result.csv',
    benchmark: 'sc_cn_calling/{sample}_{gene}_{nclones}_{seed}_benchmark.log',
    log:
        std = 'sc_cn_calling/{sample}_{gene}_{nclones}_{seed}.log',
        err = 'sc_cn_calling/{sample}_{gene}_{nclones}_{seed}.err.log',
    params:
        prefix = 'sc_cn_calling/{sample}_{gene}_{nclones}_{seed}',
        amplicon_df = config['amplicon_parameters'],
        nrestarts = config['nrestarts'],
        maxcn = config['maxcn'],
    shell:
        '''
        python mixed_NB_EM_gene_level.py \
            --sample {wildcards.sample} \
            --readcounts {input.tsv_file} \
            --amplicon {params.amplicon_df} \
            --gene {wildcards.gene} \
            --nclones {wildcards.nclones} \
            --nrestarts {params.nrestarts} \
            --seed {wildcards.seed} \
            --prefix {params.prefix} \
            --maxcn {params.maxcn} \
            1> {log.std} 2> {log.err}'''

rule cn_calling_panel_level:
    input:
        tsv_file = lambda wildcards: config['tsv_file_dict'][wildcards.sample],
    output:
        result_file = 'sc_cn_calling_panel_level/{sample}/intermediate_results/{sample}_nclones={nclones}_seed={seed}_result.csv',
    benchmark: 'sc_cn_calling_panel_level/{sample}/benchmark/{sample}_nclones={nclones}_seed={seed}_benchmark.log',
    params:
        python_script = config['scripts']['cn_calling_panel_level'],
        prefix = 'sc_cn_calling_panel_level/{sample}/intermediate_results/{sample}_nclones={nclones}_seed={seed}',
        amplicon_df = config['panel_amplicon_parameters'],
        maxcn = config['panel_maxcn'],
    log:
        std = 'sc_cn_calling_panel_level/{sample}/intermediate_results/std/{sample}_nclones={nclones}_seed={seed}.log',
        err = 'sc_cn_calling_panel_level/{sample}/intermediate_results/std/{sample}_nclones={nclones}_seed={seed}.err.log',
    conda: 
        "envs/sc_cn_calling.conda.yaml",
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2000,
        time_min = lambda wildcards, attempt: attempt * 119,
    shell:
        '''
        python {params.python_script} \
            --sample {wildcards.sample} \
            --readcounts {input.tsv_file} \
            --amplicon {params.amplicon_df} \
            --nclones {wildcards.nclones} \
            --seed {wildcards.seed} \
            --prefix {params.prefix} \
            --maxcn {params.maxcn} \
            1> {log.std} 2> {log.err}
        '''

rule gather_NB_EM_results:
# scatter by sample and nclones
    input:
        EM_result_files = expand('sc_cn_calling_panel_level/{{sample}}/intermediate_results/{{sample}}_nclones={{nclones}}_seed={seed}_result.csv', seed=panel_seeds),
    output:
        optimal_solution = 'sc_cn_calling_panel_level/{sample}/solution/{sample}_nclones={nclones}_solution-cell_assignments.csv',
        sc_amplicon_ploidy = 'sc_cn_calling_panel_level/{sample}/solution/{sample}_nclones={nclones}_solution-sc_amplicon_ploidy.csv',
    params:
        input_rc_tsv_file = lambda wildcards: config['tsv_file_dict'][wildcards.sample],
        python_script = config['scripts']['gather_NB_EM_results'],
        intermediate_results_dir = 'sc_cn_calling_panel_level/{sample}/intermediate_results',
        prefix = 'sc_cn_calling_panel_level/{sample}/solution/{sample}_nclones={nclones}',
        amplicon_parameters_f = config['panel_amplicon_parameters'],
    log:
        std = 'sc_cn_calling_panel_level/{sample}/std/gather_NB_EM_results-{sample}_nclones={nclones}.log',
        err = 'sc_cn_calling_panel_level/{sample}/std/gather_NB_EM_results-{sample}_nclones={nclones}.err.log',
    conda: 
        "envs/sc_cn_calling.conda.yaml",
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2000,
        time_min = lambda wildcards, attempt: attempt * 29,
    shell:
        '''
        python {params.python_script} \
            --sample_name {wildcards.sample} \
            --nclones {wildcards.nclones} \
            --readcounts {params.input_rc_tsv_file} \
            --amplicon_parameters {params.amplicon_parameters_f} \
            --inputs_dir {params.intermediate_results_dir} \
            --prefix {params.prefix} \
            1> {log.std} 2> {log.err}
        '''

rule add_ploidy_layers_to_h5:
    # scatter by sample
    input:
        # gather for each nclones
        sc_amplicon_ploidy_all_nclones = expand('sc_cn_calling_panel_level/{{sample}}/solution/{{sample}}_nclones={nclones}_solution-sc_amplicon_ploidy.csv', nclones=config['nclones']),
        input_H5 = lambda wildcards: config['input_H5'][wildcards.sample],
    output:
        H5_ploidy_added = 'sc_cn_calling_panel_level/{sample}/ploidy_added_H5/{sample}_m2_f.ploidy_added.h5',
    params:
        python_script = config['scripts']['add_ploidy_layers_to_h5'],
    conda: 
        "envs/mosaic-custom.yaml",
    log:
        std = 'sc_cn_calling_panel_level/{sample}/std/add_ploidy_layers_to_h5-{sample}.log',
        err = 'sc_cn_calling_panel_level/{sample}/std/add_ploidy_layers_to_h5-{sample}.err.log',
    threads: 4
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 2000,
        time_min = lambda wildcards, attempt: attempt * 29,
    shell:
        '''
        python {params.python_script} \
            --sample_name {wildcards.sample} \
            --sc_amplicon_ploidy_csvs {input.sc_amplicon_ploidy_all_nclones} \
            --input_h5 {input.input_H5} \
            --output_h5 {output.H5_ploidy_added} \
            1> {log.std} 2> {log.err}
        '''

rule analyze_EM_results:
    # scatter by sample
    input:
        # gather all nclones
        optimal_solution = expand('sc_cn_calling_panel_level/{{sample}}/solution/{{sample}}_nclones={nclones}_solution-cell_assignments.csv', nclones = config['nclones']),
    output:
        # plot1: number of iterations vs number of clones
        num_itr_vs_num_clones = 'sc_cn_calling_panel_level/{sample}/analysis/{sample}_NB_EM-num_itr_vs_num_clones.png',
        # plot2: final BIC vs number of clones
        final_BIC_vs_num_clones = 'sc_cn_calling_panel_level/{sample}/analysis/{sample}_NB_EM-final_BIC_vs_num_clones.png',
    params:
        solution_dir = 'sc_cn_calling_panel_level/{sample}/solution',
        output_dir = 'sc_cn_calling_panel_level/{sample}/analysis',
    conda:
        "envs/mosaic-custom.yaml",
    log:
        std = 'sc_cn_calling_panel_level/{sample}/std/EM_results_analysis-{sample}.log',
        err = 'sc_cn_calling_panel_level/{sample}/std/EM_results_analysis-{sample}.err.log',
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
