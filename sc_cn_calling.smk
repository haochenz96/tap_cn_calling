#
configfile: "sc_cn_calling_config.yaml"

panel_seeds=[i for i in range(config['panel_nseeds'])]

rule all:
    input:
        #expand('sc_cn_calling/{sample}_{gene}_{nclones}_{seed}_result.csv', sample=config['samples'], gene=config['genes'], nclones=config['nclones'], seed=config['seeds']),
        expand('sc_cn_calling_panel_level/{sample}_{nclones}_{seed}_result.csv', sample=config['samples'], nclones=config['nclones'], seed=panel_seeds)

rule NB_EM:
    input:
        tsv_file=lambda wildcards: config['tsv_file_dict'][wildcards.sample],
    output:
        result_file='sc_cn_calling/{sample}_{gene}_{nclones}_{seed}_result.csv',
    benchmark: 'sc_cn_calling/{sample}_{gene}_{nclones}_{seed}_benchmark.log',
    log:
        std='sc_cn_calling/{sample}_{gene}_{nclones}_{seed}.log',
        err='sc_cn_calling/{sample}_{gene}_{nclones}_{seed}.err.log',
    params:
        prefix='sc_cn_calling/{sample}_{gene}_{nclones}_{seed}',
        amplicon_df=config['amplicon_parameters'],
        nrestarts=config['nrestarts'],
        maxcn=config['maxcn'],
    shell:
        'python mixed_NB_EM_gene_level.py --sample {wildcards.sample} --readcounts {input.tsv_file} --amplicon {params.amplicon_df} --gene {wildcards.gene} '
        ' --nclones {wildcards.nclones} --nrestarts {params.nrestarts} --seed {wildcards.seed} --prefix {params.prefix} --maxcn {params.maxcn} ' 
        ' 1> {log.std} 2> {log.err}'

rule NB_panel:
    input:
        tsv_file=lambda wildcards: config['tsv_file_dict'][wildcards.sample],
    output:
        result_file='sc_cn_calling_panel_level/{sample}_{nclones}_{seed}_result.csv',
    benchmark: 'sc_cn_calling_panel_level/{sample}_{nclones}_{seed}_benchmark.log',
    log:
        std='sc_cn_calling_panel_level/{sample}_{nclones}_{seed}.log',
        err='sc_cn_calling_panel_level/{sample}_{nclones}_{seed}.err.log',
    params:
        prefix='sc_cn_calling_panel_level/{sample}_{nclones}_{seed}',
        amplicon_df=config['panel_amplicon_parameters'],
        maxcn=config['panel_maxcn'],
    shell:
        'python mixed_NB_EM_panel_level_with_normal.py --sample {wildcards.sample} --readcounts {input.tsv_file} --amplicon {params.amplicon_df} '
        ' --nclones {wildcards.nclones} --seed {wildcards.seed} --prefix {params.prefix} --maxcn {params.maxcn} '
        ' 1> {log.std} 2> {log.err}'

