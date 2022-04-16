#
configfile: "brute_force_cn_calling_config.yaml"

rule all:
    input:
        expand('brute_force_cn_calling/{sample}_{gene}_{seed}_result.csv', sample=config['samples'], gene=config['genes'], seed=config['seeds']),

rule NB_EM:
    input:
        tsv_file=lambda wildcards: config['tsv_file_dict'][wildcards.sample],
    output:
        result_file='brute_force_cn_calling/{sample}_{gene}_{seed}_result.csv',
    benchmark: 'brute_force_cn_calling/{sample}_{gene}_{seed}_benchmark.log',
    log:
        std='brute_force_cn_calling/{sample}_{gene}_{seed}.log',
        err='brute_force_cn_calling/{sample}_{gene}_{seed}.err.log',
    params:
        prefix='brute_force_cn_calling/{sample}_{gene}_{seed}',
        amplicon_df=config['amplicon_parameters'],
        nrestarts=config['nrestarts'],
        max_nclones=config['max_nclones'],
        max_cn=config['max_cn']
    shell:
        'python mixed_NB_brute_force_EM_gene_level.py --sample {wildcards.sample} --readcounts {input.tsv_file} --amplicon {params.amplicon_df} --gene {wildcards.gene} '
        ' --maxclones {params.max_nclones} --maxcn {params.max_cn} --nrestarts {params.nrestarts} --seed {wildcards.seed} --prefix {params.prefix} ' 
        ' 1> {log.std} 2> {log.err}'
