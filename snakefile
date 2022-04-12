#
configfile: "config.yaml"

rule all:
    input:
        #expand('data/{sample}_vaf.csv', sample=config['samples']),
        #expand('data/{sample}_sifit_snv_mat.txt', sample=config['samples']),
        expand('data_{paramSet}/{sample}_sifit_snv_mat.txt', sample=config['samples'], paramSet=config['paramSetList']),

rule mb_parser_default:
    output:
        vaf_dataframe="data/{sample}_sifit_snv_mat.txt",
        #vaf_dataframe="data/{sample}_vaf.csv",
    input:
        loom_file=lambda wildcards: config['loom'][wildcards.sample],
    params:
        output_prefix="data/{sample}",
    benchmark: "data/{sample}.benchmark",
    log:
        std='data/{sample}.log',
        err='data/{sample}.err.log'
    shell:
        "python mb_parser.py --loom {input.loom_file} --prefix {params.output_prefix}"
        " > {log.std} 2> {log.err}"

rule mb_parser_settings:
    output:
        vaf_dataframe="data_{paramSet}/{sample}_sifit_snv_mat.txt",
    input:
        loom_file=lambda wildcards: config['loom'][wildcards.sample],
    params:
        read_depth_threshold=lambda wildcards: config['read_depth_threshold'][wildcards.paramSet],
        vaf_threshold=lambda wildcards: config['vaf_threshold'][wildcards.paramSet],
        presence_threshold=lambda wildcards: config['presence_threshold'][wildcards.paramSet],
        output_prefix="data_{paramSet}/{sample}",
    benchmark: "data_{paramSet}/{sample}.benchmark"
    log:
        std='data_{paramSet}/{sample}.log',
        err='data_{paramSet}/{sample}.err.log'
    shell:
        "python mb_parser.py --loom {input.loom_file} --prefix {params.output_prefix} "
        " --depth {params.read_depth_threshold} --vaf {params.vaf_threshold} --presence {params.presence_threshold}"
        " > {log.std} 2> {log.err}"
