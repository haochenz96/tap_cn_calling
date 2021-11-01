#
configfile: "config.yaml"

rule all:
    input:
        expand('data/{sample}_vaf.csv', sample=config['samples']),

rule mb_parser:
    output:
        vaf_dataframe="data/{sample}_vaf.csv",
    input:
        loom_file=lambda wildcards: config['loom'][wildcards.sample],
    params:
        output_prefix="data/{sample}",
    shell:
        "python mb_parser.py --loom {input.loom_file} --prefix {params.output_prefix}"
