#
configfile: "phylogeny_config.yaml"

import numpy as np

def get_ncells(sample):
    return np.loadtxt(f'data/{sample}_scite_snv_mat.txt').shape[1]

def get_npos(sample):
    return np.loadtxt(f'data/{sample}_scite_snv_mat.txt').shape[0]

rule all:
    input:
        #expand('scite/{sample}_output_ml0.newick', sample=config['samples']),
        expand('phiscs-B/{sample}_phiscs_snv_mat.CFMatrix', sample=config['samples']),
        expand('phiscs-B_2hr/{sample}_phiscs_snv_mat.CFMatrix', sample=config['samples']),
        #expand('sifit/{sample}_sifit_snv_mat_mlTree.newick', sample=config['samples']),
        #expand('phiscs-I/{sample}_phiscs_snv_mat.CFMatrix', sample=config['samples']),

rule scite:
    output:
        newick_file="scite/{sample}_output_ml0.newick",
    input:
        snv_file="data/{sample}_scite_snv_mat.txt",
    params:
        output_prefix="scite/{sample}_output",
        ncells=lambda wildcards: get_ncells(wildcards.sample),
        npos=lambda wildcards: get_npos(wildcards.sample),
        cc=config['cc'],
        fd=config['fd'],
        ad=config['ad'],
        niter=config['niter'],
        nrestarts=config['nrestarts'],
    benchmark: "scite/{sample}.benchmark",
    log:
        std = "scite/{sample}.log",
        err = "scite/{sample}.err.log",
    shell:
        "scite/scite -i {input.snv_file} -n {params.ncells} -m {params.npos} -o {params.output_prefix}"
        " -cc {params.cc} -l {params.niter} -r {params.nrestarts} -fd {params.fd} -ad {params.ad} {params.ad} "
        " > {log.std} 2> {log.err}"

rule sifit:
    output:
        newick_file="sifit/{sample}_sifit_snv_mat_mlTree.newick",
    input:
        snv_file="data/{sample}_sifit_snv_mat.txt",
    params:
        ncells=lambda wildcards: get_ncells(wildcards.sample),
        npos=lambda wildcards: get_npos(wildcards.sample),
        fp=config['fp'],
        fn=config['fn'],
        niter=config['niter'],
        nrestarts=config['nrestarts'],
    benchmark: "sifit/{sample}.benchmark",
    log:
        std = "sifit/{sample}.log",
        err = "sifit/{sample}.err.log",
    shell:
        "java -jar ../SiFit/SiFit.jar -ipMat {input.snv_file} -m {params.ncells} -n {params.npos} -df 1 "
        " -l {params.niter} -r {params.nrestarts} -fp {params.fp} -fn {params.fn} "
        " > {log.std} 2> {log.err}; mv data/{wildcards.sample}_sifit_snv_mat_mlTree.newick sifit/."


rule phiscsB:
    output:
        cfmatrix_file="phiscs-B/{sample}_phiscs_snv_mat.CFMatrix",
    input:
        snv_file="data/{sample}_phiscs_snv_mat.txt",
    params:
        output_dir="phiscs-B/",
        fp=config['fp'],
        fn=config['fn'],
        kmax=config['kmax'],
    benchmark: "phiscs-B/{sample}.benchmark",
    log:
        std = "phiscs-B/{sample}.log",
        err = "phiscs-B/{sample}.err.log",
    shell:
        "python /n/fs/ragr-data/users/palash/PhISCS/PhISCS-B -SCFile {input.snv_file} -fn {params.fn} -fp {params.fp} -o {params.output_dir} -kmax {params.kmax} --drawTree "
        " > {log.std} 2> {log.err}"

rule phiscsB_2hr:
    output:
        cfmatrix_file="phiscs-B_2hr/{sample}_phiscs_snv_mat.CFMatrix",
    input:
        snv_file="data/{sample}_phiscs_snv_mat.txt",
    params:
        output_dir="phiscs-B_2hr/",
        fp=config['fp'],
        fn=config['fn'],
        kmax=config['kmax'],
    benchmark: "phiscs-B_2hr/{sample}.benchmark",
    log:
        std = "phiscs-B_2hr/{sample}.log",
        err = "phiscs-B_2hr/{sample}.err.log",
    shell:
        "python /n/fs/ragr-data/users/palash/PhISCS/PhISCS-B -SCFile {input.snv_file} -fn {params.fn} -fp {params.fp} -o {params.output_dir} -time 7200 -kmax {params.kmax} --drawTree "
        " > {log.std} 2> {log.err}"

rule phiscsI:
    output:
        cfmatrix_file="phiscs-I/{sample}_phiscs_snv_mat.CFMatrix",
    input:
        snv_file="data/{sample}_phiscs_snv_mat.txt",
    params:
        output_dir="phiscs-I/",
        fp=config['fp'],
        fn=config['fn'],
        kmax=config['kmax'],
    benchmark: "phiscs-I/{sample}.benchmark",
    log:
        std = "phiscs-I/{sample}.log",
        err = "phiscs-I/{sample}.err.log",
    shell:
        "python /n/fs/ragr-data/users/palash/PhISCS/PhISCS-I -SCFile {input.snv_file} -fn {params.fn} -fp {params.fp} -o {params.output_dir} -kmax {params.kmax} "
        " > {log.std} 2> {log.err}"

