
'''
Snakefile for MAJIQ

Usage:
    snakemake -s Snakefile_MAJIQ --configfile config.yaml --cores <int> --use-singularity
'''

import configparser
def parse_config_ini(config_ini):
    '''
    Parse config.ini file and return a dictionary
    '''
    config = configparser.ConfigParser()
    config.read(config_ini)
    samples_Ref = config["experiments"]["Ref"].split(",")
    samples_Alt = config["experiments"]["Alt"].split(",")
    return samples_Ref, samples_Alt

samples_Ref, samples_Alt = parse_config_ini(config["setting_file"])

workdir: config["workdir"]

rule all:
    input:
        heterogen_tsv = "heterogen/Ref-Alt.het.tsv",
        summary = "voila_modulize/summary.tsv"

rule build:
    container:
        "/rhome/naotok/bigdata/singularity/majiq_0.1.sif"
    input:
        setting_file = config["setting_file"],
        gff = config["gff"]
    output:
        splicegraph = "build/splicegraph.sql",
        majiq_ref = expand("build/{sample_ref}.majiq", sample_ref = samples_Ref),
        majiq_alt = expand("build/{sample_alt}.majiq", sample_alt = samples_Alt)
    threads:
        workflow.cores
    benchmark:
        "benchmark/build.txt"
    log:
        "log/build.log"
    shell:
        """
        majiq build \
        -c {input.setting_file} \
        {input.gff} \
        -o build \
        -j {threads} \
        > {log} 2>&1
        """

rule heterogen:
    container:
        "/rhome/naotok/bigdata/singularity/majiq_0.1.sif"
    input:
        majiq_ref = expand("build/{sample_ref}.majiq", sample_ref = samples_Ref),
        majiq_alt = expand("build/{sample_alt}.majiq", sample_alt = samples_Alt)
    output:
        voila = "heterogen/Ref-Alt.het.voila"
    params:
        majiq_ref = " ".join(["build/" + sample + ".majiq" for sample in samples_Ref]),
        majiq_alt = " ".join(["build/" + sample + ".majiq" for sample in samples_Alt])
    threads:
        workflow.cores
    benchmark:
        "benchmark/heterogen.txt"
    log:
        "log/heterogen.log"
    shell:
        """
        majiq heterogen \
        -j {threads} \
        -o heterogen \
        -n Ref Alt \
        -grp1 {params.majiq_ref} \
        -grp2 {params.majiq_alt} \
        > {log} 2>&1
        """

rule voila_tsv:
    container:
        "/rhome/naotok/bigdata/singularity/majiq_0.1.sif"
    input:
        voila = "heterogen/Ref-Alt.het.voila",
        splicegraph = "build/splicegraph.sql"
    output:
        tsv = "heterogen/Ref-Alt.het.tsv"
    threads:
        workflow.cores
    benchmark:
        "benchmark/voila_tsv.txt"
    log:
        "log/voila_tsv.log"
    shell:
        """
        voila tsv \
        -j {threads} \
        -f {output.tsv} \
        {input.splicegraph} \
        {input.voila} \
        > {log} 2>&1
        """

rule voila_modulize:
    container:
        "/rhome/naotok/bigdata/singularity/majiq_0.1.sif"
    input:
        splicegraph = "build/splicegraph.sql",
        voila = "heterogen/Ref-Alt.het.voila"
    output:
        summary = "voila_modulize/summary.tsv"
    threads:
        workflow.cores
    benchmark:
        "benchmark/voila_modulize.txt"
    log:
        "log/voila_modulize.log"
    shell:
        """
        voila modulize \
        -j {threads} \
        -d voila_modulize \
        --show-all \
        {input.splicegraph} \
        {input.voila} \
        > {log} 2>&1
        """
