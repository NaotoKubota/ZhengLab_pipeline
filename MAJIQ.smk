
'''
Snakefile for MAJIQ

Usage:
    snakemake -s Snakefile_MAJIQ --configfile config.yaml --cores <int> --use-singularity
'''

import os
from pathlib import Path
import pandas as pd
import configparser
workdir: config["workdir"]

'''experiment_table.tsv
sample	bam	group
sample_01	/path/to/sample_01/sample_01.bam	Ref
sample_02	/path/to/sample_02/sample_02.bam	Ref
sample_03	/path/to/sample_03/sample_03.bam	Ref
sample_04	/path/to/sample_04/sample_04.bam	Alt
sample_05	/path/to/sample_05/sample_05.bam	Alt
sample_06	/path/to/sample_06/sample_06.bam	Alt
'''

'''setting_file.ini
[info]
bamdirs=/path/to/sample_01,/path/to/sample_02,/path/to/sample_03,/path/to/sample_04,/path/to/sample_05,/path/to/sample_06
genome=hg19
strandness=None
[experiments]
Ref=sample_01,sample_02,sample_03
Alt=sample_04,sample_05,sample_06
'''

# Convert experiment_table.tsv to setting_file.ini
def convert_experiment_table_to_setting_file(experiment_table, genome, strandness):
    '''
    Convert experiment_table.tsv to setting_file.ini
    '''
    df = pd.read_csv(experiment_table, sep = '\t')
    bamdirs = ",".join(list(set([os.path.dirname(x) for x in df['bam'].tolist()])))
    Ref_bam = ",".join([os.path.splitext(os.path.basename(x))[0] for x in df[df['group'] == 'Ref']['bam'].tolist()])
    Alt_bam = ",".join([os.path.splitext(os.path.basename(x))[0] for x in df[df['group'] == 'Alt']['bam'].tolist()])
    with open("setting_file.ini", "w") as f:
        f.write(f"[info]\n")
        f.write(f"bamdirs={bamdirs}\n")
        f.write(f"genome={genome}\n")
        f.write(f"strandness={strandness}\n")
        f.write(f"[experiments]\n")
        f.write(f"Ref={Ref_bam}\n")
        f.write(f"Alt={Alt_bam}\n")

def parse_config_ini(config_ini):
    '''
    Parse config.ini file and return a dictionary
    '''
    config = configparser.ConfigParser()
    config.read(config_ini)
    samples_Ref = config["experiments"]["Ref"].split(",")
    samples_Alt = config["experiments"]["Alt"].split(",")
    return samples_Ref, samples_Alt

convert_experiment_table_to_setting_file(config["experiment_table"], config["genome"], config["strandness"])
samples_Ref, samples_Alt = parse_config_ini("setting_file.ini")

rule all:
    input:
        heterogen_tsv = "heterogen/Ref-Alt.het.tsv",
        summary = "voila_modulize/summary.tsv"

rule build:
    container:
        "/rhome/naotok/bigdata/singularity/majiq_0.1.sif"
    input:
        setting_file = "setting_file.ini",
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
