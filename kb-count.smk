
'''
Snakefile for gene count quantification from scRNA-seq data by kb

Usage:
    snakemake -s kb-count.smk --configfile <path to config.yaml> --cores <int> --use-singularity
'''

import sys
import pandas as pd
workdir: config["workdir"]
kb_index_dir = config["kb_index_dir"]

'''experiment_table.tsv
sample	R1	R2
1122	/rhome/naotok/bigdata/ADKnowledgePortal/MC_snRNA/fastq/1122_S1_R1_001.fastq.gz	/rhome/naotok/bigdata/ADKnowledgePortal/MC_snRNA/fastq/1122_S1_R2_001.fastq.gz
1142	/rhome/naotok/bigdata/ADKnowledgePortal/MC_snRNA/fastq/1142_S1_R1_001.fastq.gz	/rhome/naotok/bigdata/ADKnowledgePortal/MC_snRNA/fastq/1142_S1_R2_001.fastq.gz
1160	/rhome/naotok/bigdata/ADKnowledgePortal/MC_snRNA/fastq/1160_S1_R1_001.fastq.gz	/rhome/naotok/bigdata/ADKnowledgePortal/MC_snRNA/fastq/1160_S1_R2_001.fastq.gz
1258	/rhome/naotok/bigdata/ADKnowledgePortal/MC_snRNA/fastq/1258_S1_R1_001.fastq.gz	/rhome/naotok/bigdata/ADKnowledgePortal/MC_snRNA/fastq/1258_S1_R2_001.fastq.gz
14543	/rhome/naotok/bigdata/ADKnowledgePortal/MC_snRNA/fastq/14543_S1_R1_001.fastq.gz	/rhome/naotok/bigdata/ADKnowledgePortal/MC_snRNA/fastq/14543_S1_R2_001.fastq.gz
'''

# Read experiment table to get sample names and fastq files
def read_experiment_table(file):
	df = pd.read_csv(file, sep="\t")
	samples = df["sample"].tolist()
	samples_dict = {}
	for i, row in df.iterrows():
		samples_dict[str(row["sample"])] = {"R1": row["R1"], "R2": row["R2"]}
	return samples, samples_dict

samples, samples_dict = read_experiment_table(config["experiment_table"])

rule all:
	input:
		adata_h5ad = expand("kb/{sample}/counts_unfiltered/adata.h5ad", sample = samples)

rule kb_count:
	wildcard_constraints:
		sample = "|".join([re.escape(str(x)) for x in samples])
	container:
		"docker://quay.io/biocontainers/kb-python:0.28.2--pyhdfd78af_2"
	input:
		R1 = lambda wildcards: samples_dict[wildcards.sample]["R1"],
		R2 = lambda wildcards: samples_dict[wildcards.sample]["R2"]
	output:
		adata_h5ad = "kb/{sample}/counts_unfiltered/adata.h5ad"
	params:
		kb_output = "kb/{sample}"
	threads:
		workflow.cores / 4
	benchmark:
		"benchmark/kb_count_{sample}.txt"
	log:
		"log/kb_count_{sample}.log"
	shell:
		"""
		kb count \
		-i {config[kb_index_dir]}/index.idx \
		-g {config[kb_index_dir]}/t2g.txt \
		-c1 {config[kb_index_dir]}/cdna_t2c.txt \
		-c2 {config[kb_index_dir]}/intron_t2c.txt \
		-x {config[technology]} \
		-o {params.kb_output} \
		-t {threads} \
		--workflow {config[workflow]} \
		--h5ad \
		--verbose \
		{input.R1} {input.R2} >& {log}
		"""
