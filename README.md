# ZhengLab_pipeline

This repository contains the scripts for the pipelines used in Zheng Lab.

## Fetching data

- `ZhengLab_Fetch.bash`: Fetching data from the server using `ffq`

```bash
ZhengLab_Fetch.bash -i <SRA/GEO ID> -o output_dir
```

## RNA-seq

- `ZhengLab_RNAseq.bash`: Quality check and trimming by `fastp` and mapping by `STAR`

```bash
ZhengLab_RNAseq.bash -i fastq_RNAseq.txt -c config_RNAseq.txt
```

## Docker image

- `Dockerfile`: Dockerfile for the RNA-seq pipeline
