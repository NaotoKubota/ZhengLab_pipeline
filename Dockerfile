#####################################################################
# Dockerfile to build container images for ZhengLab RNA-seq pipeline
# Based on python 3.11-buster
#####################################################################

FROM python:3.11-buster

# File Author / Maintainer
LABEL maintainer="Naoto Kubota <naotok@ucr.edu>"

# Install dependencies first
RUN apt-get -qq update && \
	apt-get -qq -y install \
	make wget less vim gcc git jq

# mkdir /src
RUN mkdir /src

# Install fastp (version 0.23.2)
RUN wget http://opengene.org/fastp/fastp.0.23.2 && \
	mv fastp.0.23.2 /src/fastp && \
	chmod a+x /src/fastp

# Install STAR (version 2.7.10b)
RUN wget https://github.com/alexdobin/STAR/archive/2.7.10b.tar.gz && \
	tar -xzf 2.7.10b.tar.gz && \
	mv STAR-2.7.10b /src/STAR-2.7.10b && \
	cd /src/STAR-2.7.10b/source && \
	make STAR

# Install Samtools (version 1.16.1)
RUN wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2 && \
	tar jxf samtools-1.16.1.tar.bz2 && \
	rm -rf samtools-1.16.1.tar.bz2 && \
	mv samtools-1.16.1 /src/samtools-1.16.1 && \
	cd /src/samtools-1.16.1/ && \
	./configure --prefix=/where/to/install && \
	make && \
	make install

# Install deepTools (version 3.5.1)
RUN pip install deeptools==3.5.1

# Install MultiQC (version 1.12)
RUN pip install multiqc==1.12

# git clone ZhengLab RNA-seq pipeline
RUN git clone https://github.com/NaotoKubota/ZhengLab_pipeline.git && \
	mv ZhengLab_pipeline /src/ZhengLab_pipeline && \
	chmod a+x /src/ZhengLab_pipeline/*.bash

# Set PATH
ENV PATH $PATH:/src:/src/samtools-1.16.1:/src/STAR-2.7.10b/source:/src/ZhengLab_pipeline

# Upgrade pip
RUN pip install --upgrade pip

# Install ffq (version 0.3.0)
RUN pip install ffq==0.3.0

# Set working directory
WORKDIR /home

# bash
CMD ["bash"]
