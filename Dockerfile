#####################################################################
# Dockerfile to build container images for ZhengLab Fetch
# Based on python 3.11-buster
#####################################################################

FROM python:3.11-buster

# File Author / Maintainer
LABEL maintainer="Naoto Kubota <naotok@ucr.edu>"

# Install dependencies first
RUN apt-get -qq update && \
	apt-get -qq -y install \
	make wget less vim gcc git jq aria2

# mkdir /src
RUN mkdir /src

# git clone ZhengLab RNA-seq pipeline
RUN git clone https://github.com/NaotoKubota/ZhengLab_pipeline.git && \
	mv ZhengLab_pipeline /src/ZhengLab_pipeline && \
	chmod a+x /src/ZhengLab_pipeline/*.bash

# Set PATH
ENV PATH $PATH:/src:/src/ZhengLab_pipeline

# Upgrade pip
RUN pip install --upgrade pip

# Install ffq (version 0.3.1)
RUN pip install ffq==0.3.1

# Set working directory
WORKDIR /home

# bash
CMD ["bash"]
