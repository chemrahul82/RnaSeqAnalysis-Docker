# bam2fastq 1.1.0 (https://gsl.hudsonalpha.org/information/software/bam2fastq)
# with python 2.7-slim as base image
# Build image with:  docker build -t genomicsdocker/bam2fastq:1.1.0 .


FROM python:2.7-slim
MAINTAINER Rahul K. Das, https://github.com/chemrahul82

#working directory
RUN mkdir -p /usr/local/bin
WORKDIR /usr/local/bin

#get cutadapt & its prerequisites
RUN apt-get update \
&& apt-get -y install wget \
&& wget https://gsl.hudsonalpha.org/static/software/bam2fastq-1.1.0.tgz \
&& tar -xzf bam2fastq-1.1.0.tgz \
&& apt-get install -y build-essential zlib1g-dev\
&& cd bam2fastq-1.1.0 \ 
&& make

#copy the wrapper script to workdir
COPY bam2fastq.py /usr/local/bin/bam2fastq.py

ENV PATH=$PATH:/usr/local/bin/:/usr/local/bin/bam2fastq-1.1.0

