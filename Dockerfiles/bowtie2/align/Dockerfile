# bowtie2 2.2.9 (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
# with python 2.7-slim as base image
# Build image with:  docker build -t genomicsdocker/bowtie2:2.2.9 .


FROM python:2.7-slim
MAINTAINER Rahul K. Das, https://github.com/chemrahul82

#working directory
RUN mkdir -p /usr/local/bin
WORKDIR /usr/local/bin

#get cutadapt & its prerequisites
RUN apt-get update \
&& apt-get install -y libpython2.7-dev wget unzip dtrx \
&& wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.9/bowtie2-2.2.9-linux-x86_64.zip \
&& unzip bowtie2-2.2.9-linux-x86_64.zip 

#copy the wrapper script to workdir
COPY bowtie2.py /usr/local/bin/bowtie2.py
COPY bowtie2_build.py /usr/local/bin/bowtie2_build.py

ENV PATH=$PATH:/usr/local/bin/:/usr/local/bin/bowtie2-2.2.9

