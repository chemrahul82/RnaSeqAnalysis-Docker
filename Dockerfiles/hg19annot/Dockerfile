# hg19 annotation files required for RnaSeq workflow
# python 2.7-slim base image
# Build image with:  docker build -t genomicsdocker/hg19annot .

#FROM ubuntu:14.04
FROM python:2.7-slim
MAINTAINER Rahul K. Das, https://github.com/chemrahul82

#working directory
RUN mkdir -p /usr/local/bin
WORKDIR /usr/local/bin


#copy the annotation files to workdir
COPY refFlat /usr/local/bin/refFlat
COPY gene.gtf /usr/local/bin/gene.gtf
COPY rRNA_mask.gtf /usr/local/bin/rRNA_mask.gtf
COPY rRNA.interval /usr/local/bin/rRNA.interval
COPY xrRNA.fasta /usr/local/bin/xrRNA.fasta
COPY hg19.fasta /usr/local/bin/hg19.fasta 


#copy the wrapper script to workdir
COPY hg19_annot.py /usr/local/bin/hg19_annot.py 

ENV PATH=$PATH:/usr/local/bin/

