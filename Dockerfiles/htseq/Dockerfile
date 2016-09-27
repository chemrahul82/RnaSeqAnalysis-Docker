# HTSeq 0.6.1 (http://www-huber.embl.de/HTSeq/doc/overview.html) 
# with Ubuntu 14.04 as base image and prerequisites installed.
# Build image with:  docker build -t genomicsdocker/htseq:0.6.1 .


FROM ubuntu:14.04
MAINTAINER Rahul K. Das, https://github.com/chemrahul82

#working directory
RUN mkdir -p /usr/local/bin
WORKDIR /usr/local/bin

#get HTSeq & its prerequisites
RUN apt-get update \
	&& apt-get -y --no-install-recommends install \
	build-essential \
	zlib1g-dev \
	libpython2.7-dev \
	python-pip \
	python-numpy \
	python-matplotlib \
	wget \
	dtrx \
	&& wget https://pypi.python.org/packages/72/0f/566afae6c149762af301a19686cd5fd1876deb2b48d09546dbd5caebbb78/HTSeq-0.6.1.tar.gz \
	&& dtrx HTSeq-0.6.1.tar.gz \
	&& pip install pysam \
	&& cd HTSeq-0.6.1 \
	&& python setup.py build \
	&& python setup.py install \
	&& apt-get remove -y wget dtrx \
	&& apt-get clean 

#copy the wrapper script to workdir
COPY htseq.py /usr/local/bin/htseq.py
COPY write2json_GeneCovMetrics.py /usr/local/bin/write2json_GeneCovMetrics.py


ENV PATH=$PATH:/usr/local/bin/

