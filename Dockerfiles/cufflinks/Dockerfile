# cufflinks 2.2.1 (Trapnel et al. Nature Biotechnology, 2010, http://cole-trapnell-lab.github.io/cufflinks/) 
# with Ubuntu 14.04 as base image and prerequisites installed.
# Build image with:  docker build -t genomicsdocker/cufflinks:2.2.1 .


FROM ubuntu:14.04
MAINTAINER Rahul K. Das, https://github.com/chemrahul82

#working directory
RUN mkdir -p /usr/local/bin
WORKDIR /usr/local/bin

#get HTSeq & its prerequisites
RUN apt-get update \
	&& apt-get -y install --no-install-recommends \
	wget dtrx \
	build-essential \
	zlib1g-dev \
	libpython2.7-dev \
	&& wget http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz \
	&& dtrx cufflinks-2.2.1.Linux_x86_64.tar.gz \
	&& apt-get remove -y wget dtrx \
	&& apt-get clean 

#copy the wrapper script to workdir
COPY cufflinks.py /usr/local/bin/cufflinks.py

ENV PATH=$PATH:/usr/local/bin/:/usr/local/bin/cufflinks-2.2.1.Linux_x86_64/
CMD cufflinks 
