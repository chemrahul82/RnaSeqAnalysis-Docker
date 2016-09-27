# STAR 2.5.2a (https://github.com/alexdobin/STAR) & 
# with python 2.7-slim as base image
# Build image with:  docker build -t genomicsdocker/star:2.5.2a .


FROM python:2.7-slim
MAINTAINER Rahul K. Das, https://github.com/chemrahul82

#working directory
RUN mkdir -p /usr/local/bin
WORKDIR /usr/local/bin

#get star & its prerequisites
RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
	build-essential \
	zlib1g-dev \
	git \
	dtrx \
	wget \
	&& git clone https://github.com/alexdobin/STAR.git \
	&& apt-get remove -y git dtrx wget \
	&& apt-get clean
#&& cd STAR \
#&& make STAR

##copy the wrapper script to workdir
COPY STAR.py /usr/local/bin/STAR.py
COPY STAR_build.py /usr/local/bin/STAR.py

ENV PATH=$PATH:/usr/local/bin:/usr/local/bin/STAR/bin/Linux_x86_64
