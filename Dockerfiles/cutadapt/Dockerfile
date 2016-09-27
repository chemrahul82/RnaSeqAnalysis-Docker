# Ubuntu 14.04 with cutadapt 1.10 (http://cutadapt.readthedocs.io/en/stable/)
# and its prerequisites installed.
# Build image with:  docker build -t genomicsdocker/cutadapt:1.10 .


FROM ubuntu:14.04
MAINTAINER Rahul K. Das, https://github.com/chemrahul82

#working directory
RUN mkdir -p /usr/src/bin
WORKDIR /usr/src/bin

#get cutadapt & its prerequisites
RUN apt-get update \
&& apt-get -y install wget \
&& wget https://pypi.python.org/packages/15/a6/a05e99472b517aafd48824016f66458a31303f05256e9438ce9aec6b6bab/cutadapt-1.10.tar.gz \
&& tar -zxvf cutadapt-1.10.tar.gz \
&& apt-get install -y python-pip libpython2.7-dev \
&& cd cutadapt-1.10 \ 
&& python setup.py install \
&& apt-get clean 

#copy the wrapper script to workdir
COPY cutadapt.py /usr/src/bin/cutadapt.py

ENV PATH=$PATH:/usr/src/bin/

