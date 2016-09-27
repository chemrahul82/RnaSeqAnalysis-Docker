# samtools 0.1.19 (http://samtools.sourceforge.net)
# with python 2.7-slim as base image
# Build image with:  docker build -t genomicsdocker/samtools:0.1.19 .


FROM python:2.7-slim
MAINTAINER Rahul K. Das, https://github.com/chemrahul82

#working directory
RUN mkdir -p /usr/local/bin
WORKDIR /usr/local/bin

#install samtools-0.1.19
RUN apt-get update \ 
&& apt-get -y install samtools libpython2.7-dev 

#if want to make from zip
#&& apt-get -y install wget build-essential libpython2.7-dev zlib1g-dev libncurses5-dev \
#&& wget https://sourceforge.net/projects/samtools/files/samtools/1.3/samtools-1.3.tar.bz2 \
#&& tar jxf samtools-1.3.tar.bz2 \
#&& cd samtools-1.3 \
#&& make 

#copy the wrapper script to the workdir
COPY samtools_subsample.py /usr/local/bin/samtools_subsample.py
COPY samtools_sort.py /usr/local/bin/samtools_sort.py
COPY samtools_mapreads.py /usr/local/bin/samtools_mapreads.py
COPY samtools_depth.py /usr/local/bin/samtools_depth.py

#add the workdir to the path variable
ENV PATH=$PATH:/usr/local/bin
#ENV PATH=$PATH:/usr/local/bin:/usr/local/bin/samtools-1.3:/usr/local/bin/samtools-1.3/htslib-1.3

# default command for the container
CMD samtools 
