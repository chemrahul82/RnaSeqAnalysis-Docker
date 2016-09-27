# Picard 2.5.0 (http://broadinstitute.github.io/picard/)
# with Ubuntu 14.04 as base image and Java 8 and Python 2.7 installed.
# Build image with:  docker build -t genomicsdocker/picard:2.5.0 .


FROM ubuntu:14.04
MAINTAINER Rahul K. Das, https://github.com/chemrahul82

#working directory
RUN mkdir -p /usr/local/bin
WORKDIR /usr/local/bin

#install python and oracle java 8
RUN apt-get update \
	&& apt-get upgrade -y \
	&& apt-get install -y --no-install-recommends python python-matplotlib \
	software-properties-common \
	&& add-apt-repository ppa:webupd8team/java -y \
	&& apt-get update \
	&& echo oracle-java7-installer shared/accepted-oracle-license-v1-1 select true | /usr/bin/debconf-set-selections \
	&& apt-get -y install --no-install-recommends oracle-java8-installer \ 
	&& apt-get install -y --no-install-recommends \
	wget unzip \
	&& wget https://github.com/broadinstitute/picard/releases/download/2.5.0/picard-tools-2.5.0.zip \
	&& unzip picard-tools-2.5.0.zip \
	&& cp /usr/local/bin/picard-tools-2.5.0/picard.jar /usr/local/bin/picard.jar \
	&& apt-get autoremove \
	&& apt-get clean \
	&& apt-get autoclean 


#copy the wrapper sctipts to the work directory
COPY picard_MergeSamFiles.py /usr/local/bin/picard_MergeSamFiles.py
COPY picard_CollectAlignmentSummaryMetrics.py /usr/local/bin/picard_CollectAlignmentSummaryMetrics.py 
COPY picard_CollectRnaSeqMetrics.py /usr/local/bin/picard_CollectRnaSeqMetrics.py
COPY write2json_AlignmentMetrics.py /usr/local/bin/write2json_AlignmentMetrics.py
COPY results_RnaSeqMetrics.py /usr/local/bin/results_RnaSeqMetrics.py


#set the path
ENV PATH=$PATH:/usr/local/bin
