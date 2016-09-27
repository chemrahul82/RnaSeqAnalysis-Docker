# FastQC 0.11.5 (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
# with Ubuntu 14.04 as base image and Java 8 and Python 2.7 installed.
# Build image with:  docker build -t genomicsdocker/fastqc:0.11.5 .


FROM ubuntu:14.04
MAINTAINER Rahul K. Das, https://github.com/chemrahul82

#working directory
RUN mkdir -p /usr/local/bin
WORKDIR /usr/local/bin

#install python and oracle java 8
RUN apt-get update \
&& apt-get upgrade -y \
&& apt-get install -y python \

&& apt-get -y install software-properties-common \
&& add-apt-repository ppa:webupd8team/java -y \
&& apt-get update \
&& echo oracle-java7-installer shared/accepted-oracle-license-v1-1 select true | /usr/bin/debconf-set-selections \
&& apt-get -y install oracle-java8-installer \
&& apt-get -y install wget zip \
&& apt-get clean

#get the picard tool
RUN wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip \
&& unzip fastqc_v0.11.5.zip \
&& cd FastQC \
&& chmod 755 fastqc 


#copy the wrapper sctipts to the work directory


#set the path
ENV PATH=$PATH:/usr/local/bin/FastQC

CMD fastqc --help
