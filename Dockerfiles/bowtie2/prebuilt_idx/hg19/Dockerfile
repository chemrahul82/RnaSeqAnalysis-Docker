# package bowtie2 2.2.9 (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
# with python 2.7-slim as base image
# download hg19 genome fasta and run Bowtie2 to build genome index files
# Build image with:  docker build -t genomicsdocker/bowtie2_hg19:ion .


FROM python:2.7-slim
MAINTAINER Rahul K. Das, https://github.com/chemrahul82

#working directory
RUN mkdir -p /usr/local/bin
WORKDIR /usr/local/bin

#get bowtie2 & its prerequisites
RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
	libpython2.7-dev \
	wget unzip curl \
	&& wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.9/bowtie2-2.2.9-linux-x86_64.zip \
	&& unzip bowtie2-2.2.9-linux-x86_64.zip \
	&& apt-get clean

ENV PATH=$PATH:/usr/local/bin/:/usr/local/bin/bowtie2-2.2.9

##download hg19 reference fasta
#Chr1-22
RUN for i in `seq 1 22`; do wget -q ftp://ftp.ensembl.org/pub/release-56/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.56.dna.chromosome.$i.fa.gz \
	&& gunzip Homo_sapiens.GRCh37.56.dna.chromosome.$i.fa.gz \
	&& sed "1s/.*/>chr$i/" Homo_sapiens.GRCh37.56.dna.chromosome.$i.fa > chr$i.fa \
	&& rm Homo_sapiens.GRCh37.56.dna.chromosome.$i.fa;done \
#ChrX,Y
	&& for i in X Y; do wget -q ftp://ftp.ensembl.org/pub/release-56/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.56.dna.chromosome.$i.fa.gz \
	&& gunzip Homo_sapiens.GRCh37.56.dna.chromosome.$i.fa.gz \
	&& sed "1s/.*/>chr$i/" Homo_sapiens.GRCh37.56.dna.chromosome.$i.fa > chr$i.fa \
	&& rm Homo_sapiens.GRCh37.56.dna.chromosome.$i.fa;done \
#chrM
	&& curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=251831106&rettype=fasta&retmode=text" | sed "1s/.*/>chrM/" > chrM.fa \
#concatenate
	&& cat `ls  *.fa | sort -V | grep -i -v chrM ` chrM.fa | sed '/^$/d' > hg19.fasta \
#remove individual files
	&& rm chr*

#RUN wget -q ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/GRCh37.p13.genome.fa.gz \

#run bowtie2 to build index files
RUN mkdir -p /usr/local/bin/idx_hg19 \
	&& bowtie2-build --threads 16 hg19.fasta idx_hg19/bowtie2

#copy the wrapper script to workdir
COPY bowtie2.py /usr/local/bin/bowtie2.py
