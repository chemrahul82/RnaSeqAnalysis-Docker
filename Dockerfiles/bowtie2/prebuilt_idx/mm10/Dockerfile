# package bowtie2 2.2.9 (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
# with python 2.7-slim as base image
# download mm10 genome fasta and run Bowtie2 to build genome index files
# Build image with:  docker build -t genomicsdocker/bowtie2_mm10:ion .


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

ENV PATH=$PATH:/usr/local/bin/:/usr/local/bin/bowtie2-2.2.9

### download mm10 reference: mouse genome (GRCm38), version 4 (Ensembl 78) ###
#Chr1-19
RUN for i in `seq 1 19`; do wget -q ftp://ftp.ensembl.org/pub/release-78/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_rm.chromosome.$i.fa.gz \
&& gunzip Mus_musculus.GRCm38.dna_rm.chromosome.$i.fa.gz \
&& sed "1s/.*/>chr$i/" Mus_musculus.GRCm38.dna_rm.chromosome.$i.fa > chr$i.fa \
&& rm Mus_musculus.GRCm38.dna_rm.chromosome.$i.fa;done 

#ChrX,Y
RUN for i in X Y; do wget -q ftp://ftp.ensembl.org/pub/release-78/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_rm.chromosome.$i.fa.gz \
&& gunzip Mus_musculus.GRCm38.dna_rm.chromosome.$i.fa.gz \
&& sed "1s/.*/>chr$i/" Mus_musculus.GRCm38.dna_rm.chromosome.$i.fa > chr$i.fa \
&& rm Mus_musculus.GRCm38.dna_rm.chromosome.$i.fa;done 

#chrM
RUN curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=34538597&rettype=fasta&retmode=text" | sed "1s/.*/>chrM/" > chrM.fa

#concatenate
RUN cat `ls  *.fa | sort -V | grep -i -v chrM ` chrM.fa | sed '/^$/d' > mm10.fasta

#remove individual files
RUN rm chr*

#run bowtie2 to build index files
RUN mkdir -p /usr/local/bin/idx_mm10 \
&& bowtie2-build --threads 16 mm10.fasta idx_mm10/bowtie2

#copy the wrapper script to workdir
COPY bowtie2.py /usr/local/bin/bowtie2.py


