# Package STAR 2.5.2a (https://github.com/alexdobin/STAR) 
# with python 2.7-slim as base image
# download hg19 genome fasta and gene annotations file and
# run STAR to build genome index files
# Build image with:  docker build -t genomicsdocker/star_hg19:ion .


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
	git wget curl \
	&& git clone https://github.com/alexdobin/STAR.git \
	&& apt-get clean
#&& cd STAR \
#&& make STAR

ENV PATH=$PATH:/usr/local/bin:/usr/local/bin/STAR/bin/Linux_x86_64

#downloads the annotation of human genome (GRCh37), version 19 (Ensembl 74)
RUN wget -q ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz \
	&& gunzip gencode.v19.annotation.gtf.gz \
	&& mv gencode.v19.annotation.gtf gene.gtf 

##download hg19 reference fasta to make a facsimile of IonTorrent fasta
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

#run STAR to build index files & split large index files into chunks and put in seperate folders
RUN mkdir -p /usr/local/bin/star_idx_hg19 \
	&& STAR --runThreadN 16 --runMode genomeGenerate --genomeDir /usr/local/bin/star_idx_hg19 --genomeFastaFiles hg19.fasta --sjdbGTFfile gene.gtf --sjdbOverhang 75 

#copy the wrapper script
COPY STAR.py /usr/local/bin/STAR.py
