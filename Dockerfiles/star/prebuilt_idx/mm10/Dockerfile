# Package STAR 2.5.2a (https://github.com/alexdobin/STAR) 
# with python 2.7-slim as base image
# download mm10 genome fasta and gene annotations file and
# run STAR to build genome index files
# Build image with:  docker build -t genomicsdocker/star_mm10:ion .


FROM python:2.7-slim
MAINTAINER Rahul K. Das, https://github.com/chemrahul82

#working directory
RUN mkdir -p /usr/local/bin
WORKDIR /usr/local/bin

#get star & its prerequisites
RUN apt-get update \
&& apt-get install -y build-essential zlib1g-dev git dtrx wget curl \
&& git clone https://github.com/alexdobin/STAR.git 
#&& cd STAR \
#&& make STAR

ENV PATH=$PATH:/usr/local/bin:/usr/local/bin/STAR/bin/Linux_x86_64

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

### gene annotation file ###
##download the annotation of mouse genome (GRCm38), version 4 (Ensembl 78)
RUN wget -q ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M4/gencode.vM4.annotation.gtf.gz \
&& gunzip gencode.vM4.annotation.gtf.gz \
&& mv gencode.vM4.annotation.gtf gene.gtf

#run STAR to build index files
RUN mkdir -p /usr/local/bin/star_idx_mm10 \
&& STAR --runThreadN 16 --runMode genomeGenerate --genomeDir /usr/local/bin/star_idx_mm10 --genomeFastaFiles mm10.fasta --sjdbGTFfile gene.gtf --sjdbOverhang 75

#copy the wrapper script
COPY STAR.py /usr/local/bin/STAR.py

