# python 2.7-slim as base image
# hg19 genome fasta, gene annotations file, and other
# supplementary annotation files for RnaSeq analysis
# Build image with:  docker build -t genomicsdocker/hg19:ion .


FROM python:2.7-slim
MAINTAINER Rahul K. Das, https://github.com/chemrahul82

#working directory
RUN mkdir -p /usr/local/bin
WORKDIR /usr/local/bin

#get prerequisites
RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
	dtrx wget curl 


### hg19 reference ###
## make a facsimile of IonTorrent fasta
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


### gene annotation file ###
#download the annotation of human genome (GRCh37), version 19 (Ensembl 74)
RUN wget -q ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz \
	&& dtrx gencode.v19.annotation.gtf.gz \
	&& mv gencode.v19.annotation.gtf gene.gtf \
	&& rm gencode.v19.annotation.gtf.gz 

###make supplementary annotation files###
##use the above gene.gtf file to create the following files

##rRNA_mask.gtf
RUN grep 'gene_type "rRNA"' gene.gtf > rRNA_mask.gtf

##refFlat
RUN wget -q http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred \
	&& chmod +x gtfToGenePred \
	&& ./gtfToGenePred -ignoreGroupsWithoutExons -genePredExt gene.gtf tmp \
	&& awk 'BEGIN{FS="\t"};{print $12"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' tmp > refFlat \
	&& rm tmp 

##make rRNA.interval file
RUN if [ -f rRNA.interval ]; then rm rRNA.interval; fi \
	&& curl -s "https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes" > hg19.chrom.sizes \
	&& for i in `seq 1 22`; do egrep -w "chr$i" hg19.chrom.sizes >> temp; done \
	&& for i in X Y; do egrep -w "chr$i" hg19.chrom.sizes >> temp; done \
	&& echo "@HD VN:1.0  GO:none SO:coordinate" > rRNA.interval \
	&& awk 'BEGIN {OFS="\t"}; {print "@SQ","SN:"$1,"LN:"$2}' temp >> rRNA.interval \
	&& printf '@SQ\tSN:chrM\tLN:16569\n' >> rRNA.interval \
	&& egrep -w 'transcript_type "rRNA"|transcript_type "Mt_rRNA"' gene.gtf |awk -F'[/;]' -v OFS="\t" '{print $1,$5}' | awk 'BEGIN {OFS="\t"}; {print $1,$4,$5,$7,$12}' | uniq >> rRNA.interval \
	&& rm temp

## make xrRNA.fasta
#18S rRNA
RUN curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=225637497&rettype=fasta&retmode=text" | sed '/^$/d' > 18S_rRNA.fa \
	&& sleep 60s \
#28S rRNA
	&& curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=225637499&rettype=fasta&retmode=text" | sed '/^$/d' > 28S_rRNA.fa \
	&& sleep 60s \
#5.8S rRNA
	&& curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=142372596&rettype=fasta&retmode=text" | sed '/^$/d' > 5.8S_rRNA.fa \
	&& sleep 60s \
	&& cat 5.8S_rRNA.fa 18S_rRNA.fa 28S_rRNA.fa > xrRNA.fasta \
	&& rm 5.8S_rRNA.fa 18S_rRNA.fa 28S_rRNA.fa


#copy the wrapper script
COPY hg19_annot.py /usr/local/bin/hg19_annot.py
