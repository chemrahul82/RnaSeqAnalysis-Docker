# python 2.7-slim as base image
# mm10 genome fasta, gene annotations file and other supplementary 
# annotations files needed for RnaSeq analysis
# Build image with:  docker build -t genomicsdocker/mm10:ion .


FROM python:2.7-slim
MAINTAINER Rahul K. Das, https://github.com/chemrahul82

#working directory
RUN mkdir -p /usr/local/bin
WORKDIR /usr/local/bin

#get prerequisites
RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
	dtrx wget curl 


### mm10 reference: mouse genome (GRCm38), version 4 (Ensembl 78) ###
#Chr1-19
RUN for i in `seq 1 19`; do wget -q ftp://ftp.ensembl.org/pub/release-78/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_rm.chromosome.$i.fa.gz \
	&& gunzip Mus_musculus.GRCm38.dna_rm.chromosome.$i.fa.gz \
	&& sed "1s/.*/>chr$i/" Mus_musculus.GRCm38.dna_rm.chromosome.$i.fa > chr$i.fa \
	&& rm Mus_musculus.GRCm38.dna_rm.chromosome.$i.fa;done \ 
#ChrX,Y
	&& for i in X Y; do wget -q ftp://ftp.ensembl.org/pub/release-78/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_rm.chromosome.$i.fa.gz \
	&& gunzip Mus_musculus.GRCm38.dna_rm.chromosome.$i.fa.gz \
	&& sed "1s/.*/>chr$i/" Mus_musculus.GRCm38.dna_rm.chromosome.$i.fa > chr$i.fa \
	&& rm Mus_musculus.GRCm38.dna_rm.chromosome.$i.fa;done \
#chrM
	&& curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=34538597&rettype=fasta&retmode=text" | sed "1s/.*/>chrM/" > chrM.fa \
#concatenate
	&& cat `ls  *.fa | sort -V | grep -i -v chrM ` chrM.fa | sed '/^$/d' > mm10.fasta \
#remove individual files
	&& rm chr*


### gene annotation file ###
#download the annotation of mouse genome (GRCm38), version 4 (Ensembl 78)
RUN wget -q ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M4/gencode.vM4.annotation.gtf.gz \
	&& gunzip gencode.vM4.annotation.gtf.gz \
	&& mv gencode.vM4.annotation.gtf gene.gtf 

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
	&& curl -s "http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes" > mm10.chrom.sizes \
	&& for i in `seq 1 19`; do egrep -w "chr$i" mm10.chrom.sizes >> temp; done \
	&& for i in X Y M; do egrep -w "chr$i" mm10.chrom.sizes >> temp; done \
	&& echo "@HD VN:1.0  GO:none SO:coordinate" > rRNA.interval \
	&& awk 'BEGIN {OFS="\t"}; {print "@SQ","SN:"$1,"LN:"$2}' temp >> rRNA.interval \
	&& egrep -w 'transcript_type "rRNA"|transcript_type "Mt_rRNA"' gene.gtf |awk -F'[/;]' -v OFS="\t" '{print $1,$5}' | awk 'BEGIN {OFS="\t"}; {print $1,$4,$5,$7,$12}' | uniq >> rRNA.interval \
	&& rm temp

## make xrRNA.fasta
#18S rRNA
RUN curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=374088232&rettype=fasta&retmode=text" | sed '/^$/d' > 18S_rRNA.fa \
	&& sleep 60s \
#28S rRNA
	&& curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=120444900&rettype=fasta&retmode=text" | sed '/^$/d' > 28S_rRNA.fa \
	&& sleep 60s \
#5.8S rRNA
	&& curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=374093199&rettype=fasta&retmode=text" | sed '/^$/d' > 5.8S_rRNA.fa \
	&& sleep 60s \
#concatenate
	&& cat 5.8S_rRNA.fa 18S_rRNA.fa 28S_rRNA.fa > xrRNA.fasta \
	&& rm 5.8S_rRNA.fa 18S_rRNA.fa 28S_rRNA.fa

#copy the wrapper script
COPY mm10_annot.py /usr/local/bin/mm10_annot.py
