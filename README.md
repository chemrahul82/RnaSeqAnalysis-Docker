This app executes the NGS RnaSeq data analysis workflow/pipeline. The workflows utilizes docker, thus making this app deployable/transferrable.   
Just install docker, download the Dockerfiles, build the images, and download the wrapper script, and you will be all set to run it without the need of installing any tools, dependencies etc.   

### This version 2.0.0b is a beta version and under development. 

**It has been thoroughly tested for IonTorrent single-end RnaSeq data. Testing with Illumina TruSeq paired-end is currently undergoing. There can be bugs. Use caution if using this app for your analyses. ** 

ChangeLog:
* Added support for analyses of Paired-end reads (Single-end was supported in v0.0.1)   
* Added support for analyses of RnaSeq data from non-ion torrent platforms (ion-torrent data was supported in v0.0.1)
* Added supports for analyses with mm10 reference genome without the need of any user-supplied reference fasta & other annotation files (hg19 reference was supported in v0.0.1)
* Added supports for user-defined strand specificity, library type, average fragment length
