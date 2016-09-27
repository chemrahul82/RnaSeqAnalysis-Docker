#!/usr/bin/env python

"""wrapper that strings together all of the pieces/modules of the RnaSeqAnalysis workflow"""

__author__  = "Rahul K. Das"
__DateCreated__    = "July 26, 2016"
__email__   = "chemrahul82@gmail.com"
__status__  ="Production"


import argparse
import json
import logging
import os
import subprocess
import sys
import time
import glob



def run_command_line(system_call):
    
    # run the call and return the status
    logging.info(time.asctime( time.localtime(time.time()) ) +': Starting %s' %(system_call))
    
    try:
        output = subprocess.check_output(system_call, shell=True)
        logging.info(output)
    except subprocess.CalledProcessError as e:
        logging.exception(e)
        sys.exit(1)
    
    # capture the finishing time
    logging.info(time.asctime( time.localtime(time.time()) ) +': Completed %s' %(system_call))
    
    # return the output
    return output

           
def copy_annotations(options):
    """copy directory with hg19 annotation files from container to the host""" 
    
    system_call = 'docker run --rm -v %s:%s:ro -v %s:%s genomicsdocker/hg19annot hg19_annot.py --outdir %s' %(
            '/etc/localtime','/etc/localtime',
            options.outdir,'/workdir',
            '/workdir')

    run_command_line(system_call)
    

def bam_downsample_2fastq(outdir,inpbam,outbam,outfastq,fracread):
    """downsample input bam & then convert to fastq"""
    
    #downsample
    system_call = 'docker run --rm -v %s:%s:ro -v %s:%s:ro -v %s:%s genomicsdocker/samtools:0.1.19 samtools_subsample.py --inpbam %s --outbam %s --fracread %s' %(
            '/etc/localtime','/etc/localtime',
            os.path.dirname(inpbam),'/inpdir',
            outdir,'/workdir',
            os.path.join('/inpdir',os.path.basename(inpbam)),
            os.path.join('/workdir',os.path.basename(outbam)),
            fracread)
            
    run_command_line(system_call)

    #bam to fastq
    system_call = 'docker run --rm -v %s:%s:ro -v %s:%s genomicsdocker/bam2fastq:1.1.0 bam2fastq.py --inpbam %s --outfastq %s' %(
            '/etc/localtime','/etc/localtime',
            outdir,'/workdir',
            os.path.join('/workdir',os.path.basename(outbam)),
            os.path.join('/workdir',os.path.basename(outfastq)))

    run_command_line(system_call)

        
def bam2fastq(outdir,inpbam,outfastq):
    """convert input bam to fastq without downsampling"""

    system_call = 'docker run --rm -v %s:%s:ro -v %s:%s:ro -v %s:%s genomicsdocker/bam2fastq:1.1.0 bam2fastq.py --inpbam %s --outfastq %s' %(
            '/etc/localtime','/etc/localtime',
            os.path.dirname(inpbam),'/inpdir',
            outdir,'/workdir',
            os.path.join('/inpdir',os.path.basename(inpbam)),
            os.path.join('/workdir',os.path.basename(outfastq)))

    run_command_line(system_call)
    
        
def run_cutadapt(outdir,inpfastq,outfastq,adapter):
    """adapter trimming"""

    system_call = 'docker run --rm -v %s:%s:ro -v %s:%s:ro -v %s:%s genomicsdocker/cutadapt:1.10 cutadapt.py --adapter %s --inpfastq %s --outfastq %s' %(
            '/etc/localtime','/etc/localtime',
            os.path.dirname(inpfastq),'/inpdir',
            outdir,'/workdir',
            adapter,
            os.path.join('/inpdir',os.path.basename(inpfastq)),
            os.path.join('/workdir',os.path.basename(outfastq)))

    run_command_line(system_call)


def star_align(outdir,inpfastq,idxDir,options):
    """STAR alignment"""

    #non-hg19 reference
    if (options.reference != None) and (not options.STARidx):
        #need to implement this
        logging.info('Running STAR to build index files of user provided reference genome')
    
    #default hg19 reference; run STAR with prebuilt idx files
    elif not options.reference:
        system_call = 'docker run --rm -v %s:%s:ro -v %s:%s genomicsdocker/star_hg19:ion STAR.py --inpfastq %s --outDir %s --ncpu %s' %(
                '/etc/localtime','/etc/localtime',
                outdir,'/workdir',
                os.path.join('/workdir',os.path.basename(inpfastq)),
                '/workdir',
                options.ncpu)
            
        run_command_line(system_call)
    
    
    #user provided non-hg19 reference and index files
    else:
        system_call = 'docker run --rm -v %s:%s:ro -v %s:%s -v %s:%s:ro genomicsdocker/star:2.5.2a STAR.py --inpfastq %s --idxDir %s --outDir %s --ncpu %s' %(
                '/etc/localtime','/etc/localtime',
                outdir,'/workdir',
                os.path.dirname(idxDir),'/genomeidxDir',
                os.path.join('/workdir',os.path.basename(inpfastq)),
                '/genomeidxDir',
                '/workdir',
                options.ncpu)
            
        run_command_line(system_call)

    
def bowtie2_align(outdir,idxDir,idxBase,inpfastq,outsam,outun,ncpu,options):
    """bowtie2 alignment"""

    #non-hg19 reference
    if (options.reference != None) and (not options.bowtie2idx):
        #need to implement this
        logging.info('Running bowtie2 to build index files of user-provided reference genome and rRNA sequence')

    #default hg19 reference; run bowtie2 alignment with pre-built idx files
    elif not options.reference:
        system_call = 'docker run --rm -v %s:%s:ro -v %s:%s genomicsdocker/bowtie2_hg19:ion bowtie2.py --inpfastq %s --outsam %s --outun %s --ncpu %s' %(
                '/etc/localtime','/etc/localtime',
                outdir,'/workdir',
                os.path.join('/workdir',os.path.basename(inpfastq)),
                os.path.join('/workdir',os.path.basename(outsam)),
                os.path.join('/workdir',os.path.basename(outun)),
                ncpu)
    
        run_command_line(system_call)


    #user provided non-hg19 reference and index files
    else:
        system_call = 'docker run --rm -v %s:%s:ro -v %s:%s -v %s:%s:ro genomicsdocker/bowtie2:2.2.9 bowtie2.py --inpfastq %s --index %s --outsam %s --outun %s --ncpu %s' %(
                '/etc/localtime','/etc/localtime',
                outdir,'/workdir',
                idxDir,'/genomeidxDir',
                os.path.join('/workdir',os.path.basename(inpfastq)),
                os.path.join('/genomeidxDir',idxBase),
                os.path.join('/workdir',os.path.basename(outsam)),
                os.path.join('/workdir',os.path.basename(outun)),
                ncpu)
    
        run_command_line(system_call)


def bowtie2_build(xrRNA_fasta,idxDir):
    """bowtie2 build index of reference"""
    
    system_call = 'docker run --rm -v %s:%s:ro -v %s:%s:ro -v %s:%s genomicsdocker/bowtie2:2.2.9 bowtie2_build.py --inpfasta %s --outDir %s' %(
            '/etc/localtime','/etc/localtime',
            os.path.dirname(xrRNA_fasta),'/inpdir',
            idxDir,'/workdir',
            os.path.join('/inpdir',os.path.basename(xrRNA_fasta)),
            '/workdir')

    run_command_line(system_call)


def samtools_sort(outdir,inpsam,outbambase,ncpu):
    """
    sort alignment against genome
    """
    system_call = 'docker run --rm -v %s:%s:ro -v %s:%s genomicsdocker/samtools:0.1.19 samtools_sort.py --inpsam %s --outbambase %s --outdir %s --ncpu %s' %(
            '/etc/localtime','/etc/localtime',
            outdir,'/workdir',
            os.path.join('/workdir',os.path.basename(inpsam)),
            os.path.join('/workdir',os.path.basename(outbambase)),
            '/workdir',
            ncpu)

    run_command_line(system_call)

    
def compute_depth(outdir,inpbam,outfile):
    """compute depth at mapped positions"""

    system_call = 'docker run --rm -v %s:%s:ro -v %s:%s genomicsdocker/samtools:0.1.19 samtools_depth.py --inpbam %s --outfile %s' %(
            '/etc/localtime','/etc/localtime',
            outdir,'/workdir',
            os.path.join('/workdir',os.path.basename(inpbam)),
            os.path.join('/workdir',os.path.basename(outfile)))

    run_command_line(system_call)


def mergebam(outdir,inpbam1,inpbam2,outbam):
    """merge STAR and Bowtie2 alignments"""

    system_call = 'docker run --rm -v %s:%s:ro -v %s:%s genomicsdocker/picard:2.5.0 picard_MergeSamFiles.py --inpbam1 %s --inpbam2 %s --outbam %s' %(
            '/etc/localtime','/etc/localtime',
            outdir,'/workdir',
            os.path.join('/workdir',os.path.basename(inpbam1)),
            os.path.join('/workdir',os.path.basename(inpbam2)),
            os.path.join('/workdir',os.path.basename(outbam)))

    run_command_line(system_call)


def picard_alignMetrics(outdir,reference,inpbam,prefix):
    """generate alignment metrics"""

    #make the output directory
    if not os.path.exists(os.path.join(outdir,'out_Picard_AlignMetrics')):
        os.mkdir(os.path.join(outdir,'out_Picard_AlignMetrics'))

    system_call = 'docker run --rm -v %s:%s:ro -v %s:%s:ro -v %s:%s -v %s:%s:ro genomicsdocker/picard:2.5.0 picard_CollectAlignmentSummaryMetrics.py --inpbam %s --reference %s --output %s' %(
            '/etc/localtime','/etc/localtime',
            outdir,'/workdir',
            os.path.join(outdir,'out_Picard_AlignMetrics'),'/PicardDir',
            os.path.dirname(reference),'/referencedir',
            os.path.join('/workdir/',os.path.basename(inpbam)),
            os.path.join('/referencedir',os.path.basename(reference)),
            os.path.join('/PicardDir/',prefix+'_AlignmentSummmaryMetrics.txt'))

    run_command_line(system_call)

    # generate the results.json
    system_call = 'docker run --rm -v %s:%s:ro -v %s:%s genomicsdocker/picard:2.5.0 write2json_AlignmentMetrics.py --metricfile %s --result_json %s' %(
            '/etc/localtime','/etc/localtime',
            os.path.join(outdir,'out_Picard_AlignMetrics'),'/PicardDir',
            os.path.join('/PicardDir',prefix+'_AlignmentSummmaryMetrics.txt'),
            os.path.join('/PicardDir','temp.json'))

    run_command_line(system_call)

    # rename the json file
    os.rename(os.path.join(outdir,'out_Picard_AlignMetrics','temp.json'),os.path.join(outdir,'out_Picard_AlignMetrics','AlignmentMetrics_results.json'))
    

def picard_rnaMetrics(outdir,annotDir,inpbam,xrRNAcount,prefix,options):
    """generate the RnaSeq metrics"""

    #make the output directory
    if not os.path.exists(os.path.join(outdir,'out_Picard_RnaMetrics')):
        os.mkdir(os.path.join(outdir,'out_Picard_RnaMetrics'))

    #if user doesn't provide annotation files or reference is hg19, use hg19 annotations by default 
    if not options.reference:
        system_call = 'docker run --rm -v %s:%s:ro -v %s:%s:ro -v %s:%s -v %s:%s:ro genomicsdocker/picard:2.5.0 picard_CollectRnaSeqMetrics.py --inpbam %s --refFlat %s --rRNA_interval %s --output %s' %(
                '/etc/localtime','/etc/localtime',
                outdir,'/workdir',
                os.path.join(outdir,'out_Picard_RnaMetrics'),'/PicardDir',
                annotDir,'/annotDir',
                os.path.join('/workdir/',os.path.basename(inpbam)),
                os.path.join('/annotDir','refFlat'),
                os.path.join('/annotDir','rRNA.interval'),
                os.path.join('/PicardDir',prefix+'_RnaSeqMetrics.txt'))

        run_command_line(system_call)

    else:
        system_call = 'docker run --rm -v %s:%s:ro -v %s:%s:ro -v %s:%s -v %s:%s:ro -v %s:%s:ro genomicsdocker/picard:2.5.0 picard_CollectRnaSeqMetrics.py --inpbam %s --refFlat %s --rRNA_interval %s --output %s' %(
                '/etc/localtime','/etc/localtime',
                outdir,'/workdir',
                os.path.join(outdir,'out_Picard_RnaMetrics'),'/PicardDir',
                os.path.dirname(options.refFlat),'/refFlatdir',
                os.path.dirname(options.rRNA_interval),'/rRNAdir',
                os.path.join('/workdir/',os.path.basename(inpbam)),
                os.path.join('/refFlatdir',os.path.basename(options.refFlat)),
                os.path.join('/rRNAdir',os.path.basename(options.rRNA_interval)),
                os.path.join('/PicardDir',prefix+'_RnaSeqMetrics.txt'))

        run_command_line(system_call)

    # generate the results.json
    system_call = 'docker run --rm -v %s:%s:ro -v %s:%s:ro -v %s:%s genomicsdocker/picard:2.5.0 write2json_RnaSeqMetrics.py --metricfile %s --xrRNAbasecount %s --result_json %s' %(
            '/etc/localtime','/etc/localtime',
            options.outdir,'/workdir',
            os.path.join(options.outdir,'out_Picard_RnaMetrics'),'/PicardDir',
            os.path.join('/PicardDir',prefix+'_RnaSeqMetrics.txt'),
            os.path.join('/workdir',os.path.basename(xrRNAcount)),
            os.path.join('/PicardDir','temp.json'))

    run_command_line(system_call)

    # rename the json file
    os.rename(os.path.join(outdir,'out_Picard_RnaMetrics','temp.json'),os.path.join(outdir,'out_Picard_RnaMetrics','RnaSeqMetrics_results.json'))
       

def htseq_genecov(outdir,annotDir,inpbam,prefix,options):
    """compute the gene coverage"""
    
    #make the output directory
    if not os.path.exists(os.path.join(outdir,'out_HTSeq_GeneCov')):
        os.mkdir(os.path.join(outdir,'out_HTSeq_GeneCov'))

    #if user doesn't provide annotation files or reference is hg19, use hg19 annotations by default 
    if not options.reference:
        system_call = 'docker run --rm -v %s:%s:ro -v %s:%s:ro -v %s:%s -v %s:%s:ro genomicsdocker/htseq:0.6.1 htseq.py --gene_gtf %s --inpbam %s --output %s' %(
                '/etc/localtime','/etc/localtime',
                outdir,'/workdir',
                os.path.join(outdir,'out_HTSeq_GeneCov'),'/HtseqDir',
                annotDir,'/annotDir',
                os.path.join('/annotDir','gene.gtf'),
                os.path.join('/workdir/',os.path.basename(inpbam)),
                os.path.join('/HtseqDir/',prefix+'_GenesReadCounts.txt'))

        run_command_line(system_call)

    else:
        system_call = 'docker run --rm -v %s:%s:ro -v %s:%s:ro -v %s:%s -v %s:%s:ro genomicsdocker/htseq:0.6.1 htseq.py --gene_gtf %s --inpbam %s --output %s' %(
                '/etc/localtime','/etc/localtime',
                options.outdir,'/workdir',
                os.path.join(outdir,'out_HTSeq_GeneCov'),'/HtseqDir',
                os.path.dirname(options.gene_gtf),'/geneGTFdir',
                os.path.join('/geneGTFdir',os.path.basename(options.gene_gtf)),
                os.path.join('/workdir/',os.path.basename(inpbam)),
                os.path.join('/HtseqDir/',prefix+'_GenesReadCounts.txt'))

        run_command_line(system_call)

    covfilepath = os.path.join(outdir,'out_HTSeq_GeneCov',prefix+'_GenesReadCounts.txt')
    system_call = 'docker run --rm -v %s:%s:ro -v %s:%s genomicsdocker/htseq:0.6.1 write2json_GeneCovMetrics.py --genecovfile %s --genecovpath %s --result_json %s' %(
            '/etc/localtime','/etc/localtime',
            os.path.join(outdir,'out_HTSeq_GeneCov'),'/HtseqDir',
            os.path.join('/HtseqDir/',prefix+'_GenesReadCounts.txt'),
            covfilepath,
            os.path.join('/HtseqDir/','temp.json'))
    
    run_command_line(system_call)

    # rename the json file
    os.rename(os.path.join(outdir,'out_HTSeq_GeneCov','temp.json'), os.path.join(outdir,'out_HTSeq_GeneCov','GeneCov_results.json'))
    
    
def run_cufflinks(options,prefix,annotDir):
    """compute the transcript abundance"""

    #make the output directory
    if not os.path.exists(os.path.join(options.outdir,'out_cufflinks')):
        os.mkdir(os.path.join(options.outdir,'out_cufflinks'))
   
    #if user doesn't provide annotation files or reference is hg19, use hg19 annotations by default 
    if not options.reference:
        system_call = 'docker run --rm -v %s:%s:ro -v %s:%s:ro -v %s:%s -v %s:%s:ro genomicsdocker/cufflinks:2.2.1 cufflinks.py --inpbam %s --gene_gtf %s --mask_gtf %s --outdir %s --ncpu %s' %(
                '/etc/localtime','/etc/localtime',
                options.outdir,'/workdir',
                os.path.join(options.outdir,'out_cufflinks'),'/cufflinksDir',
                annotDir,'/annotDir',
                '/workdir/'+prefix+'.bam',
                os.path.join('/annotDir','gene.gtf'),
                os.path.join('/annotDir','rRNA_mask.gtf'),
                '/cufflinksDir',
                options.ncpu)
        
        run_command_line(system_call)

    else:
        system_call = 'docker run --rm -v %s:%s:ro -v %s:%s:ro -v %s:%s -v %s:%s:ro -v %s:%s:ro genomicsdocker/cufflinks:2.2.1 cufflinks.py --inpbam %s --gene_gtf %s --mask_gtf %s --outdir %s --ncpu %s' %(
                '/etc/localtime','/etc/localtime',
                options.outdir,'/workdir',
                os.path.join(options.outdir,'out_cufflinks'),'/cufflinksDir',
                os.path.dirname(options.gene_gtf),'/geneannot-dir',
                os.path.dirname(options.mask_gtf),'/maskannot-dir',
                '/workdir/'+prefix+'.bam',
                os.path.join('/geneannot-dir',os.path.basename(options.gene_gtf)),
                os.path.join('/maskannot-dir',os.path.basename(options.mask_gtf)),
                '/cufflinksDir',
                options.ncpu)
    
        run_command_line(system_call)


def main():
    # set up the logger
    logging.getLogger().setLevel(logging.INFO)

    # setup the option parser
    parser = argparse.ArgumentParser()

    # add the arguments for the workflow
    parser.add_argument('-i', '--inpfile', dest= 'inpfile', help = 'The input bam/fastq', metavar = '')
    parser.add_argument('-d', '--outdir', dest = 'outdir', help = 'The output directory', default = '', metavar = '')
    parser.add_argument('-f', '--fracread', dest = 'fracread', help = 'The Fraction of reads to be subsampled', type = float, default = 1.0, metavar = '')
    parser.add_argument('-a', '--trim_adapter', dest = 'trim_adapter', help = 'Whether or not adapter to be trimmed', type = bool, default = True, metavar = '') 
    parser.add_argument('-r', '--reference', dest = 'reference', help = 'The reference fasta sequence; defaults to hg19', metavar = '')
    parser.add_argument('-p', '--platform', dest = 'platform', help = 'The RnaSeq platform: iontorrent, illumina, others', default = 'iontorrent', metavar = '')
    parser.add_argument('-c', '--cleanup', dest = 'cleanup', help = 'options="alignments":clean intermediate alignments only; "all": clean index files as well', default = 'alignments', metavar = '')
    parser.add_argument('-n', '--ncpu', dest = 'ncpu', help = 'The number of CPUs to be used', type = int, default = 1, metavar = '')
    parser.add_argument('-x', '--STARidx', dest = 'STARidx', help = 'The directory of STAR prebuilt idx files', metavar = '')
    parser.add_argument('-y', '--bowtie2idx', dest = 'bowtie2idx', help = 'The path of bowtie2 prebuilt idx files (including basename) for reference genome', metavar = '')
    #parser.add_argument('-z', '--bowtie2idx_rRNA', dest = 'bowtie2idx_rRNA', help = 'The path of bowtie2 prebuilt idx files (including basename) for rRNA reference')
    parser.add_argument('-g', '--gene_gtf', dest = 'gene_gtf', help = 'The gene annotation GTF file', metavar = '')
    parser.add_argument('-m', '--mask_gtf', dest = 'mask_gtf', help = 'The annotation GTF file for regions to be ignored', metavar = '')
    parser.add_argument('-t', '--refFlat', dest = 'refFlat', help = 'The transcript annotation file in refFlat form', metavar = '')
    parser.add_argument('-R', '--rRNA_interval', dest = 'rRNA_interval', help = 'The interval file for rRNA sequence', metavar = '')

    options = parser.parse_args()

   
    #-----------------------------------------------------------------#
    #                      The workflow starts here                   #
    #-----------------------------------------------------------------#
   
    ##------- sanity checks------
    if options.platform != 'iontorrent':
        logging.info('Currently analysis of RnaSeq data from IonTorrent platform is supported only! EXITING....'.upper())
        sys.exit(1)

    if not options.inpfile:
        logging.error('You must provide the BAM/FASTQ input (--inpfile)')
        sys.exit(1)
   
    logging.info((time.asctime( time.localtime(time.time()) ) +': Execution of RnaSeq analysis workflow started....').upper())
   
    
    ##---------platform specific adapter trimming------------
    if options.platform == 'iontorrent':
        adapter = 'ATCACCGACTGCCCATAGAGAGGCTGAGAC' 

    
    ##-------------create the working/output directory if it doesn't exists---------
    if not os.path.exists(options.outdir):
        os.mkdir(options.outdir)
        logging.info('Main output directory created....'.upper())
    else:
        logging.info('Output directory already exists! not overwriting....'.upper())
    

    ##--------reference genome to be used-------------

    #set the default reference to hg19; # change the local path accordingly
    hg19ref= "/results/referenceLibrary/tmap-f3/hg19/hg19.fasta" 
    
    if not options.reference:
        logging.info('hg19 genome will be used as reference by default....'.upper())

        ##copy directory with hg19 annotation files from containers
        annotDir = os.path.join(options.outdir,'hg19-annotations')
        
        if not os.path.exists(annotDir):
            copy_annotations(options)

        reference = os.path.join(annotDir,'hg19.fasta')

        #sequence for 5.8S,16S,26S rRNA that are not in reference genome
        xrRNA_fasta = os.path.join(annotDir,'xrRNA.fasta')

        #set the aligner index dirs to blank because we will use prebuilt indices from containers
        bt2IdxDir = ''
        bt2IdxBase = ''
        starIdx = ''

        logging.info('hg19 reference and annotations files loaded....'.upper())

    else:
        reference = options.reference
        annotDir = ''
        bt2IdxDir = os.path.dirname(options.bowtie2idx)
        bt2IdxBase = os.path.basename(options.bowtie2idx)
        starIdx = options.STARidx

        logging.info('User provided reference, annotations, and index files  will be used....'.upper())


    ##filetype of input file.......
    filetype = os.path.splitext(os.path.basename(options.inpfile))[1]

    ##basename of input bam/fastq file; will be used for subsequent output files.....
    prefix = os.path.splitext(os.path.basename(options.inpfile))[0]


    if filetype == '.bam':
        ##----------------if the input is bam, will require additional processing---------------------
        logging.info(('The input file is %s....' %filetype).upper()) 
        
        if float(options.fracread) != 1.0:
            #------------------Subsample input bam & then convert to fastq----------------------
            
            logging.info('Downsampling has been requested....'.upper())
            prefix = prefix + '.downsampled'
            outbam = os.path.join(options.outdir,prefix+'.bam')
            outfastq = os.path.join(options.outdir,prefix+'.fastq')

            if not os.path.exists(outbam) and not os.path.exists(outfastq):
                bam_downsample_2fastq(options.outdir,options.inpfile,outbam,outfastq,options.fracread)
            else:
                logging.info('Downsampled bam & fastq already exist, skipping downsampling & bam2fastq....'.upper())
        
        else:
            #------------------convert input bam to fastq--------------------
            outfastq = os.path.join(options.outdir,prefix+'.fastq')
            if not os.path.exists(outfastq):
                bam2fastq(options.outdir,options.inpfile,outfastq)
            else:
                logging.info('Input bam has been already converted to fastq, skipping bam2fastq....'.upper())
            
            inpfastq = outfastq

    else:
        ##---------------if the input is fastq, the workflow starts here---------------
        logging.info(('The input file is %s, downsampling is currently not supported....' %filetype).upper())
        inpfastq = options.inpfile
    
    if options.trim_adapter == True:
        ##---------------Trim the adapters from both ends of reads & remove too short reads---------------
        oldprefix = prefix
        prefix = oldprefix + '.cutadapt'
        outfastq = os.path.join(options.outdir,prefix+'.fastq')
        
        if not os.path.exists(outfastq):
            run_cutadapt(options.outdir,inpfastq,outfastq,adapter)
        else:
            logging.info('Adapter trimmed fastq already exists, skipping cutadapt run....'.upper())


    ##--------------Run STAR: Splice junction mapping-----------------
    starInpFastq = outfastq 
    starOutSam = os.path.join(options.outdir,'Aligned.out.sam')
    starOutSortBamBase = os.path.join(options.outdir,'mappedSTAR')

    if not os.path.exists(starOutSortBamBase+'.bam'):
        star_align(options.outdir,starInpFastq,starIdx,options)

        #sort aligned sam
        samtools_sort(options.outdir,starOutSam,starOutSortBamBase,options.ncpu)


    ##------------Run bowtie2 alignment of STAR-unmapped reads:1st against ref. genome & then against rRNA ref.--------
    bt2InpFastq = os.path.join(options.outdir,'Unmapped.out.mate1')
    bt2OutSam = os.path.join(options.outdir,'unmapSTAR_alignedBowtie2.bam')
    bt2OutSortBamBase = os.path.join(options.outdir,'unmapSTAR_mappedBowtie2')
    bt2OutUnmap = os.path.join(options.outdir,'STARBowtie2_unmap.fq')

    if not os.path.exists(bt2OutSortBamBase+'.bam'):
        bowtie2_align(options.outdir,bt2IdxDir,bt2IdxBase,bt2InpFastq,bt2OutSam,bt2OutUnmap,options.ncpu,options)

        #sort aligned sam
        samtools_sort(options.outdir,bt2OutSam,bt2OutSortBamBase,options.ncpu)


    ###--------------Merge the STAR-mapped and bowtie2-remapped bam files----------
    prefix = prefix + '.STARBowtie2'
    finalbam = os.path.join(options.outdir,prefix+'.bam')

    if not os.path.exists(finalbam):
        mergebam(options.outdir,starOutSortBamBase+'.bam',bt2OutSortBamBase+'.bam',finalbam)
    else:
        logging.info('STAR+Bowtie2 merged bam already exists, skipping merging....'.upper())


    ##---------------Perform addtional alignment against 5.8S,18S,28S rRNA sequence--------------    
    if not options.reference:
        #build the index of xrRNA reference sequence
        xrRNA_idxDir = os.path.join(options.outdir,'xrRNA_idx')
        if not os.path.exists(xrRNA_idxDir):
            os.mkdir(xrRNA_idxDir)

            bowtie2_build(xrRNA_fasta,xrRNA_idxDir)

        else:
            logging.info('xrRNA index files already made! If want to rerun, delete the directory....'.upper())
        
        bt2InpFastq2 = os.path.join(options.outdir,'STARBowtie2_unmap.fq')
        bt2OutSam2 = os.path.join(options.outdir,'xrRNA_aligned.sam')
        bt2OutSortBamBase2 = os.path.join(options.outdir,'xrRNA_mapped')
        bt2OutUnmap2 = os.path.join(options.outdir,'temp_unmap.fq')


        if not os.path.exists(bt2OutSortBamBase2+'.bam'): 
            bowtie2_align(options.outdir,xrRNA_idxDir,'bowtie2',bt2InpFastq2,bt2OutSam2,bt2OutUnmap2,options.ncpu,options)
            os.remove(bt2OutUnmap2)

            #sort aligned sam
            samtools_sort(options.outdir,bt2OutSam2,bt2OutSortBamBase2,options.ncpu)
    
    
        ##------------Calculate the depth of additional aligned rRNA bases-------------
        if not os.path.exists(os.path.join(options.outdir,'xrRNA_basecounts')):
            compute_depth(options.outdir,bt2OutSortBamBase2+'.bam',os.path.join(options.outdir,'xrRNA_basecounts'))

    
    xrRNAcount = os.path.join(options.outdir,'xrRNA_basecounts')

    
    ###------------Get Alignment Metrics from Picard tool-----------
    if len(glob.glob(os.path.join(options.outdir,'out_Picard_AlignMetrics','*.json'))) == 0: 
        picard_alignMetrics(options.outdir,reference,finalbam,prefix)
    else:
        logging.info('results.json EXISTS WITHIN PICARD ALIGNMETRICS OUT DIRECTORY, SKIPPING CollectAlignmentSummaryMetrics RUN....')
    
    
    ###-------------Get RNASeq metrics from Picard tool-------------
    if len(glob.glob(os.path.join(options.outdir,'out_Picard_RnaMetrics','*.json'))) == 0:
        picard_rnaMetrics(options.outdir,annotDir,finalbam,xrRNAcount,prefix,options)
    else:
        logging.info('results.json EXISTS WITHIN PICARD RNAMETRICS OUT DIRECTORY, SKIPPING CollectRnaSeqMetrics RUN....')


    
    ###-------------Get gene coverage from HTSeq-------------
    
    #first filter the mapped reads: ideally should have used stdout from samtools run to htseq as stdin
    system_call = 'docker run --rm -v %s:%s:ro -v %s:%s genomicsdocker/samtools:0.1.19 samtools_mapreads.py --inpbam %s --outbam %s' %(
            '/etc/localtime','/etc/localtime',
            options.outdir,'/workdir',
            os.path.join('/workdir',os.path.basename(finalbam)),
            os.path.join('/workdir','mapped_final.bam'))

    run_command_line(system_call)
    
    mappedFinalBam = os.path.join(options.outdir,'mapped_final.bam')

    if len(glob.glob(os.path.join(options.outdir,'out_HTSeq_GeneCov','*.json'))) == 0:
        htseq_genecov(options.outdir,annotDir,mappedFinalBam,prefix,options)
    else:
        logging.info('results.json EXISTS WITHIN HTSEQ OUT DIRECTORY, SKIPPING HTSeq RUN....')
    
    os.remove(mappedFinalBam)


    ###-------------Run Cufflinks to get transcripts abundance------------
    #if len(glob.glob(os.path.join(options.outdir,'out_cufflinks','*.json'))) == 0:
    #    run_cufflinks(options,prefix,annotDir)   
    #else:
    #    logging.info('results.json EXISTS WITHIN CUFFLINKS OUT DIRECTORY, SKIPPING Cufflinks RUN....')

    #
    ##--------------Clean up intermediate alignment files and index files if copied--------------
    #if options.cleanup == 'alignments' or options.cleanup == 'all':
    #    if os.path.exists(os.path.join(options.outdir,'Aligned.out.sam')):
    #        os.remove(os.path.join(options.outdir,'Aligned.out.sam'))
    #    if os.path.exists(os.path.join(options.outdir,'unmappedSTAR_remapBowtie2_unsorted.sam')):
    #        os.remove(os.path.join(options.outdir,'unmappedSTAR_remapBowtie2_unsorted.sam'))
    #    if os.path.exists(os.path.join(options.outdir,'rRNA_mapped_unsorted.sam')):
    #        os.remove(os.path.join(options.outdir,'rRNA_mapped_unsorted.sam'))
    #    if os.path.exists(os.path.join(options.outdir,prefix+'_mapreads.bam')):
    #        os.remove(os.path.join(options.outdir,prefix+'_mapreads.bam'))    

    #if options.cleanup == 'all':
    #    if os.path.exists(annotDir):
    #        os.remove(annotDir)
    #    if os.path.exists(starIdxDir):
    #        os.remove(starIdxDir)
    #    if os.path.exists(bt2_genomeIdxDir):
    #        os.remove(bt2_genomeIdxDir)
    #    if os.path.exists(bt2_rRNAIdxDir):
    #        os.remove(bt2_rRNAIdxDir)

    
    logging.info((time.asctime( time.localtime(time.time()) ) +': Execution of RnaSeq analysis workflow finished').upper())




if __name__ == "__main__":
    main()


