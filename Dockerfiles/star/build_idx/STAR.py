#!/usr/bin/env python

"""
wrapper script for STAR docker container & prebuilt hg19 index files 
inside RNASEQ analysis workflow
"""

__author__='Rahul K. Das'

import argparse
import json
import logging
import os
import subprocess
import sys
import time

def run_command_line(system_call):
    logging.getLogger().setLevel(logging.INFO)
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



if __name__ == "__main__":

    # setting up the option parser
    parser = argparse.ArgumentParser()

    #arguments
    #usage:python STAR.py -i input.bam -d outdir
    parser.add_argument('-i', '--inpfastq', dest = 'inpfastq', help = 'The input fastq')
    parser.add_argument('-d', '--outDir', dest = 'outDir', help = 'The directory where STAR is run and outputs all files')
    parser.add_argument('-n', '--ncpu', dest = 'ncpu', help = 'The number of CPUs to be used', default = 1)

    options = parser.parse_args()
    
    if not options.inpfastq or not options.outDir:
        logging.error(' You must provide input fastq (--inpfastq), and output directory (--outDir)')
        sys.exit(1)
    

    #run STAR
    #usage:STAR --genomeDir directory_containing_genome_idx_file --readFilesIn input.fastq 
    #other options are similar to those used in IonTorrent RNASeqAnalysis plugin

    #directory with prebuilt genome index files inside the container
    idxDir = os.path.abspath("star_idx_hg19")

    # change directory to output/working directory; STAR outputs everything in the working directory
    os.chdir(options.outDir)

    system_call = "STAR --genomeDir %s --readFilesIn %s --outReadsUnmapped Fastx --chimSegmentMin 18 --chimScoreMin 12 --runThreadN %s" %(
            idxDir,
            options.inpfastq,
            options.ncpu)

    run_command_line(system_call)

            
