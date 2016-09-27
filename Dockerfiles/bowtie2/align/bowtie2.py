#!/usr/bin/env python

"""
wrapper script for bowtie2 docker container inside RNASEQ analysis workflow
Performs alignment against reference sequence; requires prebuilt index files
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
    #usage:python bowtie2.py -x referenceIdxPath -i input.fastq -o output.sam -u unmapped.reads -n no_cpu
    parser.add_argument('-x', '--index', dest = 'index', help = 'The path of reference idx files including basename')
    parser.add_argument('-i', '--inpfastq', dest = 'inpfastq', help = 'The input fastq')
    parser.add_argument('-o', '--outsam', dest = 'outsam', help = 'The mapped output sam')
    parser.add_argument('-u', '--outun', dest = 'outun', help = 'The unmapped output fastq')
    parser.add_argument('-n', '--ncpu', dest = 'ncpu', help = 'The number of CPUs to be used', default = 1)

    options = parser.parse_args()
    
    if not options.inpfastq or not options.outsam or not options.outun or not options.index:    
        logging.error('You must provide input fasta (--inpfastq), mapped output sam (--outsam), unmmaped fastq (--outun), prebuilt index files (--index)')
        sys.exit(1)

    #run bowtie2 alignment
    #usage:bowtie2 --local --very-sensitive-local -p ncpu -q --mm -x directory_containing_idx_file -U Unmapped.out.mate1_from_STAR_run \
    # --un sbt2_unmap.fq

    # change directory to output/working directory; run bowtie2 there; output will be in the same directory
    #os.chdir(options.outDir)

    system_call = "bowtie2 --quiet --local --very-sensitive-local -p %s -q --mm -x %s -U %s --un %s -S %s" %(
            options.ncpu,
            options.index,
            options.inpfastq,
            options.outun,
            options.outsam)

    run_command_line(system_call)

   
    #delete intermediate files
    #os.remove('Unmapped.out.mate1')
    #os.remove('STARBowtie2_unmap.fq')

