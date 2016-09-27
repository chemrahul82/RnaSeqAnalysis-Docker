#!/usr/bin/env python

"""
wrapper script for bowtie2 docker container inside RNASEQ analysis workflow
Performs alignment against prebuilt index files of mm10 genome within the container
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
    #usage:python bowtie2.py -i inp.fastq -o output.sam -u unmapped.fastq -n no_cpu 
    parser.add_argument('-n', '--ncpu', dest = 'ncpu', help = 'The number of CPUs to be used', default = 1)
    parser.add_argument('-i', '--inpfastq', dest = 'inpfastq', nargs = '+', help = 'The input fastq')
    parser.add_argument('-o', '--outsam', dest = 'outsam', help = 'The mapped output sam')
    parser.add_argument('-u', '--outun', dest = 'outun', help = 'The unmapped output fastq')
    parser.add_argument('-S', '--readtype', dest = 'readtype', help = 'The sequencing type: single end or paired end', default = 'SE')

    options = parser.parse_args()
    
    if not options.inpfastq or not options.outsam or not options.outun:
        logging.error('You must provide input fastq (--inpfastq), mapped output sam (--outsam), unmmaped fastq (--outun)')
        sys.exit(1)

    if options.readtype == 'PE' and len(options.inpfastq) < 2:
        logging.error('If you select Paired-end readtype, you must provide the two fastq files')
        sys.exit(1)

    #run bowtie2 alignment
    #usage:bowtie2 --local --very-sensitive-local -p ncpu -q --mm -x directory_containing_idx_file -U Unmapped.out.mate1_from_STAR_run \
    # --un sbt2_unmap.fq

    #directory with prebuilt genome index files inside the container
    idxDir = os.path.abspath("idx_mm10")
   
    #single-end
    if options.readtype == 'SE':
        system_call = "bowtie2 --quiet --local --very-sensitive-local -p %s -q --mm -x %s -U %s --un %s -S %s" %(
                options.ncpu,
                os.path.join(idxDir,'bowtie2'),
                options.inpfastq[0],
                options.outun,
                options.outsam)

        run_command_line(system_call)
       
    #paired-end
    if options.readtype == 'PE':
        system_call = "bowtie2 --quiet --local --very-sensitive-local -p %s -q --mm -x %s -1 %s -2 %s --un %s -S %s" %(
                options.ncpu,
                os.path.join(idxDir,'bowtie2'),
                options.inpfastq[0],
                options.inpfastq[1],
                options.outun,
                options.outsam)

        run_command_line(system_call)
     

