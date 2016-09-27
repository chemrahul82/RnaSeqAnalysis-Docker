#!/usr/bin/env python

"""
wrapper script for bowtie2 docker container inside RNASEQ analysis workflow
Builds index of a reference fasta sequence 
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
    #usage:python bowtie2_build.py -i inp.fasta -d outdir
    parser.add_argument('-i', '--inpfasta', dest = 'inpfasta', help = 'The input fasta file')
    parser.add_argument('-d', '--outDir', dest = 'outDir', help = 'The directory where index files will be generated')

    options = parser.parse_args()
    
    if not options.outDir or not options.inpfasta:
        logging.error('You must provide the output directory (--outDir) & input fasta (--inpfasta)')
        sys.exit(1)

    #run bowtie2-build
    #usage: bowtie-build $reference_fasta $outDir/bowtie2
    
    system_call = "bowtie2-build --quiet %s %s/bowtie2" %(
            options.inpfasta,
            options.outDir)

    run_command_line(system_call)

    
