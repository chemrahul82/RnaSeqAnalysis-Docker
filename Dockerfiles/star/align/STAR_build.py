#!/usr/bin/env python

"""
wrapper script for STAR docker container that builds index files 

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
    parser.add_argument('-r', '--ref_fasta', dest = 'ref_fasta', help = 'The input reference fasta')
    parser.add_argument('-d', '--outDir', dest = 'outDir', help = 'The directory where STAR will generate the index files')
    parser.add_argument('-g', '--gene_gtf', dest = 'gene_gtf', help = 'The gene annotation file')
    parser.add_argument('-n', '--ncpu', dest = 'ncpu', help = 'The number of CPUs to be used', default = 1)

    options = parser.parse_args()
    
    if not options.ref_fasta or not options.outDir:
        logging.error(' You must provide input reference fasta (--ref_fasta), and output directory (--outDir)')
        sys.exit(1)
    

    #run STAR

    # change directory to output/working directory; STAR outputs everything in the working directory
    os.chdir(options.outDir)

    system_call = "STAR --runThreadN %s --runMode genomeGenerate --genomeDir %s --genomeFastaFiles %s --sjdbGTFfile %s --sjdbOverhang 75" %(
            options.ncpu,
            options.outDir,
            options.ref_fasta,
            options.gene_gtf)

    run_command_line(system_call)

            
