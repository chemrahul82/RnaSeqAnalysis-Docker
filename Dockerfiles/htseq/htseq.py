#!/usr/bin/env python

"""
wrapper script for htseq docker container within the RNASeqAnalysis workflow
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
    #usage:python htseq.py --inpbam input.bam --gene_gtf gene_GTF_file --output gene.counts
    parser.add_argument('-i', '--inpbam', dest = 'inpbam', help = 'The input bam')
    parser.add_argument('-g', '--gene_gtf', dest = 'gene_gtf', help = 'The gene GTF file')
    parser.add_argument('-o', '--output', dest = 'output', help = 'The output file with readcounts per gene')

    options = parser.parse_args()
    
    if  not options.inpbam or not options.gene_gtf or not options.output:
        logging.error(' You must provide input bam (--inpbam), gene GTF file (--gene_gtf) and output file (--output)')
        sys.exit(1)
    

    #run HTSeq:  
    #Usage: htseq-count -q -f bam -t exon -i gene_name input.bam gene_gtf_file > gene.count
    system_call = "htseq-count -q -f bam -t exon -i gene_name %s %s > %s" %(
            options.inpbam,
            options.gene_gtf,
            options.output)
    
    
    
    run_command_line(system_call)

            
