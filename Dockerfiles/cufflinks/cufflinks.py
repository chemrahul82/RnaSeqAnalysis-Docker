#!/usr/bin/env python

"""
wrapper script for cufflinks docker container within the RNASeqAnalysis workflow
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
    #usage:python cufflinks.py --inpbam input.bam --gene_gtf gene_GTF_file --mask_gtf mask_GTF_file --outdir cufflinks_output
    parser.add_argument('-g', '--gene_gtf', dest = 'gene_gtf', help = 'The gene annotation GTF file')
    parser.add_argument('-M', '--mask_gtf', dest = 'mask_gtf', help = 'The annotation file for sequence to be ignored')
    parser.add_argument('-i', '--inpbam', dest = 'inpbam', help = 'The input bam')
    parser.add_argument('-o', '--outdir', dest = 'outdir', help = 'The output directory for cufflinks', default = '')
    parser.add_argument('-n', '--ncpu', dest = 'ncpu', help = 'The number of CPUs to be used', type = int, default = 1)
    parser.add_argument('-l', '--libtype', dest = 'libtype', help = 'The library type option for cufflinks',default = 'fr-secondstrand')
    parser.add_argument('-m', '--frag_len_mean', dest = 'frag_len_mean', help = 'The expected (mean) fragment length', type = float, default = 100)
    parser.add_argument('-R', '--readtype', dest = 'readtype', help = 'The sequencing type: single end or paired end', default = 'SE') 
    
    options = parser.parse_args()
    
    if not options.gene_gtf or not options.mask_gtf or not options.inpbam or not options.outdir:
        logging.error(' You must provide gene annotation (--gene_gtf), annotation for regions to be ignored (--mask_gtf), and input bam (--inpbam)')
        sys.exit(1)
    

    #run cufflinks: The parameters are similar to those in IonTorrent RNASeqAnalysis plugin
    if options.readtype == 'SE':
        system_call = "cufflinks -q -p %s -m %s -s 60 -G %s -M %s --library-type %s --max-bundle-length 3500000 -o %s --no-update-check %s" %(
                options.ncpu,
                options.frag_len_mean,
                options.gene_gtf,
                options.mask_gtf,
                options.libtype,
                options.outdir,
                options.inpbam)
    
    #for PE reads cufflinks can estimate the fragment length 
    elif options.readtype == 'PE':
        system_call = "cufflinks -q -p %s -G %s -M %s --library-type %s --max-bundle-length 3500000 -o %s --no-update-check %s" %(
                options.ncpu,
                options.gene_gtf,
                options.mask_gtf,
                options.libtype,
                options.outdir,
                options.inpbam)
    

    
    run_command_line(system_call)

            
