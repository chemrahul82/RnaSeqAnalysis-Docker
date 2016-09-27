#!/usr/bin/env python

"""
wrapper script for samtools docker container inside RNASEQ analysis workflow
This script takes the arguments for sorting a sam, generating a sorted bam, and performing that operation
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
    #usage:python samtools_sort.py -i input.sam -o outbam_basename -d outdir
    parser.add_argument('-i', '--inpsam', dest = 'inpsam', help = 'The input unsorted sam')
    parser.add_argument('-o', '--outbambase', dest = 'outbambase', help = 'The  basename of output sorted bam')
    parser.add_argument('-d', '--outdir', dest = 'outdir', help = 'The output directory',default = '')
    parser.add_argument('-n', '--ncpu', dest = 'ncpu', help = 'The number of CPUs to be used', default = 1)
    parser.add_argument('-t', '--sortby', dest = 'sortby', help = 'Whether to sort by leftmost coordinates or name', default = 'coordinate')

    options = parser.parse_args()
    
    if not options.outbambase or not options.inpsam or not options.outdir:
        logging.error(' You must provide input sam (--inpsam) and basename of output bam (--outbambase)')
        sys.exit(1)
    

    #run samtools
    #usage:samtools -bS input.bam | samtools sort - output 
    
    if options.sortby == 'coordinate':
        system_call = "samtools view -bS %s | samtools sort -@ %s - %s" %(
                options.inpsam,
                options.ncpu,
                os.path.join(options.outdir,options.outbambase))

        run_command_line(system_call)

    if options.sortby == 'name':
        system_call = "samtools sort -n -@ %s %s %s" %(
                options.ncpu,
                options.inpsam,
                os.path.join(options.outdir,options.outbambase))

        run_command_line(system_call)



            
