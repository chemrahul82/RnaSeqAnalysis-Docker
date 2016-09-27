#!/usr/bin/env python

"""
wrapper script for samtools docker container inside RNASEQ analysis workflow
This script takes the arguments for subsampling reads from bamfile and performing that operation
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
    #usage:python samtools_subsample.py --inpbam input.bam --fracread 0.2 --outbam output.bam
    parser.add_argument('-i', '--inpbam', dest = 'inpbam', help = 'The input bam')
    parser.add_argument('-o', '--outbam', dest = 'outbam', help = 'The output bam')
    parser.add_argument('-f', '--fracread', dest = 'fracread', help = 'Fraction of reads to be subsampled', default = '1.0')

    options = parser.parse_args()
    
    if not options.outbam or not options.inpbam:
        logging.error(' You must provide input bam (--inpbam), output bam (--outbam)')
        sys.exit(1)
    if options.outbam and options.inpbam and not options.fracread:
        logging.info('fracread not provided: switching to default of 1; not subsampling')
    

    #run samtools
    #usage:samtools -bhs fracread input.bam > output.bam 
    
    system_call = "samtools view -bhs %s %s > %s" %(
            options.fracread,
            options.inpbam,
            options.outbam)

    run_command_line(system_call)

            
