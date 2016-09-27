#!/usr/bin/env python

"""
wrapper script for bam2fastq docker container
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
    #usage:python bam2fastq.py --inpbam input.bam --outfastq output.fastq
    parser.add_argument('-o', '--outfastq', dest = 'outfastq', help = 'The output fastq')
    parser.add_argument('-i', '--inpbam', dest = 'inpbam', help = 'The input bam')

    options = parser.parse_args()
    
    if not options.outfastq or not options.inpbam:
        logging.error(' You must provide input bam (--inpbam) and output fastq (--outfastq)')
        sys.exit(1)
    

    #run bam2fastq
    #usage:bam2fastq -o out.fastq in.bam 
    system_call = "bam2fastq -o %s %s" %(
            options.outfastq,
            options.inpbam)

    run_command_line(system_call)

            
