#!/usr/bin/env python

"""
wrapper script for samtools docker container inside RNASEQ analysis workflow
This script filters the mapped reads from a bam file
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
    #usage:python samtools_mapread.py --inpbam input.bam --outbam output_mapreads.bam
    parser.add_argument('-i', '--inpbam', dest = 'inpbam', help = 'The input bam')
    parser.add_argument('-o', '--outbam', dest = 'outbam', help = 'The output bam with only mapped reads')
    options = parser.parse_args()
    
    if not options.inpbam or not options.outbam:
        logging.error(' You must provide input bam (--inpbam) and output sam (--outbam)')
        sys.exit(1)
    

    #run samtools
    #usage:samtools view -b -F4 input.bam > output.bam 
    
    system_call = "samtools view -b -F4 %s > %s" %(
            options.inpbam,
            options.outbam)

    run_command_line(system_call)

            
