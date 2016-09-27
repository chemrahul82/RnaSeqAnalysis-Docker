#!/usr/bin/env python

"""
wrapper script for samtools docker container inside RNASEQ analysis workflow
The script calculates the total depth at mapped positions and prints the number to a file
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
    #usage:python samtools_depth.py -i input.bam -o depth.txt
    parser.add_argument('-i', '--inpbam', dest = 'inpbam', help = 'The input bam')
    parser.add_argument('-o', '--outfile', dest = 'outfile', help = 'The output file to write the total depth')
    
    options = parser.parse_args()
    
    if not options.inpbam or not options.outfile:
        logging.error(' You must provide input bam (--inpbam) and output file (--outfile)')
        sys.exit(1)
    

    #run samtools
    system_call = "samtools view -b -F4 %s | samtools depth /dev/stdin | awk '{sum+=$3} END {print sum}' > %s" %(
            options.inpbam,
            options.outfile)

    run_command_line(system_call)

            
