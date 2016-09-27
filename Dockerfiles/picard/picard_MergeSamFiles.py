#!/usr/bin/env python

"""
wrapper script for picard docker container:MergeSamFiles tool
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
    #usage:python picard_MergeSamFiles.py --inpbam1 input1.bam ---inpbam2 input2.bam --outbam output.bam
    parser.add_argument('-i', '--inpbam1', dest = 'inpbam1', help = 'The input bam-1')
    parser.add_argument('-j', '--inpbam2', dest = 'inpbam2', help = 'The input bam-2')
    parser.add_argument('-o', '--outbam', dest = 'outbam', help = 'The output merged bam')

    options = parser.parse_args()
    
    #if not options.outfastq or not options.inpbam or not options.outdir:
    #    logging.error(' You must provide input fastq (-i), output fastq (-o), and output directory (-d)')
    #    sys.exit(1)
    

    #run picard
    #Usage:java -jar picard.jar MergeSamFiles I=inp1.bam I=inp2.bam O=merged.bam
    system_call = 'java -jar picard.jar MergeSamFiles ASSUME_SORTED=true MERGE_SEQUENCE_DICTIONARIES=true USE_THREADING=true I=%s I=%s O=%s' %(
            options.inpbam1,
            options.inpbam2,
            options.outbam)

    run_command_line(system_call)

            
