#!/usr/bin/env python

"""
wrapper script for picard docker container: CollectAlignmentSummaryMetrics tool
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
    #usage:python picard_CollectAlignmentSummaryMetrics.py --inpbam input.bam --reference reference.fasta --output Summary.txt
    parser.add_argument('-i', '--inpbam', dest = 'inpbam', help = 'The input bam')
    parser.add_argument('-r', '--reference', dest = 'reference', help = 'The reference file')
    parser.add_argument('-o', '--output', dest = 'output', help = 'The output alignment summary')

    options = parser.parse_args()
    
    if not options.reference or not options.inpbam or not options.output:
        logging.error(' You must provide input bam (--inpbam), output file (--output), and reference sequence file (--reference)')
        sys.exit(1)
    

    #run Picard CollectAlignmentSummaryMetrics tool
    #usage: java -jar picard.jar CollectAlignmentSummaryMetrics I=input.bam R=reference.fa O=output.txt

    system_call = 'java -jar picard.jar CollectAlignmentSummaryMetrics R=%s I=%s O=%s' %(
            options.reference,
            options.inpbam,
            options.output)

    run_command_line(system_call)




            
