#!/usr/bin/env python

"""
wrapper script for picard docker container: CollectRnaSeqMetrics tool
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

    #arguments for the wrapper script
    #usage:python picard_CollectRnaSeqMetrics.py --inpbam input.bam --refFlat gene_transcript_annotation_file --rRNA_interval rRNA.interval --output Summary.txt
    parser.add_argument('-i', '--inpbam', dest = 'inpbam', help = 'The input bam')
    parser.add_argument('-r', '--refFlat', dest = 'refFlat', help = 'The transcript annotation file')
    parser.add_argument('-R', '--rRNA_interval', dest = 'rRNA_interval', help = 'interval file with rRNA sequence location in genome')
    parser.add_argument('-o', '--output', dest = 'output', help = 'The output alignment summary')

    options = parser.parse_args()
    
    if not options.inpbam or not options.output or not options.refFlat or not options.rRNA_interval:
        logging.error(' You must provide input bam (--inpbam), output file (--output), refFlat file (--refFlat), and rRNA interval file (--rRNA_interval)')
        sys.exit(1)
    

    #run Picard CollectRnaSeqMetrics tool
    #usage: java -jar picard.jar CollectRnaSeqMetrics I=input.bam O=output.txt REF_FLAT=Gene_annotations_refFlat_form STRAND=FIRST_READ_TRANSCRIPTION_STRAND \
    #RIBOSOMAL_INTERVALS=rRNA_interval_file

    system_call = 'java -jar picard.jar CollectRnaSeqMetrics I=%s O=%s REF_FLAT=%s STRAND=%s RIBOSOMAL_INTERVALS=%s MINIMUM_LENGTH=100' %(
            options.inpbam,
            options.output,
            options.refFlat,
            'FIRST_READ_TRANSCRIPTION_STRAND',
            options.rRNA_interval)

    run_command_line(system_call)

            
