#!/usr/bin/env python

"""
wrapper script for cutadapt docker container within the RNASeqAnalysis workflow
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
    #usage:python cutadapt.py -b ADAPTERSEQUENCE -i inout.fastq -o output.fastq 
    parser.add_argument('-b', '--adapter', dest = 'adapter', help = 'The adapter sequence to be trimmed', default='')
    parser.add_argument('-o', '--outfastq', dest = 'outfastq', help = 'The output fastq after adapter trimming')
    parser.add_argument('-i', '--inpfastq', dest = 'inpfastq', help = 'The input fastq')
    #parser.add_argument('-l', '--outlog', dest = 'outlog', help = 'log file')

    options = parser.parse_args()
    
    if not options.adapter or not options.outfastq or not options.inpfastq:
        logging.error(' You must provide adapter sequence (-b), input fastq (-i), output fastq (-o)')
        sys.exit(1)
    

    #run cutadapt: filtering out very short (-m 16) reads akin to IonTorrent RnaSeqAnalysis plugin 
    #usage:cutadapt -m 1 -b ADAPTERSEQUENCE inout.fastq > output.fastq 2 > report.txt
    #system_call = "cutadapt -m 16 -b %s %s > %s 2> %s" %(
    system_call = "cutadapt -m 16 -b %s %s > %s" %(
            options.adapter,
            options.inpfastq,
            options.outfastq)

    run_command_line(system_call)

            
