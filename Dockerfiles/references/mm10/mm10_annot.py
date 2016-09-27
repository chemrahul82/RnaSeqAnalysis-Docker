#!/usr/bin/env python

"""
wrapper script for docker container that copies mm10 annotation files to the host directory inside RNASEQ analysis workflow
"""

__author__='Rahul K. Das'

import argparse
import json
import logging
import os
import subprocess
import sys
import time
import shutil

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
    parser.add_argument('-d', '--outdir', dest = 'outdir', help = 'output directory')

    options = parser.parse_args()
    
    if not options.outdir:
        logging.error(' You must provide output directory (--outdir)')
        sys.exit(1)
    
    annotDir = os.path.join(options.outdir,'mm10-annotations')
    if not os.path.exists(annotDir):
        os.mkdir(annotDir)
    
    shutil.copy('refFlat',annotDir)
    shutil.copy('gene.gtf',annotDir)
    shutil.copy('rRNA_mask.gtf',annotDir)
    shutil.copy('rRNA.interval',annotDir)
    shutil.copy('xrRNA.fasta',annotDir)
    shutil.copy('mm10.fasta',annotDir)
