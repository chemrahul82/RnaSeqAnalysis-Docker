#!/usr/bin/env python

"""
writes output metrics from HtSeq-count gene coverage Metrics to json file
"""
from __future__ import division

__author__='Rahul K. Das'

import argparse
import json
import logging
import os
import subprocess
import sys
import time
import numpy as np
from collections import OrderedDict



def write2json(geneCovOut, results_json, covFilePath):
    """
    Description: Takes the output file of HTSeq-count
    and generate a json file with metrics

    Takes: the HTSeq output file

    Returns: the json file
    """
    #initialize the primary ordered dictionary
    results = OrderedDict()

    
    #initialize the four sections as list
    results['summary_metrics'] = []
    results['additional_metrics'] = []
    results['files'] = []
    results['images'] = []

   
    # read the HtSeq out file and calculate relevant stats
    genecovf = open(geneCovOut, "r")
    counts = [int(line.split('\t')[1].strip('\n')) for line in genecovf]
    
    counts_sorted = sorted(counts[:-5], reverse=True)

    #total reads mapped to genes: last 5 lines are HTSeq summary
    map_reads = sum(counts[:-5]) 

    #median coverage across genes
    bins = [50,100,500,1000]
    median_cov  = [0]*len(bins)
    
    if int(map_reads) > 0:
        for idx,ib in enumerate(bins):
            median_cov[idx] = np.median(counts_sorted[:int(ib)])


    # --------first the "three" main metrics in summary_metrics--------
    
    # total mapped reads
    d = OrderedDict()
    d['label'] = 'READS_MAPPED2GENES'
    d['value'] = map_reads
    d['type'] = 'reads'
    results['summary_metrics'].append(d)
    
    ##reads which could not be assigned to any feature 
    d = OrderedDict()
    d['label'] = 'READS_NO_FEATURES'
    d['value'] = counts[-5]
    d['type'] = 'reads'
    results['summary_metrics'].append(d)

    # median coverage of top 100 genes
    d = OrderedDict()
    d['label'] = 'MEDIAN_COV_TOP100'
    d['value'] = median_cov[1]
    results['summary_metrics'].append(d)
    

    ##----------now the additional_metrics-------
    
    # median coverage of top 50 genes
    d = OrderedDict()
    d['label'] = 'MEDIAN_COV_TOP50'
    d['value'] = median_cov[0]
    results['additional_metrics'].append(d)
    
    # median coverage of top 500 genes
    d = OrderedDict()
    d['label'] = 'MEDIAN_COV_TOP500'
    d['value'] = median_cov[2]
    results['additional_metrics'].append(d)

    # median coverage of top 1000 genes
    d = OrderedDict()
    d['label'] = 'MEDIAN_COV_TOP1000'
    d['value'] = median_cov[3]
    results['additional_metrics'].append(d)

    ## reads which could have been assigned to more than one feature and hence were not counted for any of these
    d = OrderedDict()
    d['label'] = 'READS_AMBIGUOUS'
    d['value'] = counts[-4]
    d['type'] = 'reads'
    results['additional_metrics'].append(d)

    ## reads which were skipped due to the -a option (default=10)
    d = OrderedDict()
    d['label'] = 'READS_LT_Q10'
    d['value'] = counts[-3]
    d['type'] = 'reads'
    results['additional_metrics'].append(d)
        
    ## reads without alignments in BAM
    d = OrderedDict()
    d['label'] = 'READS_NOT_ALIGNED'
    d['value'] = counts[-2]
    d['type'] = 'reads'
    results['additional_metrics'].append(d)
     
    ## reads with more then one alignment in BAM
    d = OrderedDict()
    d['label'] = 'READS_ALIGNMENT_NOT_UNIQUE'
    d['value'] = counts[-1]
    d['type'] = 'reads'
    results['additional_metrics'].append(d)
    

    ##files
    d = OrderedDict()
    d['description'] = 'gene-count-table'
    d['file_path'] = covFilePath
    d['name'] = 'GenesReadCounts.txt'
    results['files'].append(d)


    #write into json
    with open(results_json, 'w') as jsonf:
        json.dump(results, jsonf, sort_keys=False, indent = 4)
    


if __name__ == "__main__":

    # setting up the option parser
    parser = argparse.ArgumentParser()

    #arguments
    parser.add_argument('-i', '--genecovfile', dest = 'genecovfile', help = 'The gene coverage file from HTSeq')
    parser.add_argument('-p', '--genecovpath', dest = 'genecovpath', help = 'The absolute path (string) of the gene coverage file in host system')
    parser.add_argument('-j', '--result_json', dest = 'result_json', help = 'The results.json output file')
    options = parser.parse_args()
    
    if not options.genecovfile or not options.genecovpath or not options.result_json:
        logging.error(' You must provide gene coverage file (--genecovfile), its absolute path in the host directory, and output json file (--result_json)')
        sys.exit(1)
    
    write2json(options.genecovfile, options.result_json, options.genecovpath) 
        



            
