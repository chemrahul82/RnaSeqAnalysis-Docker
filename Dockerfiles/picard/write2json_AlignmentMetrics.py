#!/usr/bin/env python

"""
writes output metrics from Picard CollectAlignmentSummaryMetrics tool to json file
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
from collections import OrderedDict



def write2json(alignMetricsOut, results_json, metFilePath, readtype):
    """
    Description: Takes the output file of Picard CollectAlignmentSummaryMetrics tool
    and generate a json file with metrics

    Takes: the Picard output file, results.json name/path, Picard metricfile absolute path in host

    Returns: the json file
    """
    #initialize the primary ordered dictionary
    results = OrderedDict()

    
    #initialize the four sections as list
    results['summary_metrics'] = []
    results['additional_metrics'] = []
    results['files'] = []
    results['images'] = []

    
    #filter out the empty lines if any
    with open(alignMetricsOut, "r") as f:
        outlines = [line for line in f if line.strip()]

    
    alignMetrics = {}
    # the last line is for the stats for SE and for both pairs in PE
    
    if options.readtype == 'SE':
        fields = [el for el in outlines[-2].strip('\n').split('\t')]
    elif options.readtype == 'PE':
        fields = [el for el in outlines[-4].strip('\n').split('\t')]
    
    values = [el for el in outlines[-1].strip('\n').split('\t')]
    
    for idx in range(len(fields)):
        alignMetrics[fields[idx]] = values[idx]

    
    # --------first the "three" main metrics in summary_metrics--------
    # no. of AQ20 reads
    d = OrderedDict()
    d['label'] = 'TOTAL_AQ20_READS'
    d['value'] = int(alignMetrics['PF_HQ_ALIGNED_READS'])
    d['type'] = 'reads'
    results['summary_metrics'].append(d)

    # no. of AQ20 bases in AQ20 reads 
    d = OrderedDict()
    d['label'] = 'TOTAL_AQ20_BASES'
    d['value'] = int(alignMetrics['PF_HQ_ALIGNED_Q20_BASES'])
    d['type'] = 'bases'
    results['summary_metrics'].append(d)

    # mean read length
    d = OrderedDict()
    d['label'] = 'MEAN_READ_LENGTH'
    d['value'] = float('%5.1f' %float(alignMetrics['MEAN_READ_LENGTH']))
    results['summary_metrics'].append(d) 

        
    #----------now the additional_metrics-------
    #total reads
    d = OrderedDict()
    d['label'] = 'TOTAL_READS'
    d['value'] = int(alignMetrics['TOTAL_READS'])
    d['type'] = 'reads'
    results['additional_metrics'].append(d)

    # % of AQ20 reads
    d = OrderedDict()
    d['label'] = 'PCT_AQ20_READS'
    d['value'] = float('%4.1f' %(100*int(alignMetrics['PF_HQ_ALIGNED_READS'])/int(alignMetrics['TOTAL_READS'])))
    results['additional_metrics'].append(d)

    # Strand balance
    d = OrderedDict()
    d['label'] = 'STRAND_BALANCE'
    d['value'] = float('%3.2f' %float(alignMetrics['STRAND_BALANCE']))
    results['additional_metrics'].append(d)
        
    # median mismatches within AQ20 reads
    d = OrderedDict()
    d['label'] = 'MEDIAN_MISMATCH_AQ20_READS'
    d['value'] = float(alignMetrics['PF_HQ_MEDIAN_MISMATCHES'])
    results['additional_metrics'].append(d)

    #mismatch rate
    d = OrderedDict()
    d['label'] = 'MISMATCH_RATE'
    d['value'] = float(alignMetrics['PF_MISMATCH_RATE'])
    results['additional_metrics'].append(d)

    #error rate within AQ20 reads
    d = OrderedDict()
    d['label'] = 'ERROR_RATE_AQ20_READS'
    d['value'] = float(alignMetrics['PF_HQ_ERROR_RATE'])
    results['additional_metrics'].append(d)

    #indel rate
    d = OrderedDict()
    d['label'] = 'INDEL_RATE'
    d['value'] = float(alignMetrics['PF_INDEL_RATE'])
    results['additional_metrics'].append(d)

    
    ##files
    d = OrderedDict()
    d['description'] = 'Alignment Summary Metrics'
    d['file_path'] = metFilePath
    d['name'] = 'AlignmentSummaryMetrics.txt'
    results['files'].append(d)
    
    #write into json
    with open(results_json, 'w') as jsonf:
        json.dump(results, jsonf, sort_keys=False, indent = 4)
    


if __name__ == "__main__":

    # setting up the option parser
    parser = argparse.ArgumentParser()

    #arguments
    parser.add_argument('-i', '--metricfile', dest = 'metricfile', help = 'The alignment metrics summary file from Picard')
    parser.add_argument('-p', '--metricpath', dest = 'metricpath', help = 'The absolute path (string) of the Alignment Summary metrics file in host')
    parser.add_argument('-j', '--result_json', dest = 'result_json', help = 'The results.json output file')
    parser.add_argument('-S', '--readtype', dest = 'readtype', help = 'The sequencing type: single end or paired end',default = 'SE')
    options = parser.parse_args()
    
    if not options.metricfile or not options.metricpath or not options.result_json:
        logging.error(' You must provide Picard output alignment summary file (--metricfile), its absolute path in the host (--metricpath), and output json file (--result_json)')
        sys.exit(1)
    
    write2json(options.metricfile, options.result_json, options.metricpath,options.readtype) 
        



            
