#!/usr/bin/env python

"""
writes output metrics from Picard CollectRnaSeqMetrics tool to json file & make figures
"""
from __future__ import division

__author__='Rahul K. Das'


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse
import json
import logging
import os
import subprocess
import sys
import time
from collections import OrderedDict


def write2json(RnaSeqMetricsOut, xrRNA_basecount, results_json, fig, metFilePath):
    """
    Description: Takes the output file of Picard CollectRnaSeqMetrics tool
    and generate a json file with metrics

    Takes: the Picard output file, file with total counts of mapped bases to xrRNA, results.json path, 
    Figure path, Picard metricfile absolute path in host

    Returns: the json file, the summary figure
    """

    #initialize the primary ordered dictionary
    results = OrderedDict()

    #initialize the four sections as list
    results['summary_metrics'] = []
    results['additional_metrics'] = []
    results['files'] = []
    results['images'] = []

   
    #get the rRNA basecount
    #some of the metrics need to be renormalized to account for rRNA mapped bases
    if os.path.exists(xrRNA_basecount):
        with open(xrRNA_basecount, 'r') as f:
            firstline = f.readline().strip('\n')
            if firstline != '':
                xrRNA_bases = int(firstline)
            else:
                xrRNA_bases = 0
    else:
        xrRNA_bases = 0


    # grab the metrics from picard output file
    rnaSeqMetrics = {}

    # the 7th & 8th lines are metrics and values, respectively
    # CAUTION: this hardcoding only works with the Picard version it was wrapped with
    with open(RnaSeqMetricsOut, "r") as f:
        for idx, line in enumerate(f):
            if idx == 6:
                fields = [el for el in line.strip('\n').split('\t')]
            elif idx == 7:
                values = [el for el in line.strip('\n').split('\t')]

    for idx in range(len(fields)):
        rnaSeqMetrics[fields[idx]] = values[idx]

    
    # --------first the "three" main metrics in summary_metrics--------
   
    tot_bases = int(rnaSeqMetrics['PF_BASES'])

    #update the aligned base to include 5.8S,18S,28S rRNA alighments
    tot_align_bases = int(rnaSeqMetrics['PF_ALIGNED_BASES']) + xrRNA_bases
    
    # Pct. of Usable bases: only this is normalized by total bases; other Pct. metrics are normalized by total aligned bases
    cod_bases = int(rnaSeqMetrics['CODING_BASES'])
    utr_bases = int(rnaSeqMetrics['UTR_BASES'])
    mRNA_bases = int(rnaSeqMetrics['CODING_BASES']) + int(rnaSeqMetrics['UTR_BASES'])

    d = OrderedDict()
    d['label'] = 'PCT_USABLE_BASES'
    d['value'] = float('%4.1f' %(100*mRNA_bases/tot_bases))
    results['summary_metrics'].append(d)

    # Pct. of mRNA bases 
    d = OrderedDict()
    d['label'] = 'PCT_mRNA_BASES'
    d['value'] = float('%4.1f' %(100*mRNA_bases/tot_align_bases))
    results['summary_metrics'].append(d)
    
    # Pct. of ribosomal bases; include 5.8S,18S,28S rRNA alignments
    rRNA_bases = int(rnaSeqMetrics['RIBOSOMAL_BASES']) + xrRNA_bases

    d = OrderedDict()
    d['label'] = 'PCT_rRNA_BASES'
    d['value'] = float('%4.1f' %(100*rRNA_bases/tot_align_bases))
    results['summary_metrics'].append(d)

   
    #----------now the additional_metrics-------
    
    #total bases
    d = OrderedDict()
    d['label'] = 'TOTAL_BASES'
    d['value'] = tot_bases
    d['type'] = 'bases'
    results['additional_metrics'].append(d) 
        
    #total aligned bases
    d = OrderedDict()
    d['label'] = 'ALIGNED_BASES'
    d['value'] = tot_align_bases
    d['type'] = 'bases'
    results['additional_metrics'].append(d) 
     
    #total coding bases
    d = OrderedDict()
    d['label'] = 'CODING_BASES'
    d['value'] = cod_bases
    d['type'] = 'bases'
    results['additional_metrics'].append(d) 
     
    #total UTR bases
    d = OrderedDict()
    d['label'] = 'UTR_BASES'
    d['value'] = utr_bases
    d['type'] = 'bases'
    results['additional_metrics'].append(d) 
     
    #total mRNA bases
    d = OrderedDict()
    d['label'] = 'mRNA_BASES'
    d['value'] = mRNA_bases
    d['type'] = 'bases'
    results['additional_metrics'].append(d) 
    
    #total intronic bases
    intron_bases = int(rnaSeqMetrics['INTRONIC_BASES'])
    d = OrderedDict()
    d['label'] = 'INTRONIC_BASES'
    d['value'] = intron_bases
    d['type'] = 'bases'
    results['additional_metrics'].append(d) 
    
    #total intergenic bases
    intergen_bases = int(rnaSeqMetrics['INTERGENIC_BASES'])
    d = OrderedDict()
    d['label'] = 'INTERGENIC_BASES'
    d['value'] = intergen_bases
    d['type'] = 'bases'
    results['additional_metrics'].append(d) 
    
    #total rRNA bases
    d = OrderedDict()
    d['label'] = 'rRNA_BASES'
    d['value'] = rRNA_bases
    d['type'] = 'bases'
    results['additional_metrics'].append(d) 
    
   
    #now the percentages
    #aligned bases
    d = OrderedDict()
    d['label'] = 'PCT_ALIGNED_BASES'
    d['value'] = float('%4.1f' %(100*tot_align_bases/tot_bases))
    results['additional_metrics'].append(d) 
     
    #coding bases
    d = OrderedDict()
    d['label'] = 'PCT_CODING_BASES'
    d['value'] = float('%4.1f' %(100*cod_bases/tot_align_bases))
    results['additional_metrics'].append(d) 
     
    #total UTR bases
    d = OrderedDict()
    d['label'] = 'UTR_BASES'
    d['value'] = float('%4.1f' %(100*utr_bases/tot_align_bases))
    results['additional_metrics'].append(d) 
     
    #intronic bases
    d = OrderedDict()
    d['label'] = 'PCT_INTRONIC_BASES'
    d['value'] = float('%4.1f' %(100*intron_bases/tot_align_bases))
    results['additional_metrics'].append(d) 

    #intergenic bases
    d = OrderedDict()
    d['label'] = 'PCT_INTERGENIC_BASES'
    d['value'] = float('%4.1f' %(100*intergen_bases/tot_align_bases))
    results['additional_metrics'].append(d) 


    ##files
    d = OrderedDict()
    d['description'] = 'RnaSeq QC metrics'
    d['file_path'] = metFilePath
    d['name'] = 'RnaSeqMetrics.txt'
    results['files'].append(d)

    
    #make plots
    #plt.figure(figsize=(20,10))
    labels = ['rRNA', 'Coding', 'UTR', 'Intronic', 'Intergenic']
    sizes = [float('%4.1f' %(100*rRNA_bases/tot_align_bases)), \
            float('%4.1f' %(100*cod_bases/tot_align_bases)), \
            float('%4.1f' %(100*utr_bases/tot_align_bases)), \
            float('%4.1f' %(100*intron_bases/tot_align_bases)), \
            float('%4.1f' %(100*intergen_bases/tot_align_bases))]

    colors = ['red','yellowgreen', 'lightgreen', 'lightcoral', 'coral']
    explode = (0.1, 0, 0, 0, 0)  # explode 1st slice
    plt.pie(sizes, explode=explode, labels=labels, colors=colors,\
            autopct='%1.1f%%', shadow=True, startangle=140)

    plt.axis('equal')
    #plt.tight_layout()
    plt.savefig(fig, format='png', dpi=300) 
    
    figpath = os.path.join(os.path.dirname(metFilePath), os.path.basename(fig))

    d = OrderedDict()
    d['description'] = 'RnaSeq QC metrics figure'
    d['file_path'] = figpath
    d['name'] = os.path.basename(fig)
    results['images'].append(d)
    
    #write into json
    with open(results_json, 'w') as jsonf:
        json.dump(results, jsonf, sort_keys=False, indent = 4)
    

if __name__ == "__main__":

    # setting up the option parser
    parser = argparse.ArgumentParser()

    #arguments
    parser.add_argument('--metricfile', dest = 'metricfile', help = 'The RnaSeq metrics summary output file from Picard')
    parser.add_argument('--xrRNAbasecount', dest = 'xrRNAbasecount', help = 'The file containing the counts for mapped rRNA bases')
    parser.add_argument('--metricpath', dest = 'metricpath', help = 'The absolute path (string) of the RnaSeqSummary metrics file in host system')
    parser.add_argument('--result_json', dest = 'result_json', help = 'The results.json output file')
    parser.add_argument('--fig', dest = 'fig', help = 'The summary figure')
    options = parser.parse_args()
    
    if not options.metricfile or not options.metricpath or not options.result_json or not options.fig or not options.xrRNAbasecount:
        logging.error(' You must provide Picard alignment metric file (--metricfile), its absolute path in the host (--metricpath), xrRNA basecount file (--xrRNAbasecount), output json (--result-json), and output figure (--fig)')
        sys.exit(1)
    
    write2json(options.metricfile, options.xrRNAbasecount, options.result_json, options.fig, options.metricpath) 
        



            
