#! /usr/bin/env python

## --- msbwtAnalyze/batchKmerCounter.py --- ##
##      Date: 18 May 2017
##      Purpose: given single probefile, query an msbwt for each kmer and record the resulting counts

import pysam
import MUSCython
import MUSCython.MultiStringBWTCython as MultiStringBWT
import glob
import pandas as pd
import os
import time
import csv
import argparse
import socket
import itertools

def gatherProbes (pf):
    probes = []
    with open(pf, 'rb') as f:
        reader = csv.reader(f)
        next(reader)
        for row in reader:
            probes.append((row[0], row[1], row[2]))
    return probes

def generate_counts(bwtfile, probe_file, out_file):
    msbwt = MultiStringBWT.loadBWT(bwtfile)
    start = time.time()
    print 'Started: {}'.format(start)
    probes = gatherProbes(probe_file)
    counts = []
    counter = 0
    for row in probes:
        counts.append(msbwt.countOccurrencesOfSeq(row[0]))
        counter += 1
        if (counter & 1023) == 0:
            print time.time() - start, counter
    finish = time.time()
    print 'Finish: {}'.format(finish-start)
    with open(out_file,'w') as outfile:
        wr = csv.writer(outfile)
        wr.writerow(['probeseq', 'vid', 'ptype', 'count'])
        wr.writerows([x+(y,) for x,y in itertools.izip(probes, counts)])
    finish = time.time()
    print 'Wrote file: {}'.format(finish-start)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('bwt_file', help='Path to BWT')
    parser.add_argument('probe_file', help='Path to probe file')
    parser.add_argument('out_file', help='Path to output file')
    args = parser.parse_args()
    generate_counts(args.bwt_file, args.probe_file, args.out_file)


