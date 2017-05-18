#! /usr/bin/env python

## --- msbwtAnalyze/sortedBatchKmerCounter.py --- ##
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

def bwtsort(atup, btup):
    aseq = atup[0]
    bseq = btup[0]
    for i in xrange(len(aseq)):
        abase = aseq[-1-i]
        bbase = bseq[-1-i]
        if (abase > bbase):
            return 1
        elif (abase < bbase):
            return -1
    else:
        if (atup[1] > btup[1]):
            return 1
        elif (atup[1] < btup[1]):
            return -1
    return 0

def gatherProbes(pf):
    probes = []
    with open(pf, 'rb') as f:
        reader = csv.reader(f)
        next(reader)
        for row in reader:
            probes.append((row[0], row[1], row[2]))
    return probes

def generate_counts(bwtfile, probe_file, out_file):
    start = time.time()
    print 'Started: {}'.format(start)
    msbwt = MultiStringBWT.loadBWT(bwtfile)
    probes = gatherProbes(probe_file)
    probesSorted = sorted(probes, cmp=bwtsort)
    counts = []
    counter = 0
    for row in probesSorted:
        counts.append(msbwt.countOccurrencesOfSeq(row[0]))
        counter += 1
        if (counter & 1023) == 0:
            print time.time() - start, counter
    finish = time.time()
    print 'Finish: {}'.format(finish-start)
    with open(out_file,'w') as outfile:
        wr = csv.writer(outfile)
        wr.writerow(['probeseq', 'vid', 'ptype','counts'])
        wr.writerows([x+(y,) for x,y in itertools.izip(probesSorted, counts)])
    finish = time.time()
    print 'Wrote file: {}'.format(finish-start)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('bwt_file', help='Path to BWT')
    parser.add_argument('probe_file', help='Path to probe file')
    parser.add_argument('out_file', help='Path to output file')
    args = parser.parse_args()
    generate_counts(args.bwt_file, args.probe_file, args.out_file)


