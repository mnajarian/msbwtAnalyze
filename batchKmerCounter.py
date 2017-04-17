#! /usr/bin/env python

## --- msbwtAnalyze/batchKmerCounter.py --- ##
##      Date: 20 March 2017
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


def generate_counts(bwtfile, probe_file, out_file):
    msbwt = MultiStringBWT.loadBWT(bwtfile)
    start = time.time()
    print 'Started: {}'.format(start)
    df = pd.read_csv(probe_file)
    counts = []
    counter = 0
    for row in df.itertuples():
        counts.append([msbwt.countOccurrencesOfSeq(row[1])])
        counter += 1
        if (counter & 1023) == 0:
            print time.time() - start, counter
    finish = time.time()
    print 'Finish: {}'.format(finish-start)
    with open(out_file,'w') as outfile:
        wr = csv.writer(outfile)
        wr.writerow(['counts'])
        wr.writerows(counts)
    finish = time.time()
    print 'Wrote file: {}'.format(finish-start)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('bwt_file', help='Path to BWT')
    parser.add_argument('probe_file', help='Path to probe file')
    parser.add_argument('out_file', help='Path to output file')
    args = parser.parse_args()
    generate_counts(args.bwt_file, args.probe_file, args.out_file)


