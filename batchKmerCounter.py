#! /usr/bin/env python

## --- msbwtAnalyze/batchKmerCounter.py --- ##
##      Date: 20 March 2017
##      Purpose: given probe files, query an msbwt for each kmer and record the resulting counts

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


def generate_counts(bwtfile, probe_files, out_dir):
    msbwt = MultiStringBWT.loadBWT(bwtfile)
    p_files = glob.glob(probe_files+'/*')
    for probe_file in p_files:
        sixmer = os.path.splitext(probe_file)[0][-6:]
        if os.path.isfile('{}/{}'.format(out_dir, os.path.basename(probe_file))):
            continue
        else:
            print probe_file
            print socket.gethostname()
            start = time.time()
            print 'Started: {}'.format(start)
            
            df = pd.read_csv(probe_file)
            counts = []
            counter = 0
            for row in df.itertuples():
                count = msbwt.countOccurrencesOfSeq(row[1])
                counts.append([count])
                counter += 1
                if (counter & 1023) == 0:
                    print time.time() - start, counter
            finish = time.time() 
            print 'Finish: {}'.format(finish-start)
            with open("{}/{}".format(out_dir, os.path.basename(probe_file)), 'w') as outfile:
                wr = csv.writer(outfile)
                wr.writerow(['counts'])
                wr.writerows(counts)
            finish = time.time()
            print 'Wrote file: {}'.format(finish-start)
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('bwt_file', help='Path to BWT')
    parser.add_argument('probe_files', help='Path to probe files')
    parser.add_argument('out_dir', help='Path to output directory')
    args = parser.parse_args()
    generate_counts(args.bwt_file, args.probe_files, args.out_dir)


