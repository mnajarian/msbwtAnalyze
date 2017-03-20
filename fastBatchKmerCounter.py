#! /usr/bin/env python

## --- msbwtAnalyze/fastBatchKmerCounter.py --- ##
##      Date: 17 March 2017
##      Purpose: given probe files, query an msbwt using a fast stack implementation 

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

def sharedPrefixLength(s1, s2):
    # Find the shared prefix length between s1 and s2
    p = 0
    for i in xrange(min(len(s1), len(s2))):
        if s1[i] == s2[i]:
            p += 1
        else:
            break
    return p

def build_stack(msbwt, stack, indices, shared_next, q):
    stack_end = len(stack[-1])
    if shared_next == stack_end:
        # Don't do anything to stack
        lo,hi = indices[-1]
        c = msbwt.countOccurrencesOfSeq(''.join(s for s in reversed(q[stack_end:])), (lo,hi))
    elif shared_next < stack_end:
        # need to cut stack
        lo,hi = indices[-1]
        if q[stack_end:] == '':
            c = hi-lo
        else:
            c = msbwt.countOccurrencesOfSeq(''.join(s for s in reversed(q[stack_end:])), (lo,hi))
        stack = stack[:-(stack_end-shared_next)]
        indices = indices[:-(stack_end-shared_next)]
    else:
        # need to add to stack
        base = stack[-1]
        for i in range(0,shared_next-stack_end):
            p_lo, p_hi = indices[-1]
            stack.append(base+q[stack_end+i])
            lo, hi = msbwt.findIndicesOfStr(q[stack_end+i], (p_lo, p_hi))
            indices.append([lo, hi])
            base += q[stack_end+i]
        lo,hi = indices[-1]
        if q[len(stack[-1]):] == '':
            c = hi-lo
        else:
            c = msbwt.countOccurrencesOfSeq(''.join(s for s in reversed(q[len(stack[-1]):])), (lo,hi))
    return stack, indices, c

def get_shared_prefixes(seqs):
    shared = []
    last = seqs[0]
    for seq in seqs[1:]:
        spl = sharedPrefixLength(seq, last)
        shared.append(spl)
        last = seq
    shared.append(0)
    return shared

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
            rev = [''.join([s for s in reversed(t)]) for t in df['probeseq']]
            # PREPROCESSING STEP: for each query, get prefix index shared with NEXT query
            shared = get_shared_prefixes(rev)
        
            sixmer_reverse = ''.join([s for s in reversed(sixmer)])
            q_len = len(df['probeseq'].iloc[0]) - 6
            stack = []
            stack.append(sixmer_reverse)
            indices = []
            lo, hi = msbwt.findIndicesOfStr(sixmer)
            indices.append([lo,hi])
            counts = []
            counter = 0
            for row in df.itertuples():
                q = ''.join([s for s in reversed(row[1])])
                shared_next = shared[counter]
                stack, indices, count = build_stack(msbwt, stack, indices, shared_next, q)
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


