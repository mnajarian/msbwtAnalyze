#! /usr/bin/env python

import csv
import argparse
import itertools
import numpy

def bwtOrder(kmer):
    basemap = {'A':0, 'C':1, 'G':2, 'T':3}
    value = 0
    for base in reversed(kmer):
        value = (value*4) + basemap[base]
    return value

def gatherProbes(pf):
    probes = []
    with open(pf, 'rb') as f:
	reader = csv.reader(f)
	next(reader)
	for row in reader:
	    probes.append((row[0], row[1], row[2]))
    return probes

def generateProbeFiles(pf, oDir):
    probes = gatherProbes(pf)
    sixmers = sorted([''.join(s) for s in itertools.product('ACGT', repeat=6)])
    for kmer in sixmers:
	probeList = []
	for p in probes:
	    if p[-6:]==s:
		probeList.append(p)
	probeInt = numpy.empty((len(probeList),), dtype='uint64')
	for i, (seq, vid, ptype) in enumerate(probeList):
    	    probeInt[i] = bwtOrder(seq)
	order = numpy.argsort(probeInt, kind='heapsort')
	fp = open(oDir+"/%s.csv" % kmer, 'w')
	fp.write("probeseq,vid,ptype\n")
	for i in order:
    	    fp.write("%s,%d,%s\n" % probeList[i])
	fp.close()


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('probeFile', help='File of probes in format: sequence, id, type')
    parser.add_argument('outDir', help='Out directory for probe files')
    args = parser.parse_args()
    generateProbeFiles(args.probeFile, args.outDir)
