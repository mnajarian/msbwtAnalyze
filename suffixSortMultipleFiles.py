import csv
import argparse
import itertools
import numpy

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

def generateProbeFiles(pf, oDir):
    probes = gatherProbes(pf)
    sixmers = sorted([''.join(s) for s in itertools.product('ACGT', repeat=6)])
    for kmer in sixmers:
	probeList = []
	for p in probes:
	    if p[0][-6:]==kmer:
		probeList.append(p)
        probeList = sorted(probeList, cmp=bwtsort)
        fp = open(oDir+"/%s.csv" % kmer, 'w')
        fp.write("probeseq,vid,ptype\n")
        for i in xrange(len(probeList)):
            fp.write("%s,%s,%s\n" % probeList[i])
        fp.close()


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('probeFile', help='File of probes in format: sequence, id, type')
    parser.add_argument('outDir', help='Out directory for probe files')
    args = parser.parse_args()
    generateProbeFiles(args.probeFile, args.outDir)
