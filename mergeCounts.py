#! /usr/bin/env python

import csv
import argparse
import glob
import pandas
import os


def mergeFiles(countFiles, probesDir, outDir):
    print countFiles
    resultList = []
    for cf in glob.glob(countFiles):
        countdf = pandas.read_csv(cf)
        probedf = pandas.read_csv(probesDir+'/'+os.path.basename(cf))
        joindf = probedf.join(countdf)
        resultList.append(joindf)
    result = pandas.concat(resultList, ignore_index=True)
    result.to_csv(outDir+'/'+os.path.basename(cf))



def merge(countsDir, probesDir, outDir):
    for countFiles in glob.glob(countsDir):
        mergeFiles(countFiles, probesDir, outDir)



if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('countsDir', help='Parent directory of all count directories')
    parser.add_argument('probesDir', help='Directory containing all probe files')
    parser.add_argument('outDir', help='Path to output directory')
    args = parser.parse_args()
    merge(args.countsDir, args.probesDir, args.outDir)
    
