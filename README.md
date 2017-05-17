# msbwtAnalyze
Scripts for analyzing next-generation sequencing data using multi-string BWTs ([msBWTs](https://github.com/holtjma/msbwt)).

- `batchKmerCounter.py`: given a single probe file, query an msbwt for each kmer and record the resulting counts
- `batchKmerCounterMultipleFiles.py`: given a set of probe files, query an msbwt for each kmer and record the resulting counts
- `fastBatchKmerCounter.py`: given a single probe file, query an msbwt using a fast stack implementation and record the resulting counts
- `fastBatchKmerCounterMultipleFiles.py`: given a set of probe files, query an msbwt using a fast stack implementation and record the resulting counts
- `suffixSort.py`: sort a given probe file by suffix
- `suffixSortMultipleFiles.py`: split up a single probe file into multiple sorted files by suffix
- `mergeCounts.py`: merge probe files and their corresponding count files
