# Instructions for using the SNP-Slice package

This repository contains the implementation of the SNP-Slice algorithm, which is a Bayesian nonparametric method to resolve multi-strain infections. You can find the motivation for this problem, a description of the algorithm, as well as our results in the Bioarxiv preprint titled **SNP-Slice Resolves Mixed Infections: Simultaneously Unveiling Strain Haplotypes and Linking Them to Hosts** (https://www.biorxiv.org/content/10.1101/2023.07.29.551098v2). 

# How to use the package. 
The structure of the directory contains:
- snpslicemain.R (the main execution file).
- inputdata/ (a directory to store input data files, named *prefix_read1.txt*, *prefix_read0.txt* and *prefix_cat.txt*.
- output/ (a directory to store output data (A,D))
- mcmcRData/ (a directory to store RData files for warm start)
- source/ (a directory containing the actual implementation of the algorithm).

When you run the algorithm, you need to specify a _prefix_ in snpslicemain.R. For example, setting `prefix <- "scenario1"' on line 21 of snpslicemain.R, the script will read scenario1_read1.txt and scenario1_read0.txt from the 'inputata' directory.
