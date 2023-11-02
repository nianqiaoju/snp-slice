# Instructions for using the SNP-Slice package

This repository contains the implementation of the SNP-Slice algorithm, which is a Bayesian nonparametric method to resolve multi-strain infections. You can find the motivation for this problem, a description of the algorithm, as well as our results in the Bioarxiv preprint titled **SNP-Slice Resolves Mixed Infections: Simultaneously Unveiling Strain Haplotypes and Linking Them to Hosts** (https://www.biorxiv.org/content/10.1101/2023.07.29.551098v2). 

## Preparing your local directory to run SNP-Slice.
The structure of the directory contains:
- `snpslicemain.R` (the main execution file).
- `inputdata/` (a directory to store input data files, named *prefix_read1.txt*, *prefix_read0.txt* and *prefix_cat.txt*.
- `output/` (a directory to store output data (A,D))
- `mcmcRData/` (a directory to store RData files for warm start)
- `source/` (a directory containing the actual implementation of the algorithm).

## Using the algorithmn.
1. First of all, specify a _prefix_ in `snpslicemain.R`.
   
   For example, setting `prefix <- "scenario1" ` on line 21 of  `snpslicemain.R`, the script will read `scenario1_read1.txt` and `scenario1_read0.txt` from the `inputata` directory.
2. Now you can run the algorithm in the command line, with, for example,

  `Rscript snpslicemain.R model=3 nmcmc=10000 alpha=2`.

3. You can also decide which model to use, by controlling the value of `model`. We recommend setting `model` in the command line instead of in the execution file. The default value is Negative Binomial model.
  This is the codebook:
- ` model <- 0` for the cat model
- ` model <- 1` for the Poisson model
- ` model <- 2` for the Binomial model
- ` model <- 3` for the Negative Binomial model.
