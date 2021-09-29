The scripts are from https://blockify.readthedocs.io/en/latest/index.html

blockify is a fast and optimal genomic peak caller for one-dimensional data (e.g. BED, qBED, CCF).

The package is built around the Bayesian blocks algorithm [SNJC13], which finds the optimal change points in time series data assuming a Poisson counting process. 
We also implement a dynamic pruning strategy which achieves linear runtime performance [KFE12]. An interactive notebook demonstrating Bayesian blocks can be found here.

While Bayesian blocks was originally developed in the astrophysics community for photon-counting experiments, we find that it has applications in genomics. 
In particular, we use it to analyze transposon calling cards experiments. 
Calling cards uses a transposase fused to a transcription factor (TF) to deposit transposons near TF binding sites. 
Bayesian blocks partitions the genome based on the local density of insertions, which in turn are used to identify peaks and candidate TF binding sites. 
We have also had success using this algorithm to perfom general-purpose genome segmentation.
