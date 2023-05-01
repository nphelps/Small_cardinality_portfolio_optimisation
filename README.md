# Small cardinality portfolio optimisation
This is an algorithm for generating efficient portfolios with small cardinality constraints. There are four variations of an algorithm created by Graham and Craven (https://www.tandfonline.com/doi/full/10.1080/01605682.2020.1718019). These four variations will be described in a paper and a link to the paper will be provided here when it is published.

The CCEF Output folder contains points along cardinality constrained efficient frontiers computed for various datasets using the algorithms. The number of sub-efficient frontiers flagged as useful by each variation of the algorithm and runtimes are also available in the corresponding folders. R code for the algorithms are in the associated folder. The files provided were designed to be run on SHARCNET, but they can be run locally.

To successfully use the code provided in this repository, it is necessary to download get_EF_lambda.R and sift_v3.R from here: https://github.com/MJCraven/SiftedQP. This repository also contains an implementation of Graham and Craven's original algorithm.
