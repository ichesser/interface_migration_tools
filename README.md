# Overview

This repository contains a set of tutorials for analysis of shuffling patterns during grain boundary (GB) migration to accompany the papers: 

1. Chesser, Ian, Brandon Runnels, and Elizabeth Holm. "A taxonomy of grain boundary migration mechanisms via displacement texture characterization." Acta Materialia 222 (2022): 117425.

2. Chesser, Ian, Elizabeth Holm, and Brandon Runnels. "Optimal transportation of grain boundaries: A forward model for predicting migration mechanisms." Acta Materialia 210 (2021): 116823.

The regularized optimal transport algorithm in paper 2 builds upon the code from the following paper: 

Peyré, Gabriel, and Marco Cuturi. "Computational optimal transport: With applications to data science." Foundations and Trends® in Machine Learning 11.5-6 (2019): 355-607.

# Tutorial Code 

**Example 1**: <br/>
Analyze displacement texture as a post-processing step of an atomistic simulation with a known initial and final state before and after migration. 

**Example 2**: <br/>
Example application of the optimal transport model to predicting migration mechanisms of the Sigma 5 (100) twist GB 

# Dependencies 
Example 1: MTEX: https://mtex-toolbox.github.io/

OVITO Basic (https://www.ovito.org/) will be convenient for visualizing input and output files.
