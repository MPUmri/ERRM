# Extended Reference Region Model

## Introduction

This repository contains code for the Extended Reference Region Model (ERRM) and the constrained ERRM (CERRM).
The fitting functions for the two models is in `./mfiles/ERRM.m` and `./mfiles/CERRM.m` respectively.

## Code

This section contains the code used in the study and is organized by the order of figures. 
The code, results, and interactive figure (if available) for each figure is listed.

Figure 1, 2, 3
- Code: `e03_1_makeSimMap.m`, `e03_2_analyzeSimMap.m`, `e03_3_collectSimMap.m`
- Results: `./dataResults/e03a...` and `./dataResults/e03b...`
- Interactive figures: [Fig1](https://cdn.rawgit.com/notZaki/ERRM-xtra/master/interactiveFigures/fig-errMap.html) [Fig2](https://cdn.rawgit.com/notZaki/ERRM-xtra/master/interactiveFigures/fig-errTRes.html) [Fig3](https://cdn.rawgit.com/notZaki/ERRM-xtra/master/interactiveFigures/fig-errTRes.html)

Figure 4, 5 & Table 1
- Code: `e04_1_doBrainSingle.m` 
- Note: This script produces raw figures which were combined with external software

Figures 6, 7
- Code: `e04_2_doBrainDownsample.m`
- Results: `./e04-downsampleGbmResults.csv`


Supplementary Figure S2
- Code: `e01_2_quickSimSingleCurve.m`

Supplementary Figure S3
- Code: `e03_3b_getKepRR_hist.m`

Supplementary Figure S4
- Code: `e04_2_doBrainDownsample.m`
