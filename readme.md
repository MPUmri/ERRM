# Extended Reference Region Model

## Introduction

This repository contains code for the Extended Reference Region Model (ERRM) and the constrained ERRM (CERRM).

These models are described in:
[Ahmed, Z., & Levesque, I. R. (2018). An extended reference region model for DCE-MRI that accounts for plasma volume. NMR in Biomedicine](https://onlinelibrary.wiley.com/doi/abs/10.1002/nbm.3924)

The fitting functions for the two models is in `./mfiles/ERRM.m` and `./mfiles/CERRM.m` respectively.

## Code

This section contains the code used in the study and is organized by the order of figures. 
The code, results, and interactive figure (if available) for each figure is listed.

**Figure 1, 2, 3**
- Code: `e03_1_makeSimMap.m`, `e03_2_analyzeSimMap.m`, `e03_3_collectSimMap.m`
- Results: `./dataResults/e03a...` and `./dataResults/e03b...`
- Interactive figures: [Fig1](https://notzaki.github.io/interactiveFigures/fig-errMap.html) [Fig2](https://notzaki.github.io/interactiveFigures/fig-errSim.html) [Fig3](https://notzaki.github.io/interactiveFigures/fig-errTRes.html)

**Figure 4, 5** & **Table 1**
- Code: `e04_1_doBrainSingle.m` 
- Note: This script produces raw figures which were combined with external software

**Figures 6, 7**
- Code: `e04_2_doBrainDownsample.m`
- Results: `./e04-downsampleGbmResults.csv`

**Supplementary Figure S2**
- Code: `e01_2_quickSimSingleCurve.m`

**Supplementary Figure S3**
- Code: `e03_3b_getKepRR_hist.m`

**Supplementary Figure S4**
- Code: `e04_2_doBrainDownsample.m`

## Note on acronyms

The manuscript discusses the following models:  
- Tofts Model (TM)
- Extended Tofts Model (ETM)
- Reference Region Model (RRM)
- Extended Reference Region Model (ERRM)
- Constrained Extended Reference Region Model (CERRM)

The code and the interactive includes three additional models/methods:
- Constrained Reference Region Model (CRRM) - [Ref](https://onlinelibrary.wiley.com/doi/abs/10.1002/mrm.26530)
- Reference Tissue Method (RTM) - briefly described in manuscript
- Extended Reference Tissue Method (ERTM)
    + Both the RTM and ERTM estimate the arterial input function from a reference tissue curve. The difference between them is that the RTM uses the input function to fit the TM to the tissue of interest, whereas the ERTM fits the ETM.
    + The RTM and ERTM have previously been explored in [Ref](http://iopscience.iop.org/article/10.1088/0031-9155/53/10/012/meta) 

Even though these latter approaches did not play a major role in the manuscript, they have been included in the code and interactive figures on the off chance that someone finds them useful.
