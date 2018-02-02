# Extended Reference Region Model

## Introduction

This repository contains code for the Extended Reference Region Model (ERRM) and the constrained ERRM (CERRM).
The fitting functions for the two models is in `./mfiles/ERRM.m` and `./mfiles/CERRM.m` respectively.

## Code

**Note:** This code is for an older revision. Will be updated soon.

This section contains the code used in the study and is organized by the order of figures. 
The code, results, and interactive figure (if available) for each figure is listed.

Figure 2
- Code: `e01_2_quickSimSingleCurve.m`

Figure 3
- Code: `e01_1_quickSim.m`

Figures 4,5,6
- Code: `e03_1_makeSimMap.m`, `e03_2_analyzeSimMap.m`, `e03_3_collectSimMap.m`
- Results: `./dataResults/e03a...` and `./dataResults/e03b...`
- Interactive figures: [Fig4](http://htmlpreview.github.io/?https://github.com/notZaki/ERRM-xtra/blob/master/interactiveFigures/fig4.html) [Fig5](http://htmlpreview.github.io/?https://github.com/notZaki/ERRM-xtra/blob/master/interactiveFigures/fig5.html) [Fig6](http://htmlpreview.github.io/?https://github.com/notZaki/ERRM-xtra/blob/master/interactiveFigures/fig6.html)
  + Note: The htmlpreview service is having a few hiccups recently and the interactive figures might not display the first time the page loads. Refreshing the page a couple of times eventually makes the page correctly. 

Figure 7,8
- Code: `e04_1_doBrainSingle.m`

Figure 9,10
- Code: `e04_2_doBrainDownsample.m`
- Results: `./e04-downsampleGbmResults.csv`
