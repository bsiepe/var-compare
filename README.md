# var-compare

## Overview

This repository contains code to reproduce the analyses in the preprint "Bayesian Estimation and Comparison of Idiographic Network Models" by Siepe & Heck (2023) (https://psyarxiv.com/uwfjc/).
Start the .Rproj **var-compare.Rproj** prior to running the scripts.

Full simulation results are too large (multiple GBs) and can be requested from the corresponding author. We try to provide intermediate simulation results where possible. 

## main folder
This folder contains the relevant main scripts to reproduce the analyses. 

All auxiliary functions used in the project can be found in **aux_funs.R**. 
The data-generating true networks can be obtained using **true-networks.Rmd**.
Simulation Study 1 can be reproduced using **simulation-1.Rmd**.
Simulation Study 2 can be reproduced using **simulation-2.Rmd**.
The empirical example can be reproduced using **empirical-example.Rmd**.
Session information can be found in **session_info.txt**.

### Revision 1

The scripts for the first revision of the paper can be found in the folder **revision-1**. 
Attempts to create tests that obtain evidence for the Null can be found under **between-compare-testing.Rmd**, 
**simulation-2-revision-diffpost-both.Rmd**, **simulation-2-revision-diffpost.Rmd**, and **simulation-2-revision-meanpost.Rmd**.

Code which tests the implementation of the model in Stan, including code for the simulation study, can be found in the repository [stan-gvar](https://github.com/bsiepe/stan-gvar). This repository also contains example code and a light wrapper to fit a Bayesian GVAR model in MPlus.


### data
Contains data-generating processes used in simulation study 1 & 2. 

### output
Contains intermediate results for simulation study 1. 





