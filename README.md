# Computing Code for the Paper "A Modified Random Survival Forest for Improving Prediction Accuracy in Case-Cohort and Generalized Case-Cohort Studies"

## Description

This repository contains computing codes for the paper "A Modified Random Survival Forest for Improving Prediction Accuracy in Case-Cohort and Generalized Case-Cohort Studies". 

## Naming Convention 

### Folders

This repo contains the following folders, corresponding to the simulation codes for reproducing Section 4 of the paper. The names and descriptions of the folders are as follows,

* *data generation* - this folder contains the code for simulating analysis datasets across 12 simulation settings considered in the Section 4 of the paper.
* *analysis* - this folder contains the analysis code for three predictive models compared in the paper: the naive, weighted, and proposed approaches.
* *installation* - this folder contains the installation of the proposed splitting rule discussed in the paper (see detailed instruction below). 
* *summary* -this folder contains the code for summarizing simulation results. Same code can be used for hyperparameter tuning in general settings.

## Installation Instruction for the Proposed Method

* *Step 1* - Download the package source version by using `download.packages("randomForestSRC", destdir = ".", type = "source")`, which will create a file like 'randomForestSRC_3.6.0.tar.gz'
* *Step 2* - Extract the source using 'untar("randomForestSRC_3.6.0.tar.gz")', This creates a folder: 'randomForestSRC/'




## References

Li, H., Zhou, H., Couper, D., & Cai, J. (2025). Super learner for survival prediction in case-cohort and generalized case-cohort studies. Biometrics, 81(4), ujaf155.

