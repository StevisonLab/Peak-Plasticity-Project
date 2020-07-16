# Peak-Plasticity-Project
This repository contains data files and scripts to reproduce the analysis for the manuscript titled: **"refining the peak time point for recombination rate plasticity in *Drosophila pseudoobscura*"**

# Citation

[![DOI](https://zenodo.org/badge/254474381.svg)](https://zenodo.org/badge/latestdoi/254474381)

# Abstract
Meiotic recombination rates vary in response to environmental factors like temperature. Variation in recombination generates an additional source for genetic variation while errors in this pathway can lead to chromosome nondisjunction. Estimating duration and sensitivity of a meiotic response to environmental perturbation requires an understanding of molecular events and well-designed experimental approaches. An earlier study suggested that the peak (most sensitive) timing of plasticity in Drosophila melanogaster occurred during the pre-meiotic interphase where DNA replication takes place in S-phase. Recently, heat stress has been shown to reveal plasticity in recombination rates in D. pseudoobscura. Here, a combination of molecular genotyping and a series of recombination analyses through visible phenotypic markers were used to determine peak plasticity timing in this species. Mutant flies were reared in either control or stress temperatures in a traditional cross design. Using mixed model analysis, the odds of crossover formation was 1.55X higher during days 7-9 (p<0.0017) and 1.41X higher on day 9 (p<0.034) in high temperature as compared to control crosses, suggesting the time period as the timing of peak plasticity. Time of peak plasticity at day 9 in D. pseudoobscura can be explained by comparison to the model organism D. melanogaster due to similar timing of key meiotic events. This comparative approach will enable future mechanistic work on the duration and the effects of temperature stress on recombination rate.  

# Data Analysis

The R code used for all statistical analysis and to produce all data figures and model tables is included in the folder ["Scripts"](https://github.com/StevisonLab/Peak-Plasticity-Project/tree/master/Scripts) as an R script divided into each section named after sequential experiments. 

The folder ["raw_data_files"](https://github.com/StevisonLab/Peak-Plasticity-Project/tree/master/raw_data_files) includes the backcrosses, raw phenotyping results, fecundity, and female count data from the mutant screens. A separate csv file with treatment information includes dates and other metadata that would be needed to validate the analysis and conclusions herein. Additionally, a processed and cleaned up data file that includes sums of males and females, as well as crossover counts across each interval per time point, per replicate vial is also included. 
As an example a tutorial from experiment 3 in the manuscript can be found in the ["Experiment3_Tutorial"](https://github.com/StevisonLab/Peak-Plasticity-Project/tree/master/Experiment3_Tutorial) folder with all of the results and figures for Experiment 3.
