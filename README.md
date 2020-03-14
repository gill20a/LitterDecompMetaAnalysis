# Litter Decomposition Meta-Analysis
## Code associated with figures and data analyses in *Gill et. al.(Year/Journal Name)* litter decomposition meta-analysis. 
This repository contains all code to replicate analyses and digures from *Gill et al. 2020. Experimental nitrogen fertilization globally accelerates, then slows leaf litter decomposition*. 

In this study, we aggregated all published litter decomposition data from N fertilization experiments and fit four decomposition models to time series of litter mass loss data. The code and datasets provided here allow you to: 
1. **Model Fitting**. Fit decomposition models (single exponential, double exponential, asymptotic exponential, and Weibull) to litter mass loss time series derived from published studies.
 - **Data**: NFert_MetaA_LitterHarvestData.csv
 - **Metadata**: LitterHarvest_Metadata.csv
 - **Code**: NFert_MetaA_DecompModelFits.R
2. **Summarize Parameters**. Evaluate model parameter fits across the dataset. Perform analyses presented in results paragraph 1 and 2; Figure 1 and S1; Table S1.
 - **Data**: NFert_MetaA_ParamSummary.csv
 - **Metadata**: ParameterSummary_Metadata.csv
 - **Code**: NFert_MetaA_ParamSummary.R  
3. **Evaluate Fertilization Responses**. Calculate parameter log response ratios and evaluate response distribution and relationship with predictor variables. Perform analyses presented in results paragraphs 3-7; Figures 2-4, Table S2-S7.
 - **Data**: NFert_MetaA_ResponseRatios.csv
 - **Metadata**: ResponseRatios_Metadata.csv
 - **Code**: NFert_MetaA_ResponseRatioAnalysis.R
 
The first line of each code file needs to be updated to match the data directory on the user's computer.  This project was built under R version3.6.1 and depends on the following packages: stats4, bbmle, readr, plyr, dplyr, Matrix, Hmisc, gglplot2, nlme, multcomp, MuMIn.
