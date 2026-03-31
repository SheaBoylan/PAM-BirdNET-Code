# PAM-BirdNET Analysis

This repository contains the R code for a dissertation 
quantifying the probability of detection and classification precision of 20 Scottish woodland birds 
when using the automated classifier BirdNET to label their vocalisation.

**Probability of detection** was classified as at least one positive ID match between the expected_species and classification.

**Classification precision** was the the proportion of correct detections, given a positive detection event.

>Note: The variable 'accuracy' used throughout the code corresponds to 
>'classification precision' as used in the dissertation text.

>Note: The variable 'masking' refers to
> 'community vocal' activity in dissertation text.

##Script order of execution:
1. Data_processing_and_probability_of_detection_models
2. Classification_precision_models + community_vocal_activity_metric
3. Community_vocal_activity_analysis
4. Confusion_matrix
5. Figures + supplementary_material
6. Detection_thresholds + correlation_analysis

## Required Packages
```r
install.packages(c("dplyr", "tidyr", "glmmTMB", "tidyverse", "car", 
                   "performance", "ggplot2", "patchwork", "scales", 
                   "DHARMa", "ggrepel"))
```

#Required data files
Analysis requires the following CSV files in order to operate:
-'Combined Detection Table.csv'
-'Diss covariates no LAI.csv'
-'LAI conversion.csv'

##Software
All analysis conducted in RStudio 2025.09.0 Build 387.
