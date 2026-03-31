PAM-BirdNET-Code description

This dissertation sought to quantify the probability of detection and classification precision of 20 Scottish woodland birds 
when using the automated classifier BirdNET to label their vocalisation.

Detectability was classified as at least one positive ID match between the expected_species and classification.

Classification precision was the the proportion of correct detections, given a positive detection event.

Scripts are listed in order of execution:
1. Data_processing_and_probability_of_detection_models
2. Classification_precision_models + community_vocal_activity_metric
3. Community_vocal_activity_analysis
4. Confusion_matrix
5. Figures + supplementary_material
6. Detection_thresholds + correlation_analysis

All analysis was conducted in R studio 2025.09.0 Build 387, with required packages listed below:
library(dplyr)
library(tidyr)
library(glmmTMB)
library(tidyverse)
library(car)
library(performance)
library(ggplot2)
library(patchwork)
library(scales)
library(DHARMa)
library(ggrepel)

Analysis requires the following CSV files in order to operate:
"Combined Detection Table.csv"
"Diss covariates no LAI.csv"
"LAI conversion.csv"
