# ALS Motor Health and Survival Analysis

This MATLAB script processes and analyzes clinical data from ALS Patients, it models the degeneration of neurons over time and visualizes the survival outcomes.

# Parts of the code:
-  Loads the csv file that contains ALS Clinical data
-  Filters and checks the ALSFRS-R and FVC scores
-  Calculates the neuron health and degeneration
-  Derives the survival time for each patient
-  Computes the baseline health score and the median deviations
-  Visualizations:
    - Histograms of baseline motor health
    -  Scatter plot of ALSFRS-R and FVC scores
    -  Deviation from median health
    -  Kaplan–Meier survival curve grouped by baseline health score

# Input Data:
-  Requires a csv file with the following columns:
    -    PatientID (Identifier for each patient)
    -    TimeMonths (time from initial observation)
    -    ALSFRS_R (ALS rating scale from 0-48)
    -    FVC (Forced vital capacity from 0-100)
    -    Status (1 = deceased, 0 = censored(alive or lost to follow up))

# Notes:
-  Max ALSFRS_R score is 48
-  Motor neuron proxy total is 10,000
-  Respiratory neuron proxy total is 8,000
-  Only scores for ALSFRS_R between 0-48 and for FVC between 0-100 are considered
-  Does not handle missing column names
-  Assumes a numeric and clean format for all input data
-  Basic error checking for present columns and the file format
