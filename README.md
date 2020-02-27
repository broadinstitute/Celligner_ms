# Celligner_ms
This repo contains code associated with our manuscript describing Celligner, our method for aligning tumor and cell line transcriptional profiles.

## Data

The data associated with this analysis will be made available shortly.

## Organization of repo

The code can be organized into config files, helper functions, preprocessing scripts, and analysis/figure generation scripts.

## Configs

global_params.R: Define global params shared across analysis scripts. Includes parameters used to run Celligner alignment and parameters used for creating plots.

## NOTES:

The functions can be run using data available with the manuscript and data from publicly available resources (primarily depmap.org and xena browser)

## Helper functions

- analysis_helpers.R: Define helper functions used throughout the analysis and creation of figures
- Celligner_helpers.R : Define helper functions used for the Celligner alignment method
- Celligner_methods.R : Functions to run the various stages and entire Celligner alignment method

## Analysis/fig-gen

There are separate scripts for each of the main and supplementary figure panels within the manuscript.
