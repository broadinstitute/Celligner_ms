# Celligner_ms
This repo contains code associated with the manuscript describing Celligner, our method for aligning tumor and cell line transcriptional profiles.

## Data

The data associated with this analysis will be made available shortly.

## Organization of repo

The code can be organized into config files, helper functions, and analysis/figure generation scripts.

NOTE: The functions can be run using data available with the manuscript and data from publicly available resources (primarily depmap.org and xena browser)

### Configs

global_params.R: Define global params shared across analysis scripts. Includes parameters used to run Celligner alignment and parameters used for creating plots.

### Helper functions

- analysis_helpers.R: Define helper functions used throughout the analysis and creation of figures
- Celligner_helpers.R : Define helper functions used for the Celligner alignment method

### Analysis/fig-gen

- Celligner_methods.R : Functions to run the various stages and entire Celligner alignment method
- There are separate scripts for each of the main and supplementary figure panels within the manuscript.
