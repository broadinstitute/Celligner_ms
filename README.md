# Celligner_ms
This repo contains code associated with the manuscript describing Celligner, our method for aligning tumor and cell line transcriptional profiles.

## Data

The data associated with this analysis is available from public data repositiories, supplementary data files associated with the manuscript (https://www.biorxiv.org/content/10.1101/2020.03.25.008342v1.supplementary-material), and in the figshare: https://figshare.com/articles/Celligner_data/11965269.

The cell line data used as input can be found at depmap.org (the file is DepMap Public 19Q4 CCLE_expression_full.csv)

The tumor data used as input is from the treehouse dataset, available here: https://xenabrowser.net/datapages/?dataset=TumorCompendium_v10_PolyA_hugo_log2tpm_58581genes_2019-07-25.tsv&host=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443

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

## Running Celligner

This method was run on macOS High Sierra v10.13.6 using RStudio version 1.2.5033 and R version 3.6.2 (2019-12-12).

### Download/clone the repo

Download or clone this repo and open as a new project in RStudio.

### Install dependencies:

Use the install_packages.R script to install the necessary R packages to run the methods.

### Download the necessary data:

Data files should be stored in the directory passed to run_Celligner(). There are 4 files needed to run Celligner, by default the files are named:
- TCGA_mat.tsv
- CCLE_mat.csv
- Celligner_info.csv
- hgnc_complete_set_7.24.2018.txt

TCGA_mat.tsv is the matrix of log2(TPM+1) expression values for the tumor samples. The file used in the paper can be download from XenaBrowser: https://xenabrowser.net/datapages/?dataset=TumorCompendium_v10_PolyA_hugo_log2tpm_58581genes_2019-07-25.tsv&host=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443 (this file should be renamed TCGA_mat.tsv to use the default naming).

CCLE_mat.csv is the matrix of log2(TPM+1) expression values for the cell line samples. The file used in the paper is the DepMap Public 19Q4 'CCLE_expression_full.csv' file, which can be dowloaded from depmap.org: https://depmap.org/portal/download (this file should be renamed CCLE_mat.csv to use the default naming).

Celligner_info.csv is a matrix of sample info, which can be downloaded from the Figshare repo here: https://figshare.com/articles/Celligner_data/11965269. This file contains the sample names for the tumors and cell lines, as well as the information such as the cancer lineage, subtype, primary vs metastatic status, and tumor purity of the samples. These features are used for plotting the data, but not for the Celligner method itself. If this file is not provided than a default matrix will be created using the row names of TCGA_mat and CCLE_mat as the sampleIDs.

hgnc_complete_set_7.24.2018.txt is a table of gene ids, and is used to convert between HGNC gene IDs and Ensembl IDs. The version of this matrix used in the paper can be downloaded from the Figshare repo here: https://figshare.com/articles/Celligner_data/11965269. This file was downloaded from HGNC, and the latest version of the file can be downloaded from here: ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt (using this version will change the genes used). 

### Running the method:

The run_Celligner() method (found in Celligner_methods.R) combines all steps of the Celligner method. It loads the data, finds differentially expressed genes, runs contrastive principal components analysis, runs mutual nearest neighbors batch correction, and creates a Seurat object containing the aligned data and a 2D UMAP projection of the aligned data. 

Running the run_Celligner() method, with default parameters takes <1 hr, while running with the exact parameters used in the manuscript (setting the global parameter fastPCA to NULL) can take ~10 hours.


### Using the output:

run_Celligner() outputs a Seurat object (named comb_obj), which is used to package the data and run dimensionality reduction methods. To learn more about Seurat, see here: https://satijalab.org/seurat/. To access various information in the Seurat object use these commands
- To get the celligner aligned output: Seurat::GetAssayData(comb_obj)
- To get the metadata: comb_obj@meta.data
- To get the coordinates for the 2D UMAP projection: Seurat::Embeddings(comb_obj, reduction ='umap')
- To use Seurat to plot the results (colored by cancer lineage): Seurat::DimPlot(comb_obj, reduction = 'umap',  group.by = 'lineage', pt.size = 0.5) + ggplot2::theme(legend.position = 'none')


### Notes:

- By default the global parameter fast_cPCA is set to 10, which reduces the time for calculating the contrastive principal components (cPCs) by estimating a calculation of only the top contrastive principal components needed for the method (by the default the methods only uses the top four cPCs, so the parameter fast_cPCA can be set to any number >=4). To recreate the exact output from the paper set the global parameter fast_cPCA to NULL so that all the contrastive principal components (cPCs) will be calculated. This is quite slow.
- If using your own data (not the data recommended above) you will need to write your own load_data method. Later methods assume that the matrix TCGA_mat is sample x gene matrix, where the rows are the tumor sample IDs and the columns are Ensembl gene IDs, the matrix CCLE_mat is sample x gene matrix, where the rows are the cell line sample IDs and the columns are Ensembl gene IDs, and that the TCGA_ann and CCLE_ann matrices output by load_data have the columns sampleID, lineage, subtype, and `Primary/Metastasis` (these columns aren't used for the method, just for plotting the results - sampleID needs to match the row names of TCGA_mat and CCLE_mat, but the other columns can be set to NA without affecting the results). 

