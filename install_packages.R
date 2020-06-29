options(repos = c("https://cran.cnr.berkeley.edu"))

cran_packages <- c('here', 'tidyverse', 'reshape2', 'plyr', 'data.table', 'Seurat',
                   'pheatmap', 'pdist', 'gridExtra', 'ggpubr', 'grDevices', 'RColorBrewer',
                   'FNN', 'ggrepel', 'ggridges', 'irlba', 'viridis')
new_cran_packages <- cran_packages[!(cran_packages %in% installed.packages()[,"Package"])]
if(length(new_cran_packages)) install.packages(new_cran_packages)

bioconductor_packages <- c('limma', 'edgeR', 'batchelor', 'BiocParallel', 'sva', 'GSEABase',
                           'piano', 'fgsea')
new_bioconductor_packages <- bioconductor_packages[!(bioconductor_packages %in% installed.packages()[,"Package"])]
if(length(new_bioconductor_packages)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install(new_bioconductor_packages)
}

