library(here)
library(magrittr)
library(tidyverse)
source(here::here('src', 'Celligner_helpers.R'))
source(here::here('src', 'analysis_helpers.R'))
source(here::here('src', 'global_params.R'))


# load tumor and cell line expression data from local directory
# TCGA_mat source: https://xenabrowser.net/datapages/?dataset=TumorCompendium_v10_PolyA_hugo_log2tpm_58581genes_2019-07-25.tsv&host=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
# CCLE_mat source: depmap.org DepMap Public 19Q4 CCLE_expression_full.csv
# hgnc_complete_set source: the version used in the paper can be found in figshare for the paper here: https://figshare.com/articles/Celligner_data/11965269
# the file was originally downloaded from here: ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt
# The annotation file, Celligner_info is available in figshare for the paper here: https://figshare.com/articles/Celligner_data/11965269
# This file contains the sampleIDs for the cell line and tumors, as well as additional info about samples (lineage, tumor purity, etc), 
# which is used for plotting the data, but not for the method itself
# if set to NULL then the annotation file is created using just the rownames from the TCGA and CCLE matrices
load_data <- function(data_dir, tumor_file = 'TCGA_mat.tsv', cell_line_file = 'CCLE_mat.csv', 
                      annotation_file = 'Celligner_info.csv', hgnc_file = "hgnc_complete_set_7.24.2018.txt") {
  hgnc.complete.set <- data.table::fread(file.path(data_dir, hgnc_file)) %>% as.data.frame()
  
  
  TCGA_mat <-  readr::read_tsv(file.path(data_dir, tumor_file)) %>% 
    as.data.frame() %>%
    tibble::column_to_rownames('Gene') %>%
    as.matrix() %>% 
    t()

  common_genes <- intersect(colnames(TCGA_mat), hgnc.complete.set$symbol)
  TCGA_mat <- TCGA_mat[,common_genes]
  hgnc.complete.set <- filter(hgnc.complete.set, symbol %in% common_genes)
  hgnc.complete.set <- hgnc.complete.set[!duplicated(hgnc.complete.set$symbol),]
  rownames(hgnc.complete.set) <- hgnc.complete.set$symbol
  hgnc.complete.set <- hgnc.complete.set[common_genes,]
  colnames(TCGA_mat) <- hgnc.complete.set$ensembl_gene_id
  
  CCLE_mat <-  readr::read_csv(file.path(data_dir, cell_line_file)) %>% 
    as.data.frame() %>%
    tibble::column_to_rownames('X1') %>%
    as.matrix()
  
  colnames(CCLE_mat) <- stringr::str_match(colnames(CCLE_mat), '\\((.+)\\)')[,2]

  if(is.null(annotation_file) | !file.exists(file.path(data_dir, annotation_file))) {
    ann <- data.frame(sampleID = c(rownames(TCGA_mat), rownames(CCLE_mat)),
                      lineage = NA,
                      subtype = NA,
                      type = c(rep('tumor', nrow(TCGA_mat)), rep('CL', nrow(CCLE_mat))))
    ann$`Primary/Metastasis` <- NA
  } else {
    ann <- data.table::fread(file.path(data_dir, annotation_file)) %>% as.data.frame()
    if('UMAP_1' %in% colnames(ann)) {
      ann <- ann %>% 
        dplyr::select(-UMAP_1)
    }
    if('UMAP_2' %in% colnames(ann)) {
      ann <- ann %>% 
        dplyr::select(-UMAP_2)
    }
    if('cluster' %in% colnames(ann)) {
      ann <- ann %>% 
        dplyr::select(-cluster)
    }
  }
  
  TCGA_ann <- dplyr::filter(ann, type=='tumor')
  CCLE_ann <- dplyr::filter(ann, type=='CL')
  
  func_genes <- dplyr::filter(hgnc.complete.set, !locus_group %in% c('non-coding RNA', 'pseudogene'))$ensembl_gene_id
  genes_used <- intersect(colnames(TCGA_mat), colnames(CCLE_mat))
  genes_used <- intersect(genes_used, func_genes)
  
  TCGA_mat <- TCGA_mat[,genes_used]
  CCLE_mat <- CCLE_mat[,genes_used]
  
  return(list(TCGA_mat = TCGA_mat, TCGA_ann = TCGA_ann, CCLE_mat = CCLE_mat, CCLE_ann = CCLE_ann))
}

# calculate the average gene expression and variance
# hgnc.complete.set source: the version used in the paper can be found in the figshare for the paper here: https://figshare.com/articles/Celligner_data/11965269
# the file was originally downloaded from here: ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt
calc_gene_stats <- function(dat, data_dir, hgnc_file = "hgnc_complete_set_7.24.2018.txt") {
  common_genes <- intersect(colnames(dat$TCGA_mat), colnames(dat$CCLE_mat))
  
  hgnc.complete.set <- data.table::fread(file.path(data_dir, hgnc_file)) %>% as.data.frame()
  hgnc.complete.set <- hgnc.complete.set %>% 
    dplyr::select(Gene = ensembl_gene_id, Symbol = symbol) %>%
    filter(Gene %in% common_genes)
  hgnc.complete.set <- hgnc.complete.set[!duplicated(hgnc.complete.set$Gene),]
  rownames(hgnc.complete.set) <- hgnc.complete.set$Gene
  hgnc.complete.set <- hgnc.complete.set[common_genes,]  
  
  gene_stats <- data.frame(
    Tumor_SD = apply(dat$TCGA_mat, 2, sd, na.rm=T),
    CCLE_SD = apply(dat$CCLE_mat, 2, sd, na.rm=T),
    Tumor_mean = colMeans(dat$TCGA_mat, na.rm=T),
    CCLE_mean = colMeans(dat$CCLE_mat, na.rm=T),
    Gene = common_genes,
    stringsAsFactors = F) %>% 
  dplyr::mutate(max_SD = pmax(Tumor_SD, CCLE_SD, na.rm=T)) #add avg and max SD per gene
  
  gene_stats <- left_join(hgnc.complete.set, gene_stats, by = "Gene")
  
  return(gene_stats)

  }

# create Seurat object of expression data and annotations
# and run dimensionality reduction
create_Seurat_object <- function(exp_mat, ann, type = NULL) {
  seu_obj <- Seurat::CreateSeuratObject(t(exp_mat),
                                         min.cells = 0,
                                         min.features = 0,
                                         meta.data = ann %>%
                                           magrittr::set_rownames(ann$sampleID))
  if(!is.null(type)) {
    seu_obj@meta.data$type <- type
  }
  # mean center the data, important for PCA
  seu_obj <- Seurat::ScaleData(seu_obj, features = rownames(Seurat::GetAssayData(seu_obj)), do.scale = F)
  
  seu_obj %<>% Seurat::RunPCA(assay='RNA',
                               features = rownames(Seurat::GetAssayData(seu_obj)),
                               npcs = global$n_PC_dims, verbose = F)
  
  seu_obj %<>% Seurat::RunUMAP(assay = 'RNA', dims = 1:global$n_PC_dims,
                                reduction = 'pca',
                                n.neighbors = global$umap_n_neighbors,
                                min.dist =  global$umap_min_dist,
                                metric = global$distance_metric, verbose=F)
  
  return(seu_obj)
}

# Cluster expression using Seurat clustering function
cluster_data <- function(seu_obj) {
  seu_obj <- Seurat::FindNeighbors(seu_obj, reduction = 'pca',
                                    dims = 1:global$n_PC_dims,
                                    k.param = 20, 
                                    force.recalc = TRUE,
                                    verbose = FALSE)
  
  seu_obj %<>% Seurat::FindClusters(reduction = 'pca',
                                     resolution = global$mod_clust_res)
  
  seu_obj@meta.data$cluster <- seu_obj@meta.data$seurat_clusters
  
  return(seu_obj)
  
  }

# Find genes which are differentially expressed between
# clusters in the data
find_differentially_expressed_genes <- function(seu_obj) {
  n_clusts <- nlevels(seu_obj@meta.data$seurat_clusters)
  if (n_clusts > 2) {
    cur_DE_genes <- run_lm_stats_limma_group(
      t(Seurat::GetAssayData(seu_obj, assay='RNA', slot='scale.data')),
      seu_obj@meta.data %>% dplyr::select(seurat_clusters),
      limma_trend = TRUE) %>%
      dplyr::select(Gene, gene_stat = F_stat)
  } else if (n_clusts == 2) {
    cur_DE_genes <- run_lm_stats_limma(t(Seurat::GetAssayData(seu_obj, assay='RNA', slot='scale.data')),
                                               seu_obj@meta.data$cluster,
                                               limma_trend = TRUE) %>%
      dplyr::mutate(gene_stat = abs(t_stat)) %>%
      dplyr::select(Gene, gene_stat)
  } else {
    cur_DE_genes <- data.frame(Gene = colnames(seu_obj), gene_stat = NA)
  }
  
  return(cur_DE_genes)
  
}


# run contrastive principal components analysis
# set pc_dims to a value >= 4 to run fast cPCA by just calculating the top contrastive principle components 
run_cPCA <- function(TCGA_obj, CCLE_obj, pc_dims = NULL) {
  cov_diff_eig <- run_cPCA_analysis(t(Seurat::GetAssayData(TCGA_obj, assay='RNA', slot='scale.data')), 
                                    t(Seurat::GetAssayData(CCLE_obj, assay='RNA', slot='scale.data')), 
                                    TCGA_obj@meta.data, CCLE_obj@meta.data, pc_dims=pc_dims)
 return(cov_diff_eig) 
}

# run mutual nearest neighbors batch correction
run_MNN <- function(CCLE_cor, TCGA_cor,  k1 = global$mnn_k_tumor, k2 = global$mnn_k_CL, ndist = global$mnn_ndist, 
                    subset_genes) {
  mnn_res <- modified_mnnCorrect(CCLE_cor, TCGA_cor, k1 = k1, k2 = k2, ndist = ndist, 
                             subset_genes = subset_genes)
  
  return(mnn_res)
}

# calculate the correlation between tumors and cell lines in the Celligner_aligned_data
calc_tumor_CL_cor <- function(Celligner_aligned_data, Celligner_info) {
  tumors_samples <- dplyr::filter(Celligner_info, type=='tumor')$sampleID
  cl_samples <- dplyr::filter(Celligner_info, type=='CL')$sampleID
  tumor_CL_cor <- cor(t(Celligner_aligned_data[tumor_samples,]), t(Celligner_aligned_data[cl_samples,]),
                      use='pairwise')
  
  
  return(tumor_CL_cor)
}

# run all Celligner methods
# outputs a Seurat object of the Celligner-aligned data, 2D UMAP embedding, and combined
# tumor and cell line annotations
run_Celligner <- function(data_dir) {
  dat <- load_data(data_dir)
  gene_stats <- calc_gene_stats(dat, data_dir)
  
  comb_ann <- rbind(
    dat$TCGA_ann %>% dplyr::select(sampleID, lineage, subtype, `Primary/Metastasis`) %>%
      dplyr::mutate(type = 'tumor'),
    dat$CCLE_ann %>% dplyr::select(sampleID, lineage, subtype, `Primary/Metastasis`) %>%
      dplyr::mutate(type = 'CL')
  )
  
  TCGA_obj <- create_Seurat_object(dat$TCGA_mat, dat$TCGA_ann, type='tumor')
  CCLE_obj <- create_Seurat_object(dat$CCLE_mat, dat$CCLE_ann, type='CL')
  
  TCGA_obj <- cluster_data(TCGA_obj)
  CCLE_obj <- cluster_data(CCLE_obj)
  
  tumor_DE_genes <- find_differentially_expressed_genes(TCGA_obj)
  CL_DE_genes <- find_differentially_expressed_genes(CCLE_obj)
  
  DE_genes <- full_join(tumor_DE_genes, CL_DE_genes, by = 'Gene', suffix = c('_tumor', '_CL')) %>%
    mutate(
      tumor_rank = dplyr::dense_rank(-gene_stat_tumor),
      CL_rank = dplyr::dense_rank(-gene_stat_CL),
      best_rank = pmin(tumor_rank, CL_rank, na.rm=T)) %>%
    dplyr::left_join(gene_stats, by = 'Gene')
  
  # take genes that are ranked in the top 1000 from either dataset, used for finding mutual nearest neighbors
  DE_gene_set <- DE_genes %>%
    dplyr::filter(best_rank < global$top_DE_genes_per) %>%
    .[['Gene']]
  
  
  cov_diff_eig <- run_cPCA(TCGA_obj, CCLE_obj, global$fast_cPCA)
  
  if(is.null(global$fast_cPCA)) {
    cur_vecs <- cov_diff_eig$vectors[, global$remove_cPCA_dims, drop = FALSE]
  } else {
    cur_vecs <- cov_diff_eig$rotation[, global$remove_cPCA_dims, drop = FALSE]
  }
  
  rownames(cur_vecs) <- colnames(dat$TCGA_mat)
  TCGA_cor <- resid(lm(t(dat$TCGA_mat) ~ 0 + cur_vecs)) %>% t()
  CCLE_cor <- resid(lm(t(dat$CCLE_mat) ~ 0 + cur_vecs)) %>% t()
  
  mnn_res <- run_MNN(CCLE_cor, TCGA_cor,  k1 = global$mnn_k_tumor, k2 = global$mnn_k_CL, ndist = global$mnn_ndist, 
                      subset_genes = DE_gene_set)
  
  combined_mat <- rbind(mnn_res$corrected, CCLE_cor)
  
  comb_obj <- create_Seurat_object(combined_mat, comb_ann)
  comb_obj <- cluster_data(comb_obj)
  
  return(comb_obj)
}







