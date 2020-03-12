library(here)
library(magrittr)
library(tidyverse)
source(here::here('src', 'analysis_helpers.R'))
source(here::here('src', 'global_params.R'))

# used as input for Supplementary Figure 1a
run_COMBAT_correction <- function(TCGA_mat, CCLE_mat, alignment) {
  library(sva)
  dat <- rbind(TCGA_mat, as.matrix(CCLE_mat))
  cov_df <- alignment[match(alignment$sampleID_CCLE_Name, rownames(dat)), c('type', 'lineage', 'sampleID_CCLE_Name')] %>% 
    dplyr::filter(!is.na(lineage), sampleID_CCLE_Name %in% rownames(dat))
  dat <- dat[cov_df$sampleID_CCLE_Name, , drop = FALSE]
  mod = model.matrix(~1, data=cov_df)
  
  zero_var_genes <- which(apply(t(dat), 1, var)==0)
  combat_output <- sva::ComBat(t(dat[,-zero_var_genes]), batch = cov_df$type, mod = mod)
  
  return(combat_output)
}

# Supplementary Figure 1a
plot_COMBAT_corrected <- function(combat_output, alignment) {
  
  combat_obj <-  Seurat::CreateSeuratObject(combat_output,
                                            min.cells = 0,
                                            min.features = 0,
                                            meta.data = alignment %>%
                                              magrittr::set_rownames(alignment$sampleID))
  combat_obj <- Seurat::ScaleData(combat_obj, 
                                  features = rownames(Seurat::GetAssayData(combat_obj)),
                                  do.scale = F)
  
  # PCA
  combat_obj %<>% Seurat::RunPCA(assay='RNA',
                                 features = rownames(Seurat::GetAssayData(combat_obj)),
                                 npcs = global$n_PC_dims, verbose = F)
  
  combat_obj %<>% Seurat::RunUMAP(assay = 'RNA', dims = 1:global$n_PC_dims,
                                  reduction = 'pca',
                                  n.neighbors = global$umap_n_neighbors,
                                  min.dist = global$umap_min_dist,
                                  metric = global$distance_metric, verbose=F)
  
  combat_res <- Seurat::Embeddings(combat_obj, reduction = global$reduction.use) %>%
    as.data.frame() %>%
    magrittr::set_colnames(c('UMAP_1', 'UMAP_2')) %>%
    cbind.data.frame(alignment[,c('sampleID_CCLE_Name','lineage', 'subtype', 'type')])
  
  combat_plot <- ggplot2::ggplot(combat_res, 
                                 ggplot2::aes(UMAP_1, UMAP_2, fill=lineage, size=type)) +
    ggplot2::geom_point(alpha=0.6, pch=21, color='white') +
    ggplot2::geom_point(data = dplyr::filter(alignment, type=='CL'), color='#333333', pch=21, alpha=0.6) +
    ggplot2::scale_size_manual(values=c(`CL`=1.5, `tumor`=0.75)) + 
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = 'none',
          text=ggplot2::element_text(size=6)) +
    ggplot2::scale_fill_manual(values=tissue_colors)
  
  return(combat_plot)
}

# Supplementary Figure 1b
plot_TCGA_clusters <- function(TCGA_ann) {
  TCGA_cluster_avgs <- TCGA_ann %>% 
    dplyr::group_by(cluster) %>% 
    dplyr::summarise(UMAP_1 = median(UMAP_1, na.rm=T),
                     UMAP_2 = median(UMAP_2, na.rm=T))
  TCGA_cluster_lab_aes <- ggplot2::geom_text(data = TCGA_cluster_avgs,
                                             mapping = ggplot2::aes(x = UMAP_1, y = UMAP_2, label = cluster),
                                             size = 3, 
                                             color="#000000")
  
  tcga_orig_cluster_plot <- ggplot2::ggplot(TCGA_ann, ggplot2::aes(UMAP_1, UMAP_2, color=as.factor(cluster))) +
    ggplot2::geom_point(alpha=0.6) +
    ggplot2::theme_classic() +
    ggplot2::xlab('UMAP 1') +
    ggplot2::ylab("UMAP 2") +
    ggplot2::theme(legend.position = 'none', 
          text=ggplot2::element_text(size=6)) +
    TCGA_cluster_lab_aes
  
  return(tcga_orig_cluster_plot)
}

# Supplementary Figure 1c
plot_CCLE_clusters <- function(CCLE_ann) {
  CCLE_cluster_avgs <- CCLE_ann %>% 
    dplyr::group_by(cluster) %>% 
    dplyr::summarise(UMAP_1 = median(UMAP_1, na.rm=T),
                     UMAP_2 = median(UMAP_2, na.rm=T))
  
  CCLE_cluster_lab_aes <- ggplot2::geom_text(data = CCLE_cluster_avgs,
                                             mapping = ggplot2::aes(x = UMAP_1, y = UMAP_2, label = cluster),
                                             size = 3,
                                             color="#000000")
  
  CCLE_orig_cluster_plot <- ggplot2::ggplot(CCLE_ann, ggplot2::aes(UMAP_1, UMAP_2, color=as.factor(cluster))) +
    ggplot2::geom_point(alpha=0.6) +
    ggplot2::theme_classic() +
    ggplot2::xlab('UMAP 1') +
    ggplot2::ylab("UMAP 2") +
    ggplot2::theme(legend.position = 'none',
          text=ggplot2::element_text(size=6)) +
    CCLE_cluster_lab_aes

  return(CCLE_orig_cluster_plot)
}

# Supplementary Figure 1d
plot_aligned_clusters <- function(alignment) {
  cluster_avgs <- alignment %>% 
    dplyr::group_by(cluster) %>% 
    dplyr::summarise(UMAP_1 = median(UMAP_1, na.rm=T),
                     UMAP_2 = median(UMAP_2, na.rm=T))
  
  cluster_lab_aes <- ggplot2::geom_text(data = cluster_avgs, 
                                        mapping = ggplot2::aes(x = UMAP_1, y = UMAP_2, label = cluster),
                                        size = 3, 
                                        color="#000000")
  
  aligned_cluster_plot <- ggplot2::ggplot(alignment, 
                                          ggplot2::aes(UMAP_1, UMAP_2, 
                                     fill=as.factor(cluster),
                                     size=type, color=type)) +
    ggplot2::geom_point(alpha=0.6, pch=21) +
    ggplot2::scale_color_manual(values=c(`CL`='black', `tumor`='white')) +
    ggplot2::xlab('UMAP 1') +
    ggplot2::ylab("UMAP 2") +
    ggplot2::scale_size_manual(values=c(`CL`=1, `tumor`=0.75)) +
    ggplot2::theme_classic() + 
    ggplot2::theme(legend.position = 'none',
          text=ggplot2::element_text(size=6)) +
    cluster_lab_aes
  
  return(aligned_cluster_plot)
}
