library(here)
library(magrittr)
library(tidyverse)
source(here::here('src', 'analysis_helpers.R'))
source(here::here('src', 'global_params.R'))

# Supplementary Figure 3a
# re-ran run_Celligner without running the run_MNN step
# to produce the input to this method
plot_cPCA_only_output <- function(cPCA_only) {
  cPCA_only_plot <- ggplot2::ggplot(cPCA_only, ggplot2::aes(UMAP_1, UMAP_2, fill=lineage, size=type, color =type)) +
    ggplot2::geom_point(alpha=0.6, pch=21) +
    ggplot2::scale_color_manual(values=c(`CL`='black', `tumor`='white')) +
    ggplot2::scale_size_manual(values=c(`CL`=1, `tumor`=0.75)) +
    ggplot2::theme_classic() +
    ggplot2::xlab("UMAP 1") +
    ggplot2::ylab("UMAP 2") +
    ggplot2::theme(legend.position = 'bottom', 
          text=ggplot2::element_text(size=8), 
          legend.margin =ggplot2::margin(0,0,0,0)) +
    ggplot2::guides(fill=FALSE, color=FALSE, size=FALSE) +
    ggplot2::scale_fill_manual(values=tissue_colors)
  
  return(cPCA_only_plot)
  
}

# Supplementary Figure 3b
# mnn_vect is the mnn correction vectors created while running modified_mnnCorrect 
MNN_heatmap <- function(mnn_vect, gene_stats, TCGA_ann, filename) {
  # heatmap of correction vectors
  avg_vect <- data.frame(diff_mean = colMeans(mnn_vect), diff_sd = apply(mnn_vect, 2, sd, na.rm=T),
                         Gene = colnames(mnn_vect)) %>% 
    dplyr::left_join(gene_stats, by = 'Gene')
  
  
  top_mnn_genes <- avg_vect %>% 
    dplyr::arrange(dplyr::desc(abs(diff_mean))) %>% 
    head(50)
  
  tissue_order <- TCGA_ann$lineage
  names(tissue_order) <- TCGA_ann$sampleID
  tissue_order <- tissue_order[order(tissue_order)]
  tissue_annot <- as.data.frame(tissue_order)
  rownames(tissue_annot) <- names(tissue_order)
  colnames(tissue_annot) <- 'lineage'
  pheatmap::pheatmap(mnn_vect[names(tissue_order), top_mnn_genes$Gene],
           show_rownames = FALSE, 
           labels_col = top_mnn_genes$Symbol,
           cluster_rows = F,
           na_col= heatmap_params$na_color, 
           fontsize = heatmap_params$title_font_size,
           fontsize_col = 4,
           fontsize_row = heatmap_params$row_font_size, 
           annotation_row = tissue_annot,
           annotation_legend = F,
           width = 4, 
           height = 3.2, 
           fontface = heatmap_params$font_face,
           angle_col=90, 
           filename = filename)
}

# Supplementary Figure 3c
# gsc data used is MSigDB v6.2 downloaded from www.gsea-msigdb.org/gsea/msigdb
# mnn_vect is the mnn correction vectors created while running modified_mnnCorrect 
MNN_GSEA_table <- function(mnn_vect, gene_stats, filename) {
  avg_vect <- data.frame(diff_mean = colMeans(mnn_vect), diff_sd = apply(mnn_vect, 2, sd, na.rm=T),
                         Gene = colnames(mnn_vect)) %>% 
    dplyr::left_join(gene_stats, by = 'Gene')
  gene_stat <- avg_vect$diff_mean %>% set_names(avg_vect$Symbol)
  
  MNN_GSEA <- run_fGSEA(gsc = gsc_data$GO_biological_process,
                              gene_stat = gene_stat,
                              nperm = 1e5,
                              perm_type = 'gene') %>% dplyr::arrange(dplyr::desc(abs(NES))) %>% 
    dplyr::select(-leadingEdge)
  
  MNN_GSEA_data <- rbind.data.frame(MNN_GSEA %>% dplyr::arrange(NES) %>% head(5), MNN_GSEA %>% dplyr::arrange(dplyr::desc(NES)) %>% head(5))
  MNN_GSEA_data$pathway <-sub("^([^_]*_[^_]*_[^_]*_[^_]*_[^_]*).*", "\\1", MNN_GSEA_data$pathway)
  MNN_GSEA_data$pathway <- factor(MNN_GSEA_data$pathway, levels = MNN_GSEA_data$pathway) 
  MNN_GSEA_data$`higher variance in` <- c(rep('tumors', 5), rep('cell lines',5))
  
  MNN_GSEA_data$`adjusted pval` <- signif(MNN_GSEA_data$padj, 3)
  MNN_GSEA_data$NES <- signif(MNN_GSEA_data$NES, 3)
  
  png(filename, height=3, width=4.5, units='in', res=200)
  grid.table(MNN_GSEA_data[,c('pathway', 'adjusted pval', 'NES', 'higher variance in')], theme=table_theme, rows=NULL)
  dev.off()
  
}

# Supplementary Figure 3d
# re-ran run_Celligner without running the cPCA step
# to produce the input to this method
plot_MNN_only_output <- function(MNN_only) {
  
  MNN_only_plot <- ggplot2::ggplot(MNN_only, 
                                   ggplot2::aes(UMAP_1, UMAP_2, fill=lineage, size=type, color = type)) +
    ggplot2::geom_point(alpha=0.6, pch=21) +
    ggplot2::scale_color_manual(values=c(`CL`='black', `tumor`='white')) + 
    ggplot2::scale_size_manual(values=c(`CL`=1, `tumor`=0.75)) +
    ggplot2::theme_classic() + 
    ggplot2::xlab("UMAP 1") + 
    ggplot2::ylab("UMAP 2") + 
    ggplot2::theme(legend.position = 'none',
                   text=ggplot2::element_text(size=8)) +
    ggplot2::scale_fill_manual(values=tissue_colors)
  
  return(MNN_only_plot)
}

# Supplementary Figure 3e
# mnn_pairs is part of the output of run_MNN
# re-ran run_Celligner without running the cPCA step
# to produce the input to this method
num_MNN_purity <- function(mnn_pairs, mnn_only_pairs, TCGA_ann) {
  # effect of doing cPCA correction on number of mnn
  df  <- TCGA_ann %>% 
    dplyr::mutate(`MNN with cPCA` = sampleID %in% mnn_pairs$targ_ID,
           `MNN, no cPCA` = sampleID %in% mnn_only_pairs$targ_ID)
  
  sample_purity <- rbind(df %>% dplyr::filter(!is.na(purity) & `MNN, no cPCA`==T) %>%
                           dplyr::select(purity) %>% dplyr::mutate(variable='MNN, no cPCA'), 
                         df %>% dplyr::filter(!is.na(purity) & `MNN with cPCA`==T) %>%
                           dplyr::select(purity) %>%
                           dplyr::mutate(variable='MNN with cPCA'))
  
    
  cPCA_num_mnn_plot <- ggpubr::ggdensity(sample_purity, x='purity', y='..count..',
                                         fill='variable', 
                                         xlab="(ABSOLUTE) estimated tumor purity", 
                                           ylab = "number of tumor mutual nearest neighbors") +
    ggplot2::theme_classic() +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
          axis.title.x = ggplot2::element_text(size=8),
          axis.title.y = ggplot2::element_text(size=8, angle=90),
          text = ggplot2::element_text(size=8),
          legend.position = 'bottom',
          legend.margin =ggplot2::margin(0,0,0,0)) +
    ggplot2::guides(fill=ggplot2::guide_legend(title=""))
  
  return(cPCA_num_mnn_plot)
    
}







