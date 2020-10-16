library(here)
library(magrittr)
library(tidyverse)
source(here::here('src', 'analysis_helpers.R'))
source(here::here('src', 'global_params.R'))

# Supplementary Figure 5a
CL_SUSA_plot <- function(CCLE_ann) {
  CCLE_ann$name <- gsub("_.*", "", CCLE_ann$sampleID_CCLE_Name)
  CCLE_ann$tissue_name <-  gsub("_", " ", CCLE_ann$lineage)
  Cl_only_testis <- ggplot2::ggplot(CCLE_ann, ggplot2::aes(UMAP_1, UMAP_2)) + 
    ggplot2::geom_point(size=0.5, alpha=0.3, ggplot2::aes(color = lineage)) +
    ggplot2::geom_point(data=dplyr::filter(CCLE_ann, sampleID_CCLE_Name %in% c("SUSA_TESTIS")),
               ggplot2::aes(UMAP_1, UMAP_2, fill=lineage), pch=21, color='black', size = 3) +  
    ggrepel::geom_text_repel(data=dplyr::filter(CCLE_ann, sampleID_CCLE_Name %in% c("SUSA_TESTIS")),
                             ggplot2::aes(UMAP_1, UMAP_2, label=name), color='black', size=3) + 
    ggplot2::xlab("UMAP 1") + 
    ggplot2::ylab("UMAP 2") +
    ggplot2::theme_classic() + 
    ggplot2::theme(legend.position = 'none',
          text=ggplot2::element_text(size=6),
          axis.text=ggplot2::element_text(size=6),
          plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::ggtitle("Cell line only embedding")

  return(Cl_only_testis)
}

# Supplementary Figure 5b
CL_SUSA_highlight <- function(CCLE_ann) {
  CCLE_ann$name <- gsub("_.*", "", CCLE_ann$sampleID_CCLE_Name)
  CCLE_ann$tissue_name <-  gsub("_", " ", CCLE_ann$lineage)

  Cl_only_testis_highlight <- ggplot2::ggplot(dplyr::filter(CCLE_ann, UMAP_1 < 32.5 & UMAP_1 > 29 & UMAP_2 < 26.5 & UMAP_2 > 25),
                                              ggplot2::aes(UMAP_1, UMAP_2)) + 
    ggplot2::geom_point(size=1.5, alpha=0.7, ggplot2::aes(color = tissue_name)) +
    ggplot2::geom_point(data=dplyr::filter(CCLE_ann, sampleID_CCLE_Name %in% c("SUSA_TESTIS")),
                        ggplot2::aes(UMAP_1, UMAP_2, fill=tissue_name), pch=21, color='black', size = 3) + 
    ggrepel::geom_text_repel(data=dplyr::filter(CCLE_ann, sampleID_CCLE_Name %in% c("SUSA_TESTIS")),
                             ggplot2::aes(UMAP_1, UMAP_2, label=name), color='black', size=3) +
    ggplot2::xlab("UMAP 1") + 
    ggplot2::ylab("UMAP 2") +
    ggplot2::theme_classic() + 
    ggplot2::theme(legend.position = 'right', 
          text=ggplot2::element_text(size=6),
          axis.text=ggplot2::element_text(size=6),
          plot.title = ggplot2::element_text(hjust = 0.5),
          legend.margin=ggplot2::margin(0,0.3,0,0)) + 
    ggplot2::guides(fill=FALSE, color=ggplot2::guide_legend(title="lineage", ncol=2))
  
  return(Cl_only_testis_highlight)
}

# Supplementary Figure 5c
Celligner_alignment_testis <- function(alignment) {
  alignment$name <- gsub("_.*", "", alignment$sampleID_CCLE_Name)
  
  alignment_testis <- ggplot2::ggplot(alignment, ggplot2::aes(UMAP_1, UMAP_2, shape = type)) + 
    ggplot2::geom_point(size=0.5, alpha=0.5, color = 'gray80') +
    ggplot2::geom_point(data = dplyr::filter(alignment, lineage %in% c('germ_cell') & type=='tumor'),
               size=1, alpha=0.5, ggplot2::aes(color = lineage))  +
    ggplot2::geom_point(data=dplyr::filter(alignment, sampleID_CCLE_Name %in% c("SUSA_TESTIS")),
                        ggplot2::aes(UMAP_1, UMAP_2, fill=lineage), color='black', size = 3, show.legend = F) + 
    ggplot2::scale_shape_manual(values=c(`CL`=21, `tumor`=16)) +
    ggrepel::geom_text_repel(data=dplyr::filter(alignment, sampleID_CCLE_Name %in% c("SUSA_TESTIS")),
                    aes(UMAP_1, UMAP_2, label=name), color='black', size=3) +
    ggplot2::xlab("UMAP 1") + 
    ggplot2::ylab("UMAP 2") +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = 'none', 
          text=ggplot2::element_text(size=6),
          axis.text=ggplot2::element_text(size=6),
          plot.title = ggplot2::element_text(hjust = 0.5)) + 
    ggplot2::ggtitle("Cell line and tumor Celligner-alignment")
  
  return(alignment_testis)
}

# Supplementary Figure 5d
Celligner_alignment_testis_highlight <- function(alignment) {
  alignment$name <- gsub("_.*", "", alignment$sampleID_CCLE_Name)
  
  alignment_testis_subtype <- ggplot2::ggplot(dplyr::filter(alignment, lineage == 'germ_cell' & UMAP_1 <7),
                                              ggplot2::aes(UMAP_1, UMAP_2, fill=subtype)) + 
    ggplot2::geom_point(alpha=0.7, pch=21, color='white', size=1.5) +
    ggplot2::geom_point(data=dplyr::filter(alignment, sampleID_CCLE_Name %in% c("SUSA_TESTIS")), 
                        ggplot2::aes(UMAP_1, UMAP_2, fill=lineage), color='black', size = 3, show.legend = F, pch=21) + 
    ggrepel::geom_text_repel(data=dplyr::filter(alignment, sampleID_CCLE_Name %in% c("SUSA_TESTIS")),
                             ggplot2::aes(UMAP_1, UMAP_2, label=name), color='black', size=3) +
    ggplot2::xlab("UMAP 1") + 
    ggplot2::ylab("UMAP 2") + 
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = 'bottom',
          text=ggplot2::element_text(size=8),
          axis.text=ggplot2::element_text(size=6), 
          axis.title=ggplot2::element_text(size=6),
          legend.margin =ggplot2::margin(0,0,0,0)) +  
    ggplot2::scale_fill_manual(breaks = c("non-seminoma", "seminoma"), 
                      values = c("#F8766D", "#7CAE00", "#00BFC4", "#777777")) 

     
     return(alignment_testis_subtype)                                                                                                                                                                             

}                                    


# Supplementary Figure 5e
plot_misanntations_melanoma <- function(alignment) {
  melanoma_cluster <- dplyr::filter(alignment, lineage=='skin')$cluster %>%
    table() %>%
    as.data.frame() %>%
    dplyr::arrange(dplyr::desc(Freq)) %>%
    head(2) %>% 
    .[['.']] %>%
    as.character()
  
  melanoma_cluster_data <- dplyr::filter(alignment, cluster %in% melanoma_cluster & UMAP_2 < -8)
  melanoma_cluster_data$name <- gsub("_.*", "", melanoma_cluster_data$sampleID_CCLE_Name)
  melanoma_cluster_data$tissue_name <- gsub("_", " ", melanoma_cluster_data$lineage)
  
  misannotation_plot <- ggplot2::ggplot(melanoma_cluster_data, ggplot2::aes(UMAP_1, UMAP_2, size=type, color=type)) +
    ggplot2::geom_point(alpha=0.5, pch=21, fill='gray70') +
    ggplot2::geom_point(data = dplyr::filter(melanoma_cluster_data, type=='CL' & lineage != 'skin'),
                        ggplot2::aes(fill=tissue_name), color='black', alpha=0.7, pch=21) +
    ggrepel::geom_label_repel(data = dplyr::filter(melanoma_cluster_data, type=='CL' & lineage != 'skin'),
                              ggplot2::aes(label=name)) +
    ggplot2::scale_color_manual(values=c(`CL`='black', `tumor`='white'), guide='none') +
    ggplot2::scale_size_manual(values=c(`CL`=2.5, `tumor`=1.5), guide='none') +
    ggplot2::xlab("UMAP 1") + 
    ggplot2::ylab("UMAP 2") +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = 'bottom',
          text=ggplot2::element_text(size=7),
          legend.margin =ggplot2::margin(0,0,0,0),
          legend.box.margin=ggplot2::margin(-10,-30,-10,-10),
          axis.text =ggplot2::element_text(size = 6)) +
    ggplot2::guides(size=FALSE, fill=ggplot2::guide_legend(nrow=2,byrow=TRUE,title='lineage')) 
  
  
  return(misannotation_plot)
  
}

# Supplementary Figure 5f
ALL_subtype_expression <- function(TCGA_mat, alignment, filename) {
  t_cell_cluster <- c(43)
  b_cell_cluster <- c(16,26,59)
  
  ALL_tumors <- dplyr::filter(alignment, subtype=='acute lymphoblastic leukemia' &
                         cluster %in% c(t_cell_cluster, b_cell_cluster))$sampleID
  
  ALL_marker_exp <- TCGA_mat[ALL_tumors, ALL_marker_genes$ensembl]
  colnames(ALL_marker_exp) <- ALL_marker_genes$gene
  
  row_ann <- ifelse(ALL_tumors %in% dplyr::filter(alignment, cluster %in% b_cell_cluster)$sampleID,
                    "B-cell", "T-cell") %>% as.data.frame()
  colnames(row_ann) <- 'cell line cluster'
  rownames(row_ann) <- ALL_tumors
  
  col_ann <- ALL_marker_genes$type %>% as.data.frame()
  colnames(col_ann) <- 'marker gene'
  rownames(col_ann) <- ALL_marker_genes$gene
  
  pheatmap::pheatmap(ALL_marker_exp,
           show_rownames = F,
           annotation_row = row_ann,
           annotation_col = col_ann,
           filename = filename, 
           main="ALL subtype marker expression",
           fontsize = heatmap_params$title_font_size, 
           fontsize_col = heatmap_params$column_font_size,
           fontsize_row = 4,
           width = 3.8,
           height = 3.4,
           fontface = heatmap_params$font_face,
           angle_col=90)

}





