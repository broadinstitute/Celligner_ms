library(here)
library(magrittr)
library(tidyverse)
source(here::here('src', 'analysis_helpers.R'))
source(here::here('src', 'global_params.R'))

# Supplementary Figure 6a
# re-ran run_Celligner, modifying the global param remove_cPCA_dims to:
# c(), c(1), c(1,2), c(1,2,3), c(1,2,3,4), c(1,2,3,4,5), and c(1,2,3,4,5,6)
# and observed the total number of mutual nearest neighbor pairs identified 
# to produce the input to this method
number_of_cPCs_plot <- function() {
  num_cPCs <- c(0, 1, 2, 3, 4, 5, 6)
  num_MNN <- c(8436,8879, 10074, 10220, 10415, 10248, 10273)
  num_cPCs_used <- cbind.data.frame(num_cPCs, num_MNN)
  num_cPCs_used$param_type <- ifelse(num_cPCs_used$num_cPCs==4, 'original parameter', 'modified parameter')
  
  cPCs_MNN <- ggplot2::ggplot(num_cPCs_used, ggplot2::aes(num_cPCs, num_MNN, shape=param_type, group= interaction('param_type'))) +
    ggplot2::geom_line(color='gray80') + 
    ggplot2::geom_point() +
    ggplot2::xlab('number of cPCs used') +
    ggplot2::ylab('number of MNN pairs') +
    ggplot2::scale_shape_manual(values=c(`modified parameter`=16, `original parameter`=8)) +
    ggplot2::guides(shape=ggplot2::guide_legend(title="")) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position='bottom', 
          text = ggplot2::element_text(size=8),
          axis.title = ggplot2::element_text(size=6), 
          axis.text = ggplot2::element_text(size=6),
          legend.margin =ggplot2::margin(0,0,0,0),
          legend.box.margin=ggplot2::margin(-10,-30,-10,-30))
  
  return(cPCs_MNN)
  
}

# Supplementary Figure 6b
# re-ran run_Celligner, modifying the global param mnn_k_CL to:
# 3, 4, 5, 7, and 10 to produce the input to this method
# Then ran cell_line_tumor_class and compared the proportion of agreement
# between the cell line tumor classes calculated at mnn_k_CL = 5 to the 
# cell line tumor classes calculated at the other parameter values.
modified_k_CL_plot <- function() {
  k_CL_agreement <- c(0.917534, 0.9439552, 1.0000000, 0.9495596, 0.9367494)
  k_CL_names <- c(3, 4, 5, 7, 10)
  k_CL_type <- c('modified parameter', 'modified parameter', 'original parameter', 'modified parameter', 'modified parameter' )
  k_CL_df <- cbind.data.frame(k_CL_agreement, k_CL_names, k_CL_type)
  k_CL_agreement <- ggplot2::ggplot(k_CL_df, 
                                    ggplot2::aes(k_CL_names, k_CL_agreement, shape= k_CL_type, group = interaction('k_CL_type'))) + 
    ggplot2::geom_line(col='gray80') + 
    ggplot2::geom_point() + 
    ggplot2::xlab("number of mutual nearest neighbors in cell line data") +
    ggplot2::ylab("agreement with original cell line tumor type classifications") +
    ggplot2::guides(shape=ggplot2::guide_legend(title="")) +
    ggplot2::scale_shape_manual(values=c(`modified parameter`=16, `original parameter`=8)) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position='bottom', 
          text = ggplot2::element_text(size=8),
          axis.title = ggplot2::element_text(size=6), 
          axis.text = ggplot2::element_text(size=6),
          legend.margin =ggplot2::margin(0,0,0,0), 
          legend.box.margin=ggplot2::margin(-10,-30,-10,-30))
  
  return(k_CL_agreement)
}

# Supplementary Figure 6c
# re-ran run_Celligner, modifying the global param mnn_k_tumor to:
# 40, 45, 50, 55, and 60 to produce the input to this method. 
# Then ran cell_line_tumor_class and compared the proportion of agreement
# between the cell line tumor classes calculated at mnn_k_tumor = 50 to the 
# cell line tumor classes calculated at the other parameter values.
modified_k_tumor_plot <- function() {
  k_tumor_agreement <- c(0.9543635, 0.9615693, 1.0000000, 0.9751801, 0.9639712)
  k_tumor_names <- c(40, 45, 50, 55, 60)
  k_tumor_type <- c('modified parameter', 'modified parameter', 'original parameter', 'modified parameter', 'modified parameter' )
  k_tumor_df <- cbind.data.frame(k_tumor_agreement, k_tumor_names, k_tumor_type)
  k_tumor_agreement <- ggplot2::ggplot(k_tumor_df, 
                                       ggplot2::aes(k_tumor_names, k_tumor_agreement, shape= k_tumor_type, group = interaction('k_tumor_type'))) +
    ggplot2::geom_line(color='gray80') +
    ggplot2::geom_point() +
    ggplot2::xlab("number of mutual nearest neighbors in tumor data") +
    ggplot2::ylab("agreement with original cell line tumor type classifications") +
    ggplot2::guides(shape=ggplot2::guide_legend(title="")) +
    ggplot2::scale_shape_manual(values=c(`modified parameter`=16, `original parameter`=8))+
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position='bottom',
          text = ggplot2::element_text(size=8),
          axis.title = ggplot2::element_text(size=6),
          axis.text = ggplot2::element_text(size=6), 
          legend.margin =ggplot2::margin(0,0,0,0), 
          legend.box.margin=ggplot2::margin(-10,-30,-10,-30))
  
  return(k_tumor_agreement)
  
}

# Supplementary Figure 6d
# re-ran run_Celligner, removing all tumors annotated as skin
# to produce the input to this method
remove_skin_tumors_plot <- function(rm_skin_tumors_alignment) {
  rm_skin_tumors <- ggplot2::ggplot(rm_skin_tumors_alignment, 
                                    ggplot2::aes(UMAP_1, UMAP_2, fill=lineage=='skin', color=type, size=type)) +
    ggplot2::geom_point(alpha=0.7, pch=21) +
    ggplot2::scale_color_manual(values=c(`CL`='black', `tumor`='white')) +
    ggplot2::scale_size_manual(values=c(`CL`=1.5, `tumor`=0.75)) +
    ggplot2::theme_classic() + 
    ggplot2::xlab("UMAP 1") + 
    ggplot2::ylab("UMAP 2") +
    ggplot2::theme(legend.position = 'right', 
          text=ggplot2::element_text(size=6),
          axis.text=ggplot2::element_text(size=6)) +
    ggplot2::labs(fill="skin")
  
  return(rm_skin_tumors)
}

# Supplementary Figure 6d
remove_skin_tumors_highlight <- function(rm_skin_tumors_alignment) {
  rm_skin_tumors_alignment$lineage <- gsub("_", " ", rm_skin_tumors_alignment$lineage)
  rm_skin_tumors_highlight <- ggplot2::ggplot(dplyr::filter(rm_skin_tumors_alignment,
                                            UMAP_1 < 10 & UMAP_1 > 6 & UMAP_2 < -4.8 & UMAP_2 > -8),
                                            ggplot2::aes(UMAP_1, UMAP_2, fill=lineage, color=type, size=type)) +
    ggplot2::geom_point(alpha=0.7, pch=21) +
    ggplot2::scale_color_manual(values=c(`CL`='black', `tumor`='white')) +
    ggplot2::scale_size_manual(values=c(`CL`=2, `tumor`=1.5)) +
    ggplot2::theme_classic() +
    ggplot2::xlab("UMAP 1") +
    ggplot2::ylab("UMAP 2") +
    ggplot2::theme(legend.position = 'right',
          text=ggplot2::element_text(size=6),
          axis.text=ggplot2::element_text(size=6)) +
    ggplot2::labs(fill="lineage") +
    ggplot2::guides(fill=ggplot2::guide_legend(ncol=2))
  
  return(rm_skin_tumors_highlight)
}

# Supplementary Figure 6e
# re-ran run_Celligner, removing all cell lines annotated as skin
# to produce the input to this method
remove_skin_CLs_plot <- function(rm_skin_CLs_alignment) {
  rm_skin_CLs <- ggplot2::ggplot(rm_skin_CLs_alignment, 
                                 ggplot2::aes(UMAP_1, UMAP_2, fill=lineage=='skin', color=type, size=type)) + 
    ggplot2::geom_point(alpha=0.7, pch=21) +
    ggplot2::scale_color_manual(values=c(`CL`='black', `tumor`='white')) +
    ggplot2::scale_size_manual(values=c(`CL`=1.5, `tumor`=0.75)) +
    ggplot2::theme_classic() + 
    ggplot2::xlab("UMAP 1") + 
    ggplot2::ylab("UMAP 2") +
    ggplot2::theme(legend.position = 'right',
          text=ggplot2::element_text(size=6),
          axis.text=ggplot2::element_text(size=6)) +
    ggplot2::labs(fill="skin")
  
  return(rm_skin_CLs)
}

# Supplementary Figure 6e
remove_skin_CLs_highlight <- function(rm_skin_CLs_alignment) {
  rm_skin_CLs_alignment$lineage <- gsub("_", " ", rm_skin_CLs_alignment$lineage)
  rm_skin_CLs_highlight <- ggplot2::ggplot(dplyr::filter(rm_skin_CLs_alignment, 
                                         UMAP_1 < 11 & UMAP_1 > 6 & UMAP_2 < 5 & UMAP_2 > 1),
                                         ggplot2::aes(UMAP_1, UMAP_2, fill=lineage, color=type, size=type)) +
    ggplot2::geom_point(alpha=0.7, pch=21) +
    ggplot2::scale_color_manual(values=c(`CL`='black', `tumor`='white')) +
    ggplot2::scale_size_manual(values=c(`CL`=2, `tumor`=1.5)) +
    ggplot2::theme_classic() + 
    ggplot2::xlab("UMAP 1") +
    ggplot2::ylab("UMAP 2") +
    ggplot2::theme(legend.position = 'right',
          text=ggplot2::element_text(size=6), 
          axis.text=ggplot2::element_text(size=6)) +
    ggplot2::labs(fill="lineage") +
    ggplot2::guides(fill=ggplot2::guide_legend(ncol=2))
  
  return(rm_skin_CLs_highlight)
}
