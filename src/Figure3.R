library(here)
library(magrittr)
library(tidyverse)
source(here::here('src', 'analysis_helpers.R'))
source(here::here('src', 'global_params.R'))


# used as input for 3a
# tumor_CL_cor included in the figshare
get_cell_line_tumor_class <- function(tumor_CL_cor, alignment) {
  cl_tumor_classes <- apply(tumor_CL_cor, 2, function(x) cell_line_tumor_class(x, tumor_CL_cor, alignment)) %>% 
    as.character()
  names(cl_tumor_classes) <- colnames(tumor_CL_cor)
  
  return(cl_tumor_classes)
}

# figure 3a
cell_line_tumor_class_plot <- function(cl_tumor_classes, alignment, tumor_CL_cor, filename) {
  cl_tissue_type <- dplyr::filter(alignment, type=='CL')
 #cl_tissue_type[grep('rhabdomyosarcoma', cl_tissue_type$subtype),'tissue'] <- 'rhabdomyosarcoma'
  rownames(cl_tissue_type) <- cl_tissue_type$sampleID
  classification_freq <- table(cl_tumor_classes, cl_tissue_type[colnames(tumor_CL_cor),'lineage']) %>% as.data.frame()
  classification_freq <- reshape2::dcast(classification_freq, cl_tumor_classes ~ Var2, value.var = 'Freq') %>%
    tibble::column_to_rownames('cl_tumor_classes')
  print(setdiff(intersect(unique(dplyr::filter(alignment, type=='CL')$lineage),
                          unique(dplyr::filter(alignment, type=='tumor')$lineage)), 
                rownames(classification_freq)))
  
  esophagus_tumor <- rep(0, ncol(classification_freq))
  classification_freq <- rbind(classification_freq,`esophagus`= esophagus_tumor) 
  common_types <- intersect(rownames(classification_freq), colnames(classification_freq))
  
  prop_agree <- sum(diag(as.matrix(classification_freq[common_types, common_types])))/sum(as.matrix(classification_freq[common_types, common_types]))
  
  for(i in 1:ncol(classification_freq)) {
    classification_freq[,i] <- classification_freq[,i]/sum(classification_freq[,i])
  }
  
  
  agreement <- diag(as.matrix(classification_freq[common_types, common_types]))
  agreement_CL <- agreement
  names(agreement_CL) <- common_types
  agreement_tumor <- agreement
  names(agreement_tumor) <- common_types
  
  agreement_CL <- base::sort(agreement_CL, decreasing=T)
  agreement_tumor <- base::sort(agreement_tumor, decreasing=T)
  
  
  
  
  classification_freq <- classification_freq[names(agreement_tumor), names(agreement_CL)]
  rownames(classification_freq) <- gsub("_", " ", rownames(classification_freq))
  colnames(classification_freq) <- gsub("_", " ", colnames(classification_freq))
  
  pheatmap::pheatmap(classification_freq, 
           border_color = heatmap_params$square_border_color, 
           na_col= heatmap_params$na_color, 
           cluster_rows = F, 
           cluster_cols = F, 
           main="", 
           fontsize = heatmap_params$title_font_size,
           fontsize_col = heatmap_params$column_font_size,
           fontsize_row = heatmap_params$row_font_size,
           width = 3.5,
           height = 3,
           fontface = heatmap_params$font_face,
           angle_col=90, 
           filename = filename,
           color= heatmap_params$color_vector)
  
}


# figure 3b
# input files included in Supplementary data
cell_line_tumor_distance_distribution <- function(alignment, tumor_CL_cor) {
  alignment$compare_types <- alignment$lineage
  alignment[which(alignment$subtype=='Ewing sarcoma'),'compare_types'] <- 'Ewing sarcoma'
  alignment[which(alignment$subtype=='osteosarcoma'),'compare_types'] <- 'osteosarcoma'
  alignment[grep('rhabdomyosarcoma', alignment$subtype),'compare_types'] <- 'rhabdomyosarcoma'
  alignment[which(alignment$subtype=='uveal melanoma'),'compare_types'] <- 'uveal melanoma'
  alignment[which(alignment$subtype=='acute myeloid leukemia'),'compare_types'] <- 'acute myeloid leukemia'
  alignment[grep('acute lymphoblastic leukemia', alignment$subtype),'compare_types'] <- 'acute lymphoblastic leukemia'
  alignment[dplyr::filter(alignment, lineage=='breast')$sampleID[grep('basal', dplyr::filter(alignment, lineage=='breast')$subtype)],'compare_types'] <- 'basal breast'
  alignment[dplyr::filter(alignment, lineage=='breast')$sampleID[grep('luminal', dplyr::filter(alignment, lineage=='breast')$subtype)],'compare_types'] <- 'luminal breast'
  
  
  common_cancer_types <- intersect(dplyr::filter(alignment, type=='tumor')$compare_types, 
                                   dplyr::filter(alignment, type=='CL')$compare_types)
  common_cancer_types <- setdiff(common_cancer_types, c('eye', 'bone', 'breast', 'blood'))
  
  tumor_names <- character()
  CL_names <- character()
  dist_list <- numeric()
  tissue_types <- character()
  for(cancer in common_cancer_types) {
    cur_tumors <- dplyr::filter(alignment, type=='tumor' & compare_types==cancer)$sampleID
    cur_CLs <- dplyr::filter(alignment, type=='CL' & compare_types==cancer)$sampleID
    cur_dist <- reshape2::melt(as.matrix(tumor_CL_cor[cur_tumors, cur_CLs]))
    tumor_names <- c(tumor_names, as.character(cur_dist$Var1))
    CL_names <- c(CL_names, as.character(cur_dist$Var2))
    dist_list <- c(dist_list, cur_dist$value)
    tissue_types <- c(tissue_types, rep(cancer, nrow(cur_dist)))
    
  }
  
  dist_df <- cbind.data.frame(tumor_names, CL_names, dist_list, tissue_types)
  dist_df$tissue_types <- gsub("_", " ", dist_df$tissue_types)
  mean_dist <- aggregate(dist_df$dist_list, list(dist_df$tissue_types), 
                         FUN = quantile, probs = 0.25) %>% dplyr::arrange(desc(x))
  
  mean_dist$Group.1 <- rev(mean_dist$Group.1)
  dist_df$tissue_types <- factor(dist_df$tissue_types, levels = mean_dist$Group.1)
  
  tumor_dist_spread <- ggplot2::ggplot(dplyr::filter(dist_df, tissue_types != 'all'),
                                       ggplot2::aes(x = dist_list, y = tissue_types, fill = tissue_types)) +
    ggridges::geom_density_ridges(alpha=0.8) +
    ggridges::theme_ridges() +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none", text=ggplot2::element_text(size=6),
          axis.text = ggplot2::element_text(size=6)) +
    ggplot2::xlab("correlation between cell lines and tumors") +
    ggplot2::ylab('cancer type')
  
  return(tumor_dist_spread)
  
}






