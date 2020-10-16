library(here)
library(magrittr)
library(tidyverse)
source(here::here('src', 'analysis_helpers.R'))
source(here::here('src', 'Figure3.R'))

# Parse correlation matrix from Yu et al.
# to match cell line and tumor samples between Yu et al and Celligner
# yu_tumor_CL_cor is taken from Supplementary Data 4 from 
# Yu et al, Nature Communications 2019
parse_yu_et_al_cor <- function(yu_tumor_CL_cor, Celligner_info) {
  colnames(yu_tumor_CL_cor) <- gsub("_.*", "", colnames(yu_tumor_CL_cor))
  yu_sample_info <- filter(Celligner_info,
                           gsub("_.*", "", sampleID_CCLE_Name) %in% colnames(yu_tumor_CL_cor))
  yu_sample_info <- yu_sample_info[-grep("TT_OESOPHAGUS", yu_sample_info$sampleID_CCLE_Name),]
  yu_sample_info <- yu_sample_info[-grep("RESISTANT", yu_sample_info$sampleID_CCLE_Name),]
  
  rownames(yu_sample_info) <- gsub("_.*", "", yu_sample_info$sampleID_CCLE_Name)
  colnames(yu_tumor_CL_cor) <- yu_sample_info[colnames(yu_tumor_CL_cor),'sampleID']
  
  yu_tumor_CL_cor <- yu_tumor_CL_cor[intersect(Celligner_info$sampleID, rownames(yu_tumor_CL_cor)),]
  yu_tumor_CL_cor <- as.matrix(yu_tumor_CL_cor)
  return(yu_tumor_CL_cor)
} 

# Get nearest neighbor classifications for Yu et al method
get_yu_et_al_classifications <- function(yu_tumor_CL_cor, Yu_tumor_ann) {
  Yu_tumor_ann$sampleID <- substr(Yu_tumor_ann$sample,1,nchar(Yu_tumor_ann$sample)-1)
  Yu_tumor_ann$lineage <- Yu_tumor_ann$Disease
  Yu_tumor_ann[which(Yu_tumor_ann$lineage=='COADREAD'),'lineage'] <- "COAD_READ"
  
  yu_tumor_CL_classes <- get_cell_line_tumor_class(yu_tumor_CL_cor, Yu_tumor_ann)

  return(yu_tumor_CL_classes)
}

# Get updated Celligner tumor/CL correlation using subset of data from
# Yu et al 
get_Celligner_updated_cor <- function(Celligner_aligned_data, yu_tumor_CL_cor) {
  common_CLs <- colnames(yu_tumor_CL_cor)
  common_tumors <- rownames(yu_tumor_CL_cor)
  
  Celligner_cor <- cor(Celligner_aligned_data[,common_tumors], 
                       Celligner_aligned_data[,common_CLs], use='pairwise')
  
  return(Celligner_cor)
  
}

# Get updated Celligner classifications using updated correlation matrix
# Yu_tumor_ann is taken from Supplementary Data 1 from 
# Yu et al, Nature Communications 2019
get_Celligner_updated_classes <- function(Celligner_cor, Yu_tumor_ann) {
  Yu_tumor_ann$sampleID <- substr(Yu_tumor_ann$sample,1,nchar(Yu_tumor_ann$sample)-1)
  Yu_tumor_ann$lineage <- Yu_tumor_ann$Disease
  Yu_tumor_ann[which(Yu_tumor_ann$lineage=='COADREAD'),'lineage'] <- "COAD_READ"
  
  Celligner_tumor_CL_classes <- get_cell_line_tumor_class(Celligner_cor, Yu_tumor_ann)
  
  return(Celligner_tumor_CL_classes)
}

# Compare classification accuracy between the Yu et al method and Celligner
# CL_annotations is taken from Supplementary Data 2 from 
# Yu et al, Nature Communications 2019
get_classification_accuracy <- function(CL_annotations, Celligner_info,
                                        yu_tumor_CL_classes, Celligner_tumor_CL_classes) {
  
  CL_annotations$CCLE <- gsub("_.*", "", CL_annotations$CCLE_name)
  yu_sample_info <- filter(Celligner_info,
                           gsub("_.*", "", sampleID_CCLE_Name) %in% CL_annotations$CCLE)
  yu_sample_info <- yu_sample_info[-grep("TT_OESOPHAGUS", yu_sample_info$sampleID_CCLE_Name),]
  yu_sample_info <- yu_sample_info[-grep("RESISTANT", yu_sample_info$sampleID_CCLE_Name),]
  
  rownames(yu_sample_info) <- gsub("_.*", "", yu_sample_info$sampleID_CCLE_Name)
  
  CL_annotations$sampleID <- yu_sample_info[CL_annotations$CCLE,'sampleID']
  rownames(CL_annotations) <- CL_annotations$sampleID
  
  CL_annotations[names(yu_tumor_CL_classes), 'Yu_classification'] <- yu_tumor_CL_classes 
  CL_annotations[names(Celligner_tumor_CL_classes), 'Celligner_classification'] <- Celligner_tumor_CL_classes  
  CL_annotations[which(CL_annotations$disease == 'COADREAD'), 'disease'] <- "COAD_READ"
  
  print(length(which(CL_annotations$disease == CL_annotations$Yu_classification))/nrow(CL_annotations))
  print(length(which(CL_annotations$disease == CL_annotations$Celligner_classification))/nrow(CL_annotations))
  
  return(CL_annotations)
  
}

# Calculate median correlation per cancer type (using Yu et al categories)
median_cor <- function(cor_mat, Yu_tumor_ann) {
  Yu_tumor_ann$sampleID <- substr(Yu_tumor_ann$sample,1,nchar(Yu_tumor_ann$sample)-1)
  Yu_tumor_ann[which(Yu_tumor_ann$Disease=='COADREAD'),'Disease'] <- "COAD_READ"
  yu_tumor_types <- unique(Yu_tumor_ann$Disease)
  
  med_cor_mat <- matrix(data=NA, nrow=length(yu_tumor_types), ncol=ncol(cor_mat))
  rownames(med_cor_mat) <- yu_tumor_types
  colnames(med_cor_mat) <- colnames(cor_mat)
  for(tcga in yu_tumor_types) {
    cur_tumors <- intersect(filter(Yu_tumor_ann, Disease  == tcga)$sampleID, rownames(cor_mat))
    cur_cor_med <- apply(cor_mat[cur_tumors,], 2, function(x) median(x))
    med_cor_mat[tcga,] <- cur_cor_med
  }
  
  return(med_cor_mat)
  
}

# Supplementary Figure 6
create_comparison_plots <- function(cur_cancer_type, yu_tumor_CL_cor, Celligner_cor,
                                    Yu_tumor_ann, CL_annotations) {
  rownames(CL_annotations) <- CL_annotations$sampleID
  yu_median_cor <- median_cor(yu_tumor_CL_cor, Yu_tumor_ann)
  Celligner_median_cor <- median_cor(Celligner_cor, Yu_tumor_ann)
  
  cur_samples <- filter(CL_annotations, disease == cur_cancer_type)$sampleID
  cur_yu <- yu_median_cor[cur_cancer_type, cur_samples] %>% t()
  cur_yu <- cur_yu[1,]
  
  cur_celligner <- Celligner_median_cor[cur_cancer_type, cur_samples] %>% t()
  cur_celligner <- cur_celligner[1,]
  
  cur_yu_celligner_compare <- data.frame(yu_et_al = cur_yu, 
                                         celligner = cur_celligner, 
                                         sample = cur_samples,
                                         cur_class = ifelse(CL_annotations[cur_samples,]$Celligner_classification == cur_cancer_type,
                                                            TRUE, FALSE))
  
  cur_comparison_plot <- ggplot(cur_yu_celligner_compare, aes(celligner, yu_et_al)) + 
    geom_point(aes(color=cur_class)) + 
    theme_classic() +
    xlab('Celligner') + ylab("Yu et al") +
    scale_color_manual(values=c(`TRUE`="#00BFC4", `FALSE`= "#F8766D")) +
    ggpubr::stat_cor(label.y.npc = 'top', size=2.8, color='black') +
    ggplot2::geom_smooth(method = 'lm', color='black') + ggtitle(cur_cancer_type) +
    theme(legend.position = 'none', text = element_text(size=6), 
          axis.text = element_text(size=4), axis.title = element_text(size=6))
  
  return(cur_comparison_plot)
  
}







