library(here)
library(magrittr)
library(tidyverse)
source(here::here('src', 'analysis_helpers.R'))
source(here::here('src', 'Figure3.R'))

# Parse CCL_sampleTab (part of Supplemental Table1)
# to match cell line and tumor samples between CancerCellNet and Celligner
# CCL_sampleTab is taken from the first tab of Supplemental Table 1 of 
# Peng et al, BioRxiv 2020
parse_CancerCellNet <- function(CCL_sampleTab, Celligner_info) {
  
  CCL_sampleTab$CCLE_Name <- CCL_sampleTab$cancer_Category
  CCL_sampleTab$CCLE_Name <- gsub(" ", "", CCL_sampleTab$CCLE_Name)
  CCL_sampleTab$CCLE_Name <- gsub("-", "", CCL_sampleTab$CCLE_Name)
  CCL_sampleTab$CCLE_Name <- gsub("\\.", "", CCL_sampleTab$CCLE_Name)
  CCL_sampleTab$CCLE_Name <- gsub("\\(", "", CCL_sampleTab$CCLE_Name)
  CCL_sampleTab$CCLE_Name <- gsub("\\)", "", CCL_sampleTab$CCLE_Name)
  CCL_sampleTab$CCLE_Name <- gsub("\\/", "", CCL_sampleTab$CCLE_Name)
  CCL_sampleTab$CCLE_Name <- gsub(",", "", CCL_sampleTab$CCLE_Name)
  CCL_sampleTab$CCLE_Name <- gsub(":", "", CCL_sampleTab$CCLE_Name)
  CCL_sampleTab$CCLE_Name <- gsub("_.*", "", CCL_sampleTab$CCLE_Name)
  CCL_sampleTab[which(CCL_sampleTab$cancer_Category=='TT'),'CCLE_Name'] <- 'TDOTT'
  CCL_sampleTab$CCLE_Name <- toupper(CCL_sampleTab$CCLE_Name)
  
  CCL_samples <- filter(Celligner_info, 
                        gsub("_.*", "", sampleID_CCLE_Name) %in% CCL_sampleTab$CCLE_Name)
  CCL_samples[which(CCL_samples$sampleID_CCLE_Name=='TT_OESOPHAGUS'),'sampleID_CCLE_Name'] <- "TDOTT_OESOPHAGUS"
  CCL_samples <- CCL_samples[-grep("RESISTANT", CCL_samples$sampleID_CCLE_Name),]

  rownames(CCL_sampleTab) <- CCL_sampleTab$CCLE_Name
  CCL_sampleTab[gsub("_.*", "", CCL_samples$sampleID_CCLE_Name),'sampleID'] <- CCL_samples$sampleID
  
  return(CCL_sampleTab)
} 

# Classify each cell line using CCN probabilities
# CCL_Classification is taken from the second tab of Supplemental Table 1 of 
# Peng et al, BioRxiv 2020
get_CCN_classifications <- function(CCL_Classification, CCL_sampleTab) {
  colnames(CCL_Classification) <- CCL_sampleTab$sampleID
  
  # get classification as tumor type (not unknown) using max probability output by models
  CCN_tumor_CL_classes <- apply(CCL_Classification[1:22,], 2, 
                                function(x) rownames(CCL_Classification)[which.max(x)])
  return(CCN_tumor_CL_classes)
}

# Get updated Celligner tumor/CL correlation using subset of data from
# CancerCellNet. Tumor samples used by CancerCellNet were not provided, so 
# tumor samples are approximated using TCGA data that fits with the categories
# used by CancerCellNet
get_Celligner_updated_cor <- function(Celligner_aligned_data, Celligner_info) {
  common_CLs <- filter(Celligner_info, type=='CL' & !is.na(CancerCellNet_annotation))$sampleID
  common_tumors <- filter(Celligner_info, type=='tumor' & !is.na(CancerCellNet_annotation))$sampleID
  
  Celligner_cor <- cor(Celligner_aligned_data[,common_tumors], 
                       Celligner_aligned_data[,common_CLs], use='pairwise')
  
  return(Celligner_cor)
  
}

# Get updated Celligner classifications using updated correlation matrix
get_Celligner_updated_classes <- function(Celligner_cor, Celligner_info) {
  CCN_sample_info <- filter(Celligner_info, !is.na(CancerCellNet_annotation))
  CCN_sample_info$lineage <- CCN_sample_info$CancerCellNet_annotation
  Celligner_tumor_CL_classes <- get_cell_line_tumor_class(Celligner_cor, CCN_sample_info)
  
  return(Celligner_tumor_CL_classes)
}

# Compare classification accuracy between CancerCellNet and Celligner
get_CCN_classification_accuracy <- function(Celligner_info,
                                        CCN_tumor_CL_classes, Celligner_tumor_CL_classes) {
  
  compare_classes <- filter(Celligner_info, type=='CL' & !is.na(CancerCellNet_annotation))
  rownames(compare_classes) <- compare_classes$sampleID
   
  compare_classes[names(CCN_tumor_CL_classes),'CCN_classification'] <- CCN_tumor_CL_classes
  compare_classes[names(Celligner_tumor_CL_classes),'Celligner_classification'] <- Celligner_tumor_CL_classes
  
  print(length(which(compare_classes$CancerCellNet_annotation == compare_classes$CCN_classification))/nrow(compare_classes))
  print(length(which(compare_classes$CancerCellNet_annotation == compare_classes$Celligner_classification))/nrow(compare_classes))
  
  
}

# Calculate median correlation per cancer type (using CCN categories)
Celligner_median_cor <- function(Celligner_cor, Celligner_info) {
  CCN_tumor_types <- unique(filter(Celligner_info, type=='tumor' & !is.na(CancerCellNet_annotation))$CancerCellNet_annotation)
  
  Celligner_med_cor <- matrix(data=NA, nrow=length(CCN_tumor_types), ncol=ncol(Celligner_cor))
  rownames(Celligner_med_cor) <- CCN_tumor_types
  colnames(Celligner_med_cor) <- colnames(Celligner_cor)
  for(tcga in CCN_tumor_types) {
    cur_tumors <- filter(Celligner_info, type=='tumor' & CancerCellNet_annotation  == tcga)$sampleID
    cur_cor_med <- apply(Celligner_cor[cur_tumors,], 2, function(x) median(x))
    Celligner_med_cor[tcga,] <- cur_cor_med
  }
  
  return(Celligner_med_cor)
  
}

# Supplementary Figure 7
create_CCN_comparison_plots <- function(cur_cancer_type, CCL_Classification, Celligner_med_cor,
                                    Celligner_info, Celligner_tumor_CL_classes, CCL_sampleTab) {
    colnames(CCL_Classification) <- CCL_sampleTab$sampleID
  
  
    rownames(Celligner_info) <- Celligner_info$sampleID
    Celligner_info$Celligner_class <- NA
    Celligner_info[names(Celligner_tumor_CL_classes),]$Celligner_class <- Celligner_tumor_CL_classes
    
    cur_CLs <- filter(Celligner_info, type=='CL' & CancerCellNet_annotation == cur_cancer_type)$sampleID
    cur_ccn <- CCL_Classification[cur_cancer_type, cur_CLs] %>% t()
    cur_ccn <- cur_ccn[1,]
    
    cur_celligner <- Celligner_med_cor[cur_cancer_type, cur_CLs] %>% t()
    cur_celligner <- cur_celligner[1,]
    
    cur_ccn_celligner_compare <- data.frame(cancercellnet = cur_ccn, 
                                            celligner = cur_celligner, 
                                            sample = cur_CLs,
                                            cur_class = ifelse(Celligner_info[cur_CLs,]$Celligner_class == cur_cancer_type,
                                                               TRUE, FALSE))
    CCN_Celligner_comparison_plot <- ggplot(cur_ccn_celligner_compare, aes(celligner, cancercellnet)) + 
      geom_point(aes(color=cur_class)) + 
      theme_classic() +
      xlab('Celligner') + ylab("CancerCellNet") +
      geom_hline(aes(yintercept = 0.3), color = 'gray30', linetype = "dashed") +
      scale_color_manual(values=c(`TRUE`="#00BFC4", `FALSE`= "#F8766D")) +
      ggpubr::stat_cor(label.y.npc = 'top', size=3, color='black') +
      ggplot2::geom_smooth(method = 'lm', color='black') + ggtitle(cur_cancer_type) +
      theme(legend.position = 'none', text = element_text(size=6), 
            axis.text = element_text(size=4), axis.title = element_text(size=6))
    
    return(CCN_Celligner_comparison_plot)
}







