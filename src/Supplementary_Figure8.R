library(here)
library(magrittr)
library(tidyverse)
source(here::here('src', 'analysis_helpers.R'))
source(here::here('src', 'Figure3.R'))

CFE_jaccard_simmilarity <- function(CFEs, Celligner_info) {
  cfe_bin <- CFEs %>% 
    dplyr::mutate(value=1) %>%
    reshape2::acast(sampleID~`CFE nodes`, value.var="value",fun.aggregate = function(x){return(ifelse(any(!is.na(x)), 1, 0))})
  
  cls <- intersect(filter(Celligner_info, type=='CL')$sampleID, rownames(cfe_bin))
  tmrs <- intersect(filter(Celligner_info, type=='tumor')$sampleID, rownames(cfe_bin))
  
  CFE_jaccard <- matrix(data=0, nrow = length(tmrs), ncol=length(cls))
  for(i in 1:nrow(cfe_bin[cls,])) {
    CFE_jaccard[,i] <- calc_jaccard(cfe_bin[cls[i],], cfe_bin[tmrs,])
  }
  colnames(CFE_jaccard) <- cls
  rownames(CFE_jaccard) <- tmrs
  
  return(CFE_jaccard)

}

calc_jaccard <- function(x, tumor_cfe) {
  cur_jac <- apply(tumor_cfe, 1, function(tmr) jaccard::jaccard(tmr, x))
  return(cur_jac)
}

get_CFE_classifications <- function(CFE_jaccard, Celligner_info) {
  
  CFE_tumor_CL_classes <- get_cell_line_tumor_class(CFE_jaccard, Celligner_info)
  
  return(CFE_tumor_CL_classes)
}

get_updated_Celligner_cor <- function(Celligner_aligned_data, CFEs) {
  cls <- intersect(filter(Celligner_info, type=='CL')$sampleID, CFEs$sampleID)
  tmrs <- intersect(filter(Celligner_info, type=='tumor')$sampleID, CFEs$sampleID)
  
  
  Celligner_cor <- cor(Celligner_aligned_data[,tmrs], 
                       Celligner_aligned_data[,cls], use='pairwise')
  
  return(Celligner_cor)
  
}

get_updated_Celligner_CFE_classifications <- function(Celligner_cor, Celligner_info) {
  
  Celligner_tumor_CL_classes <- get_cell_line_tumor_class(Celligner_cor, Celligner_info)
  
  return(Celligner_tumor_CL_classes)
}

get_CFE_class_accuracy <- function(CFE_tumor_CL_classes, Celligner_tumor_CL_classes, Celligner_info, CFEs) {
  cfe_sample_info <- filter(Celligner_info, sampleID %in% names(CFE_tumor_CL_classes))
  
  rownames(cfe_sample_info) <- cfe_sample_info$sampleID
  
  cfe_sample_info[names(CFE_tumor_CL_classes),'CFE_class'] <- CFE_tumor_CL_classes
  cfe_sample_info[names(Celligner_tumor_CL_classes),'Celligner_class'] <- Celligner_tumor_CL_classes
  
  tumor_types <- unique(filter(Celligner_info, sampleID %in% CFEs$sampleID & type=='tumor')$lineage)
  cfe_sample_info <- filter(cfe_sample_info, lineage %in% tumor_types)
  
  print(length(which(cfe_sample_info$CFE_class == cfe_sample_info$lineage))/nrow(cfe_sample_info))
  print(length(which(cfe_sample_info$Celligner_class == cfe_sample_info$lineage))/nrow(cfe_sample_info))
  
  
}


create_CFE_Celligner_plot <- function(cur_cancer_type, CFE_jaccard, Celligner_cor, Celligner_tumor_CL_classes, Celligner_info) {
  cfe_sample_info <- filter(Celligner_info, sampleID %in% names(Celligner_tumor_CL_classes))
  
  rownames(cfe_sample_info) <- cfe_sample_info$sampleID
  cfe_sample_info[names(Celligner_tumor_CL_classes),'Celligner_class'] <- Celligner_tumor_CL_classes

  cfe_sample_info <- filter(cfe_sample_info, lineage == cur_cancer_type)
  

  cur_tumors <- filter(Celligner_info, sampleID %in% rownames(CFE_jaccard) & lineage == cur_cancer_type)$sampleID
  
  CFE_median <- apply(CFE_jaccard[cur_tumors, cfe_sample_info$sampleID], 2, median)
  Celligner_median <- apply(Celligner_cor[cur_tumors, cfe_sample_info$sampleID], 2, median)
  
  cfe_sample_info$CFE_median <- CFE_median
  cfe_sample_info$Celligner_median <- Celligner_median
  cfe_sample_info$classification <- cfe_sample_info$Celligner_class == cfe_sample_info$lineage
  
  cur_plot <- ggplot(cfe_sample_info, aes(Celligner_median, CFE_median)) +
    geom_point(aes(color=classification)) + 
    theme_classic() +
    xlab('Celligner') + ylab("CFE") +
    scale_color_manual(values=c(`TRUE`="#00BFC4", `FALSE`= "#F8766D")) +
    ggpubr::stat_cor(label.y.npc = 'top', size=3, color='black') +
    ggplot2::geom_smooth(method = 'lm', color='black') + ggtitle(gsub("_", " ", cur_cancer_type)) +
    theme(legend.position = 'none', text = element_text(size=6), 
          axis.text = element_text(size=4), axis.title = element_text(size=6))
  
  return(cur_plot)
}



