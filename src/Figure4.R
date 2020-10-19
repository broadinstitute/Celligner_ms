library(here)
library(magrittr)
library(tidyverse)
source(here::here('src', 'analysis_helpers.R'))
source(here::here('src', 'global_params.R'))

# use the Celligner_info file for the alignment matrix used within these methods

# figure 4a
make_select_subtypes_plot <- function(alignment) {
  select_subtypes <- c("Kidney Cancer", "Skin Cancer", "Breast Cancer", 'Myeloma', "Lymphoma", "Leukemia")
  
  select_subtypes_overall_plot <- ggplot2::ggplot(alignment,
                                                  ggplot2::aes(UMAP_1, UMAP_2, fill=disease %in% select_subtypes, size=type)) +
    ggplot2::geom_point(alpha=0.7, pch=21, color='white') +
    ggplot2::geom_point(data = dplyr::filter(alignment, type=='CL' & disease %in% select_subtypes),
                        ggplot2::aes(fill=disease %in% select_subtypes), color='black', pch=21, alpha=0.7) +
    ggplot2::scale_size_manual(values=c(`CL`=1.5, `tumor`=0.75)) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = 'none',
          text=ggplot2::element_text(size=7),
          axis.text=ggplot2::element_text(size=7)) +
    ggplot2::scale_fill_manual(values=c(`TRUE`='#EE4B33', `FALSE`='#999999')) + 
    ggplot2::guides(color=FALSE, fill=FALSE) +
    ggplot2::xlab("UMAP 1") +
    ggplot2::ylab("UMAP 2")
  
  return(select_subtypes_overall_plot)
  
}

# figure 4b
make_breast_subtypes_plot <- function(alignment) {
  # consolidate subtype annotations
  breast_subtypes <- c("luminal", "basal", 'HER2-enriched', 'basal A', 'basal B')
  breast_data <- dplyr::filter(alignment, disease=='Breast Cancer')
  rownames(breast_data) <- breast_data$sampleID
  breast_data[grep("luminal", breast_data$subtype), "subtype"] <- 'luminal'
  breast_data[grep("HER2", breast_data$subtype), "subtype"] <- 'HER2-enriched'
  breast_data[which(!(breast_data$subtype %in% breast_subtypes)), 'subtype'] <- NA
  
  breast_data$subtype_shape <- 'other'
  breast_data[which(breast_data$subtype == 'basal A'), 'subtype_shape'] <- 'A'
  breast_data[which(breast_data$subtype == 'basal B'), 'subtype_shape'] <- 'B'
  
  breast_data[grep("basal", breast_data$subtype), "subtype"] <- 'basal'
 
  
  breast_subtypes_plot <- ggplot2::ggplot(breast_data, ggplot2::aes(UMAP_1, UMAP_2, fill=subtype, color = type, size=type)) +
    ggplot2::geom_point(data = filter(breast_data, !subtype_shape %in% c('A', 'B')), alpha=0.7, pch = 21) +
    ggplot2::geom_point(data = filter(breast_data, subtype_shape == 'A'), pch = 22, show.legend = F, alpha=0.7) +
    ggplot2::geom_point(data = filter(breast_data, subtype_shape == 'B'), pch = 25, show.legend = F, alpha=0.7) +
    ggplot2::scale_size_manual(values=c(`CL`=1.5, `tumor`=1)) +
    ggplot2::scale_color_manual(values=c(`CL`='black', `tumor`='white')) +
    ggplot2::theme_classic() + 
    ggplot2::theme(legend.position = c(0.85,0.6),
                            text=ggplot2::element_text(size=9),
          axis.text = ggplot2::element_text(size=7), axis.title = ggplot2::element_text(size=7)) +
    ggplot2::guides(size=FALSE, color=FALSE, shape=FALSE) +
    ggplot2::guides(fill=ggplot2::guide_legend(title="")) +
    ggplot2::xlab("UMAP 1") +
    ggplot2::ylab("UMAP 2")
  
  return(breast_subtypes_plot)
}

# figure 4c
make_heme_types_plot <- function(alignment) {
  heme_types <- c('Myeloma', "Leukemia", "Lymphoma")
  heme_data <- dplyr::filter(alignment, disease %in% heme_types)
  rownames(heme_data) <- heme_data$sampleID

  # consolidate subtype annotations
  t_cell_tumor <- dplyr::filter(alignment, subtype=='acute lymphoblastic leukemia' & UMAP_2 < 3)$sampleID
  cur_subtype_focus <- heme_data$sampleID[which(heme_data$subtype %in% c("acute lymphoblastic leukemia", 'acute myeloid leukemia', 
                                                                         "acute lymphoblastic leukemia, b-cell", "acute lymphoblastic leukemia, t-cell"))]
  
  all_target <- grep("TARGET", dplyr::filter(heme_data, subtype=='acute lymphoblastic leukemia')$sampleID, value=T)
  heme_data[all_target,'subtype'] <- paste0('ALL', ', b-cell')
  heme_data[which(heme_data$subtype == 'acute lymphoblastic leukemia'),'subtype'] <- 'ALL, unknown'
  heme_data[which(heme_data$subtype == 'acute lymphoblastic leukemia, t-cell'),'subtype'] <- 'ALL, t-cell'
  heme_data[which(heme_data$subtype == 'acute lymphoblastic leukemia, b-cell'),'subtype'] <- 'ALL, b-cell'
  heme_data[which(heme_data$subtype == 'acute myeloid leukemia'),'subtype'] <- 'AML'
  
  heme_subtypes_plot <- ggplot2::ggplot(heme_data, ggplot2::aes(UMAP_1, UMAP_2, size=type)) +
    ggplot2::geom_point(data = dplyr::filter(heme_data, type=='tumor' & !sampleID %in% cur_subtype_focus),
                        ggplot2::aes(fill=disease),alpha=0.7,  color='white', pch=21) +
    ggplot2::geom_point(data = dplyr::filter(heme_data, type=='tumor' & sampleID %in% cur_subtype_focus),
                        ggplot2::aes(fill=subtype), alpha=0.7, color='white', pch=21) +
    ggplot2::geom_point(data = dplyr::filter(heme_data, type=='CL' & sampleID %in% cur_subtype_focus), 
               aes(fill=subtype), color='black', alpha=0.7, pch=21) +
    ggplot2::geom_point(data = dplyr::filter(heme_data, type=='CL' & !sampleID %in% cur_subtype_focus & UMAP_1 < 4), 
                        ggplot2::aes(fill=disease), color='black', alpha=0.7, pch=21) +
    ggplot2::geom_point(data = dplyr::filter(heme_data, type=='tumor' & sampleID %in% cur_subtype_focus & sampleID %in% t_cell_tumor),
                        ggplot2::aes(fill=subtype), alpha=0.5, color='white', pch=21) +
    ggplot2::scale_size_manual(values=c(`CL`=1.5, `tumor`=1)) +
    ggplot2::theme_classic() + 
    ggplot2::theme(legend.position = c(0.91,0.55), 
          text=ggplot2::element_text(size=9),
          axis.text = ggplot2::element_text(size=7),
          axis.title = ggplot2::element_text(size=7),
          plot.margin = ggplot2::unit(c(0.2,1.3,0.2,0.2), "cm")) +
    guides(size=FALSE, fill=ggplot2::guide_legend(title="")) +
    ggplot2::scale_fill_manual(values=c(`ALL, unknown`='grey30',
                               `ALL, b-cell`='#F8766D',
                               `AML`='#C77CFF',
                               `ALL, t-cell`='#FFB625',
                               `Leukemia`='#3CB44B',
                               `Lymphoma`='#00BFC4',
                               `Myeloma`='#FABEBE')) +
    ggplot2::xlab("UMAP 1") +
    ggplot2::ylab("UMAP 2")
  
  return(heme_subtypes_plot)

}

# figure 4d
make_plasma_cell_plot <- function(alignment) {
  heme_types <- c('Myeloma', "Leukemia", "Lymphoma")
  heme_data <- dplyr::filter(alignment, disease %in% heme_types)
  rownames(heme_data) <- heme_data$sampleID
  
  
  
  plasma_cell_highlight <- ggplot2::ggplot(dplyr::filter(heme_data, UMAP_1 > -2 & UMAP_1 < 3 & UMAP_2 < 2 & UMAP_2 > -1),
                                           ggplot2::aes(UMAP_1, UMAP_2, size=type, fill=disease, color=type), alpha=0.7) +
    ggplot2::geom_point(alpha=0.7, pch=21) + 
    ggplot2::geom_point(data = dplyr::filter(heme_data, type=='CL' & UMAP_1 > 0.8 & UMAP_1 < 1.4 & UMAP_2 < 1.4 & UMAP_2 > 0.6),
               alpha=0.5, pch=21) +
    ggplot2::guides(size=FALSE, shape=FALSE) +
    ggplot2::scale_size_manual(values=c(`CL`=1.5, `tumor`=1)) +
    ggplot2::scale_color_manual(values=c(`CL`='black', `tumor`='white')) + 
    ggplot2::theme_classic() + 
    ggplot2::theme(legend.position = c(1.1, 0.4),
          text=ggplot2::element_text(size=9),
          axis.text = ggplot2::element_text(size = 7),
          axis.title = ggplot2::element_text(size = 7),
          plot.margin = ggplot2::unit(c(0.2,1.6,0.2,0.2), 'cm')) +
    ggplot2::guides(size=FALSE, color=FALSE, fill=ggplot2::guide_legend(title="")) +
    ggplot2::scale_fill_manual(values=c(`Leukemia`='#3CB44B', `Lymphoma`='#00BFC4', `Myeloma`='#FABEBE')) +
    ggplot2::xlab("UMAP 1") +
    ggplot2::ylab("UMAP 2")
  
  return(plasma_cell_highlight)
  
}

# figure 4e
make_skin_subtypes_plot <- function(alignment) {
  skin_subtype_data <- dplyr::filter(alignment, disease == 'Skin Cancer' & UMAP_2 <1 & UMAP_1 > 1)
  skin_subtypes <- c('transitory', 'melanocytic', "undifferentiated", "neural crest-like")
  skin_subtype_data[which(!(skin_subtype_data$subtype %in% skin_subtypes)), 'subtype'] <- NA
  
  skin_subtypes_plot <- ggplot2::ggplot(skin_subtype_data, ggplot2::aes(UMAP_1, UMAP_2, fill=subtype, size=type, color=type)) +
    ggplot2::geom_point(alpha=0.7, pch=21) +
    ggplot2::geom_point(data = dplyr::filter(skin_subtype_data, type=='CL' & !is.na(subtype)),
                        ggplot2::aes(fill=subtype), color='black', alpha=0.6, pch=21) +
    ggplot2::geom_point(data = filter(skin_subtype_data, type=='CL' & is.na(subtype)),
                        ggplot2::aes(fill=subtype), color='black', alpha=0.2, pch=21) +
    ggplot2::scale_color_manual(values=c(`CL`='black', `tumor`='white'), guide='none') +
    ggplot2::scale_size_manual(values=c(`CL`=1.5, `tumor`=1), guide='none') +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = c(0.75, 0.6),
          text=ggplot2::element_text(size=9),
          legend.margin =ggplot2::margin(0,0,0,0),
          legend.box.margin=ggplot2::margin(-10,-30,-10,-10),
          axis.text = ggplot2::element_text(size = 7),
          axis.title = ggplot2::element_text(size = 7)) +
    ggplot2::guides(size=FALSE, fill=ggplot2::guide_legend(title="")) +
    ggplot2::xlab("UMAP 1") + 
    ggplot2::ylab("UMAP 2")
  
  return(skin_subtypes_plot)
  
}

# figure 4e inset
make_skin_subtypes_highlight <- function(alignment) {
  skin_subtype_data <- dplyr::filter(alignment, disease == 'Skin Cancer' & UMAP_2 <1 & UMAP_1 > 1)
  skin_subtypes <- c('transitory', 'melanocytic', "undifferentiated", "neural crest-like")
  skin_subtype_data[which(!(skin_subtype_data$subtype %in% skin_subtypes)), 'subtype'] <- NA
  
  skin_subtypes_highlight <- ggplot2::ggplot(dplyr::filter(skin_subtype_data, UMAP_2 < -9 & UMAP_1 < 14),
                                             ggplot2::aes(UMAP_1, UMAP_2, fill=subtype, size=type, color=type), alpha=0.6) +
    ggplot2::geom_point(alpha=0.7, pch=21) +
    ggplot2::geom_point(data = dplyr::filter(skin_subtype_data, type=='CL' & !is.na(subtype) & UMAP_2 < -9 & UMAP_1 < 14),
                        ggplot2::aes(fill=subtype), color='black', alpha=0.6, pch=21) +
    ggplot2::scale_color_manual(values=c(`CL`='black', `tumor`='white'), guide='none') +
    ggplot2::scale_size_manual(values=c(`CL`=1.5, `tumor`=1), guide='none') +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = 'none', 
          axis.title.x=ggplot2::element_blank(),
          axis.text.x=ggplot2::element_blank(),
          axis.ticks.x=ggplot2::element_blank(), 
          axis.title.y=ggplot2::element_blank(),
          axis.text.y=ggplot2::element_blank(),
          axis.ticks.y=ggplot2::element_blank())
  
  return(skin_subtypes_highlight)
}


# figure 4f
make_kidney_subtypes_plot <- function(alignment) {
  kidney_subtype_data <- dplyr::filter(alignment, disease=='Kidney Cancer')
  kidney_subtypes <- c('wilms tumor', 'papillary renal cell carcinoma', "kidney clear cell carcinoma", "kidney chromophobe")
  kidney_subtype_data[which(!(kidney_subtype_data$subtype %in% kidney_subtypes)), 'subtype'] <- NA
  
  kidney_subtypes_plot <- ggplot2::ggplot(kidney_subtype_data, 
                                          ggplot2::aes(UMAP_1, UMAP_2, fill=stringr::str_wrap(subtype, 10), size=type, color=type)) +
    ggplot2::geom_point(alpha=0.7, pch=21) +
    ggplot2::geom_point(data = dplyr::filter(kidney_subtype_data, type=='CL'),
               color='black', alpha=0.5, fill = 'gray80', pch=21) +
    ggplot2::scale_color_manual(values=c(`CL`='black', `tumor`='white'), guide='none') +
    ggplot2::scale_size_manual(values=c(`CL`=1.5, `tumor`=1), guide='none') +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = c(1.1, 0.7),
          text=ggplot2::element_text(size=9),
          axis.text = ggplot2::element_text(size = 7),
          axis.title = ggplot2::element_text(size = 7),
          plot.margin = ggplot2::unit(c(0.2,1.8,0.2,0.2), 'cm')) + 
    ggplot2::labs(fill = "") +
    ggplot2::guides(size=FALSE) +
    ggplot2::xlab("UMAP 1") +
    ggplot2::ylab("UMAP 2")
  
  return(kidney_subtypes_plot)
}




