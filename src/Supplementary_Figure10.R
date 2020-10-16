
# Supplementary Figure 10a
# Sampling year is taken from file model_list_20200702.csv taken
# from Sanger cell model passports (https://cellmodelpassports.sanger.ac.uk/downloads)
cell_line_sampling_year <- function(Celligner_info) {

  sampling_year <- ggplot2::ggplot(dplyr::filter(Celligner_info, !is.na(sampling_year)),
                                   ggplot2::aes(undifferentiated_cluster, sampling_year, color=undifferentiated_cluster)) +
    ggplot2::geom_boxplot() + 
    ggplot2::xlab('undifferentiated') + 
    ggplot2::ylab('sampling year') +
    ggpubr::stat_compare_means(method = 'wilcox.test', label = "p.format", size=3, paired = F) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = 'none', text=ggplot2::element_text(size=8))
  
  return(sampling_year)
  
}

# Supplementary Figure 10b
undifferentiated_age <- function(Celligner_info) {
  CL_undifferentiated_age <- filter(Celligner_info, type=='CL') %>%
    dplyr::select('sampleID', 'undifferentiated_cluster', 'age')
  age_of_patient <- ggplot(CL_undifferentiated_age, 
                           aes(undifferentiated_cluster, age, color=undifferentiated_cluster)) +
    geom_boxplot() + 
    xlab('undifferentiated') +
    ggpubr::stat_compare_means(method = 'wilcox.test', label = "p.format", size=3) +
    theme_classic() +
    theme(legend.position = 'none', text=element_text(size=8))
  
  return(age_of_patient)
}

# Supplementary Figure 10c
undifferentiated_sex <- function(Celligner_info) {
  clinical_features <- filter(Celligner_info, type=='CL' & sex != 'Unknown') %>% 
    dplyr::select('sampleID', 'sex', 'undifferentiated_cluster', 'lineage')

  clinical_features_prop <- clinical_features %>%
    group_by(undifferentiated_cluster, sex) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n))
  
  
  undifferentiated_sex <- ggplot(clinical_features_prop, aes(sex, freq, fill = undifferentiated_cluster)) +
    geom_col(position='dodge') +
    ylab('proportion') +
    theme_classic() +
    guides(fill=guide_legend(title="undifferentiated")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), text=element_text(size=8))
  
  return(undifferentiated_sex)
}

# Supplementary Figure 10d
undifferentiated_metastatic <- function(Celligner_info) {
  met_prop <- filter(Celligner_info, type=='CL' & `Primary/Metastasis` != "Unknown") %>%
    group_by(undifferentiated_cluster, `Primary/Metastasis`) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n))
  
  
  undifferentiated_met <- ggplot(met_prop, aes(`Primary/Metastasis`, freq, fill = undifferentiated_cluster)) +
    geom_col(position='dodge') +
    theme_classic() +
    ylab("proportion") + xlab('type') +
    guides(fill=guide_legend(title="undifferentiated")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), text=element_text(size=8))
  
  return(undifferentiated_metastatic)
}