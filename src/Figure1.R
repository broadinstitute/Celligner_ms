library(here)
library(magrittr)
library(tidyverse)
source(here::here('src', 'analysis_helpers.R'))
source(here::here('src', 'global_params.R'))

# figure 1a
# TCGA_mat source: https://xenabrowser.net/datapages/?dataset=TumorCompendium_v10_PolyA_hugo_log2tpm_58581genes_2019-07-25.tsv&host=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
# CCLE_mat source: depmap.org DepMap Public 19Q4 CCLE_expression_full.csv
# using genes that are common to both datasets, except genes labeled as 'non-coding RNA' or 'pseudogene' by HGNC
plot_uncorrected_data <- function(CCLE_mat, TCGA_mat) {
  
  comb_ann <- cbind.data.frame(`sampleID` = c(rownames(TCGA_mat), rownames(CCLE_mat)),
                               `type` = c(rep('tumor', nrow(TCGA_mat)), rep("CL", nrow(CCLE_mat))))
  
  # using Seurat object to run cPCA and UMAP
  original_combined_obj <-  Seurat::CreateSeuratObject(t(rbind(TCGA_mat,
                                                               CCLE_mat)),
                                                       min.cells = 0,
                                                       min.features = 0,
                                                       meta.data = comb_ann %>%
                                                         magrittr::set_rownames(comb_ann$sampleID))
  
  original_combined_obj <- Seurat::ScaleData(original_combined_obj, 
                                             features = rownames(Seurat::GetAssayData(original_combined_obj)), 
                                             do.scale = F)
  
  original_combined_obj %<>% Seurat::RunPCA(assay='RNA',
                                            features = rownames(Seurat::GetAssayData(original_combined_obj)),
                                            npcs = global$n_PC_dims,
                                            verbose = F)
  
  original_combined_obj %<>% Seurat::RunUMAP(assay = 'RNA', dims = 1:global$n_PC_dims,
                                             reduction = 'pca',
                                             n.neighbors = global$umap_n_neighbors,
                                             min.dist = global$umap_min_dist,
                                             metric = global$distance_metric,
                                             verbose=F)
  
  uncorrected_alignment <- Seurat::Embeddings(original_combined_obj, reduction = 'umap') %>%
    as.data.frame() %>%
    set_colnames(c('UMAP_1', 'UMAP_2')) %>%
    rownames_to_column(var = 'sampleID') %>%
    left_join(comb_ann, by = 'sampleID')
  
  uncorrected_combined_type_plot <- ggplot2::ggplot(uncorrected_alignment, ggplot2::aes(UMAP_1, UMAP_2)) +
    ggplot2::geom_point(data = filter(uncorrected_alignment, type=='tumor'), alpha=0.6, size=0.5, pch=21, color='white', aes(fill=type)) +
    ggplot2::geom_point(data = filter(uncorrected_alignment, type=='CL'), alpha=0.6, size=0.6, pch=3, aes(color=type), stroke=0.5) +
    ggplot2::scale_color_manual(values=c(CL="#F8766D")) +
    ggplot2::scale_fill_manual(values=c(tumor="#00BFC4")) +
    ggplot2::xlab('UMAP 1') + ggplot2::ylab("UMAP 2") +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position='bottom',
          text = ggplot2::element_text(size=8),
          axis.text = ggplot2::element_text(size=6),
          axis.title = ggplot2::element_text(size=8),
          legend.margin =ggplot2::margin(0,0,0,0), 
          legend.box.margin=ggplot2::margin(-10,-30,-10,-30),
          axis.line = ggplot2::element_line(size = .3))
  
  return(uncorrected_combined_type_plot)
  
}

# figure 1c and Supplementary figure 2a
# cPCA_values included in the figshare: https://figshare.com/articles/Celligner_data/11965269
make_cPC_eigenspectrum <- function(cPCA_values) {
  cPCA_eigenvalues_plot <- data.frame(eval =cPCA_values$cPCA_values,
                                      rank = seq(length(cPCA_values$cPCA_values))) %>%
    ggplot2::ggplot(ggplot2::aes(rank, eval)) +
    ggplot2::geom_point(size=1) +
    ggplot2::geom_hline(yintercept = 0, linetype = 'dashed') +
    ggplot2::ylab('Eigenvalue') + 
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text=ggplot2::element_text(size=6), 
          axis.title.x = ggplot2::element_text(size=6), 
          axis.title.y = ggplot2::element_text(size=6),
          axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
  
  return(cPCA_eigenvalues_plot)
}



# used as input for 1d and 1e
# input files: 
# TCGA_mat source: https://xenabrowser.net/datapages/?dataset=TumorCompendium_v10_PolyA_hugo_log2tpm_58581genes_2019-07-25.tsv&host=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
# CCLE_mat source: depmap.org DepMap Public 19Q4 CCLE_expression_full.csv
get_current_eigenvector <- function(TCGA_mat, CCLE_mat, alignment, cPCA_vectors, cur_eig_vec = 1) {
  tcga_ccle <- rbind(TCGA_mat, CCLE_mat)
  cur_eig <- cPCA_vectors[,cur_eig_vec] %>% set_names(rownames(cPCA_vectors))
  cur_eig_proj <- get_signature_proj(cur_eig, tcga_ccle, alignment[rownames(tcga_ccle),])
  flip_sign <- choose_eig_sign(cur_eig_proj)
  cur_eig <- cur_eig * flip_sign
  
  return(cur_eig)
}

# figure 1d
make_cPC_purity_plot <- function(TCGA_mat, TCGA_ann, cur_eig) {
  dproj <- lm(t(TCGA_mat) ~ cur_eig)$coefficients[2,]
  df <- TCGA_ann %>% 
    dplyr::inner_join(data.frame(proj = dproj, sampleID = rownames(TCGA_mat)), by = 'sampleID')
  
  cPCA_purity_cor_plot <- ggplot2::ggplot(df %>% dplyr::filter(!is.na(purity)), 
                                          ggplot2::aes(proj, purity)) + 
    ggplot2::geom_point(size=0.5) + 
    ggplot2::ylab('Tumor purity estimate') + 
    ggplot2::xlab('Eigenvector projection') +
    ggpubr::stat_cor(label.y.npc = 'bottom', size=2) +
    ggplot2::geom_smooth(method = 'lm') +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text=ggplot2::element_text(size=6),
          axis.title.x = ggplot2::element_text(size=6),
          axis.title.y = ggplot2::element_text(size=6)) +
    xlab("cPC2 projection")
  
  return(cPCA_purity_cor_plot)
  
}

# used to run analysis for 1e and Supplementary figures 2b-e
# gsc data used is MSigDB v6.2 downloaded from www.gsea-msigdb.org/gsea/msigdb
run_cPCA_GSEA <- function(cur_eig, gene_stats, gsc_data) {
  df <- data.frame(vec = cur_eig, Gene = rownames(cPCA_vectors)) %>% 
    dplyr::left_join(gene_stats, by = "Gene")
  gene_stat <- with(df %>% dplyr::filter(!is.na(Symbol)), vec %>% set_names(Symbol))
  cPCA_GSEA <- run_fGSEA(gsc = gsc_data$GO_biological_process,
                               gene_stat = gene_stat,
                               nperm = 1e5,
                               perm_type = 'gene') %>% 
    dplyr::arrange(dplyr::desc(NES)) %>% 
    dplyr::select(-leadingEdge)
  
  return(cPCA_GSEA)
}

# GSEA table in figure 1e
make_cPC_GSEA_table <- function(cPCA_GSEA, filename) {
  cPCA_GSEA_data <- rbind.data.frame(cPCA_GSEA %>% dplyr::arrange(dplyr::desc(NES)) %>% head(5))
  cPCA_GSEA_data$pathway <- factor(cPCA_GSEA_data$pathway, levels = cPCA_GSEA_data$pathway) 
  cPCA_GSEA_data$pathway <- sub("^([^_]*_[^_]*_[^_]*_[^_]*_[^_]*).*", "\\1", cPCA_GSEA_data$pathway)
  cPCA_GSEA_data$pathway <- factor(cPCA_GSEA_data$pathway, levels = unique(cPCA_GSEA_data$pathway))
  cPCA_GSEA_data <- as.data.frame(cPCA_GSEA_data)
  
  cPCA_GSEA_data$`adjusted pval` <- signif(cPCA_GSEA_data$padj, 3)
  cPCA_GSEA_data$NES <- signif(cPCA_GSEA_data$NES, 3)
  
  pdf(filename, height=0.3*nrow(cPCA_GSEA_data), width=3.6)
  grid.table(cPCA_GSEA_data[,c('pathway', 'adjusted pval', 'NES')], theme=table_theme, rows=NULL)
  dev.off()
  
}
