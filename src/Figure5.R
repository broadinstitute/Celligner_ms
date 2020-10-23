library(here)
library(magrittr)
library(tidyverse)
source(here::here('src', 'analysis_helpers.R'))
source(here::here('src', 'global_params.R'))

# use the Celligner_info file for the alignment matrix used within these methods

# figure 5a
plot_undifferentiated_cluster <- function(alignment, cluster_selected = c(12,28)) {
  type_in_cluster <- rep('tumor_out', nrow(alignment))
  type_in_cluster[which(alignment$type == 'tumor' & alignment$cluster %in% cluster_selected)] <- 'tumor_in'
  type_in_cluster[which(alignment$type == 'CL' & alignment$cluster %in% cluster_selected)] <- 'CL_in'
  type_in_cluster[which(alignment$type == 'CL' & !alignment$cluster %in% cluster_selected)] <- 'CL_out'
  alignment$type_in_cluster <- type_in_cluster
  
  undifferentiated_cluster_plot <- ggplot2::ggplot(alignment, ggplot2::aes(UMAP_1, UMAP_2, fill=cluster  %in% cluster_selected , size=type, alpha=type_in_cluster)) +
    ggplot2::geom_point(pch=21, color='white') +
    ggplot2::geom_point(data = dplyr::filter(alignment, type=='CL'), color='black', pch=21) +
    ggplot2::scale_size_manual(values=c(`CL`=1.5, `tumor`=0.75)) + 
    ggplot2::scale_alpha_manual(values=c(`tumor_in`=0.8, `CL_in`=0.7, `tumor_out`=0.5, `CL_out`=0.2)) + 
    ggplot2::scale_fill_manual(values=c(`TRUE`='#EE4B33', `FALSE`='#999999')) +
    ggplot2::theme_classic() + 
    ggplot2::theme(legend.position = "none", 
          text=ggplot2::element_text(size=6),
          axis.text=ggplot2::element_text(size=6)) +
    ggplot2::xlab("UMAP 1") +
    ggplot2::ylab("UMAP 2")
  
  return(undifferentiated_cluster_plot)
}

# figure 5b
undifferentiated_cluster_CL_composition <- function(alignment, cluster_selected=c(12,28)) {
  CL_count_data <- dplyr::filter(alignment, type=='CL')
  CL_count_data$lineage <- gsub("engineered_","", CL_count_data$lineage)
  CL_count_data$lineage <- factor(CL_count_data$lineage, levels = unique(CL_count_data$lineage))
  CL_cluster_composition_data <- dplyr::filter(CL_count_data, cluster %in% cluster_selected)$lineage %>%
    table() %>%
    as.data.frame()
  CL_counts <- dplyr::filter(CL_count_data, !cluster %in% cluster_selected)$lineage %>% 
    table() %>%
    as.data.frame()
  colnames(CL_cluster_composition_data) <- c("lineage", "in cluster")
  colnames(CL_counts) <- c("lineage", 'other')
  CL_cluster_composition_data <- merge(CL_cluster_composition_data, CL_counts, by='lineage')
  CL_cluster_composition_data <- reshape2::melt(CL_cluster_composition_data)
  CL_cluster_composition_data$lineage <- gsub("_", " ", CL_cluster_composition_data$lineage)
  CL_cluster_composition_data$variable <- factor(CL_cluster_composition_data$variable, levels = c('other', "in cluster"))
  
  CL_cluster_composition_plot <- ggplot2::ggplot(data=CL_cluster_composition_data, ggplot2::aes(x=lineage, y=value, fill=variable)) +
    ggplot2::geom_bar(stat="identity") +
    ggplot2::theme_classic() + 
    ggplot2::theme(legend.position = 'left',
          text = ggplot2::element_text(size=8),
          axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size=8), 
          axis.title.x = ggplot2::element_text(vjust=10)) +
    ggplot2::scale_fill_manual(values=c(`in cluster`='#EE4B33', `other`= '#BFBFBF')) + 
    ggplot2::ylab("number of samples") +
    ggplot2::scale_y_continuous(expand = c(0, 0)) + 
    ggplot2::labs(fill='')
  
  return(CL_cluster_composition_plot)
}

# figure 5c
# The CRISPR data used can be found on depmap.org. The data used is DepMap Public 19Q4 Achilles_gene_effect.csv
SOX10_dependency_expression_plot <- function(CCLE_mat, CRISPR, alignment) {
  int_cell_lines <- intersect(rownames(CCLE_mat), rownames(CRISPR))
  cur_cell_lines <- dplyr::filter(alignment, type=='CL' & lineage == 'skin')$sampleID

  mel_cluster <- dplyr::filter(alignment, lineage=='skin')$cluster %>% 
    table() %>%
    as.data.frame() %>%
    dplyr::arrange(dplyr::desc(Freq)) %>%
    head(2) %>% 
    .[['.']] %>%
    as.character()
  
  sox10_dep_exp <- cbind.data.frame(CCLE_mat[int_cell_lines,  "ENSG00000100146"], 
                                  CRISPR[int_cell_lines, 'SOX10'], 
                                  int_cell_lines %in% cur_cell_lines)
  
  colnames(sox10_dep_exp) <- c("expression", "gene_effect", "cur_lineage")
  sox10_dep_exp$`skin clusters` <- "other samples"
  sox10_dep_exp[dplyr::filter(alignment, sampleID %in% int_cell_lines & cluster %in% mel_cluster)$sampleID,'skin clusters'] <- "melanoma"
  sox10_dep_exp[dplyr::filter(alignment, cluster %in% c(12,28) & lineage=='skin' & sampleID %in% int_cell_lines)$sampleID,'skin clusters'] <- "undifferentiated"
  sox10_dep_exp[dplyr::filter(alignment, !cluster %in% c(12,28, mel_cluster) & lineage=='skin' & sampleID %in% int_cell_lines)$sampleID,'skin clusters'] <- "other"

  sox10_dep_exp$`skin clusters` <- as.factor(sox10_dep_exp$`skin clusters`)

  sox10_dep_exp_plot <- ggplot2::ggplot(sox10_dep_exp, ggplot2::aes(expression, gene_effect)) +
    ggplot2::geom_point(alpha=0.5, fill = '#BFBFBF', color='white', pch=21, size=1.5) +
    ggplot2::geom_point(data=dplyr::filter(sox10_dep_exp, `skin clusters`!="other samples"), 
                        ggplot2::aes(fill=`skin clusters`), alpha=0.7, pch=21, color='gray10', size=1.75) +
    ggplot2::theme_classic() +
    ggplot2::xlab("SOX10 expression") +
    ggplot2::ylab("SOX10 dependency") +
    ggplot2::theme(legend.position = 'bottom',
          text=ggplot2::element_text(size=7),
          legend.margin=ggplot2::margin(0,0,0,0),
          legend.box.margin=ggplot2::margin(-10,-10,-10,-10),
          legend.justification = "left",
          legend.spacing.y = ggplot2::unit(.0001, 'cm'),
          plot.margin = ggplot2::unit(c(0.2,1,0.2,0.2), "cm")) +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow=2, byrow=TRUE, title="skin \nclusters"))

  return(sox10_dep_exp_plot)
}

# figure 5d
# The CRISPR data used can be found on depmap.org. The data used is DepMap Public 19Q4 Achilles_gene_effect.csv
dependency_expression_plot <- function(CCLE_mat, CRISPR, cur_gene_symbol = 'HNF4A', alignment,
                                             cur_gene_ensembl = 'ENSG00000101076', cur_type = 'liver') {
  
  int_cell_lines <- intersect(rownames(CCLE_mat), rownames(CRISPR))
  
  
  other_dep_exp <- cbind.data.frame(CCLE_mat[int_cell_lines,cur_gene_ensembl],
                                    CRISPR[int_cell_lines, cur_gene_symbol],
                                    ifelse(int_cell_lines %in% dplyr::filter(alignment, cluster %in% c(12,28) &
                                                                               lineage == cur_type)$sampleID,
                                           'undifferentiated', cur_type),
                                    int_cell_lines %in% dplyr::filter(alignment, lineage == cur_type)$sampleID)
  
  colnames(other_dep_exp) <- c("expression", "gene_effect", "clusters", "cur_lineage")
  
  dep_exp_plot <- ggplot2::ggplot(other_dep_exp, ggplot2::aes(expression, gene_effect)) +
    ggplot2::geom_point(alpha=0.5, fill = '#BFBFBF', color='white', size=1.5, pch=21) +
    ggplot2::geom_point(data=dplyr::filter(other_dep_exp, cur_lineage==T), 
                        ggplot2::aes(fill=`clusters`), alpha=0.7, pch=21, color='gray10', size=1.75) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = 'bottom',
          text=ggplot2::element_text(size=7),
          legend.margin=ggplot2::margin(0,0,0,0), 
          legend.box.margin=ggplot2::margin(-10,-10,-10,-10),
          legend.justification = "left", legend.spacing.y = ggplot2::unit(.0001, 'cm')) +
    ggplot2::xlab(paste(cur_gene_symbol, "expression")) +
    ggplot2::ylab(paste(cur_gene_symbol, "dependency")) +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow=2, byrow=TRUE))
  
  return(dep_exp_plot)
  
}

# analysis used as input for figure 5e
# CCLE_count data source is depmap.org, DepMap Public 19Q4 CCLE_RNAseq_reads.csv
# genes were subset to the genes used elsewhere (a set of 19,188 genes 
# that are common to both datasets and are not non-coding RNA or pseudogenes)
undifferentiated_DE_analysis <- function(aligment, clusters = c(12,28), CCLE_count, gene_stats, CCLE_mat) {
  colnames(CCLE_count) <- stringr::str_match(colnames(CCLE_count), '\\((.+)\\)')[,2]
  CCLE_count <- CCLE_count[,colnames(CCLE_mat)]
  rownames(gene_stats) <- gene_stats$Gene
  
  vec <- (dplyr::filter(alignment, type=='CL')$cluster %in% clusters) %>%
    set_names(dplyr::filter(alignment, type=='CL')$sampleID)
  
  covars <- dplyr::filter(alignment, type=='CL') %>% 
    dplyr::select(lineage)
  rownames(covars) <- rownames(dplyr::filter(alignment, type=='CL'))
  
  cell_line_info <- dplyr::filter(alignment, type=='CL')
  gene_info <- gene_stats
  rownames(gene_info) <- gene_info$Gene
  gene_info <- gene_info[colnames(CCLE_count),]
  
  d <- edgeR::DGEList(counts= t(CCLE_count), group=vec, samples=cell_line_info, genes=gene_info)
  cpm <- edgeR::cpm(d)
  keep.exprs <- rowSums(cpm>1)>=10 
  d <- d[keep.exprs,, keep.lib.sizes=FALSE] 
  d <- edgeR::calcNormFactors(d, method='TMM')
  lcpm <- edgeR::cpm(d, log=TRUE)
  
  design <- model.matrix(~0 + group + lineage, d$samples)
  colnames(design) <- gsub("lineage", "", colnames(design))
  contr.matrix <- limma::makeContrasts(isUndifferentiated = groupTRUE-groupFALSE, levels=colnames(design))
  vfit <- limma::lmFit(lcpm, design)
  vfit <- limma::contrasts.fit(vfit, contrasts=contr.matrix)
  efit <- limma::eBayes(vfit, trend = T)
  de_table <- limma::topTable(efit, coef='isUndifferentiated', n=Inf)
  de_table$Gene <- rownames(de_table)
  de_table <- merge(de_table, gene_stats)
  
  return(de_table)
  
}

# figure 5e
undifferentiated_DE_plot <- function(de_table) {
  de_table$EMT <- ifelse(de_table$Symbol %in% EMT_genes, 'EMT gene', 'NA')
  de_table$`EMT genes` <- ifelse(de_table$logFC > 0, 'up regulated', 'down regulated')
  
  undiff_de_analysis_plot <- ggplot2::ggplot(de_table, ggplot2::aes(logFC, -log10(P.Value))) + 
    ggplot2::geom_point(fill = '#BFBFBF', color='white', pch=21, stroke=0.2, alpha=0.5) + 
    ggplot2::geom_point(data=dplyr::filter(de_table, EMT != 'NA'), ggplot2::aes(logFC, -log10(P.Value), fill=`EMT genes`, color=`EMT genes`), pch=21, color='white', size=1.75) +
    ggplot2::scale_fill_manual(values=c(`down regulated`='#1F295C', `up regulated`='#90191C')) +
    ggplot2::scale_color_manual(values=c(`down regulated`='#1F295C', `up regulated`='#90191C')) +
    ggrepel::geom_text_repel(data = dplyr::filter(de_table, EMT != 'NA') %>% 
                               dplyr::arrange(dplyr::desc(abs(logFC))) %>% 
                               head(5), 
                             ggplot2::aes(label = Symbol, color=`EMT genes`), show.legend = FALSE, size=3) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = 'bottom',
          text=ggplot2::element_text(size=8),
          axis.text = ggplot2::element_text(size=6),
          legend.margin=ggplot2::margin(0,0,0,0),
          legend.box.margin=ggplot2::margin(-9,-9,-9,-9)) +
    ggplot2::ylab("-log10(p.value)")
  
  return(undiff_de_analysis_plot)
  
}

# figure 5f
# gsc data used is MSigDB v6.2 downloaded from www.gsea-msigdb.org/gsea/msigdb
GSEA_on_undifferentiated_DE <- function(de_table, gsc_data, filepath) {
  gene_stat <- with(de_table %>% dplyr::filter(!is.na(Symbol)), logFC %>% set_names(Symbol))
  cur_GSEA <- run_fGSEA(gsc = gsc_data$hallmark,
                              gene_stat = gene_stat,
                              nperm = 1e5,
                              perm_type = 'gene') %>%
    dplyr::arrange(dplyr::desc(abs(NES))) %>%
    dplyr::select(-leadingEdge)
  
  GSEA_data <- rbind.data.frame(cur_GSEA %>% dplyr::arrange(dplyr::desc(NES)) %>% head(5))
  GSEA_data$pathway <- gsub('HALLMARK_', "", GSEA_data$pathway)
  GSEA_data$pathway <- factor(GSEA_data$pathway, levels = GSEA_data$pathway) 
  
  GSEA_data$`adjusted pval` <- signif(GSEA_data$padj, 3)
  GSEA_data$NES <- signif(GSEA_data$NES, 3)
  
  pdf(filepath, height=2, width=3.5)
  grid.table(GSEA_data[,c('pathway', 'adjusted pval', 'NES')], theme=table_theme, rows=NULL)
  dev.off()
}

# analysis used as input for 5g
# secondary.screen.replicate.collapsed.logfold.change, secondary.screen.replicate.collapsed.treatment.info,
# and sample.info files can be found on depmap.org
run_differential_drug_sensitivities <- function(alignment,
                                                secondary.screen.replicate.collapsed.logfold.change, 
                                                secondary.screen.replicate.collapsed.treatment.info,
                                                sample.info, clusters = c(12,28)) {
  # switch to CCLE names
  new_rownames <- plyr::llply(rownames(secondary.screen.replicate.collapsed.logfold.change), function(x) {
    sample.info$CCLE_Name[sample.info$DepMap_ID == x]
  }) %>% 
    as.character()
  
  rownames(secondary.screen.replicate.collapsed.logfold.change) <- new_rownames
  
  # take dose-level data at 2.5uM
  prism_metadata_one_dose <- secondary.screen.replicate.collapsed.treatment.info[
    (secondary.screen.replicate.collapsed.treatment.info$dose == 2.5),
    ]
  
  # get shared cell lines for celligner and prism, and format names
  prism_shared_cell_lines <- rownames(secondary.screen.replicate.collapsed.logfold.change)[rownames(secondary.screen.replicate.collapsed.logfold.change) %in% alignment$sampleID_CCLE_Name]
  prism_data <- secondary.screen.replicate.collapsed.logfold.change[prism_shared_cell_lines, prism_metadata_one_dose$column_name] %>% t()
  prism_data_with_names <- merge(prism_data, secondary.screen.replicate.collapsed.treatment.info[,c("column_name", "name")], by.x = "row.names", by.y = "column_name")
  prism_data_with_names %<>% tibble::column_to_rownames("Row.names")
  
  # collapse by pert
  prism_data_collapsed <- prism_data_with_names %>% 
    dplyr::group_by(name) %>% 
    dplyr::summarise_all(mean, na.rm = T)
  prism_data_formatted <- prism_data_collapsed %>%
    tibble::column_to_rownames("name") %>% 
    t() %>% 
    as.data.frame() %>% 
    as.matrix()
  
  # set NAs and infs to 0
  k <- which(is.na(prism_data_formatted), arr.ind=TRUE)
  prism_data_formatted[k] <- 0
  
  k <- which(is.infinite(prism_data_formatted %>% as.matrix()), arr.ind=TRUE)
  prism_data_formatted[k] <- rowMeans(prism_data_formatted*is.finite(prism_data_formatted), na.rm = T)[k[,1]]
  
  still_na_rows <- names(which(rowSums(is.na(prism_data_formatted)) > 0))
  prism_data_formatted_finally <- prism_data_formatted[!(rownames(prism_data_formatted) %in% still_na_rows),]
  
  # compounds and cell lines to test
  cpds_to_test <- colnames(prism_data_formatted_finally)
  shared_cell_lines <- intersect(rownames(prism_data_formatted_finally), alignment$sampleID_CCLE_Name)
  
  rownames(alignment) <- alignment$sampleID_CCLE_Name
  is_undifferentiated <- alignment[shared_cell_lines,'cluster'] %in% clusters
  # fit linear model
  lmstats_out_chemical_prism <- run_lm_stats_limma(mat = prism_data_formatted_finally[shared_cell_lines, cpds_to_test], 
                                                         vec = is_undifferentiated, 
                                                         covars = alignment[match(shared_cell_lines,alignment$sampleID_CCLE_Name),
                                                                            'lineage'])
  
 return(lmstats_out_chemical_prism)
}

# figure 5g
# repurposing.export can be found on clue.io/data labeled as
# Repurposing related drug annotations
plot_differential_drug_sensitivities <- function(lmstats_out_chemical_prism,
                                                 repurposing.export) {
# add compound metadata
  lmstats_out_chemical_prism_with_annots <- merge(lmstats_out_chemical_prism, repurposing.export[,c("Name", "MOA", "Target")], by.x = "Gene", by.y = "Name")

  lmstats_out_chemical_prism_with_annots$moa_trimmed <- plyr::llply(lmstats_out_chemical_prism_with_annots$MOA, function(x) {
    if (grepl("EGFR inhibitor", x)) {
      return("EGFR inhibitor")
    } else if (grepl("tubulin polymerization inhibitor", x)) {
      return("Tubulin polymerization inhibitor")
    } else {
      return(NA)
    }
  }) %>% as.character()
  
  undifferentiated_differential_drug <- lmstats_out_chemical_prism_with_annots %>% 
    ggplot2::ggplot(aes(EffectSize, -log10(p.value))) + 
    ggplot2::geom_point(fill = "#BFBFBF", color = 'white', pch=21, stroke=0.2, alpha=0.5) + 
    ggplot2::geom_point(data = dplyr::filter(lmstats_out_chemical_prism_with_annots, moa_trimmed != 'NA'), 
               aes(EffectSize, -log10(p.value), fill = moa_trimmed), size=1.75, color='white', pch=21) +
    ggplot2::scale_fill_manual(values=c(`EGFR inhibitor`="#90191C", `Tubulin polymerization inhibitor`="#1F295C")) + 
    ggplot2::scale_color_manual(values=c(`EGFR inhibitor`="#90191C", `Tubulin polymerization inhibitor`="#1F295C")) + 
    ggrepel::geom_text_repel(data = rbind(lmstats_out_chemical_prism_with_annots[grepl("EGFR", lmstats_out_chemical_prism_with_annots$MOA),] %>%
                                            dplyr::arrange(-EffectSize) %>% head(3),
                                          lmstats_out_chemical_prism_with_annots[grepl("tubulin polymerization inhibitor", lmstats_out_chemical_prism_with_annots$MOA),] %>%
                                            dplyr::arrange(EffectSize) %>% head(3)),
                             ggplot2::aes(label = Gene, color=moa_trimmed), show.legend = FALSE, size=3) + 
    ggplot2::theme_classic() +
    ggplot2::labs(fill = "Drug class") + 
    ggplot2::guides(fill = ggplot2::guide_legend(nrow=2,byrow=TRUE), alpha=FALSE) +
    ggplot2::theme(legend.position = "bottom",
          text=ggplot2::element_text(size=8),
          axis.text = ggplot2::element_text(size=6),
          axis.title = ggplot2::element_text(size=6),  
          legend.margin =ggplot2::margin(0,0,0,0),
          legend.box.margin=ggplot2::margin(-10,-10,-10,-10))
  
  return(undifferentiated_differential_drug)
}

# analysis used as input for 5h
# The CRISPR data used can be found on depmap.org. The data used is DepMap Public 19Q4 Achilles_gene_effect.csv
run_differential_dependency_analysis <- function(CRISPR, alignment, clusters = c(12,28)) {
  vec <- (dplyr::filter(alignment, type=='CL')$cluster %in% clusters) %>%
    set_names(dplyr::filter(alignment, type=='CL')$sampleID)
  covars <- dplyr::filter(alignment, type=='CL') %>% 
    dplyr::select(lineage)
  rownames(covars) <- dplyr::filter(alignment, type=='CL')$sampleID

  common_CLs <- intersect(rownames(CRISPR), names(vec))
  DE_res <- run_lm_stats_limma(CRISPR[common_CLs,],
                                     vec[common_CLs],
                                     covars = covars[common_CLs,,drop=F], 
                                     limma_trend = FALSE) 
  
  return(DE_res)
}

# figure 5h
plot_differentiated_dependencies <- function(DE_res) {
  more_dependent <- DE_res %>% dplyr::arrange(EffectSize) %>% head(4) %>% .[['Gene']]
  less_dependent <- 'EGFR'
  DE_res$label <- 'NA'
  DE_res[which(DE_res$Gene %in% more_dependent),'label'] <- 'more dependent'
  DE_res[which(DE_res$Gene %in% less_dependent),'label'] <- 'less dependent'
  DE_res$label <- factor(DE_res$label, levels = c('more dependent', 'less dependent', 'NA'))
  
  CL_cluster_dif_dep_plot <- ggplot2::ggplot(DE_res, ggplot2::aes(EffectSize, -log10(p.value))) + 
    ggplot2::geom_point(fill ='#BFBFBF', color='white', pch=21, stroke=0.1) + 
    ggplot2::geom_point(data=dplyr::filter(DE_res, label != 'NA'), ggplot2::aes(EffectSize, -log10(p.value), fill = label), color='white', size=1.75, pch=21) +
    ggrepel::geom_text_repel(data = dplyr::filter(DE_res, 
                                           label != 'NA'),
                             ggplot2::aes(label = Gene, color=label),
                             size = 3, show.legend = FALSE) + 
    ggplot2::scale_color_manual(values=c(`more dependent`='#1F295C', `less dependent`='#90191C')) +
    ggplot2::scale_fill_manual(values=c(`more dependent`='#1F295C', `less dependent`='#90191C')) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = 'bottom',
          text=ggplot2::element_text(size=8),
          axis.text = ggplot2::element_text(size=6),
          axis.title = ggplot2::element_text(size=6),
          legend.title = ggplot2::element_blank(),
          legend.margin=ggplot2::margin(0,0,0,0),
          legend.box.margin=ggplot2::margin(-10,-10,-10,-10))
  
  return(CL_cluster_dif_dep_plot)
  
}
