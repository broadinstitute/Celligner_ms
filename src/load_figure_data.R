library(magrittr)
# TCGA_mat source: https://xenabrowser.net/datapages/?dataset=TumorCompendium_v10_PolyA_hugo_log2tpm_58581genes_2019-07-25.tsv&host=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
load_TCGA_mat <- function(data_dir, tumor_file='TCGA_mat.tsv') {
  TCGA_mat <-  readr::read_tsv(file.path(data_dir, tumor_file)) %>% 
    as.data.frame() %>%
    tibble::column_to_rownames('Gene') %>%
    as.matrix() %>% 
    t()
  
  return(TCGA_mat)
}

# CCLE_mat source: depmap.org DepMap Public 19Q4 CCLE_expression_full.csv
load_CCLE_mat <- function(data_dir, cell_line_file = 'CCLE_mat.csv') {
  CCLE_mat <-  readr::read_csv(file.path(data_dir, cell_line_file)) %>% 
    as.data.frame() %>%
    tibble::column_to_rownames('X1') %>%
    as.matrix()
  
  colnames(CCLE_mat) <- stringr::str_match(colnames(CCLE_mat), '\\((.+)\\)')[,2]
  
  return(CCLE_mat)
  
}

# Celligner_info file available on figshare: https://figshare.com/articles/Celligner_data/11965269
load_alignment <- function(data_dir, filename = 'Celligner_info.csv') {
  alignment <- data.table::fread(file.path(data_dir, filename)) %>%
    as.data.frame()
  
  rownames(alignment) <- alignment$sampleID
  
  return(alignment)
}

load_CCLE_ann <- function(data_dir, filename = 'Celligner_info.csv') {
  CCLE_ann <- data.table::fread(file.path(data_dir, filename)) %>%
    as.data.frame()
  
  CCLE_ann <- dplyr::filter(CCLE_ann, type=='CL') %>%
    dplyr::select(-UMAP_1, -UMAP_2, -cluster, -uncorrected_tumor_UMAP_1, -uncorrected_tumor_UMAP_2, -uncorrected_tumor_cluster) %>%
    dplyr::rename(
      UMAP_1 = uncorrected_CL_UMAP_1,
      UMAP_2 = uncorrected_CL_UMAP_2,
      cluster = uncorrected_CL_cluster
    ) 
  rownames(CCLE_ann) <- CCLE_ann$sampleID
  
  return(CCLE_ann)

}

load_TCGA_ann <- function(data_dir, filename = 'Celligner_info.csv') {
  TCGA_ann <- data.table::fread(file.path(data_dir, filename)) %>%
    as.data.frame()
  
  TCGA_ann <- dplyr::filter(TCGA_ann, type=='tumor') %>%
    dplyr::select(-UMAP_1, -UMAP_2, -cluster, -uncorrected_CL_UMAP_1, -uncorrected_CL_UMAP_2, -uncorrected_CL_cluster) %>%
    dplyr::rename(
      UMAP_1 = uncorrected_tumor_UMAP_1,
      UMAP_2 = uncorrected_tumor_UMAP_2,
      cluster = uncorrected_tumor_cluster
    ) 
  rownames(TCGA_ann) <- TCGA_ann$sampleID
  
  return(TCGA_ann)
  
}

load_cPCA_values <- function(data_dir, cPCA_values = 'cPCA_values.csv') {
  cPCA_values <- read_csv(file.path(data_dir, cPCA_values))
  
  return(cPCA_values)
}

load_cPCA_vectors <- function(data_dir, cPCs = 'cPCs.csv') {
  cPCA_vectors <- readr::read_csv(file.path(data_dir, cPCs)) %>%
    as.data.frame() %>%
    tibble::column_to_rownames('X1')
  
  return(cPCA_vectors)
}
  
load_tumor_CL_cor <- function(data_dir, cor_file = 'tumor_CL_cor.csv') {
  tumor_CL_cor <- readr::read_csv(file.path(data_dir, cor_file)) %>%
    tibble::column_to_rownames('X1')
  
  return(tumor_CL_cor)
}

load_data <- function(data_dir, tumor_file = 'TCGA_mat.tsv', cell_line_file = 'CCLE_mat.csv', 
                      annotation_file = 'Celligner_info.csv', hgnc_file = "hgnc_complete_set_7.24.2018.txt") {
  
  TCGA_mat <- load_TCGA_mat(data_dir, tumor_file)
  CCLE_mat <- load_CCLE_mat(data_dir, cell_line_file)
  hgnc.complete.set <- data.table::fread(file.path(data_dir, hgnc_file)) %>% as.data.frame()
  
  common_genes <- intersect(colnames(TCGA_mat), hgnc.complete.set$symbol)
  TCGA_mat <- TCGA_mat[,common_genes]
  hgnc.complete.set <- dplyr::filter(hgnc.complete.set, symbol %in% common_genes)
  hgnc.complete.set <- hgnc.complete.set[!duplicated(hgnc.complete.set$symbol),]
  rownames(hgnc.complete.set) <- hgnc.complete.set$symbol
  hgnc.complete.set <- hgnc.complete.set[common_genes,]
  colnames(TCGA_mat) <- hgnc.complete.set$ensembl_gene_id
  
  
  if(is.null(annotation_file) | !file.exists(file.path(data_dir, annotation_file))) {
    ann <- data.frame(sampleID = c(rownames(TCGA_mat), rownames(CCLE_mat)),
                      lineage = NA,
                      subtype = NA,
                      type = c(rep('tumor', nrow(TCGA_mat)), rep('CL', nrow(CCLE_mat))))
    ann$`Primary/Metastasis` <- NA
  } else {
    ann <- data.table::fread(file.path(data_dir, annotation_file)) %>% as.data.frame()
    if('UMAP_1' %in% colnames(ann)) {
      ann <- ann %>% 
        dplyr::select(-UMAP_1)
    }
    if('UMAP_2' %in% colnames(ann)) {
      ann <- ann %>% 
        dplyr::select(-UMAP_2)
    }
    if('cluster' %in% colnames(ann)) {
      ann <- ann %>% 
        dplyr::select(-cluster)
    }
  }
  
  TCGA_ann <- dplyr::filter(ann, type=='tumor')
  CCLE_ann <- dplyr::filter(ann, type=='CL')
  
  func_genes <- dplyr::filter(hgnc.complete.set, !locus_group %in% c('non-coding RNA', 'pseudogene'))$ensembl_gene_id
  genes_used <- intersect(colnames(TCGA_mat), colnames(CCLE_mat))
  genes_used <- intersect(genes_used, func_genes)
  
  TCGA_mat <- TCGA_mat[,genes_used]
  CCLE_mat <- CCLE_mat[,genes_used]
  
  return(list(TCGA_mat = TCGA_mat, TCGA_ann = TCGA_ann, CCLE_mat = CCLE_mat, CCLE_ann = CCLE_ann))
}


