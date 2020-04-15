library(here)
library(magrittr)
library(tidyverse)
source(here::here('src', 'analysis_helpers.R'))
source(here::here('src', 'global_params.R'))

# Differentially expressed genes --------------------------------------------------------------------

# Estimate linear-model stats for a matrix of data with respect to a group of phenotype variables
# using limma with empirical Bayes moderated F-stats for p-values
run_lm_stats_limma_group <- function (mat, phenos, covars = NULL, weights = NULL, target_type = "Gene", 
          limma_trend = FALSE) 
{
  require(limma)
  require(magrittr)
  require(tibble)
  require(plyr)
  require(dplyr)
  udata <- rownames(mat) %>% intersect(rownames(phenos))
  if (!is.null(covars)) {
    udata %<>% intersect(rownames(covars))
  }
  form <- as.formula(paste("~", paste0(colnames(phenos), collapse = " + ")))
  design <- model.matrix(form, data = phenos[udata, , drop = F])
  if (!is.null(covars)) {
    covars <- data.frame(covars)
    form <- as.formula(paste("~", paste0(colnames(covars), 
                                         collapse = " + ")))
    Cdesign <- model.matrix(form, data = covars[udata, , 
                                                drop = F])
    Cdesign <- Cdesign[, setdiff(colnames(Cdesign), "(Intercept)"), 
                       drop = FALSE]
    stopifnot(length(intersect(colnames(Cdesign), colnames(design))) == 
                0)
    design %<>% cbind(Cdesign)
  }
  if (!is.null(weights)) {
    if (is.matrix(weights)) {
      weights <- t(weights[udata, ])
    }
    else {
      weights <- weights[udata]
    }
  }
  design <- design[, colSums(design) > 2, drop = FALSE]
  targ_coefs <- setdiff(colnames(design), "(Intercept)")
  fit <- limma::lmFit(t(mat[udata, ]), design, weights = weights)
  fit <- limma::eBayes(fit, trend = limma_trend)
  targ_coef <- which(colnames(design) %in% targ_coefs)
  results <- limma::topTable(fit, coef = targ_coef, number = Inf, 
                             sort.by = "F", genelist = colnames(mat))
  results %<>% tibble::rownames_to_column(var = target_type)
  results %<>% magrittr::set_colnames(revalue(colnames(.), c(AveExpr = "Avg", 
                                                   F = "F_stat", P.Value = "p.value", adj.P.Val = "q.value"))) %>% 
    na.omit() %>% dplyr::select(-ProbeID)
  return(results)
}

# cPCA --------------------------------------------------------------------

# run contrastive principal components analysis, first removing average cluster expression, to
# estimate the average intra-cluster covariance
# if pc_dims = NULL, all cPCs are calculated. Faster cPCA can be run by setting pc_dims to a
# value >=4 and approximating just those cPCs 
run_cPCA_analysis <- function(TCGA_dat, CCLE_dat, tumor_cluster_df, CL_cluster_df, pc_dims=NULL) {
  tumor_clust_avgs <- get_cluster_averages(TCGA_dat, tumor_cluster_df)
  CL_clust_avgs <- get_cluster_averages(CCLE_dat, CL_cluster_df)
  
  TCGA_subtype_ms <- TCGA_dat - tumor_clust_avgs[tumor_cluster_df$seurat_clusters,]
  CCLE_subtype_ms <- CCLE_dat - CL_clust_avgs[CL_cluster_df$seurat_clusters,]
  
  TCGA_cov <- cov(TCGA_subtype_ms)
  CCLE_cov <- cov(CCLE_subtype_ms)
  
  if(!is.null(pc_dims)) {
    cov_diff_eig <- irlba::prcomp_irlba(TCGA_cov - CCLE_cov, n = pc_dims)
  } else {
    cov_diff_eig <- eigen(TCGA_cov - CCLE_cov)
  }
  return(cov_diff_eig)
}

# calculate the average expression per cluster
get_cluster_averages <- function(mat, cluster_df) {
  n_clusts <- nlevels(cluster_df$seurat_clusters)
  clust_avgs <- matrix(NA, nrow = n_clusts, ncol = ncol(mat)) %>% 
    magrittr::set_colnames(colnames(mat)) %>% 
    magrittr::set_rownames(levels(cluster_df$seurat_clusters))
  for (ii in levels(cluster_df$seurat_clusters)) {
    clust_avgs[ii,] <- colMeans(mat[cluster_df$seurat_clusters == ii,], na.rm=T)
  }
  return(clust_avgs)
}

# MNN --------------------------------------------------------------------

# Modification of the scran::fastMNN (https://github.com/MarioniLab/scran)
# Allows for separate k values per dataset, and simplifies some of the IO and doesn't use PCA reduction
modified_mnnCorrect <- function(ref_mat, targ_mat, k1 = 20, k2 = 20, 
                            ndist = 3, subset_genes = NULL) {
  if (is.null(subset_genes)) {
    subset_genes <- colnames(ref_mat)
  }
  
  sets <- batchelor::findMutualNN(ref_mat[, subset_genes], 
                                  targ_mat[, subset_genes], 
                                  k1 = k2, k2 = k1, 
                                  BPPARAM = BiocParallel::SerialParam())
  mnn_pairs <- as.data.frame(sets) %>% 
    dplyr::mutate(ref_ID = rownames(ref_mat)[first],
           targ_ID = rownames(targ_mat)[second],
           pair = seq(nrow(.))) %>% 
    dplyr::select(-first, -second)
  
  # Estimate the overall batch vector.
  ave.out <- .average_correction(ref_mat, sets$first, targ_mat, sets$second)
  overall.batch <- colMeans(ave.out$averaged)
  
  #remove variation along the overall batch vector
  ref_mat <- .center_along_batch_vector(ref_mat, overall.batch)
  targ_mat <- .center_along_batch_vector(targ_mat, overall.batch)
  
  # Recompute correction vectors and apply them.
  re.ave.out <- .average_correction(ref_mat, sets$first, targ_mat, sets$second)
  targ_mat <- .tricube_weighted_correction(targ_mat, re.ave.out$averaged, re.ave.out$second, k=k2, ndist=ndist, subset_genes, BPPARAM=BiocParallel::SerialParam())
  
  final <- list(corrected = targ_mat, 
                pairs = mnn_pairs)
  return(final)
}

# Copied from dev version of scran (2018-10-28) with slight modifications as noted
#https://github.com/MarioniLab/scran
.average_correction <- function(refdata, mnn1, curdata, mnn2)
  # Computes correction vectors for each MNN pair, and then
  # averages them for each MNN-involved cell in the second batch.
{
  corvec <- refdata[mnn1,,drop=FALSE] - curdata[mnn2,,drop=FALSE]
  corvec <- rowsum(corvec, mnn2)
  npairs <- table(mnn2)
  stopifnot(identical(names(npairs), rownames(corvec)))
  corvec <- unname(corvec)/as.vector(npairs)
  list(averaged=corvec, second=as.integer(names(npairs)))
}


.center_along_batch_vector <- function(mat, batch.vec) 
  # Projecting along the batch vector, and shifting all cells to the center _within_ each batch.
  # This removes any variation along the overall batch vector within each matrix.
{
  batch.vec <- batch.vec/sqrt(sum(batch.vec^2))
  batch.loc <- as.vector(mat %*% batch.vec)
  central.loc <- mean(batch.loc)
  mat <- mat + outer(central.loc - batch.loc, batch.vec, FUN="*")
  return(mat)
}

#' @importFrom BiocNeighbors queryKNN
#' @importFrom BiocParallel SerialParam
.tricube_weighted_correction <- function(curdata, correction, in.mnn, k=20, ndist=3, subset_genes, BNPARAM=NULL, BPPARAM=BiocParallel::SerialParam())
  # Computing tricube-weighted correction vectors for individual cells,
  # using the nearest neighbouring cells _involved in MNN pairs_.
  # Modified to use FNN rather than queryKNN for nearest neighbor finding
{
  cur.uniq <- curdata[in.mnn,,drop=FALSE]
  safe.k <- min(k, nrow(cur.uniq))
  # closest <- queryKNN(query=curdata, X=cur.uniq, k=safe.k, BNPARAM=BNPARAM, BPPARAM=BPPARAM)
  closest <- FNN::get.knnx(cur.uniq[, subset_genes], query=curdata[, subset_genes], k=safe.k)
  # weighted.correction <- .compute_tricube_average(correction, closest$index, closest$distance, ndist=ndist)
  weighted.correction <- .compute_tricube_average(correction, closest$nn.index, closest$nn.dist, ndist=ndist)
  curdata + weighted.correction
}

.compute_tricube_average <- function(vals, indices, distances, bandwidth=NULL, ndist=3) 
  # Centralized function to compute tricube averages.
  # Bandwidth is set at 'ndist' times the median distance, if not specified.
{
  if (is.null(bandwidth)) {
    middle <- ceiling(ncol(indices)/2L)
    mid.dist <- distances[,middle]
    bandwidth <- mid.dist * ndist
  }
  bandwidth <- pmax(1e-8, bandwidth)
  
  rel.dist <- distances/bandwidth
  rel.dist[rel.dist > 1] <- 1 # don't use pmin(), as this destroys dimensions.
  tricube <- (1 - rel.dist^3)^3
  weight <- tricube/rowSums(tricube)
  
  output <- 0
  for (kdx in seq_len(ncol(indices))) {
    output <- output + vals[indices[,kdx],,drop=FALSE] * weight[,kdx]
  }
  
  if (is.null(dim(output))) {
    matrix(0, nrow(vals), ncol(vals))
  } else {
    output
  }
}

