library(here)
library(magrittr)
library(tidyverse)
source(here::here('src', 'global_params.R'))

get_signature_proj <- function(eig_vec, mat, ann = NULL) {
  proj <- lm(t(mat) ~ eig_vec)$coefficients[2,]
  dproj_df <- data.frame(proj = proj,
                         sampleID = rownames(mat)) 
  if (!is.null(ann)) {
    dproj_df %<>% dplyr::left_join(ann, by = 'sampleID')
  }
  return(dproj_df)
}

choose_eig_sign <- function(df, positive_for = 'tumor') {
  if (positive_for != 'tumor') {
    flip_sign <- -1
  } else {
    flip_sign <- 1
  }
  greater_type <- df %>% 
    dplyr::filter(!is.na(type)) %>%
    dplyr::group_by(type) %>% 
    dplyr::summarise(avg_proj = mean(proj)) %>% 
    dplyr::arrange(dplyr::desc(avg_proj)) %>% 
    head(1) %>% 
    .[['type']]
  if (greater_type == 'CL') {
    flip_sign  <- -flip_sign
  }
  return(flip_sign)
}

# calculate the tumor classification for a cell line
# by identifying the most common tumor type within its 25
# nearest tumor neighbors
cell_line_tumor_class <- function(x, dist_mat, ann_mat, k=25, decreasing = T) {
  names(x) <- rownames(dist_mat)
  x <- sort(x, decreasing = decreasing)
  tumor_type <- dplyr::filter(ann_mat, sampleID %in% names(x[1:k]))$lineage %>%
    table() %>%
    as.data.frame() %>%
    dplyr::arrange(dplyr::desc(Freq)) %>%
    head(1) %>%
    .[['.']] %>% as.character()
  
  return(tumor_type)
}

# run GSEA using the fgsea package, with either gene or label permutation
run_fGSEA <- function (gsc, X = NULL, y = NULL, perm_type = "label", nperm = 10000, 
          min_set_size = 3, max_set_size = 300, stat_type = "t_stat", 
          stat_trans = "none", gseaParam = 1, nproc = 0, gene_stat = NULL) 
{
  if (any(grepl("package:piano", search()))) {
    detach("package:piano", unload = TRUE)
  }
  library(fgsea)
  stopifnot(perm_type %in% c("label", "gene"))
  if (perm_type == "label") {
    stopifnot(is.matrix(X) & !is.null(y))
    print(sprintf("Running sample-permutation testing with %d perms", 
                  nperm))
    used_samples <- which(!is.na(y))
    used_genes <- which(apply(X[used_samples, ], 2, var, 
                              na.rm = T) > 0)
    fgseaRes <- fgsea::fgsea(pathways = GSEABase:::geneIds(gsc), mat = t(X[used_samples, 
                                                                 used_genes]), labels = y[used_samples], minSize = min_set_size, 
                              maxSize = max_set_size, nperm = nperm, gseaParam = gseaParam, 
                              nproc = nproc)
  }
  else if (perm_type == "gene") {
    print(sprintf("Running gene-permutation testing with %d perms", 
                  nperm))
    fgseaRes <- fgsea::fgsea(pathways = GSEABase:::geneIds(gsc), stats = gene_stat, 
                             minSize = min_set_size, maxSize = max_set_size, 
                             nperm = nperm, gseaParam = gseaParam, nproc = nproc)
  }
  return(fgseaRes)
}

# Estimate linear-model stats for a matrix of data using limma with empirical Bayes moderated t-stats for p-values
run_lm_stats_limma <- function (mat, vec, covars = NULL, weights = NULL, target_type = "Gene", 
          limma_trend = FALSE) 
{
  require(limma)
  require(magrittr)
  require(tibble)
  require(plyr)
  require(dplyr)
  udata <- which(!is.na(vec))
  if (!is.numeric(vec)) {
    pred <- factor(vec[udata])
    stopifnot(length(levels(pred)) == 2)
    n_out <- colSums(!is.na(mat[udata[pred == levels(pred)[1]], 
                                , drop = F]))
    n_in <- colSums(!is.na(mat[udata[pred == levels(pred)[2]], 
                               , drop = F]))
    min_samples <- pmin(n_out, n_in) %>% set_names(colnames(mat))
  }
  else {
    pred <- vec[udata]
    min_samples <- colSums(!is.na(mat[udata, ]))
  }
  if (length(unique(pred)) <= 1) {
    return(NULL)
  }
  if (!is.null(covars)) {
    if (!is.data.frame(covars)) {
      covars <- data.frame(covars)
    }
    combined <- covars[udata, , drop = FALSE]
    combined[["pred"]] <- pred
    form <- as.formula(paste("~", paste0(colnames(combined), 
                                         collapse = " + ")))
    design <- model.matrix(form, combined)
    design <- design[, colSums(design) != 0, drop = FALSE]
  }
  else {
    design <- model.matrix(~pred)
  }
  if (!is.null(weights)) {
    if (is.matrix(weights)) {
      weights <- t(weights[udata, ])
    }
    else {
      weights <- weights[udata]
    }
  }
  fit <- limma::lmFit(t(mat[udata, ]), design, weights = weights)
  fit <- limma::eBayes(fit, trend = limma_trend)
  targ_coef <- grep("pred", colnames(design), value = TRUE)
  results <- limma::topTable(fit, coef = targ_coef, number = Inf)
  if (colnames(results)[1] == "ID") {
    colnames(results)[1] <- target_type
  }
  else {
    results %<>% tibble::rownames_to_column(var = target_type)
  }
  results$min_samples <- min_samples[results[[target_type]]]
  two_to_one_sided <- function(two_sided_p, stat, test_dir) {
    one_sided_p <- two_sided_p/2
    if (test_dir == "right") {
      one_sided_p[stat < 0] <- 1 - one_sided_p[stat < 
                                                 0]
    }
    else {
      one_sided_p[stat > 0] <- 1 - one_sided_p[stat > 
                                                 0]
    }
    return(one_sided_p)
  }
  results %<>% magrittr::set_colnames(revalue(colnames(.), c(logFC = "EffectSize", 
                                                   AveExpr = "Avg", t = "t_stat", B = "log_odds", P.Value = "p.value", 
                                                   adj.P.Val = "q.value", min_samples = "min_samples"))) %>% 
    na.omit()
  results %<>% dplyr::mutate(p.left = two_to_one_sided(p.value, 
                                                       EffectSize, "left"), p.right = two_to_one_sided(p.value, 
                                                                                                       EffectSize, "right"), q.left = p.adjust(p.left, method = "BH"), 
                             q.right = p.adjust(p.right, method = "BH"))
  return(results)
}

