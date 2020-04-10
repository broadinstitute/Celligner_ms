# Parameters
global <- list(
  n_genes = 'all', # set to 'all' to use all protein coding genes found in both datasets 
  umap_n_neighbors = 10, # num nearest neighbors used to create UMAP plot
  umap_min_dist = 0.5, # min distance used to create UMAP plot
  mnn_k_CL = 5, # number of nearest neighbors of tumors in the cell line data
  mnn_k_tumor = 50, # number of nearest neighbors of cell lines in the tumor data
  top_DE_genes_per = 1000, # differentially expressed genes with a rank better than this is in the cell line or tumor data
  # are used to identify mutual nearest neighbors in the MNN alignment step
  remove_cPCA_dims = c(1,2,3,4), # which cPCA dimensions to regress out of the data 
  distance_metric = 'euclidean', # distance metric used for the UMAP projection
  mod_clust_res = 5, # resolution parameter used for clustering the data
  mnn_ndist = 3, # ndist parameter used for MNN
  n_PC_dims = 70, # number of PCs to use for dimensionality reduction
  reduction.use = 'umap', # 2D projection used for plotting
  fast_cPCA = 10 # to run fast cPCA (approximate the cPCA eigenvectors instead of calculating all) set this to a value >= 4
)

table_theme <- gridExtra::ttheme_minimal(
  core=list(bg_params = list(fill = rep("white",3), col='black'),
            fg_params=list(fontsize=6)),
  colhead = list(fg_params=list(col="black", fontsize=7),
                 bg_params=list(fill="white", col='black')))

tissue_colors <- c(`central_nervous_system`= "#f5899e",`engineered_central_nervous_system` = "#f5899e",
                   `teratoma` = "#f5899e",
                   `bone` = "#9f55bb",   
                   `pancreas` = "#b644dc", 
                   `soft_tissue` = "#5fdb69",
                   `skin` = "#6c55e2",    
                   `liver` = "#9c5e2b",
                   `blood` = "#da45bb",
                   `lymphocyte`=  "#abd23f",
                   `peripheral_nervous_system` = "#73e03d",
                   `ovary` = "#56e79d",`engineered_ovary` = "#56e79d",
                   `adrenal` = "#e13978",  `adrenal_cortex` = "#e13978",
                   `upper_aerodigestive` = "#5da134",
                   `kidney` = "#1f8fff",`engineered_kidney` = "#1f8fff",
                   `gastric` = "#dfbc3a",
                   `eye` = "#349077",
                   `nasopharynx` = "#a9e082",
                   `nerve` = "#c44c90",
                   `unknown` = "#999999",
                   `cervix` = "#5ab172",
                   `thyroid` = "#d74829",
                   `lung` = "#51d5e0",`engineered_lung` = "#51d5e0",
                   `rhabdoid` = "#d04850",
                   `germ_cell` = "#75dfbb",   `embryo` = "#75dfbb",
                   `colorectal` = "#96568e",
                   `endocrine` = "#d1d684",
                   `bile_duct` = "#c091e3",                        
                   `pineal` = "#949031",
                   `thymus` = "#659fd9",
                   `mesothelioma` = "#dc882d",
                   `prostate` = "#3870c9", `engineered_prostate` = "#3870c9",
                   `uterus` = "#e491c1",
                   `breast` = "#45a132",`engineered_breast` = "#45a132",
                   `urinary_tract` = "#e08571",
                   `esophagus` = "#6a6c2c",
                   `fibroblast` = "#d8ab6a",
                   `plasma_cell` = "#e6c241")   

heatmap_params <- list(
  na_color = '#666666',
  title_font_size = 6,
  row_font_size = 6,
  column_font_size = 6,
  font_face = "plain",
  square_border_color = "white",
  color_palette = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11, 'YlOrRd')), space='Lab'),
  color_vector = c('#e3e3e3', rev(grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11, 'YlOrRd')), space='Lab')(100)))
)

EMT_genes <- c("ZEB1", "LIX1L", "VIM", "AXL", 'MMP2', "ANTXR2", "C3ORF21", "FN1", "NRP1", "TGFB1",
               "GALNT5", "PPARG", "HNMT", "CARD6", "RBPMS", "TNFRSF21", "TMEM45B", "MPP7", "SSH3",
               "MTAC2D1", "MUC1", "EPPK1", "SHROOM3", "EPN3", "PRSS22", "AP1M2", "SH3YL1", "KLC3",
               "SERINC2", "EVPL", "FXYD3", "CLDN4", "CRB3", "LRRC54", "MAPK13", "EPPK1", "FALNT3",
               "STAP2", "AP1M2", "DSP", "ELMO3", "KRTCAP3", "MAL2", "F11R", "GPR110", "GPR56",
               "KRT19", "GRHL1", "BSPRY", "C1ORF116", "S100A14", "SPINT2", "ANKRD22", "ST14", "GRHL2",
               "PRR5", "TJP3", "TACSTD2", "CH3", "C1ORD172", "CDS1", "MPZL2","INADL", "EPN3", "RBM35A", 
               "TMC4", "ITGB6", "TMEM125", "EPHA1", "CDS1", "ENPP5", "RAB25", "PRSS8", "TMEM30B", 
               "CLDN7", "RBM35A", "TACSTD1", "CDS1", "SCNN1A", "CDH1")

ALL_marker_genes <- cbind.data.frame(`gene` = c("CD3G", "CD3D", "CD3E", "CD8A", "CD19", "CD22", "CD2",  "CD10"),
                                    `ensembl` = c("ENSG00000160654", "ENSG00000167286", "ENSG00000198851",
                                                  "ENSG00000153563", "ENSG00000177455", "ENSG00000012124",
                                                  "ENSG00000116824", "ENSG00000196549"),
                                    `type` = c("T-cell", "T-cell", "T-cell", "T-cell", "B-cell", "B-cell", "T-cell",
                                               "B-cell"))
        
                                     
                                     
                                     
                                     
