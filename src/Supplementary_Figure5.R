
# revised supplementary plot 5a
fibroblast_cluster_plot <- function(Celligner_info, cur_lineage = 'fibroblast') {
  
  fibrolast_cl_plot <- ggplot2::ggplot(Celligner_info, 
                                   ggplot2::aes(UMAP_1, UMAP_2, fill=lineage==cur_lineage, size=type, color=type)) +
    ggplot2::geom_point(alpha=0.7, pch=21) +
    ggplot2::scale_size_manual(values=c(`CL`=1.5, `tumor`=1)) +
    ggplot2::scale_color_manual(values=c(`CL`='black', `tumor`='white')) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = 'right', text=ggplot2::element_text(size=6), axis.text = ggplot2::element_text(size=6)) +
    ggplot2::guides(size=FALSE, fill=ggplot2::guide_legend(title=cur_lineage), color=FALSE) +
    ggplot2::xlab("UMAP 1") +
    ggplot2::ylab("UMAP 2")
  
  return(fibroblast_cl_plot)
  
}

# revised supplementary plot 5b
fibroblast_highlight_plot <- function(Celligner_info) {
  fibrolast_cls_highlight <- ggplot2::ggplot(filter(Celligner_info, UMAP_1 >6.5 & UMAP_1 < 9 & UMAP_2 > -3 & UMAP_2 < -1.5),
                                             ggplot2::aes(UMAP_1, UMAP_2, fill=lineage=='fibroblast', size=type, color=type)) +
    ggplot2::geom_point(alpha=0.7, pch=21) +
    ggplot2::scale_size_manual(values=c(`CL`=1.5, `tumor`=1)) +
    ggplot2::scale_color_manual(values=c(`CL`='black', `tumor`='white')) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = 'none', text=ggplot2::element_text(size=6), axis.text = ggplot2::element_text(size=6)) +
    ggplot2::guides(size=FALSE, fill=ggplot2::guide_legend(title="fibroblast"), color=FALSE) +
    ggplot2::xlab("UMAP 1") +
    ggplot2::ylab("UMAP 2")
  
  return(fibroblast_cls_highlight)

}
