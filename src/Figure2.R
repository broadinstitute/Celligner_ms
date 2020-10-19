library(here)
library(magrittr)
library(tidyverse)
source(here::here('src', 'analysis_helpers.R'))
source(here::here('src', 'global_params.R'))

# figure 2
Celligner_alignment_plot <- function(alignment) {
  Celligner_alignment <- ggplot2::ggplot(alignment, 
                                         ggplot2::aes(UMAP_1, UMAP_2, fill=lineage, size=type, color = type)) +
    ggplot2::geom_point(pch=21, alpha=0.7)  +
    ggplot2::scale_color_manual(values=c(`CL`='black', `tumor`='white')) +
    ggplot2::scale_size_manual(values=c(`CL`=1, `tumor`=0.75)) +
    ggplot2::theme_classic() + 
    ggplot2::theme(legend.position = 'bottom', 
          text=ggplot2::element_text(size=8),
          legend.margin =ggplot2::margin(0,0,0,0)) +
    ggplot2::guides(fill=FALSE, color=FALSE) +
    ggplot2::scale_fill_manual(values=tissue_colors) +
    ggplot2::xlab("UMAP 1") +
    ggplot2::ylab("UMAP 2")
  
  return(Celligner_alignment)
}
