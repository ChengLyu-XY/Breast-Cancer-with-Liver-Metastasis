
plot_ked_umap <- function(se_obj,reduction = 'harmony_umap',n = 200, h = 0.2, quench = 0.015){
  
  library(tibble)
  library(cowplot)
  library(dplyr)
  library(tidyr)
  library(stringr)
  
  rdbu <- rev(RColorBrewer::brewer.pal(7, "RdBu"))
  myfeatures <- c("harmonyumap1", "harmonyumap2", "sample", "nCount_RNA", "nFeature_RNA", "percent.mt", "cell_types","tissue",'seurat_clusters')
  
  
  se_obj$harmonyumap1 = se_obj[[reduction]]@cell.embeddings[,1]
  se_obj$harmonyumap2 = se_obj[[reduction]]@cell.embeddings[,2]
  
  plot_data <- as_tibble(FetchData(se_obj, c(myfeatures)))%>% 
    mutate(umapscaled_1 = scales::rescale(harmonyumap1),
           umapscaled_2 = scales::rescale(harmonyumap2))
  
  kernel_tbl <- kde2d_contrast(
    data_tbl = plot_data,
    x = harmonyumap1, y = harmonyumap2, kernel_group = tissue,
    kernel_subsets_a = c("BC"), kernel_subsets_b = c("BCLM"),
    n = n, h = h, quench = quench
  ) %>%
    filter(value.x > 0.00001 & value.y > 0.00001)
  
  plot_data$cell_types = as.factor(plot_data$cell_types)
  
  plot_data_label <- plot_data %>% 
    group_by(cell_types) %>% 
    summarise(harmonyumap1 = median(harmonyumap1),
              harmonyumap2 = median(harmonyumap2),
              umapscaled_1 = median(umapscaled_1),
              umapscaled_2 = median(umapscaled_2)) %>% 
    mutate(cluster_number = as.numeric(cell_types))
  
  kernel_umap <- ggplot() +
    ggrastr::geom_point_rast(aes(umapscaled_1, umapscaled_2), color = "grey90", size = 0.01, alpha = 0.02,
                             data = plot_data,
                             raster.dpi = 50, raster.width = 10, raster.height = 10) +
    ggrastr::geom_point_rast(aes(x, y, color = value_scaled), data = kernel_tbl, size = 0.01, 
                             raster.dpi = 50, raster.width = 10, raster.height = 10) +
    geom_point(aes(umapscaled_1, umapscaled_2), color = "white",
               data = plot_data_label, alpha = 0.5, size = 6) +
    geom_text(aes(umapscaled_1, umapscaled_2, label = cluster_number), color = "black",
              data = plot_data_label) +
    scale_color_gradientn(colours = rdbu,
                          values = scales::rescale(c(min(kernel_tbl$value_scaled), 0,
                                                     max(kernel_tbl$value_scaled)))) +
    guides(color = guide_colorbar(label.position = "right",
                                  title.position = "top",
                                  title.hjust = 0,
                                  title.vjust = 1,
                                  direction = "vertical")) +
    theme(aspect.ratio = 1,
          legend.key.height = unit(0.05, "npc"),
          legend.key.width = unit(0.03, "npc"),
          legend.position = c(-0.1, 0.95),
          legend.justification = c("left", "top"),
          plot.title = element_text(hjust = 0.5, vjust = 0.5, face = "plain", size = 18)) + 
    labs(color = "Enrichment\nin BC") +
    remove_xaxis +
    remove_yaxis
  
  add_umap_coord(kernel_umap)
}
