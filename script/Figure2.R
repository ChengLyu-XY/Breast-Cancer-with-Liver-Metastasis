library(Seurat);library(dplyr);library(ggplot2)
library(data.table);library(ggpubr);library(ArchR);library(RColorBrewer)
source('umap_kernel.R')
source('plot_utils.R')

colour_bk2 <- c(colorRampPalette(c("#3381b2","#9ACF9F"))(15),
               colorRampPalette(c("#9ACF9F","#F8F6BC"))(15),
               colorRampPalette(c("#F8F6BC","#F19A5A"))(15),
               colorRampPalette(c("#F19A5A","#C93E4C"))(16))

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

##### Figure2 B #####

plot_ked_umap(bc_T_nonk_fil)

##### Figure2 C #####

bc_tcr_data = bc_T_nonk_fil@meta.data[,c("Patient_cloneID","Tissue","Patient","Frac.Patient","Freq.Patient","have_tcr","sample")]
bc_tcr_data = bc_tcr_data[bc_tcr_data$have_tcr == 'yes',]
bc_tcr_data = data.table::data.table(bc_tcr_data)
clone_num_plot = bc_tcr_data[,.(count=.N),by=.(Tissue,Patient_cloneID,Patient,sample)]
allcount_table <- table(bc_tcr_data$sample) %>% data.table()
colnames(allcount_table) <- c('sample','all_count')

t1 <- clone_num_plot[Tissue=='BC',]
t1 <- t1[,.(all_count=sum(count),Patient_cloneID,Tissue,count,Patient),by=.(sample)]
t1[,frequency:=count/all_count]
colnames(t1) <- paste0(colnames(t1),'_bc')

t2 <- clone_num_plot[Tissue=='BCLM',]
t2 <- t2[,.(all_count=sum(count),Patient_cloneID,Tissue,count,Patient),by=.(sample)]
t2[,frequency:=count/all_count]
colnames(t2) <- paste0(colnames(t2),'_bclm')

plot_data_1 <- full_join(t1,t2,by=c('Patient_cloneID_bc'='Patient_cloneID_bclm'))
print(plot_data_1[!is.na(frequency_bc),frequency_bc] %>% min())
print(plot_data_1[!is.na(frequency_bclm),frequency_bclm] %>% min())
plot_data_1[is.na(frequency_bc),frequency_bc:=0.0001] 
plot_data_1[is.na(frequency_bclm),frequency_bclm:=0.0001]
plot_data_1[is.na(count_bc),count_bc:=0]
plot_data_1[is.na(count_bclm),count_bclm:=0]
plot_data_1[is.na(Patient_bclm),'Patient_bclm'] = plot_data_1[is.na(Patient_bclm),'Patient_bc']
plot_data_1[is.na(Patient_bc),'Patient_bc'] = plot_data_1[is.na(Patient_bc),'Patient_bclm']
plot_data_1$Tissue_bc = 'BC'
plot_data_1$Tissue_bclm = 'BCLM'

plot_data_1[count_bc==0,]
tmp <- left_join(data.table(sample=paste0(plot_data_1[count_bc==0,Patient_bclm],'_BC')),allcount_table,by=c('sample'='sample'))
plot_data_1[count_bc==0,'all_count_bc'] <- tmp$all_count

plot_data_1[count_bclm==0,]
tmp <- left_join(data.table(sample=paste0(plot_data_1[count_bclm==0,Patient_bc],'_BCLM')),allcount_table,by=c('sample'='sample'))
plot_data_1[count_bclm==0,'all_count_bclm'] <- tmp$all_count

plot_data_1$fisher_p <- 0
for (i in 1:dim(plot_data_1)[1]) {
  tmp <- plot_data_1[i,]
  tmp1 <- matrix(c(tmp$count_bc,tmp$all_count_bc-tmp$count_bc,tmp$count_bclm,tmp$all_count_bclm-tmp$count_bclm),nrow=2)
  tmp_res <- fisher.test(tmp1)
  plot_data_1[i,fisher_p:=tmp_res$p.value]
  
}

plot_data_1[,label:='Persistent']
plot_data_1[fisher_p<0.05&frequency_bc>frequency_bclm,label:='Contract']
plot_data_1[fisher_p<0.05&frequency_bc<frequency_bclm,label:='Expand']
plot_data_1[fisher_p<0.05&count_bc==0&count_bclm>0,label:='Novel']
table(plot_data_1$label)
plot_data_1[count_bc==0&count_bclm>0,]

ggplot(plot_data_1[!plot_data_1$Patient_bc %in% c('CASE5','CASE7'),],aes(x=frequency_bc,y=frequency_bclm))+
  geom_point(aes(shape=label,fill=label,color=label),stroke =0.2)+
  geom_abline(intercept = 0, slope = 1, color = "gray", linetype = "dashed") +
  scale_shape_manual(values = c('contract'=25,'Novel'=24,'expand'=24,"persistent"=21))+
  scale_fill_manual(values = c('Novel'='#00bc12','expand'='#AACE80','contract'='#ffd3d3','persistent'='#dadada'))+
  scale_color_manual(values = c('Novel'='black','expand'='black','contract'='black','persistent'='#b7b6b5'))+
  scale_x_continuous(trans = "log10", breaks = c(0.0001, 0.001, 0.01, 0.1,0.5), labels = c("1e-04", "0.001", "0.01", "0.1",'0.5'))+
  scale_y_continuous(trans = "log10", breaks = c(0.0001, 0.001, 0.01, 0.1,0.5), labels = c("1e-04", "0.001", "0.01", "0.1",'0.5'))+
  theme_classic()+
  ggtitle('')+ylab('Clone Frequency BC Liver Metastasis')+xlab('Clone frequency BC primary')+
  theme(axis.text.x=element_text(size = 8),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.text.y = element_text(size = 8),
        axis.line = element_line(size = 0.2),
        axis.ticks = element_line(size = 0.2),
        plot.title = element_text(hjust=0.5,size = 10),
        # legend.position = 'none',
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.key.size = unit(0.1, "inches"),
        legend.key.height = unit(0.1, "inches"),
        legend.key.width = unit(0.1, "inches")
  )

##### Figure 2H #####

tcr_multi = unique(bc_T_nonk_fil@meta.data[bc_T_nonk_fil$Freq.Patient >3,c("Patient_cloneID","T_CellType","Tissue")])
tcr_multi = na.omit(tcr_multi)
tcr_multi = tcr_multi[tcr_multi$T_CellType!='CD8_gdT',]

tcr_share_prop = matrix(0,ncol = 14,nrow = 14)
colnames(tcr_share_prop) = rownames(tcr_share_prop) = unique(tcr_multi$T_CellType)
for(i in rownames(tcr_share_prop)){
  for (j in colnames(tcr_share_prop)) {
    tmp2 = tcr_multi[(tcr_multi$T_CellType == i & tcr_multi$Tissue == 'BC') | (tcr_multi$T_CellType == j & tcr_multi$Tissue == 'BCLM'),]
    tmp3 = sum(tcr_multi$T_CellType == i & tcr_multi$Tissue == 'BC')
    tcr_share_prop[i,j] = length(intersect(tmp2$Patient_cloneID[tmp2$Tissue == 'BC'],tmp2$Patient_cloneID[tmp2$Tissue == 'BCLM']))/tmp3
  }
}

pheatmap(tcr_share_prop,cluster_rows = F,cluster_cols = F,color = colour_bk2,border_color = NA)

