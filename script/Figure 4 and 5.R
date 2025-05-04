library(Seurat);library(readxl);library(dplyr);library(ggplot2);library(pheatmap);library(harmony)
library(data.table);library(ggpubr);library(ArchR);library(DoubletFinder);library(RColorBrewer)

##### Figure 4D #####

for(i in unique(bc_epi$patient)){
  infercnv_count = GetAssayData(subset(bc_infercnv,subset = patient == i | cellType =='Fibroblast'), slot = 'counts')
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix = infercnv_count,
                                      annotations_file = data.frame(infercnv_anno[colnames(infercnv_count),],row.names = colnames(infercnv_count)),
                                      gene_order_file = infercnv_genes,
                                      min_max_counts_per_cell = c(100, +Inf),
                                      ref_group_names = c("Fibroblast"))
  infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff = 0.1,
                               out_dir = paste0('infercnv_all/',i), 
                               cluster_by_groups = F,
                               analysis_mode = 'samples',
                               HMM = F,
                               denoise = TRUE,
                               num_threads = 10,
                               no_plot = T,
                               #leiden_resolution_per_chr = 0.05,
                               #tumor_subcluster_pval = 0.01,
                               up_to_step=17
  )
  gc()
}

all_cnv_dir = list.dirs('infercnv_all',recursive = F)
all_sample_paired_cnv = lapply(all_cnv_dir,function(d){
  return(readRDS(file.path('infercnv_all',d,'run.final.infercnv_obj')))
})
names(all_sample_paired_cnv) = all_cnv_dir

group_list = lapply(all_sample_paired_cnv,function(cnv){
  cutree(cnv@tumor_subclusters$hc$all_observations,k = 2)
})
names(group_list) = all_cnv_dir

infercnv_obj_clu_list = lapply(names(all_sample_paired_cnv),function(n){
  aggregate(.~group,data = data.frame(t(all_sample_paired_cnv[[n]]@expr.data[,group_list.names[[n]]]),group = group_list[[n]]),median)
})

infercnv_obj_clu_mtx = do.call('rbind',infercnv_obj_clu_list)

Heatmap(data.matrix(infercnv_obj_clu_mtx),
        cluster_rows = F,cluster_columns = F,show_column_names = F,show_row_names = T,
        column_split = factor(gene_order$chr, paste("chr",1:22,sep = "")),
        column_gap = unit(0.01, "mm"),
        heatmap_legend_param = list(title = "Modified expression",
                                    direction = "vertical",
                                    title_position = "leftcenter-rot",
                                    at=c(0.7,1,1.3),legend_height = unit(3, "cm")),
        top_annotation = top_anno,#left_annotation = left_anno,
        row_title = NULL,column_title = NULL,raster_by_magick = T)

##### Figure 4E #####

pam50 <- read.table(file ="E:/Demo/PAM50/data/pam50_centroids.txt")
pam50 <- data.table('geneid'=rownames(pam50),pam50)
count <- as.matrix(bc_epi_orig@assays$RNA@counts)

count1 <- data.table('geneid'=rownames(count1),count1)
all_data <- inner_join(pam50,count1)  
pam50.1 <- all_data[,2:6]
all_data1 <- all_data[,-c(2:6)]

#Calculate euclidean distance
result <- data.frame()
for(i in 2:dim(all_data1)[2]){
    exp_cell <- all_data1[,i,with=F]
    base <- dist(t(cbind(pam50.1,exp_cell)),method = 'euclidean')
    data1 <- as.matrix(base)
    data2 <- data1[6,-6]
    lable <- names(data2)[which (data2==min(data2))]
    cellid <- names(exp_cell)
    result <- rbind(result, data.frame('cellid'=cellid,'label'=lable))
}

ggplot(data, aes(Sample,value, fill=condition))+
  scale_y_continuous(labels=scales::percent)+
  geom_bar(stat='identity',position='fill',
           width=0.75)+
  scale_fill_brewer(palette = "Set1")+
  theme(axis.text.x=element_text(angle=45,
                                 size=8.5,vjust = 1, hjust=1),
        strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid"),
        strip.text.x = element_text(size=8, face = "bold"),
        axis.text = element_text(size=17.5), 
        axis.title = element_text(size=25, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.key.size = unit(0.3, 'cm'), #change legend key size
        legend.key.height = unit(0.3, 'cm'), #change legend key height
        legend.key.width = unit(0.3, 'cm'), #change legend key width
        legend.title = element_text(size=15), #change legend title font size
        legend.text = element_text(size=10) #change legend text font size
  ) +
  scale_x_discrete(limits=as.character(c("CASE1BC" , "CASE1LM", 
                                        "CASE2BC",  "CASE2LM",  "CASE3BC",  "CASE3LM" , "CASE4BC" , "CASE4LM" , "CASE5BC" , "CASE5LM" , "CASE6BC" , "CASE6LM", 
                                        "CASE7BC" , "CASE7LM" , "CASE8BC" , "CASE8LM",  "CASE9BC" , "CASE9LM", "CASE10BC" ,"CASE10LM" ,"CASE11BC", "CASE11LM", 
                                        "CASE12BC" ,"CASE12LM" ,"CASE13BC" ,"CASE13LM")))+ 
  labs(x = 'Euclidean',y = 'Percentage') +
  theme(axis.title =element_text(size = 12),
        axis.text =element_text(size = 8, color = 'black'))
##### Figure 5C #####

library(scMetabolism)
bc_epi_metab = sc.metabolism.Seurat(bc_epi,method = 'AUCell',metabolism.type = 'KEGG')
bc_epi@meta.data = data.frame(bc_epi@meta.data,bc_epi_metab[,sel_path])

ggplot(bc_epi@meta.data, aes(x = cellType_2, y = Glycolysis, fill = tissue)) +
  introdataviz::geom_split_violin(alpha = .7, trim = FALSE,scale = 'width') +
  stat_summary(fun.data = "mean_se", geom = "pointrange", show.legend = F,
      position = position_dodge(.175)) +
  scale_y_continuous(name = "VEGFA") +
  scale_fill_manual(values = tissue_col2, name = "Tissue")+
  theme_pubr()

##### Figure 5D #####

ggplot(epi_point_df[epi_point_df$tissue == 'BC',],aes(x = PROM1,y = CD44)) +
  geom_point(mapping = aes(color = patient))+
  scale_fill_gradientn(colors = viridis::plasma(7))
ggplot(epi_point_df[epi_point_df$tissue == 'BC',],aes(x = PROM1+ runif(33135, min = 0.005, max = 0.01),
                                                      y = CD44+ runif(33135, min = 0.005, max = 0.01))) +
  geom_hdr(method = "kde",n=50,probs = c(0.99,0.95,0.9,0.85,0.8,0.75,0.5))+
  scale_fill_gradientn(colors = viridis::plasma(7))

ggplot(epi_point_df[epi_point_df$tissue == 'BCLM',],aes(x = PROM1,y = CD44)) +
  geom_point(mapping = aes(color = patient))+
  scale_fill_gradientn(colors = viridis::plasma(7))
ggplot(epi_point_df[epi_point_df$tissue == 'BCLM',],aes(x = PROM1+runif(33307, min = 0.005, max = 0.01), 
                                                        y = CD44+runif(33307, min = 0.005, max = 0.01))) +
  geom_hdr(method = "freqpoly",n=20,probs = c(0.999,0.99,0.95,0.9,0.85,0.8,0.75,0.5,0.25))+
  scale_fill_gradientn(colors = viridis::plasma(7))

##### The End #####
