library(Seurat);library(dplyr);library(ggplot2);library(pheatmap)
library(data.table);library(ggpubr);library(ArchR);library(RColorBrewer)

colour_bk <- c(colorRampPalette(c("#4575b4","#e0f3f8"))(35),
               colorRampPalette(c("#e0f3f8","#ffffbf"))(10),
               colorRampPalette(c("#ffffbf","#d73027"))(36))
bk <- bk <- c(seq(-2,2,by=0.05))

only_draw_pair_boxplot <- function(plot_df,compare_list,test.method = 'wilcox.test',
                                   x,y,fill = NA,color = NA,group = NA,y_lab,pt.size = 1,alpha = 0.9,fill_col){  
  ggplot(data=plot_df, aes(x = plot_df[,x], y = plot_df[,y])) +
    geom_boxplot(alpha =0.9,size=1,outlier.shape = NA, mapping = aes(fill = plot_df[,fill]),)+
    #scale_fill_manual(limits=c("BC","BCLM"), values=c("#D72729","#2B69B2"))+
    stat_compare_means(method = test.method, paired = T, comparisons=compare_list)+
    geom_jitter(size=pt.size, shape=16,aes(group=plot_df[,group]),alpha = alpha,position = position_dodge(0))+
    scale_color_manual(values = 'black')+
    scale_fill_manual(values = fill_col)+
    geom_line(aes(group = plot_df[,group]), color = 'grey40', lwd = 0.3,position = position_dodge(0),alpha = alpha)+
    theme_classic() +
    theme(panel.grid =element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 17,colour = "black"),
          axis.title = element_text(size = 20,colour = "black"),
          legend.text = element_text(size = 15,colour = "black"),
          legend.title = element_text(size = 17,colour = "black"),
          legend.position = 'top')+
    labs(y = y_lab)
}

calc_celltype_roe <- function(se_obj, cluster, tissue){
  ct_roe = list()
  for(i in unique(se_obj@meta.data[cluster,])){
    ct_roe_tmp = calTissueDist(se_obj@meta.data[se_obj@meta.data[cluster,] == i,],colname.cluster = "cellType_1",colname.tissue = "tissue")
    ct_roe_tmp = reshape2::melt(ct_roe_tmp)
    colnames(ct_roe_tmp) = c('cell_types','tissue','Roe')
    ct_roe_tmp$Patient = i
    ct_roe = append(ct_roe,list(ct_roe_tmp))
  }
  return(do.call(rbind,ct_roe))
}
##### Figure1 D #####
cellType_1_roe = calc_celltype_roe(all_cells_fil,'cellType_1','tissue')
top_anno = unique(all_cells_fil@meta.data[,c('subtype','tissue','orig.ident')])
rownames(top_anno) = top_anno$orig.ident; top_anno$orig.ident = NULL
top_anno = top_anno[order(top_anno$tissue,top_anno$subtype),]
pheatmap(ct_roe[,rownames(top_anno)],scale = 'row',cluster_cols = F,color = colour_bk,breaks = bk,show_colnames = F)

##### Figure1 E-H #####
bc_all_plot_df = table()
