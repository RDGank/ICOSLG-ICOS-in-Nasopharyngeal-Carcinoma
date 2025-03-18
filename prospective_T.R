library(Seurat)
library(SCP)
library(stringr)
library(ggplot2)
library(patchwork)
load(file = 'sce_prospective_T.RData')
DotPlot(sce,features = CD8_markers_list2,assay = 'RNA')
DotPlot(sce,features = c('CD3D','CD3E','NKG7','GNLY'),assay = 'RNA')
save(sce,file = 'sce_prospective_T.RData')
sce$celltype<-sce$seurat_clusters
sce$celltype<-str_replace(sce$celltype,'14','Treg')
sce$celltype<-str_replace(sce$celltype,'13','T_Naive')
sce$celltype<-str_replace(sce$celltype,'10|16','NK')
sce$celltype<-str_replace(sce$celltype,'11|12|15|18|17','T_Memory')
sce$celltype<-str_replace(sce$celltype,'','')
sce$celltype<-str_replace(sce$celltype,'','')
sce$celltype<-str_replace(sce$celltype,'','')
sce$celltype<-str_replace(sce$celltype,'','')
sce$celltype<-str_replace(sce$celltype,'','')
sce$celltype<-str_replace(sce$celltype,'','')
DimPlot(sce,group.by = 'celltype',label = T)
sce$celltype<-str_replace(sce$celltype,'2','Treg')
sce$celltype<-str_replace(sce$celltype,'5','Tfh')
sce$celltype<-str_replace(sce$celltype,'0|7|8|9','T_Naive')
sce$celltype<-str_replace(sce$celltype,'1|3|4','T_Memory')
sce$celltype<-str_replace(sce$celltype,'6','Teff')


markers_list =list(
  T_cell=c("CD3D","CD3E","CD4","CD8A","CD8B"),
  NK=c('NKG7','GNLY'),
  Treg=c("TNFRSF4","BATF","TNFRSF18","FOXP3","IL2RA","IKZF2"),
  T_Naive=c("CCR7","SELL","CD5"),
  T_Memory=c("GZMK","EOMES","ITM2C" ),
  Tfh=c("CXCR5","BCL6","ICA1","TOX","TOX2","IL6ST"),
  Teff=c("TBX21","FCGR3A","FGFBP2")#滤泡辅助性T细胞
  #ILC=c("TNFRSF25","KRT81","LST1","AREG","LTB","CD69")
) 

marker<-unlist(markers_list)
###CD8T
marker<-as.data.frame(marker)
marker$type<-rownames(marker)
marker$type<-str_replace(marker$type,'1|2|3|4|5|6|7|8|9','')


p1<-CellDimPlot(sce,group.by = 'celltype',label = T,split.by = 'Group',bg_color = 'white'
            ,show_stat = F,theme_use = 'theme_blank')

ht4 <- GroupHeatmap(sce,
                    features = marker$marker, feature_split = marker$type, group.by = "celltype",
                    heatmap_palette = "YlOrRd",
                    #cell_annotation = c("Phase", "G2M_score", "Neurod2"), cell_annotation_palette = c("Dark2", "Paired", "Paired"),
                    #cell_annotation_params = list(height = unit(10, "mm")),
                    #feature_annotation = c("TF", "CSPA"),
                    #feature_annotation_palcolor = list(c("gold", "steelblue"), c("forestgreen")),
                    add_dot = TRUE, add_bg = TRUE, nlabel = 0, show_row_names = TRUE
)

p2<-ht4$plot
py<-(p1/px)
p3/p4
p3<-CellStatPlot(sce,stat.by = 'Group',group.by = 'celltype',label = T)
p4<-CellStatPlot(sce,stat.by = 'celltype',group.by = 'Group',label.size = 3,label = T)
px<-p3+p4+plot_layout(design = 'AAAB')
pdf(file = 'prospective_T1.pdf',width = 7,height =3 )
p1
dev.off()
pdf(file = 'prospective_T2.pdf',width = 10,height =4.5 )
px
dev.off()
pdf(file = 'prospective_T3.pdf',width = 8,height =6.5 )
p2
dev.off()
p4<-FeatureStatPlot(sce,stat.by  = c('ICOS'),group.by = 'celltype',
                    split.by = 'Group',add_point = T,
                    pt.size = 0.0000001,stack = T,comparisons = T)
pdf(file = 'prospective_4T.pdf',width = 17,height =3 )
p4
dev.off()
