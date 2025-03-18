
DimPlot(sce,label = T)
DotPlot(sce,features = myeloids_markers_list2,assay = 'RNA')
DotPlot(sce,features = gastric_cancer_markers,assay = 'RNA')
DotPlot(sce,features = Bcels_markers_list,assay = 'RNA')
DotPlot(sce,features = c('CD3D','CD3E','NKG7','GNLY'),assay = 'RNA')
DotPlot(sce,features = stromal_markers,assay = 'RNA')

sceepican<-sce[,sce$seurat_clusters%in%c(22,29,33,38,43,44,42)]

save(sceepican,file = 'sceepi_prospect.RData')

DotPlot(sce,features = 'ICOSLG',assay = 'RNA',group.by = 'celltype',split.by = 'Group')
FeatureDimPlot(sce,features = 'ICOSLG')
library(stringr)
sce$celltype<-sce$seurat_clusters
sce$celltype<-str_replace(sce$celltype,'39','Mast')
sce$celltype<-str_replace(sce$celltype,'32','DC')
sce$celltype<-str_replace(sce$celltype,'15|35|36','Macrophage')
sce$celltype<-str_replace(sce$celltype,'20|41|47','pDC')
sce$celltype<-str_replace(sce$celltype,'22|33|38|43|44|42','Cancer cells')
sce$celltype<-str_replace(sce$celltype,'29','Epithelial cells')
sce$celltype<-str_replace(sce$celltype,'14|30|31','Plasma')
sce$celltype<-str_replace(sce$celltype,'17|21|25|27|28|34|45','B')
sce$celltype<-str_replace(sce$celltype,'10|11|12|13|16|18|19|21|23|24|25|26|37|40|46|47','T/NK')
sce$celltype<-str_replace(sce$celltype,'44','Stromal')

DimPlot(sce,group.by = 'celltype',label = T)
sce$celltype<-str_replace(sce$celltype,'2|4|6','B')
sce$celltype<-str_replace(sce$celltype,'0|1|3|5|7|8|9','T/NK')
sce$celltype<-str_replace(sce$celltype,'','')

sce$Group<-str_replace(sce$Group,'Pr','Pre-treatment')
sce$Group<-str_replace(sce$Group,'Po','Post-treatment')

DimPlot(sce,group.by = 'Group')


save(sce,file = 'sce_prospective.RData')
sce<-sce[,sce$celltype%in%c('T/NK')]
library(harmony)
HAR_rm_batch<-function(project,split.by){
  print(dim(project))
  project <- NormalizeData(project, normalization.method = "LogNormalize",scale.factor = 1e4)%>% 
    FindVariableFeatures()%>%
    ScaleData()%>%
    RunPCA(features = VariableFeatures(object = project))%>%
    RunHarmony(split.by)%>%
    RunUMAP(dims = 1:15,reduction = "harmony")%>%
    RunTSNE(dims = 1:15,reduction = "harmony",chech_duplicates=F)%>%
    FindNeighbors(reduction = "harmony",dims = 1:15) 
  
  for (res in c(0.5,0.8,1)) {
    project=FindClusters(project,resolution = res, algorithm = 1)
  }
  return(project)
}
#eg

sce<-HAR_rm_batch(sce,'orig.ident')
save(sce,file = 'sce_prospective_T.RData')
DimPlot(sce)

load(file = 'sce_prospective.RData')
p1<-CellDimPlot(sce,group.by = 'celltype',label = T,split.by = 'Group',bg_color = 'white'
                ,show_stat = F,theme_use = 'theme_blank',label.size = 2)
p2<-CellStatPlot(sce,stat.by = 'Group',group.by = 'celltype',label = T)

markers_list =list(
  T_NK=c("CD3D","CD3E","NKG7",'GNLY'),
  #NK=c('NKG7','GNLY'),
  B = c('MS4A1','SDC1','CD27','CD38','CD19', 'CD79A'),
  Plasma= c ( 'IGHA1','IGHG1'), 
  Mast=c('TPSB2','TPSAB1'),
  Epi_Cancer=c('EPCAM', 'KRT18', 'MUC1'),
  pDC = c("CLEC4C","IRF7","TCF4","GZMB"),
  Macrophages = c("APOC1","HLA-DRB5","C1QA","C1QB"),
  DC = c("CCL19","LAMP3","IDO1","IDO2","LAD1","FSCN1","CCR7","LY75","CCL22")
  
  
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

p3<-ht4$plot

p4<-FeatureStatPlot(sce,stat.by  = c('ICOSLG','ICOS'),group.by = 'celltype',
                split.by = 'Group',add_point = T,
                pt.size = 0.0000001,stack = T,comparisons = T)


ht4 <- GroupHeatmap(sce,
                    features =c('ICOSLG','ICOS'), group.by = "celltype",
                    heatmap_palette = "YlOrRd",split.by = 'Group',
                    #cell_annotation = c("Phase", "G2M_score", "Neurod2"), cell_annotation_palette = c("Dark2", "Paired", "Paired"),
                    #cell_annotation_params = list(height = unit(10, "mm")),
                    #feature_annotation = c("TF", "CSPA"),
                    #feature_annotation_palcolor = list(c("gold", "steelblue"), c("forestgreen")),
                    add_dot = TRUE, add_bg = TRUE, nlabel = 0, show_row_names = TRUE
)

p5<-ht4$plot


pdf(file = 'prospective_1.pdf',width = 7.5,height =3.5 )
p1
dev.off()
pdf(file = 'prospective_2.pdf',width = 7,height =3 )
p2
dev.off()
pdf(file = 'prospective_3.pdf',width = 8.5,height =8 )
p3
dev.off()
pdf(file = 'prospective_4.pdf',width = 17,height =5 )
p4
dev.off()
pdf(file = 'prospective_5.pdf',width = 10,height =3 )
p5
dev.off()

