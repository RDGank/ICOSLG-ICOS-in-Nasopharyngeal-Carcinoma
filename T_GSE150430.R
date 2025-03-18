sce<-sce[,sce$celltype%in%c('T/NK')]
library(Seurat)
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
DimPlot(sce,label =T )

k<-!duplicated(cd4_and_cd8T_markers_list)

DotPlot(sce,features = c(CD4_markers_list,CD8_markers_list2),assay = 'RNA')
DotPlot(sce,features = c(CD4_markers_list),assay = 'RNA')
DotPlot(sce,features = c(CD8_markers_list2),assay = 'RNA')

library(stringr)
sce$celltype<-sce$seurat_clusters
sce$celltype<-str_replace(sce$celltype,'13','T_Memory')
sce$celltype<-str_replace(sce$celltype,'12','Treg')
sce$celltype<-str_replace(sce$celltype,'11','T_Naive')
sce$celltype<-str_replace(sce$celltype,'10','T_Memory')
sce$celltype<-str_replace(sce$celltype,'1|9','T_Naive')
sce$celltype<-str_replace(sce$celltype,'8','NKT')
sce$celltype<-str_replace(sce$celltype,'2|4|7','Treg')
sce$celltype<-str_replace(sce$celltype,'5','Tfh')
sce$celltype<-str_replace(sce$celltype,'0|3|6','T_Memory')
sce$celltype<-str_replace(sce$celltype,'','')
sce$celltype<-str_replace(sce$celltype,'','')
sce$celltype<-str_replace(sce$celltype,'','')
sce$celltype<-str_replace(sce$celltype,'','')
sce$celltype<-str_replace(sce$celltype,'','')
sce$celltype<-str_replace(sce$celltype,'','')
sce$celltype<-str_replace(sce$celltype,'','')

DimPlot(sce,group.by = 'celltype')

FeatureStatPlot(sce,stat.by = 'ICOS',bg.by = 'celltype',group.by = 'celltype',add_point = T,pt.size = 0.0001)
save(sce,file = 'sce_GSE150430_Tcell.RData')

