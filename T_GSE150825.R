load(file = 'sce_GSE150825_Tcell.RData')
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

sce<-FindClusters(sce,resolution = 0.8)
DimPlot(sce,label = T)
DotPlot(sce,features = c('MS4A1','CD4','CD8A','GNLY','NKG7','FOXP3'),assay = 'RNA')
DotPlot(sce,features = c(CD4_markers_list,CD8_markers_list2),assay = 'RNA')
DotPlot(sce,features = c(CD4_markers_list),assay = 'RNA')
DotPlot(sce,features = c(CD8_markers_list2),assay = 'RNA')

library(stringr)
sce$celltype<-sce$seurat_clusters
sce$celltype<-str_replace(sce$celltype,'13','T_Memory')
sce$celltype<-str_replace(sce$celltype,'10|12','Treg')
sce$celltype<-str_replace(sce$celltype,'15','T_Naive')
sce$celltype<-str_replace(sce$celltype,'11','Teff')
sce$celltype<-str_replace(sce$celltype,'14','NK')
DimPlot(sce,group.by = 'celltype')
sce$celltype<-str_replace(sce$celltype,'0|5|6','T_Naive')
sce$celltype<-str_replace(sce$celltype,'','')
sce$celltype<-str_replace(sce$celltype,'2','Treg')
sce$celltype<-str_replace(sce$celltype,'3','Tfh')
sce$celltype<-str_replace(sce$celltype,'1|4|8','T_Memory')
sce$celltype<-str_replace(sce$celltype,'7|9','T_Memory')
sce$celltype<-str_replace(sce$celltype,'','')
sce$celltype<-str_replace(sce$celltype,'','')
sce$celltype<-str_replace(sce$celltype,'','')
sce$celltype<-str_replace(sce$celltype,'','')
sce$celltype<-str_replace(sce$celltype,'','')
sce$celltype<-str_replace(sce$celltype,'','')

DimPlot(sce,group.by = 'celltype')
sce[['RNA']]<-as(object = sce[['RNA']],Class = 'Assay')
FeatureStatPlot(sce,stat.by = 'ICOS',bg.by = 'celltype',group.by = 'celltype',add_point = T,pt.size = 0.0001)
save(sce,file = 'sce_GSE150825_Tcell.RData')

