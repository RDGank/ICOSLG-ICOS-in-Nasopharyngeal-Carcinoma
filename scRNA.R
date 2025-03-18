library(Seurat)
library(harmony)
load(file = 'sce_NPC.RData')
data<-Read10X(data.dir = 'GSE120926/')
sce<-CreateSeuratObject(data,min.cells = 3,min.features = 200,meta.data = meta)
meta<-read.table(file = 'GSE120926/GSE120926_NPC_10X_cell-barcode-process.txt',header = T,row.names = 1)
sce<-AddMetaData(sce,metadata = meta)
sce[['RNA']]<-as(object = sce[['RNA']],Class = 'Assay')
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
sce<-HAR_rm_batch(sce,'library')


FeatureStatPlot(sce,stat.by = c('ICOSLG','CD70','EPCAM','MS4A1'),stack = T,add_point = T,group.by = 'seurat_clusters'
                ,flip = T,pt.size = 0.001)

library(SCP)
GroupHeatmap(sce,features =c('ICOSLG','CD70','EPCAM','MS4A1','FOXP3','ICOS'),add_dot = T,group.by = 'seurat_clusters' )

save(sce,file = 'sce_NPC.RData')

DotPlot(sce,features = c('EPCAM','ICOSLG'),assay = 'RNA')
load(file = 'sce_GSE120.RData')
CellStatPlot(sce,group.by = 'celltype',stat.by = 'type')
DimPlot(sce,label = T)
DotPlot(sce,features =stromal_markers ,assay = 'RNA')
table(sce$orig.ident)
DimPlot(sce,label = T)
sce<-FindClusters(sce,resolution = 0.8)
library(stringr)
sce$celltype<-sce$seurat_clusters
sce$celltype<-str_replace(sce$celltype,'11|25|26|18|20','Tumor cells')
sce$celltype<-str_replace(sce$celltype,'24','Epithelial cells')
sce$celltype<-str_replace(sce$celltype,'15|21|14','Stromal cells')
sce$celltype<-str_replace(sce$celltype,'10','B')
sce$celltype<-str_replace(sce$celltype,'22','Plasma')
sce$celltype<-str_replace(sce$celltype,'12|16','T/NK')
sce$celltype<-str_replace(sce$celltype,'17','pDC')
sce$celltype<-str_replace(sce$celltype,'13','DC')
sce$celltype<-str_replace(sce$celltype,'23|26','Macrophage')
sce$celltype<-str_replace(sce$celltype,'19','Mast')
DimPlot(sce,group.by = 'celltype',label = T)                          
sce$celltype<-str_replace(sce$celltype,'3','Tumor cells')
sce$celltype<-str_replace(sce$celltype,'','Epithelial cells')
sce$celltype<-str_replace(sce$celltype,'','Stromal cells')
sce$celltype<-str_replace(sce$celltype,'0|5|9','B')
sce$celltype<-str_replace(sce$celltype,'7','Plasma')
sce$celltype<-str_replace(sce$celltype,'1|2|4|6','T/NK')
sce$celltype<-str_replace(sce$celltype,'','pDC')
sce$celltype<-str_replace(sce$celltype,'','DC')
sce$celltype<-str_replace(sce$celltype,'8','Macrophage')
DimPlot(sce,group.by = 'celltype',label = T)
DotPlot(sce,features = 'ICOSLG',group.by = 'celltype')
VlnPlot(sce,features = 'ICOSLG',group.by = 'celltype',assay = 'RNA')

save(sce,file = 'sce_GSE120.RData')
