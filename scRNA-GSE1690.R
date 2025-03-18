

filenames<-list.files('GSE160/')
filenames1<-filenames[1]
filenames10<-filenames[10]
setwd('GSE160/')
for (i in filenames) {
  data<-read.csv(file = i,header = T,row.names = 1)
  sce1<-CreateSeuratObject(data,min.cells = 3,min.features = 200)
  if(i==filenames1) sce<-sce1
  else sce<-merge(sce,sce1)
  if(i==filenames10) return(sce)
}

sce<-JoinLayers(sce,assay = 'RNA')
sce<-HAR_rm_batch(sce,'orig.ident')
DotPlot(sce,features = c('EPCAM','TUBB2B','ICOSLG','MS4A1','CD70'),assay = 'RNA')

setwd('E:/NPCsc/pipeline')
save(sce,file = 'sce_GSE160.RData')
load(file = 'sce_GSE160.RData')

VlnPlot(sce,features = c('EPCAM','ICOSLG'),assay = 'RNA')

CellStatPlot(sce,group.by = 'celltype',stat.by = 'Group')
DimPlot(sce,label = T)
DotPlot(sce,features =stromal_markers ,assay = 'RNA')
table(sce$orig.ident)
library(stringr)
sce$celltype<-sce$seurat_clusters
sce$celltype<-str_replace(sce$celltype,'13|15','Tumor cells')
sce$celltype<-str_replace(sce$celltype,'','Epithelial cells')
sce$celltype<-str_replace(sce$celltype,'20','Stromal cells')
sce$celltype<-str_replace(sce$celltype,'18','Mast')
sce$celltype<-str_replace(sce$celltype,'16','Plasma')
sce$celltype<-str_replace(sce$celltype,'10|11|12','T/NK')
sce$celltype<-str_replace(sce$celltype,'19','pDC')
sce$celltype<-str_replace(sce$celltype,'17','DC')
sce$celltype<-str_replace(sce$celltype,'14','Macrophage')
DimPlot(sce,label = T,group.by = 'celltype')
sce$celltype<-str_replace(sce$celltype,'','Tumor cells')
sce$celltype<-str_replace(sce$celltype,'','Epithelial cells')
sce$celltype<-str_replace(sce$celltype,'','Stromal cells')
sce$celltype<-str_replace(sce$celltype,'2|6','B')
sce$celltype<-str_replace(sce$celltype,'','Plasma')
sce$celltype<-str_replace(sce$celltype,'0|1|3|4|5|7|8|9','T/NK')
sce$celltype<-str_replace(sce$celltype,'','pDC')
sce$celltype<-str_replace(sce$celltype,'','DC')
sce$celltype<-str_replace(sce$celltype,'','Macrophage')

DimPlot(sce,group.by = 'celltype',label = T)
DotPlot(sce,features = 'ICOSLG',group.by = 'celltype')
VlnPlot(sce,features = 'ICOSLG',group.by = 'celltype',assay = 'RNA')
save(sce,file = 'sce_GSE160.RData')
