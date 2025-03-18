library(Seurat)
library(SCP)
data<-read.table('E:/NPCsc/pipeline/GSE150430_npc_scRNA_hg19_processed_data.txt/npc_scRNA_hg19_processed_data.txt',header=T,row.names = 1)
sce<-CreateSeuratObject(data,min.cells = 3,min.features = 200)

save(sce,file = 'sce_GSE150430.RData')

load(file = 'sce_GSE150430.RData')
table(sce$orig.ident)

library(stringr)
sce$Group<-sce$orig.ident
sce$Group<-str_replace(sce$Group,'N','HC')
sce$Group<-str_replace(sce$Group,'P','NPC')
sce$Group<-str_replace(sce$Group,'0','')
sce$Group<-str_replace(sce$Group,'1','')
sce$Group<-str_replace(sce$Group,'2','')
sce$Group<-str_replace(sce$Group,'3','')
sce$Group<-str_replace(sce$Group,'4','')
sce$Group<-str_replace(sce$Group,'5','')
sce$Group<-str_replace(sce$Group,'6','')
sce$Group<-str_replace(sce$Group,'7','')
sce$Group<-str_replace(sce$Group,'8|9','')
table(sce$Group)

sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")

VlnPlot(sce,features = c('percent.mt','nFeature_RNA'))

sce<-HAR_rm_batch(sce,split.by = 'orig.ident')
DimPlot(sce)
sce<-FindClusters(sce,resolution =0.8)
FeaturePlot(sce,features = c('EPCAM','ICOSLG'))
CellStatPlot(sce,stat.by = 'Group',group.by = 'celltype')
DotPlot(sce,features = last_markers,assay = 'RNA')
DotPlot(sce,features =Tcells_markers ,assay = 'RNA')


library(stringr)
sce$celltype<-sce$seurat_clusters
sce$celltype<-str_replace(sce$celltype,'14|16|18|20|21|22','Tumor cells')
sce$celltype<-str_replace(sce$celltype,'12','Epithelial cells')
sce$celltype<-str_replace(sce$celltype,'','Stromal cells')
sce$celltype<-str_replace(sce$celltype,'11|17|19','B')
sce$celltype<-str_replace(sce$celltype,'','Plasma')
sce$celltype<-str_replace(sce$celltype,'','T/NK')
sce$celltype<-str_replace(sce$celltype,'15','pDC')
sce$celltype<-str_replace(sce$celltype,'10|13','DC')
sce$celltype<-str_replace(sce$celltype,'18','Macrophage')
sce$celltype<-str_replace(sce$celltype,'','')
sce$celltype<-str_replace(sce$celltype,'5|9','Tumor cells')
sce$celltype<-str_replace(sce$celltype,'','Epithelial cells')
sce$celltype<-str_replace(sce$celltype,'','Stromal cells')
sce$celltype<-str_replace(sce$celltype,'3|0','B')
sce$celltype<-str_replace(sce$celltype,'7','Plasma')
sce$celltype<-str_replace(sce$celltype,'1|2|4|8','T/NK')
sce$celltype<-str_replace(sce$celltype,'','pDC')
sce$celltype<-str_replace(sce$celltype,'','DC')
sce$celltype<-str_replace(sce$celltype,'6','Macrophage')

DimPlot(sce,group.by = 'celltype',label = T)
DotPlot(sce,features = 'ICOSLG',assay = 'RNA',group.by = 'celltype')

save(sce,file = 'sce_GSE150430.RData')

