library(Seurat)
library(infercnv)
library(stringr)
#查看数据
seurat_object <- readRDS('seurat_object.rds')
levels(seurat_object)
#筛选分析需要的细胞类型
seurat_object <- subset(seurat_object, idents=c('Epithelial cells', 'Myeloid cells', 'T cells'))
#抽样，仅用于该教程
seurat_object <- subset(seurat_object, downsample=200)

setwd("E:/NPCsc/pipeline/")
load(file = 'sce_GSE150430.RData')
load(file = 'sceepi_prospect.RData')
table(sce$celltype) 
seurat_object<-sce
sce<-merge(sce,seurat_object)
sce<-JoinLayers(sce)
sce<-sce[,sce$celltype%in%c('Epithelial cells')]
sce<-sce[,sce$Group%in%c('HC')]
DimPlot(sce)
counts <- GetAssayData(sce, slot = 'counts')
sce$Group<-sce$type
library(stringr)
sce$Group<-str_replace(sce$Group,'npc','NPC')
sce$Group<-str_replace(sce$Group,'non','HC')

Idents(sce)<-'Group'
anno <- data.frame(Idents(sce))
gene_order <- 'hg38_gencode_v27.txt'

x<-read.table('hg38_gencode_v27.txt')
k<-!duplicated(x$V1)
table(k)
x<-x[k,]
#write.table(x,file = 'hg38_gencode_v27.txt')
infercnv_obj = CreateInfercnvObject(raw_counts_matrix = counts,
                                     annotations_file = anno,
                                     delim="\t",
                                     gene_order_file = gene_order,
                                     min_max_counts_per_cell = c(100, +Inf),
                                     ref_group_names = c("HC"))

setwd(dir = 'inferCNV_GSE160')
infercnv_obj = infercnv::run(infercnv_obj,
                              cutoff = 0.1,
                              out_dir = ".", 
                              cluster_by_groups = F,
                              k_obs_groups = 8,
                              HMM = FALSE,
                              denoise = TRUE,
                              num_threads = 8
                                ) 


