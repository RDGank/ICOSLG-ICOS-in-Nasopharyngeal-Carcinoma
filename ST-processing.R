
library(Seurat)
spe2 = Read10X("E:/NPCsc/pipeline/GSE206245_RAW/ST5/out/filtered_feature_bc_matrix/")
image2 <- Read10X_Image(image.dir = file.path("E:/NPCsc/pipeline/GSE206245_RAW/ST5/out", 
                                              "spatial"), filter.matrix = TRUE)
spe2 <- CreateSeuratObject(counts = spe2, assay = "Spatial")

image2 <- image2[Cells(x = spe2)]
DefaultAssay(spe2 = image2) <- "Spatial"
spe2[["slice1"]] <- image2
SpatialFeaturePlot(spe2, features = "nFeature_Spatial")

GBM4<-spe2
GBM4 <- SCTransform(GBM4, assay = "Spatial", verbose = FALSE)
GBM4 <- RunPCA(GBM4, assay = "SCT", verbose = FALSE) 

##数据聚类
GBM4<- FindNeighbors(GBM4, reduction = "pca", dims = 1:10)
GBM4 <- FindClusters(GBM4, verbose = FALSE,resolution = 0.4)
p1<-SpatialPlot(GBM4, label = TRUE, label.size = 5)

#UMAP降维
GBM4 <- RunUMAP(GBM4, reduction = "pca", dims = 1:10)
p2 <- DimPlot(GBM4, reduction = "umap", label = TRUE)
p1+p2

ST5<-GBM4


spe2 = Read10X("E:/NPCsc/pipeline/GSE206245_RAW/ST6/out/filtered_feature_bc_matrix/")
image2 <- Read10X_Image(image.dir = file.path("E:/NPCsc/pipeline/GSE206245_RAW/ST6/out", 
                                              "spatial"), filter.matrix = TRUE)
spe2 <- CreateSeuratObject(counts = spe2, assay = "Spatial")

image2 <- image2[Cells(x = spe2)]
DefaultAssay(spe2 = image2) <- "Spatial"
spe2[["slice1"]] <- image2
SpatialFeaturePlot(spe2, features = "nFeature_Spatial")

GBM4<-spe2
GBM4 <- SCTransform(GBM4, assay = "Spatial", verbose = FALSE)
GBM4 <- RunPCA(GBM4, assay = "SCT", verbose = FALSE) 

##数据聚类
GBM4<- FindNeighbors(GBM4, reduction = "pca", dims = 1:10)
GBM4 <- FindClusters(GBM4, verbose = FALSE,resolution = 0.4)
p1<-SpatialPlot(GBM4, label = TRUE, label.size = 5)

#UMAP降维
GBM4 <- RunUMAP(GBM4, reduction = "pca", dims = 1:10)
p2 <- DimPlot(GBM4, reduction = "umap", label = TRUE)
p1+p2

ST6<-GBM4

spe2 = Read10X("E:/NPCsc/pipeline/GSE206245_RAW/ST8/out/filtered_feature_bc_matrix/")
image2 <- Read10X_Image(image.dir = file.path("E:/NPCsc/pipeline/GSE206245_RAW/ST8/out", 
                                              "spatial"), filter.matrix = TRUE)
spe2 <- CreateSeuratObject(counts = spe2, assay = "Spatial")

image2 <- image2[Cells(x = spe2)]
DefaultAssay(spe2 = image2) <- "Spatial"
spe2[["slice1"]] <- image2
SpatialFeaturePlot(spe2, features = "nFeature_Spatial")

GBM4<-spe2
GBM4 <- SCTransform(GBM4, assay = "Spatial", verbose = FALSE)
GBM4 <- RunPCA(GBM4, assay = "SCT", verbose = FALSE) 

##数据聚类
GBM4<- FindNeighbors(GBM4, reduction = "pca", dims = 1:10)
GBM4 <- FindClusters(GBM4, verbose = FALSE,resolution = 0.4)
p1<-SpatialPlot(GBM4, label = TRUE, label.size = 5)

#UMAP降维
GBM4 <- RunUMAP(GBM4, reduction = "pca", dims = 1:10)
p2 <- DimPlot(GBM4, reduction = "umap", label = TRUE)
p1+p2

ST8<-GBM4

spe2 = Read10X("E:/NPCsc/pipeline/GSE206245_RAW/ST7/out/filtered_feature_bc_matrix/")
image2 <- Read10X_Image(image.dir = file.path("E:/NPCsc/pipeline/GSE206245_RAW/ST7/out", 
                                              "spatial"), filter.matrix = TRUE)
spe2 <- CreateSeuratObject(counts = spe2, assay = "Spatial")

image2 <- image2[Cells(x = spe2)]
DefaultAssay(spe2 = image2) <- "Spatial"
spe2[["slice1"]] <- image2
SpatialFeaturePlot(spe2, features = "nFeature_Spatial")

GBM4<-spe2
GBM4 <- SCTransform(GBM4, assay = "Spatial", verbose = FALSE)
GBM4 <- RunPCA(GBM4, assay = "SCT", verbose = FALSE) 

##数据聚类
GBM4<- FindNeighbors(GBM4, reduction = "pca", dims = 1:10)
GBM4 <- FindClusters(GBM4, verbose = FALSE,resolution = 0.4)
p1<-SpatialPlot(GBM4, label = TRUE, label.size = 5)

#UMAP降维
GBM4 <- RunUMAP(GBM4, reduction = "pca", dims = 1:10)
p2 <- DimPlot(GBM4, reduction = "umap", label = TRUE)
p1+p2
ST7<-GBM4

spe2 = Read10X("E:/NPCsc/pipeline/GSE206245_RAW/ST9/out/filtered_feature_bc_matrix/")
image2 <- Read10X_Image(image.dir = file.path("E:/NPCsc/pipeline/GSE206245_RAW/ST9/out", 
                                              "spatial"), filter.matrix = TRUE)
spe2 <- CreateSeuratObject(counts = spe2, assay = "Spatial")

image2 <- image2[Cells(x = spe2)]
DefaultAssay(spe2 = image2) <- "Spatial"
spe2[["slice1"]] <- image2
SpatialFeaturePlot(spe2, features = "nFeature_Spatial")

GBM4<-spe2
GBM4 <- SCTransform(GBM4, assay = "Spatial", verbose = FALSE)
GBM4 <- RunPCA(GBM4, assay = "SCT", verbose = FALSE) 

##数据聚类
GBM4<- FindNeighbors(GBM4, reduction = "pca", dims = 1:10)
GBM4 <- FindClusters(GBM4, verbose = FALSE,resolution = 0.4)
p1<-SpatialPlot(GBM4, label = TRUE, label.size = 5)

#UMAP降维
GBM4 <- RunUMAP(GBM4, reduction = "pca", dims = 1:10)
p2 <- DimPlot(GBM4, reduction = "umap", label = TRUE)
p1+p2

ST9<-GBM4


spe2 = Read10X("E:/NPCsc/pipeline/GSE206245_RAW/ST10/out/filtered_feature_bc_matrix/")
image2 <- Read10X_Image(image.dir = file.path("E:/NPCsc/pipeline/GSE206245_RAW/ST10/out", 
                                              "spatial"), filter.matrix = TRUE)
spe2 <- CreateSeuratObject(counts = spe2, assay = "Spatial")

image2 <- image2[Cells(x = spe2)]
DefaultAssay(spe2 = image2) <- "Spatial"
spe2[["slice1"]] <- image2
SpatialFeaturePlot(spe2, features = "nFeature_Spatial")

GBM4<-spe2
GBM4 <- SCTransform(GBM4, assay = "Spatial", verbose = FALSE)
GBM4 <- RunPCA(GBM4, assay = "SCT", verbose = FALSE) 

##数据聚类
GBM4<- FindNeighbors(GBM4, reduction = "pca", dims = 1:10)
GBM4 <- FindClusters(GBM4, verbose = FALSE,resolution = 0.4)
p1<-SpatialPlot(GBM4, label = TRUE, label.size = 5)

#UMAP降维
GBM4 <- RunUMAP(GBM4, reduction = "pca", dims = 1:10)
  p2 <- DimPlot(GBM4, reduction = "umap", label = TRUE)
  p1+p2

ST10<-GBM4


spe2 = Read10X("E:/NPCsc/pipeline/GSE206245_RAW/ST14/out/filtered_feature_bc_matrix/")
image2 <- Read10X_Image(image.dir = file.path("E:/NPCsc/pipeline/GSE206245_RAW/ST14/out", 
                                              "spatial"), filter.matrix = TRUE)
spe2 <- CreateSeuratObject(counts = spe2, assay = "Spatial")

image2 <- image2[Cells(x = spe2)]
DefaultAssay(spe2 = image2) <- "Spatial"
spe2[["slice1"]] <- image2
SpatialFeaturePlot(spe2, features = "nFeature_Spatial")

GBM4<-spe2
GBM4 <- SCTransform(GBM4, assay = "Spatial", verbose = FALSE)
GBM4 <- RunPCA(GBM4, assay = "SCT", verbose = FALSE) 

##数据聚类
GBM4<- FindNeighbors(GBM4, reduction = "pca", dims = 1:10)
GBM4 <- FindClusters(GBM4, verbose = FALSE,resolution = 0.4)
p1<-SpatialPlot(GBM4, label = TRUE, label.size = 5)

#UMAP降维
GBM4 <- RunUMAP(GBM4, reduction = "pca", dims = 1:10)
p2 <- DimPlot(GBM4, reduction = "umap", label = TRUE)
p1+p2

ST14<-GBM4

spe2 = Read10X("E:/NPCsc/pipeline/GSE206245_RAW/ST11/out/filtered_feature_bc_matrix/")
image2 <- Read10X_Image(image.dir = file.path("E:/NPCsc/pipeline/GSE206245_RAW/ST11/out", 
                                              "spatial"), filter.matrix = TRUE)
spe2 <- CreateSeuratObject(counts = spe2, assay = "Spatial")

image2 <- image2[Cells(x = spe2)]
DefaultAssay(spe2 = image2) <- "Spatial"
spe2[["slice1"]] <- image2
SpatialFeaturePlot(spe2, features = "nFeature_Spatial")

GBM4<-spe2
GBM4 <- SCTransform(GBM4, assay = "Spatial", verbose = FALSE)
GBM4 <- RunPCA(GBM4, assay = "SCT", verbose = FALSE) 

##数据聚类
GBM4<- FindNeighbors(GBM4, reduction = "pca", dims = 1:10)
GBM4 <- FindClusters(GBM4, verbose = FALSE,resolution = 0.4)
p1<-SpatialPlot(GBM4, label = TRUE, label.size = 5)

#UMAP降维
GBM4 <- RunUMAP(GBM4, reduction = "pca", dims = 1:10)
p2 <- DimPlot(GBM4, reduction = "umap", label = TRUE)
p1+p2

ST11<-GBM4


spe2 = Read10X("E:/NPCsc/pipeline/GSE206245_RAW/ST12/out/filtered_feature_bc_matrix/")
image2 <- Read10X_Image(image.dir = file.path("E:/NPCsc/pipeline/GSE206245_RAW/ST12/out", 
                                              "spatial"), filter.matrix = TRUE)
spe2 <- CreateSeuratObject(counts = spe2, assay = "Spatial")

image2 <- image2[Cells(x = spe2)]
DefaultAssay(spe2 = image2) <- "Spatial"
spe2[["slice1"]] <- image2
SpatialFeaturePlot(spe2, features = "nFeature_Spatial")

GBM4<-spe2
GBM4 <- SCTransform(GBM4, assay = "Spatial", verbose = FALSE)
GBM4 <- RunPCA(GBM4, assay = "SCT", verbose = FALSE) 

##数据聚类
GBM4<- FindNeighbors(GBM4, reduction = "pca", dims = 1:10)
GBM4 <- FindClusters(GBM4, verbose = FALSE,resolution = 0.4)
p1<-SpatialPlot(GBM4, label = TRUE, label.size = 5)

#UMAP降维
GBM4 <- RunUMAP(GBM4, reduction = "pca", dims = 1:10)
p2 <- DimPlot(GBM4, reduction = "umap", label = TRUE)
p1+p2

ST12<-GBM4

#brain<-JoinLayers(brain,assay = 'Spatial')

#ST.FeaturePlot(GBM4, features = c("CXCL13",'ITGA2','SPP1'))

spe2 = Read10X("E:/NPCsc/pipeline/GSE206245_RAW/ST16/out/filtered_feature_bc_matrix/")
image2 <- Read10X_Image(image.dir = file.path("E:/NPCsc/pipeline/GSE206245_RAW/ST16/out", 
                                              "spatial"), filter.matrix = TRUE)
spe2 <- CreateSeuratObject(counts = spe2, assay = "Spatial")

image2 <- image2[Cells(x = spe2)]
DefaultAssay(spe2 = image2) <- "Spatial"
spe2[["slice1"]] <- image2
SpatialFeaturePlot(spe2, features = "nFeature_Spatial")

GBM4<-spe2
GBM4 <- SCTransform(GBM4, assay = "Spatial", verbose = FALSE)
GBM4 <- RunPCA(GBM4, assay = "SCT", verbose = FALSE) 

##数据聚类
GBM4<- FindNeighbors(GBM4, reduction = "pca", dims = 1:10)
GBM4 <- FindClusters(GBM4, verbose = FALSE,resolution = 0.4)
p1<-SpatialPlot(GBM4, label = TRUE, label.size = 5)

#UMAP降维
GBM4 <- RunUMAP(GBM4, reduction = "pca", dims = 1:10)
p2 <- DimPlot(GBM4, reduction = "umap", label = TRUE)
p1+p2

ST16<-GBM4

spe2 = Read10X("E:/NPCsc/pipeline/GSE206245_RAW/ST17/out/filtered_feature_bc_matrix/")
image2 <- Read10X_Image(image.dir = file.path("E:/NPCsc/pipeline/GSE206245_RAW/ST17/out", 
                                              "spatial"), filter.matrix = TRUE)
spe2 <- CreateSeuratObject(counts = spe2, assay = "Spatial")

image2 <- image2[Cells(x = spe2)]
DefaultAssay(spe2 = image2) <- "Spatial"
spe2[["slice1"]] <- image2
SpatialFeaturePlot(spe2, features = "nFeature_Spatial")

GBM4<-spe2
GBM4 <- SCTransform(GBM4, assay = "Spatial", verbose = FALSE)
GBM4 <- RunPCA(GBM4, assay = "SCT", verbose = FALSE) 

##数据聚类
GBM4<- FindNeighbors(GBM4, reduction = "pca", dims = 1:10)
GBM4 <- FindClusters(GBM4, verbose = FALSE,resolution = 0.4)
p1<-SpatialPlot(GBM4, label = TRUE, label.size = 5)

#UMAP降维
GBM4 <- RunUMAP(GBM4, reduction = "pca", dims = 1:10)
p2 <- DimPlot(GBM4, reduction = "umap", label = TRUE)
p1+p2

ST17<-GBM4
spe2 = Read10X("E:/NPCsc/pipeline/GSE206245_RAW/ST18/out/filtered_feature_bc_matrix/")
image2 <- Read10X_Image(image.dir = file.path("E:/NPCsc/pipeline/GSE206245_RAW/ST18/out", 
                                              "spatial"), filter.matrix = TRUE)
spe2 <- CreateSeuratObject(counts = spe2, assay = "Spatial")

image2 <- image2[Cells(x = spe2)]
DefaultAssay(spe2 = image2) <- "Spatial"
spe2[["slice1"]] <- image2
SpatialFeaturePlot(spe2, features = "nFeature_Spatial")

GBM4<-spe2
GBM4 <- SCTransform(GBM4, assay = "Spatial", verbose = FALSE)
GBM4 <- RunPCA(GBM4, assay = "SCT", verbose = FALSE) 

##数据聚类
GBM4<- FindNeighbors(GBM4, reduction = "pca", dims = 1:10)
GBM4 <- FindClusters(GBM4, verbose = FALSE,resolution = 0.4)
p1<-SpatialPlot(GBM4, label = TRUE, label.size = 5)

#UMAP降维
GBM4 <- RunUMAP(GBM4, reduction = "pca", dims = 1:10)
p2 <- DimPlot(GBM4, reduction = "umap", label = TRUE)
p1+p2

ST18<-GBM4

spe2 = Read10X("E:/NPCsc/pipeline/GSE206245_RAW/ST19/out/filtered_feature_bc_matrix/")
image2 <- Read10X_Image(image.dir = file.path("E:/NPCsc/pipeline/GSE206245_RAW/ST19/out", 
                                              "spatial"), filter.matrix = TRUE)
spe2 <- CreateSeuratObject(counts = spe2, assay = "Spatial")

image2 <- image2[Cells(x = spe2)]
DefaultAssay(spe2 = image2) <- "Spatial"
spe2[["slice1"]] <- image2
SpatialFeaturePlot(spe2, features = "nFeature_Spatial")

GBM4<-spe2
GBM4 <- SCTransform(GBM4, assay = "Spatial", verbose = FALSE)
GBM4 <- RunPCA(GBM4, assay = "SCT", verbose = FALSE) 

##数据聚类
GBM4<- FindNeighbors(GBM4, reduction = "pca", dims = 1:10)
GBM4 <- FindClusters(GBM4, verbose = FALSE,resolution = 0.4)
p1<-SpatialPlot(GBM4, label = TRUE, label.size = 5)

#UMAP降维
GBM4 <- RunUMAP(GBM4, reduction = "pca", dims = 1:10)
p2 <- DimPlot(GBM4, reduction = "umap", label = TRUE)
p1+p2

ST19<-GBM4

library(scCustomize)

save(T3,T7,file = 'S-ST-lung-1-23.RData')
