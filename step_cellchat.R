# devtools::install_github("sqjin/CellChat")
rm(list = ls())
load("sce.rds")
library(CellChat)
library(Seurat)
library(Seurat)
load(file = 'cellchat_GSE150430.RData')
options(stringsAsFactors = FALSE)
load(file = 'sce_GSE150430.RData')
sce1<-sce
table(sce1$celltype)
sce1<-sce1[,sce1$celltype%in%c('B','Cancer cells','pDC','Epithelial cells')]
load(file = 'sce_GSE150430_Tcell.RData')
sce<-merge(sce1,sce)
Idents(sce)<-'celltype'
#sce<-JoinLayers(sce,assay = 'RNA')
data.input <- GetAssayData(sce, assay = "RNA", slot = "data") # normalized data matrix
labels <- Idents(sce)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")

cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "group") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
# set the used database in the object
cellchat@DB <- CellChatDB.use
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
#future::plan("multiprocess", workers = 4) # do scerallel
cellchat <- identifyOverExpressedGenes(cellchat,do.fast = F)
cellchat <- identifyOverExpressedInteractions(cellchat)
# project gene expression data onto PPI network (optional)
#cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
netVisual_bubble(cellchat, 
                 sources.use = ,
                 signaling = 'ICOS',
                 targets.use = , 
                 remove.isolate = FALSE)
save(sce,cellchat,file = 'cellchat_GSE150430.RData')
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
dev.off()
netVisual_circle(cellchat@net$count, 
                 vertex.weight = groupSize, 
                 weight.scale = T, 
                 label.edge= T, 
                 title.name = "Number of interactions")


netVisual_circle(cellchat@net$weight, 
                 vertex.weight = groupSize,
                 weight.scale = T, 
                 label.edge= F, 
                 title.name = "Interaction weights/strength")

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
 
pdf(file = 'our-heat.pdf',width = 9.7,height = 13.5)
netAnalysis_signalingRole_heatmap(cellchat, pattern = "all",height = 28,width = 20)
dev.off()

netVisual_heatmap(cellchat)
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

library(tidyverse)
df.net <- subsetCommunication(cellchat,thresh = 0.05)
df<-df.net%>%filter(source=="TSK")
df<-df.net%>%filter(target=="Fibroblast")

dev.off()
pathways.show <- c("MIF") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling =pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
#> Note: The first link end is drawn out of sector 'Inflam. FIB'.


dev.off()
unique(labels)
pathway.show<-c("GALECTIN","VCAM"  ,   "LIGHT" ,    "FASLG", "CD40" ,  "PECAM1" ,  "NCAM"  ,"CD137"    )
netVisual_bubble(cellchat, 
                 sources.use = ,
                 signaling = 'ICOS',
                 targets.use = , 
                 remove.isolate = FALSE)

pdf(file = '8n-cellchat3.pdf',width = 11,height = 6)

PATH<-c("CD86")
netVisual_bubble(cellchat, 
                 sources.use =c(9), 
                 targets.use =, 
                 #signaling = PATH,
                 remove.isolate = FALSE)
dev.off()

PATH<-c(#"MHC-I",
        #"MHC-II",
        #"MIF",
        "GALECTIN",
        "CD99",
        "APP",
        #"MK",
        #"COLLAGEN",
        #"LAMININ",
        #"ITGB2",
        "CLEC",
        "CD45",
        "CD22",
        #"ICAM",
        "SELPLG",
        "LCK",
        "ANNEXIN",
       # "CCL",
        "FN1",
        "ADGRE5",
        #"PTN",
        "BAFF",
        "IL16",
        #"CXCL",
        "SPP1",
        "RESISTIN",
        "JAM",
        "SEMA4",
        "PARs",
        "VCAM",
        "TNF",
        "CD48",
        "ALCAM",
        "CD6",
        "VISFATIN",
        "CDH1",
        "IFN-II",
        "CD86",
        "CD70",
        "GRN",
        "CD23",
        "COMPLEMENT",
        "PECAM1",
        "EGF",
        "SELL",
        "LT",
        "THBS",
        "CDH",
        "CEACAM",
        "ICOS",
        "CD40",
        "MPZ",
        "DESMOSOME",
        "CD46",
        "APRIL",
        "THY1",
        "SN",
        "GAS",
        "BTLA",
        "CSF",
        "NOTCH",
        "BAG",
        #"EPHA",
        "LIGHT",
        "SEMA7",
        "KIT",
        "WNT",
        "IL1",
        "CX3C",
        "GDF",
        "CD80",
        "AGRN",
        "TENASCIN",
        "FASLG",
        "HSPG",
        "EPHB",
        "ANGPTL",
        "NRG",
        "TWEAK",
        "SEMA3",
        "FLT3",
        "IL2",
        "CADM",
        "CD137",
        "NT",
        "NECTIN",
        "TIGIT",
        "FGF",
        "IL10",
        "SEMA6",
        "SEMA5",
        "SAA",
        "OSM",
        "PD-L1",
        "PDL2",
        "CHEMERIN",
        "IL6",
        "PERIOSTIN",
        "BMP",
        "ANGPT",
        "CALCR",
        "ncWNT",
        "OCLN",
        "PDGF",
        "VEGI",
        "LIFR",
        "NCAM",
        "CD30",
        "ACTIVIN")

netVisual_chord_gene(cellchat, 
                     #sources.use  , 
                     #targets.use, 
                     lab.cex =0.1,
                     legend.pos.y = 20,)



StackedVlnPlot(sce,features = c("EPCAM","PTPRC"))

load("position.Rds")
load("stRNA.Rds")
embed_umap2=data.frame(UMAP_1=position_sub_sub$x,
                       UMAP_2=position_sub_sub$y,
                       row.names = rownames(position_sub_sub))

stRNA@reductions$umap@cell.embeddings=as.matrix(embed_umap2)
DimPlot(stRNA)
FeaturePlot(stRNA,"CD99")

saveRDS(cellchat, file = "cellchat.rds")



pdf(file = 'Tchat.pdf',height = 7,width = 14)
p1<-netVisual_heatmap(cellchat, color.heatmap = 'Reds',)
p2<-netVisual_heatmap(cellchat, color.heatmap = 'Reds',measure = 'weight')
p1+p2
dev.off()
