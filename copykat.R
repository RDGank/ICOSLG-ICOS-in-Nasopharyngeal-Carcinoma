library(devtools)
install_github("navinlabcode/copykat")
library(copykat)
library(Seurat)
library(stringr)
load(file = 'sce_GSE160.RData')
table(sce$celltype)
sce$celltype<-str_replace(sce$celltype,'Tumor','Cancer')
scRNA <- subset(sce,celltype==c("Epithelial cells","Cancer cells"))
scRNA <- subset(scRNA,downsample = 1000)

exp.rawdata <- as.matrix(scRNA@assays$RNA@counts)

copykat.test <- copykat(rawmat=exp.rawdata, 
                        id.type="S", # S是symbol的含义
                        ngene.chr=5, # 规定每条染色体中至少5个基因来计算DNA拷贝数(可调整)
                        win.size=25, # 每个片段至少25个基因
                        KS.cut=0.1, # 增加KS.cut会降低敏感度，通常范围在0.05-0.15
                        sam.name="241016",  # 自行命名
                        distance="euclidean", 
                        norm.cell.names="",
                        output.seg="FLASE", 
                        plot.genes="TRUE", 
                        genome="hg20",
                        n.cores=1)

