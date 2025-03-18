library(ktplots)
library(SingleCellExperiment)
library(Seurat)
library(ggplot2)
load(file = 'E:/NM30/scp/cellchat_prospective_post.RData')
pvals <- read.delim("cellphone_prospective_post/out/pvalues.txt", check.names = FALSE)
means <- read.delim("cellphone_prospective_post/out/means.txt", check.names = FALSE)
decon = read.delim("cellphone_prospective_post/out/deconvoluted.txt")

load(file = 'cellchat_GSE150430.RData')
pvals <- read.delim("cellphone_GSE150430/out/pvalues.txt", check.names = FALSE)
means <- read.delim("cellphone_GSE150430/out/means.txt", check.names = FALSE)
decon = read.delim("cellphone_GSE150430/out/deconvoluted.txt")
sce  <- as.SingleCellExperiment(sce)

load(file = 'cellchat_GSE150430.RData')
pvals <- read.delim("cellphone_GSE150430/out/pvalues.txt", check.names = FALSE)
means <- read.delim("cellphone_GSE150430/out/means.txt", check.names = FALSE)
decon = read.delim("cellphone_GSE150430/out/deconvoluted.txt")
sce  <- as.SingleCellExperiment(sce)

c('#330066','#336699','#66CC66','#FFCC33')

p1<-plot_cpdb_heatmap(pvals = pvals, cellheight = 20,
                  cellwidth = 20, low_col = "#336699",
                  mid_col = '#66CC66',
                  high_col = '#FFCC33',)

p2<-plot_cpdb(
  scdata = sce,min_interaction_score = 1,
  highlight_col = 'red',keep_significant_only = T,
  max_size=8,
  max_highlight_size=10,
  #cell_type1 = 'B',
  #cell_type1 = 'pDC',
  #cell_type1 = 'Cancer cells',
  cell_type1 = "B|Cancer cells|pDC",
  #cell_type2 = c('Treg','Tfh'), # this means all cell-types
  cell_type2 = 'Treg|Tfh',
  celltype_key = "celltype",
  means = means,
  pvals = pvals,
  #gene_family = 'Treg',
  genes = c("ICOS",'TGFB'),
  title = "Results of GSE150430 (CellphoneDB)",
)+theme(axis.text  = element_text(size = 10, color = 'black'))

library(patchwork)
library(ggplot2)
library(ggplotify)
p1<-as.ggplot(p1)
p1+p2

pdf(file = 'cellphone_GSE150430.pdf',height =5 ,width =14.5 )
p1+p2
dev.off()


load(file = 'cellchat_GSE120926.RData')
pvals <- read.delim("cellphone_GSE120926/out/pvalues.txt", check.names = FALSE)
means <- read.delim("cellphone_GSE120926/out/means.txt", check.names = FALSE)
decon = read.delim("cellphone_GSE120926/out/deconvoluted.txt")
sce  <- as.SingleCellExperiment(sce)

c('#330066','#336699','#66CC66','#FFCC33')

p1<-plot_cpdb_heatmap(pvals = pvals, cellheight = 20,
                      cellwidth = 20, low_col = "#336699",
                      mid_col = '#66CC66',
                      high_col = '#FFCC33',)

p2<-plot_cpdb(
  scdata = sce,min_interaction_score = 1,
  highlight_col = 'red',keep_significant_only = T,
  max_size=8,
  max_highlight_size=10,
  #cell_type1 = 'B',
  #cell_type1 = 'pDC',
  #cell_type1 = 'Cancer cells',
  cell_type1 = "B|Cancer cells|pDC",
  #cell_type2 = c('Treg','Tfh'), # this means all cell-types
  cell_type2 = 'Treg|Tfh',
  celltype_key = "celltype",
  means = means,
  pvals = pvals,
  #gene_family = 'Treg',
  genes = c("ICOS",'TGFB'),
  title = "Results of GSE120926 (CellphoneDB)",
)+theme(axis.text  = element_text(size = 10, color = 'black'))

library(patchwork)
library(ggplot2)
library(ggplotify)
p1<-as.ggplot(p1)
p1+p2

pdf(file = 'cellphone_GSE120926.pdf',height =6 ,width =14.5 )
p1+p2
dev.off()


load(file = 'cellchat_GSE162025.RData')
pvals <- read.delim("cellphone_GSE162025/out/pvalues.txt", check.names = FALSE)
means <- read.delim("cellphone_GSE162025/out/means.txt", check.names = FALSE)
decon = read.delim("cellphone_GSE162025/out/deconvoluted.txt")
sce  <- as.SingleCellExperiment(sce)

c('#330066','#336699','#66CC66','#FFCC33')

p1<-plot_cpdb_heatmap(pvals = pvals, cellheight = 20,
                      cellwidth = 20, low_col = "#336699",
                      mid_col = '#66CC66',
                      high_col = '#FFCC33',)

p2<-plot_cpdb(
  scdata = sce,min_interaction_score = 1,
  highlight_col = 'red',keep_significant_only = T,
  max_size=8,
  max_highlight_size=10,
  #cell_type1 = 'B',
  #cell_type1 = 'pDC',
  #cell_type1 = 'Cancer cells',
  cell_type1 = "B|Cancer cells|pDC",
  #cell_type2 = c('Treg','Tfh'), # this means all cell-types
  cell_type2 = 'Treg|Tfh',
  celltype_key = "celltype",
  means = means,
  pvals = pvals,
  #gene_family = 'Treg',
  genes = c("ICOS",'TGFB'),
  title = "Results of GSE162025 (CellphoneDB)",
)+theme(axis.text  = element_text(size = 10, color = 'black'))

library(patchwork)
library(ggplot2)
library(ggplotify)
p1<-as.ggplot(p1)
p1+p2

pdf(file = 'cellphone_GSE162025.pdf',height =5 ,width =14.5 )
p1+p2
dev.off()


load(file = 'E:/NM30/scp/cellchat_prospective_post.RData')
pvals <- read.delim("cellphone_prospective_post/out/pvalues.txt", check.names = FALSE)
means <- read.delim("cellphone_prospective_post/out/means.txt", check.names = FALSE)
decon = read.delim("cellphone_prospective_post/out/deconvoluted.txt")
sce  <- as.SingleCellExperiment(sce)

c('#330066','#336699','#66CC66','#FFCC33')

p1<-plot_cpdb_heatmap(pvals = pvals, cellheight = 20,
                      cellwidth = 20, low_col = "#336699",
                      mid_col = '#66CC66',
                      high_col = '#FFCC33',)

p2<-plot_cpdb(
  scdata = sce,min_interaction_score = 1,
  highlight_col = 'red',keep_significant_only = T,
  max_size=8,
  max_highlight_size=10,
  #cell_type1 = 'B',
  #cell_type1 = 'pDC',
  #cell_type1 = 'Cancer cells',
  cell_type1 = "B|Cancer cells|pDC",
  #cell_type2 = c('Treg','Tfh'), # this means all cell-types
  cell_type2 = 'Treg|Tfh',
  celltype_key = "celltype",
  means = means,
  pvals = pvals,
  #gene_family = 'Treg',
  genes = c("ICOS",'TGFB'),
  title = "Results of Post-treatment (CellphoneDB)",
)+theme(axis.text  = element_text(size = 10, color = 'black'))

library(patchwork)
library(ggplot2)
library(ggplotify)
p1<-as.ggplot(p1)
p1+p2

pdf(file = 'cellphone_Post_treatment.pdf',height =5.5 ,width =14.5 )
p1+p2
dev.off()


load(file = 'E:/NM30/scp/cellchat_prospective_pre.RData')
pvals <- read.delim("cellphone_prospective_pre/out/pvalues.txt", check.names = FALSE)
means <- read.delim("cellphone_prospective_pre/out/means.txt", check.names = FALSE)
decon = read.delim("cellphone_prospective_pre/out/deconvoluted.txt")
sce  <- as.SingleCellExperiment(sce)

c('#330066','#336699','#66CC66','#FFCC33')

p1<-plot_cpdb_heatmap(pvals = pvals, cellheight = 20,
                      cellwidth = 20, low_col = "#336699",
                      mid_col = '#66CC66',
                      high_col = '#FFCC33',)

p2<-plot_cpdb(
  scdata = sce,min_interaction_score = 1,
  highlight_col = 'red',keep_significant_only = T,
  max_size=8,
  max_highlight_size=10,
  #cell_type1 = 'B',
  #cell_type1 = 'pDC',
  #cell_type1 = 'Cancer cells',
  cell_type1 = "B|Cancer cells|pDC",
  #cell_type2 = c('Treg','Tfh'), # this means all cell-types
  cell_type2 = 'Treg|Tfh',
  celltype_key = "celltype",
  means = means,
  pvals = pvals,
  #gene_family = 'Treg',
  genes = c("ICOS",'TGFB'),
  title = "Results of Pre-treatment (CellphoneDB)",
)+theme(axis.text  = element_text(size = 10, color = 'black'))

library(patchwork)
library(ggplot2)
library(ggplotify)
p1<-as.ggplot(p1)
p1+p2

pdf(file = 'cellphone_Pre_treatment.pdf',height =5.5 ,width =14.5 )
p1+p2
dev.off()



library(CellChat)
load(file = 'cellchat_GSE120926.RData')
netVisual_bubble(cellchat, 
                 sources.use = c('Cancer cells','B','pDC'),
                 #signaling = 'ICOS',
                 targets.use = c('Treg','Tfh'), 
                 remove.isolate = FALSE)

netVisual_heatmap(cellchat, color.heatmap = 'Reds',measure = 'weight')



load(file = 'cellchat_GSE150430.RData')
p1<-netVisual_bubble(cellchat, 
                 sources.use = c('Cancer cells','B','pDC'),
                 signaling = 'ICOS',
                title.name = 'Results of GSE150430 (CellChatDB)',
                 targets.use = c('Treg','Tfh'), 
                 remove.isolate = FALSE,angle.x = 45)#+element_text(angle = 45,vjust = 0.5,hjust = 0.5)

pdf(file = 'cellchat_GSE150430.pdf',height = 2.5,width = 5)
p1
dev.off()




load(file = 'cellchat_GSE150825.RData')
p1<-netVisual_bubble(cellchat, 
                     sources.use = c('Cancer cells','B','pDC'),
                     signaling = 'ICOS',
                     title.name = 'Results of Post-treatment (CellChatDB)',
                     targets.use = c('Treg','Tfh'), 
                     remove.isolate = FALSE,angle.x = 45)#+element_text(angle = 45,vjust = 0.5,hjust = 0.5)

pdf(file = 'cellchat_GSE150825.pdf',height = 2.5,width = 5)
p1
dev.off()

load(file = 'cellchat_GSE162025.RData')
p1<-netVisual_bubble(cellchat, 
                     sources.use = c('Cancer cells','B','pDC'),
                     signaling = 'ICOS',
                     title.name = 'Results of GSE162025 (CellChatDB)',
                     targets.use = c('Treg','Tfh'), 
                     remove.isolate = FALSE,angle.x = 45)#+element_text(angle = 45,vjust = 0.5,hjust = 0.5)

pdf(file = 'cellchat_GSE162025.pdf',height = 2.5,width = 5)
p1
dev.off()


load(file = 'E:/NM30/scp/cellchat_prospective_pre.RData')
p1<-netVisual_bubble(cellchat, 
                     sources.use = c('Cancer cells','B','pDC'),
                     signaling = 'ICOS',
                     title.name = 'Results of Pre-treatment (CellChatDB)',
                     targets.use = c('Treg','Tfh'), 
                     remove.isolate = FALSE,angle.x = 45)#+element_text(angle = 45,vjust = 0.5,hjust = 0.5)

pdf(file = 'cellchat_pre.pdf',height = 2.5,width = 5)
p1
dev.off()


load(file = 'E:/NM30/scp/cellchat_prospective_post.RData')
p1<-netVisual_bubble(cellchat, 
                     sources.use = c('Cancer cells','B','pDC'),
                     signaling = 'ICOS',
                     title.name = 'Results of Post-treatment (CellChatDB)',
                     targets.use = c('Treg','Tfh'), 
                     remove.isolate = FALSE,angle.x = 45)#+element_text(angle = 45,vjust = 0.5,hjust = 0.5)

pdf(file = 'cellchat_post.pdf',height = 2.5,width = 5)
p1
dev.off()
