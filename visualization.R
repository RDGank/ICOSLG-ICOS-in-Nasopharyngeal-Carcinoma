
library(SCP)
library(Seurat)
load(file = 'sce_GSE150430.RData')
library(stringr)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(dplyr)
library(tidyr)
load(file = 'sce_GSE150430.RData')
sce[['RNA']]<-as(object = sce[['RNA']],Class = 'Assay')
sce$celltype<-str_replace(sce$celltype,'Tumor','Cancer')
p1<-CellDimPlot(sce,theme_use = 'theme_blank',
            show_stat = F,label.size = 3,label = T,group.by = 'celltype',
            title = 'GSE150430')
pdf(file = 'p1_GSE150430.pdf',height = 3.5,width = 4.5)
p1
dev.off()
#p2<-CellStatPlot(sce,stat.by = 'celltype',flip = T,group.by = 'Group',label = T)
#p2
#p2 percentage

#p3拷贝数变异,其中有一组用copykat做(因为没有正常上皮)


p4<-FeatureStatPlot(sce,stat.by =  c('ICOSLG'),stack = T, bg.by = 'celltype',add_point = T,pt.size = 0.000001,
                    group.by = 'celltype',flip = T,title = 'GSE150430')+labs(x=NULL,y=NULL)
#p5<-DotPlot(sce,features = 'ICOSLG',group.by = 'celltype',assay = 'RNA')
p5 <- DotPlot(sce, c('ICOSLG' ),group.by = 'celltype',assay = 'RNA' ) + 
  #coord_flip() + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 0, hjust = 0.5,vjust=0.5))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33')) #
p5
pdf(file = 'p4_GSE150430.pdf',height =8 ,width =6)
p4/p5
dev.off()

#save(sce,file = 'sce_GSE150430.RData')

load(file = 'sce_GSE150430_Tcell.RData')
sce[['RNA']]<-as(object = sce[['RNA']],Class = 'Assay')
p1<-CellDimPlot(sce,theme_use = 'theme_blank',
                show_stat = F,label.size = 3,label = T,group.by = 'celltype',
                title = 'GSE150430')
pdf(file = 'p1_T_GSE150430.pdf',height = 3.5,width = 4.5)
p1
dev.off()
#p2<-CellStatPlot(sce,stat.by = 'celltype',flip = T,group.by = 'Group',label = T)
#p2
#p2 percentage

#p3拷贝数变异,其中有一组用copykat做(因为没有正常上皮)


p4<-FeatureStatPlot(sce,stat.by =  c('ICOS'),stack = T, bg.by = 'celltype',add_point = T,pt.size = 0.000001,
                    group.by = 'celltype',flip = T,title = 'GSE150430')+labs(x=NULL,y=NULL)
#p5<-DotPlot(sce,features = 'ICOSLG',group.by = 'celltype',assay = 'RNA')
p5 <- DotPlot(sce, c('ICOS' ),group.by = 'celltype',assay = 'RNA' ) + 
  #coord_flip() + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 0, hjust = 0.5,vjust=0.5))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33')) #
p5
pdf(file = 'p4_T_GSE150430.pdf',height =8 ,width =6)
p4/p5
dev.off()










load(file = 'sce_GSE120926.RData')
sce[['RNA']]<-as(object = sce[['RNA']],Class = 'Assay')
sce$celltype<-str_replace(sce$celltype,'Tumor','Cancer')
p1<-CellDimPlot(sce,theme_use = 'theme_blank',
                show_stat = F,label.size = 3,label = T,group.by = 'celltype',
                title = 'GSE120926')
pdf(file = 'p1_GSE120926.pdf',height = 3.5,width = 4.5)
p1
dev.off()
#p2<-CellStatPlot(sce,stat.by = 'celltype',flip = T,group.by = 'Group',label = T)
#p2
#p2 percentage

#p3拷贝数变异,其中有一组用copykat做(因为没有正常上皮)


p4<-FeatureStatPlot(sce,stat.by =  c('ICOSLG'),stack = T, bg.by = 'celltype',add_point = T,pt.size = 0.000001,
                    group.by = 'celltype',flip = T,title = 'GSE120926')+labs(x=NULL,y=NULL)
#p5<-DotPlot(sce,features = 'ICOSLG',group.by = 'celltype',assay = 'RNA')
p5 <- DotPlot(sce, c('ICOSLG' ),group.by = 'celltype',assay = 'RNA' ) + 
  #coord_flip() + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 0, hjust = 0.5,vjust=0.5))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33')) #
p5
pdf(file = 'p4_GSE120926.pdf',height =8 ,width =6)
p4/p5
dev.off()

#save(sce,file = 'sce_GSE120926.RData')

load(file = 'sce_GSE120926_Tcell.RData')
sce[['RNA']]<-as(object = sce[['RNA']],Class = 'Assay')
p1<-CellDimPlot(sce,theme_use = 'theme_blank',
                show_stat = F,label.size = 3,label = T,group.by = 'celltype',
                title = 'GSE120926')
pdf(file = 'p1_T_GSE120926.pdf',height = 3.5,width = 4.5)
p1
dev.off()
#p2<-CellStatPlot(sce,stat.by = 'celltype',flip = T,group.by = 'Group',label = T)
#p2
#p2 percentage

#p3拷贝数变异,其中有一组用copykat做(因为没有正常上皮)


p4<-FeatureStatPlot(sce,stat.by =  c('ICOS'),stack = T, bg.by = 'celltype',add_point = T,pt.size = 0.000001,
                    group.by = 'celltype',flip = T,title = 'GSE120926')+labs(x=NULL,y=NULL)
#p5<-DotPlot(sce,features = 'ICOSLG',group.by = 'celltype',assay = 'RNA')
p5 <- DotPlot(sce, c('ICOS' ),group.by = 'celltype',assay = 'RNA' ) + 
  #coord_flip() + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 0, hjust = 0.5,vjust=0.5))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33')) #
p5
pdf(file = 'p4_T_GSE120926.pdf',height =8 ,width =6)
p4/p5
dev.off()





load(file = 'sce_GSE150825.RData')
sce[['RNA']]<-as(object = sce[['RNA']],Class = 'Assay')
sce$celltype<-str_replace(sce$celltype,'Tumor','Cancer')
p1<-CellDimPlot(sce,theme_use = 'theme_blank',
                show_stat = F,label.size = 3,label = T,group.by = 'celltype',
                title = 'GSE150825')
pdf(file = 'p1_GSE150825.pdf',height = 3.5,width = 4.5)
p1
dev.off()
#p2<-CellStatPlot(sce,stat.by = 'celltype',flip = T,group.by = 'Group',label = T)
#p2
#p2 percentage

#p3拷贝数变异,其中有一组用copykat做(因为没有正常上皮)


p4<-FeatureStatPlot(sce,stat.by =  c('ICOSLG'),stack = T, bg.by = 'celltype',add_point = T,pt.size = 0.000001,
                    group.by = 'celltype',flip = T,title = 'GSE150825')+labs(x=NULL,y=NULL)
#p5<-DotPlot(sce,features = 'ICOSLG',group.by = 'celltype',assay = 'RNA')
p5 <- DotPlot(sce, c('ICOSLG' ),group.by = 'celltype',assay = 'RNA' ) + 
  #coord_flip() + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 0, hjust = 0.5,vjust=0.5))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33')) #
p5
pdf(file = 'p4_GSE150825.pdf',height =8 ,width =6)
p4/p5
dev.off()

#save(sce,file = 'sce_GSE150825.RData')

load(file = 'sce_GSE150825_Tcell.RData')
sce[['RNA']]<-as(object = sce[['RNA']],Class = 'Assay')
p1<-CellDimPlot(sce,theme_use = 'theme_blank',
                show_stat = F,label.size = 3,label = T,group.by = 'celltype',
                title = 'GSE150825')
pdf(file = 'p1_T_GSE150825.pdf',height = 3.5,width = 4.5)
p1
dev.off()
#p2<-CellStatPlot(sce,stat.by = 'celltype',flip = T,group.by = 'Group',label = T)
#p2
#p2 percentage

#p3拷贝数变异,其中有一组用copykat做(因为没有正常上皮)


p4<-FeatureStatPlot(sce,stat.by =  c('ICOS'),stack = T, bg.by = 'celltype',add_point = T,pt.size = 0.000001,
                    group.by = 'celltype',flip = T,title = 'GSE150825')+labs(x=NULL,y=NULL)
#p5<-DotPlot(sce,features = 'ICOSLG',group.by = 'celltype',assay = 'RNA')
p5 <- DotPlot(sce, c('ICOS' ),group.by = 'celltype',assay = 'RNA' ) + 
  #coord_flip() + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 0, hjust = 0.5,vjust=0.5))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33')) #
p5
pdf(file = 'p4_T_GSE150825.pdf',height =8 ,width =6)
p4/p5
dev.off()






load(file = 'sce_GSE162025.RData')
sce[['RNA']]<-as(object = sce[['RNA']],Class = 'Assay')
sce$celltype<-str_replace(sce$celltype,'Tumor','Cancer')
p1<-CellDimPlot(sce,theme_use = 'theme_blank',
                show_stat = F,label.size = 3,label = T,group.by = 'celltype',
                title = 'GSE162025')
pdf(file = 'p1_GSE162025.pdf',height = 3.5,width = 4.5)
p1
dev.off()
#p2<-CellStatPlot(sce,stat.by = 'celltype',flip = T,group.by = 'Group',label = T)
#p2
#p2 percentage

#p3拷贝数变异,其中有一组用copykat做(因为没有正常上皮)


p4<-FeatureStatPlot(sce,stat.by =  c('ICOSLG'),stack = T, bg.by = 'celltype',add_point = T,pt.size = 0.000001,
                    group.by = 'celltype',flip = T,title = 'GSE162025')+labs(x=NULL,y=NULL)
#p5<-DotPlot(sce,features = 'ICOSLG',group.by = 'celltype',assay = 'RNA')
p5 <- DotPlot(sce, c('ICOSLG' ),group.by = 'celltype',assay = 'RNA' ) + 
  #coord_flip() + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 0, hjust = 0.5,vjust=0.5))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33')) #
p5
pdf(file = 'p4_GSE162025.pdf',height =8 ,width =6)
p4/p5
dev.off()

#save(sce,file = 'sce_GSE162025.RData')

load(file = 'sce_GSE162025_Tcell.RData')
sce[['RNA']]<-as(object = sce[['RNA']],Class = 'Assay')
p1<-CellDimPlot(sce,theme_use = 'theme_blank',
                show_stat = F,label.size = 3,label = T,group.by = 'celltype',
                title = 'GSE162025')
pdf(file = 'p1_T_GSE162025.pdf',height = 3.5,width = 4.5)
p1
dev.off()
#p2<-CellStatPlot(sce,stat.by = 'celltype',flip = T,group.by = 'Group',label = T)
#p2
#p2 percentage

#p3拷贝数变异,其中有一组用copykat做(因为没有正常上皮)


p4<-FeatureStatPlot(sce,stat.by =  c('ICOS'),stack = T, bg.by = 'celltype',add_point = T,pt.size = 0.000001,
                    group.by = 'celltype',flip = T,title = 'GSE162025')+labs(x=NULL,y=NULL)
#p5<-DotPlot(sce,features = 'ICOSLG',group.by = 'celltype',assay = 'RNA')
p5 <- DotPlot(sce, c('ICOS' ),group.by = 'celltype',assay = 'RNA' ) + 
  #coord_flip() + #翻转
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 0, hjust = 0.5,vjust=0.5))+ #轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33')) #
p5
pdf(file = 'p4_T_GSE162025.pdf',height =8 ,width =6)
p4/p5
dev.off()

