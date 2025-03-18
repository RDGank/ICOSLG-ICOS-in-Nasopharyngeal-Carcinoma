p1<-SpatialPlot(T4, label = TRUE, label.size = 5,pt.size.factor = 4)
p2 <- DimPlot(T4, reduction = "umap", label = TRUE)
p1+p2
p3<-SpatialFeaturePlot(T4,features = c('CD68','SPP1','EPCAM','ITGA2'),pt.size.factor = 3,ncol = 4)
#p3<-
genes<-c('CD68','SPP1','EPCAM','ITGA2')
genes1<-c('CD68','SPP1')
genes2<-c('EPCAM','ITGA2')

T4<-AddModuleScore(T4,features = as.list(genes),name = 'Col-NSCLC(ITGA2)_Mac(SPP1)')
T4<-AddModuleScore(T4,features = as.list(genes1),name = 'Abundance_Mac-SPP')
T4<-AddModuleScore(T4,features = as.list(genes2),name = 'Abundance_NSCLC(ITGA2)_')
p4<-SpatialFeaturePlot(T4,features = c('Col-NSCLC(ITGA2)_Mac(SPP1)1'),pt.size.factor = 3)

metadata<-T4@meta.data
cor.test(metadata$`Abundance_Mac-SPP1`,metadata$`Abundance_NSCLC(ITGA2)_1`,method = 'pearson',use="complete.obs")
#Mal_TUBB2B_abundance <- metadata$Treg_TNFRSF4_abundance1
#ICOSLG_abundance<- metadata$ICOS_1
#plot(Mal_TUBB2B_abundance, ICOSLG_abundance, main = "r=0.17,p-value p-value = 1.465e-10",
 #    xlab = "Treg_TNFRSF4_abundance", ylab = "ICOS_abundance",
  #   pch = 18,frame = F)
# 添加回归线
#abline(lm(y ~ x, data = metadata), col = "red")
#pdf(file = 'cor-TregFOX-CTLA4.pdf',height = 3.5,width = 4)
#dev.off()

ggplot(metadata, aes(x=metadata$`Abundance_Mac-SPP1`, y=metadata$`Abundance_NSCLC(ITGA2)_1`)) + 
  geom_point(color="#A1CB4E",size=1)+ geom_smooth(method = 'lm', se = F, color = "#4DBBD5FF")+theme_bw()+
  stat_cor(data=metadata, method = "pearson")+
  theme_bw()+
  labs(x='Treg(FOXP3+) abundance',y='CTLA4')
px
p0<-SpatialFeaturePlot(T4, features = "nFeature_Spatial",pt.size.factor = 4)
px<-p0+p1+p4+plot_layout(design = 'ABC')
py<-(p2+p3)+plot_layout(design = 'ABBBB')
pdf(file = 'col-Mal_macro-ST-T4.pdf',height = 6,width = 12)
px/py
dev.off()

