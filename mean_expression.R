library(Seurat)
library(stringr)
load(file = 'sce_GSE150430.RData')
ICOSLG<-AverageExpression(object = sce,features = 'ICOSLG',group.by = 'celltype')
ICOSLG<-unlist(ICOSLG)
ICOSLG<-as.data.frame(ICOSLG)
ICOSLG<-t(ICOSLG)
ICOSLG<-as.data.frame(ICOSLG)
colnames(ICOSLG)<-'GSE150430'
ICOSLG_GSE150430<-ICOSLG

load(file = 'sce_GSE150825.RData')
sce$celltype<-str_replace(sce$celltype,'Tumor','Cancer')
ICOSLG<-AverageExpression(object = sce,features = 'ICOSLG',group.by = 'celltype')
ICOSLG<-unlist(ICOSLG)
ICOSLG<-as.data.frame(ICOSLG)
ICOSLG<-t(ICOSLG)
ICOSLG<-as.data.frame(ICOSLG)
colnames(ICOSLG)<-'GSE150825'
ICOSLG_GSE150825<-ICOSLG

load(file = 'sce_GSE162025.RData')
sce$celltype<-str_replace(sce$celltype,'Tumor','Cancer')
ICOSLG<-AverageExpression(object = sce,features = 'ICOSLG',group.by = 'celltype')
ICOSLG<-unlist(ICOSLG)
ICOSLG<-as.data.frame(ICOSLG)
ICOSLG<-t(ICOSLG)
ICOSLG<-as.data.frame(ICOSLG)
colnames(ICOSLG)<-'GSE162025'
ICOSLG_GSE162025<-ICOSLG


load(file = 'sce_GSE120926.RData')
sce$celltype<-str_replace(sce$celltype,'Tumor','Cancer')
ICOSLG<-AverageExpression(object = sce,features = 'ICOSLG',group.by = 'celltype')
ICOSLG<-unlist(ICOSLG)
ICOSLG<-as.data.frame(ICOSLG)
ICOSLG<-t(ICOSLG)
ICOSLG<-as.data.frame(ICOSLG)
colnames(ICOSLG)<-'GSE120926'
ICOSLG_GSE120926<-ICOSLG

load(file = 'E:/NM30/scp/sce_prospective.RData')
sce1<-sce
sce<-sce1[,sce1$Group%in%c('Pre-treatment')]
sce$celltype<-str_replace(sce$celltype,'Tumor','Cancer')
ICOSLG<-AverageExpression(object = sce,features = 'ICOSLG',group.by = 'celltype')
ICOSLG<-unlist(ICOSLG)
ICOSLG<-as.data.frame(ICOSLG)
ICOSLG<-t(ICOSLG)
ICOSLG<-as.data.frame(ICOSLG)
colnames(ICOSLG)<-'Pre-treatment'
ICOSLG_Pre_treatment<-ICOSLG

sce<-sce1[,sce1$Group%in%c('Post-treatment')]
sce$celltype<-str_replace(sce$celltype,'Tumor','Cancer')
ICOSLG<-AverageExpression(object = sce,features = 'ICOSLG',group.by = 'celltype')
ICOSLG<-unlist(ICOSLG)
ICOSLG<-as.data.frame(ICOSLG)
ICOSLG<-t(ICOSLG)
ICOSLG<-as.data.frame(ICOSLG)
colnames(ICOSLG)<-'Post-treatment'
ICOSLG_Post_treatment<-ICOSLG
ICOSLG<-inner_join(ICOSLG_GSE120926,ICOSLG_GSE162025)

ICOSLG<-cbind(ICOSLG_Pre_treatment,ICOSLG_Post_treatment)
k<-rownames(ICOSLG)%in%rownames(ICOSLG_GSE150430)
table(k)
ICOSLG<-ICOSLG[k,]
k<-rownames(ICOSLG_GSE150430)%in%rownames(ICOSLG)
table(k)
#ICOSLG<-ICOSLG[k,]
ICOSLG<-cbind(ICOSLG_GSE150430,ICOSLG)
ICOSLG<-cbind(ICOSLG_GSE150825,ICOSLG)
#ICOSLG<-cbind(ICOSLG_Pre_treatment,ICOSLG_Post_treatment)
k<-rownames(ICOSLG)%in%rownames(ICOSLG_GSE162025)
table(k)
ICOSLG<-ICOSLG[k,]
k<-rownames(ICOSLG_GSE162025)%in%rownames(ICOSLG)
table(k)
ICOSLG_GSE162025<-ICOSLG_GSE162025[k,]
ICOSLG<-cbind(ICOSLG_GSE162025,ICOSLG)

k<-rownames(ICOSLG)%in%rownames(ICOSLG_GSE120926)
table(k)
ICOSLG<-ICOSLG[k,]
k<-rownames(ICOSLG_GSE120926)%in%rownames(ICOSLG)
table(k)
ICOSLG_GSE120926<-ICOSLG_GSE120926[k,]
ICOSLG<-merge(ICOSLG_GSE120926,ICOSLG)


ICOSLG<-merge(ICOSLG_Pre_treatment,ICOSLG_Post_treatment)

ICOSLG_Pre_treatment$celltype<-rownames(ICOSLG_Pre_treatment)
ICOSLG_Post_treatment$celltype<-rownames(ICOSLG_Post_treatment)
ICOSLG_GSE150430$celltype<-rownames(ICOSLG_GSE150430)
ICOSLG_GSE120926$celltype<-rownames(ICOSLG_GSE120926)
ICOSLG_GSE150825$celltype<-rownames(ICOSLG_GSE150825)
ICOSLG_GSE162025$celltype<-rownames(ICOSLG_GSE162025)


ICOSLG<-inner_join(ICOSLG_Pre_treatment,ICOSLG_Post_treatment)
ICOSLG<-full_join(ICOSLG,ICOSLG_GSE150430)
ICOSLG<-full_join(ICOSLG,ICOSLG_GSE150825)
ICOSLG<-full_join(ICOSLG,ICOSLG_GSE162025)
ICOSLG<-full_join(ICOSLG,ICOSLG_GSE120926)
rownames(ICOSLG)<-ICOSLG$celltype
ICOSLG<-ICOSLG[,-2]
