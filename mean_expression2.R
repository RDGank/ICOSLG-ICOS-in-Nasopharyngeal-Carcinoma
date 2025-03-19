load(file = 'sce_GSE150430_Tcell.RData')
ICOS<-AverageExpression(object = sce,features = 'ICOS',group.by = 'celltype')
ICOS<-unlist(ICOS)
ICOS<-as.data.frame(ICOS)
ICOS<-t(ICOS)
ICOS<-as.data.frame(ICOS)
colnames(ICOS)<-'GSE150430'
ICOS_GSE150430<-ICOS

load(file = 'sce_GSE150825_Tcell.RData')
sce$celltype<-str_replace(sce$celltype,'Tumor','Cancer')
ICOS<-AverageExpression(object = sce,features = 'ICOS',group.by = 'celltype')
ICOS<-unlist(ICOS)
ICOS<-as.data.frame(ICOS)
ICOS<-t(ICOS)
ICOS<-as.data.frame(ICOS)
colnames(ICOS)<-'GSE150825'
ICOS_GSE150825<-ICOS

load(file = 'sce_GSE162025_Tcell.RData')
sce$celltype<-str_replace(sce$celltype,'Tumor','Cancer')
ICOS<-AverageExpression(object = sce,features = 'ICOS',group.by = 'celltype')
ICOS<-unlist(ICOS)
ICOS<-as.data.frame(ICOS)
ICOS<-t(ICOS)
ICOS<-as.data.frame(ICOS)
colnames(ICOS)<-'GSE162025'
ICOS_GSE162025<-ICOS


load(file = 'sce_GSE120926_Tcell.RData')
sce$celltype<-str_replace(sce$celltype,'Tumor','Cancer')
ICOS<-AverageExpression(object = sce,features = 'ICOS',group.by = 'celltype')
ICOS<-unlist(ICOS)
ICOS<-as.data.frame(ICOS)
ICOS<-t(ICOS)
ICOS<-as.data.frame(ICOS)
colnames(ICOS)<-'GSE120926'
ICOS_GSE120926<-ICOS

load(file = 'E:/NM30/scp/sce_prospective_T.RData')
sce1<-sce
sce<-sce1[,sce1$Group%in%c('Pre-treatment')]
sce$celltype<-str_replace(sce$celltype,'Tumor','Cancer')
ICOS<-AverageExpression(object = sce,features = 'ICOS',group.by = 'celltype')
ICOS<-unlist(ICOS)
ICOS<-as.data.frame(ICOS)
ICOS<-t(ICOS)
ICOS<-as.data.frame(ICOS)
colnames(ICOS)<-'Pre-treatment'
ICOS_Pre_treatment<-ICOS

sce<-sce1[,sce1$Group%in%c('Post-treatment')]
sce$celltype<-str_replace(sce$celltype,'Tumor','Cancer')
ICOS<-AverageExpression(object = sce,features = 'ICOS',group.by = 'celltype')
ICOS<-unlist(ICOS)
ICOS<-as.data.frame(ICOS)
ICOS<-t(ICOS)
ICOS<-as.data.frame(ICOS)
colnames(ICOS)<-'Post-treatment'
ICOS_Post_treatment<-ICOS

ICOS_Pre_treatment$celltype<-rownames(ICOS_Pre_treatment)
ICOS_Post_treatment$celltype<-rownames(ICOS_Post_treatment)
ICOS_GSE150430$celltype<-rownames(ICOS_GSE150430)
ICOS_GSE120926$celltype<-rownames(ICOS_GSE120926)
ICOS_GSE150825$celltype<-rownames(ICOS_GSE150825)
ICOS_GSE162025$celltype<-rownames(ICOS_GSE162025)


ICOS<-inner_join(ICOS_Pre_treatment,ICOS_Post_treatment)
ICOS<-full_join(ICOS,ICOS_GSE150430)
ICOS<-full_join(ICOS,ICOS_GSE150825)
ICOS<-full_join(ICOS,ICOS_GSE162025)
ICOS<-full_join(ICOS,ICOS_GSE120926)
rownames(ICOS)<-ICOS$celltype
ICOS<-ICOS[,-2]
