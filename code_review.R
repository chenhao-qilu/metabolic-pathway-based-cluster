####TCGA EXPRESSION####
setwd("E:/R/ESO3/DEG")
install.packages("rio")
library(rio)
counts2 = import(file = 'data source.xlsx', sheet = 1) 
counts<-counts2[!duplicated(counts2$Hugo_Symbol),]
counts<-counts[-1,]
rownames(counts)<-counts[,1]
counts<-counts[,-c(1:2)]
counts<-log2(counts+1)
counts<-2^counts-1
counts<-round(counts)
library(DMwR2)
counts <- knnImputation(counts,k=10)
write.table(counts,"countsknnESCC.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
counts1 <- read.table("countsknn.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
colnames(counts) <- substr(colnames(counts),1,12)

####GEO EXPRESSION####
setwd("E:/R/ESO3/GEO")
counts3 = import(file = 'data source.xlsx',sheet = 3) 
counts3<-counts3%>%column_to_rownames("Gene Symbol")
surv3<- import(file = 'data source.xlsx', sheet = 4)
surv3<-surv3[which(surv3$Sample_Type=="Tumor"),]
comsample3<-intersect(surv3$Accession,colnames(counts3))
counts3<-counts3[,comsample3]
counts4 = import(file = 'data source.xlsx', sheet = 5)
counts4<-counts4%>%column_to_rownames("Gene Symbol")
write.table(counts3,"countsgeo.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

####DEG####
setwd("E:/R/ESO3/DEG")
library(limma)
group$group<-paste0("MP",group$group)
group<-group%>%column_to_rownames("sample")
identical(rownames(group),colnames(counts2))
group<-as.character(group$group)
design<-model.matrix(~0+factor(group))
colnames(design)=levels(factor(group))
rownames(design)=colnames(counts2)
contrast.matrix<-makeContrasts(MPS1 - MPS2,levels=design)
##step1
fit<-lmFit(counts2,design)
##step2
fit2<-contrasts.fit(fit,contrast.matrix)
fit2<-eBayes(fit2)
##step3
tempOutput=topTable(fit2,coef=1,n=Inf,adjust='fdr')
nrDEG3=na.omit(tempOutput)
nrDEG3<-nrDEG3[which(nrDEG3$adj.P.Val <= 0.001),]
union<-union(rownames(nrDEG1),rownames(nrDEG2))
uninter<-intersect(union,rownames(nrDEG3))
write.table(uninter,"uni.txt",sep = "\t")
write.table(uninter,"uninter.txt",sep = "\t",row.names = T,col.names = NA,quote = F)


####GO####
setwd("E:/R/ESO3/INTERACT")
library(tidyverse)
library("BiocManager")
library(org.Hs.eg.db)
library(clusterProfiler)
DEG<-read.table("cytoscape.txt",sep = "\t",header = F)


genelist <- bitr(DEG$V2, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
DEG <- inner_join(DEG,genelist,by=c("V2"="SYMBOL"))

ego <- enrichGO(gene = DEG$ENTREZID,
                OrgDb = org.Hs.eg.db, 
                ont = "all",
                pAdjustMethod = "BH",
                minGSSize = 1,
                pvalueCutoff =0.05, 
                qvalueCutoff =0.05,
                readable = TRUE)

ego_res <- ego@result
save(ego,ego_res,file = "GO_SPP1_DEG.Rdata")

barplot(ego, showCategory = 20,color = "pvalue")

dotplot(ego, showCategory = 20)

barplot(ego, drop = TRUE, showCategory =10,split="ONTOLOGY") + 
  facet_grid(ONTOLOGY~., scale='free')
dotplot(ego,showCategory = 10,split="ONTOLOGY") + 
  facet_grid(ONTOLOGY~., scale='free')

####GSVA####
setwd("E:/R/ESO3/GSVA")
library(GSVA)
library(limma)
library(GSEABase)
library(org.Hs.eg.db)
library(clusterProfiler)
hallset=getGmt('h.all.v7.5.1.symbols.gmt')
reactset=getGmt('Reactome_Metabolism.gmt')
keggset=getGmt('c2.cp.kegg.v7.5.1.symbols.gmt')
counts<-counts3
keggEs1=gsva(as.matrix(counts2),reactset,abs.ranking=TRUE,method="gsva")
keggEs2=gsva(as.matrix(counts3),reactset,abs.ranking=TRUE,method="gsva")
keggEs3=gsva(as.matrix(counts4),reactset,abs.ranking=TRUE,method="gsva")

write.table(keggEs1,"reactomegeossgsea.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

keggEs4=gsva(as.matrix(counts2),keggset,abs.ranking=TRUE,method="gsva")
keggEs5=gsva(as.matrix(counts3),keggset,abs.ranking=TRUE,method="gsva")
keggEs7=gsva(as.matrix(counts3),keggset,abs.ranking=TRUE,method="gsva")
keggEs8=gsva(as.matrix(counts4),keggset,abs.ranking=TRUE,method="gsva")
keggEs9<-gsva(as.matrix(counts2),hallset,abs.ranking=TRUE,method="gsva")
keggEs10<-gsva(as.matrix(counts3),hallset,abs.ranking=TRUE,method="gsva")


####GSEA####
setwd("E:/R/ESO3/GSEA")
library(limma)
group$group<-paste0("MPC",group$group)
group<-group%>%column_to_rownames("sample")
group<-as.character(group$group)
design<-model.matrix(~0+factor(group))
colnames(design)=levels(factor(group))
rownames(design)=colnames(counts3)
contrast.matrix<-makeContrasts(S3 - S1,levels=design)

##step1
fit<-lmFit(counts3,design)
##step2
fit2<-contrasts.fit(fit,contrast.matrix)
fit2<-eBayes(fit2)
##step3
tempOutput=topTable(fit2,coef=1,n=Inf,adjust='fdr')
nrDEG=na.omit(tempOutput)
nrDEG<-nrDEG[which(nrDEG$adj.P.Val<=0.001),]
write.table(nrDEG,"S1S2",sep = "\t",row.names = T,col.names = NA,quote = F)
library(tidyverse)
library("BiocManager")
library(org.Hs.eg.db)
library(clusterProfiler)
kegmt <- read.gmt("c2.cp.kegg.v7.5.1.symbols.gmt")
geneList = nrDEG[,1]
names(geneList) = as.character(rownames(nrDEG))
head(geneList)

names(geneList) = as.character(rownames(tempOutput))


geneList = sort(geneList, decreasing = T)
head(geneList)

KEGG<-GSEA(geneList,TERM2GENE = kegmt,by="fgsea") 
KEGG_result_df <- as.data.frame(KEGG)


write.table(KEGG_result_df,file="GSEA_MSigDb_keggk3_result.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
save(KEGG,KEGG_result_df,file = "GSEA_deg_SPP1.rda")

library(enrichplot)
gseaplot2(KEGG,22,color="")
gseaplot2(KEGG,3,rel_heights=c(1, .2, .6),pvalue_table = T)
gseaplot2(KEGG, geneSetID = c(1,2,4,6,8), subplots = (1:3))
gseaplot2(KEGG, geneSetID = c(2,6,8,10,13), subplots = (1:2),base_size = 18)
gseaplot2(KEGG,1,color="red",pvalue_table = T, title="", base_size = 10,ES_geom = "line")

####NMF####
setwd("E:/R/ESO3/NMF")
library(NMF)
library(doMPI)
res1 <- nmf(keggEs1,2:7,"lee",nrun=10)
plot(res1)
V.random1 <- randomize(keggEs1)
res.random1 <- nmf(V.random1, 2:7,"lee",nrun=10)
plot(res1, res.random1)
res4<- nmf(keggEs1,3,"lee",nrun=10)
coefmap(res4)
consensusmap(res4)

group <- predict(res4)
group<- as.data.frame(group)
rownames(group)<-colnames(keggEs1)
group$group <- factor(group$group,levels=c(1,2,3))
save(group,file = 'grouptcgaleek3.rda')
write.csv(group,'group.csv')


####HEATMAP####
setwd("E:/R/ESO3/PHEATMAP")
library(tidyverse)
library(pheatmap)
library(rio)

reactomees<-keggEs1
reactomees<-reactomees%>%t()%>%as.data.frame()
Pvaluekw<-c(rep(0,ncol(reactomees)))
design<-data.frame(rownames(reactomees),group$group)
for(i in 1:ncol(reactomees)){
  ab<-as.numeric(reactomees[1:nrow(reactomees),i])
  b<-design$group.group
  aa<-data.frame(ab,b)
  y1=kruskal.test(ab~b,data=aa)
  Pvaluekw[i]<-y1$p.value
}
reactomees<-t(reactomees)
reactomees<-as.data.frame(reactomees)
reactomees$P<-Pvaluekw
reactomees<-reactomees[which(reactomees$P < 0.05),]
reactomees<-reactomees[,-ncol(reactomees)]

group <- group %>% 
  rownames_to_column("sample")
annotation <- group %>% arrange(group) %>% column_to_rownames("sample")
a <- group %>% arrange(group) %>% mutate(sample=substring(.$sample,1,12))
b <- t(reactomees) %>% 
  as.data.frame() %>% 
  rownames_to_column("sample") %>% 
  mutate(sample=substring(.$sample,1,12))
c <- inner_join(a,b,"sample") %>% .[,-2] %>% column_to_rownames("sample") %>% t(.)
bk <- c(seq(-2,-0.01,by=0.01),seq(0,2,by=0.01))

pheatmap(c,annotation = annotation,
         cluster_cols = F,fontsize=5,fontsize_row=9,
         scale="row",show_colnames=F,
         fontsize_col=5,
         color = colorRampPalette(c("navy", "white", "red"))(100),
         legend_breaks=seq(-2.5,2.5,1),
         breaks=bk)
ann_colors = list(
  group= c( S1="#2a9d8c",S2= "#e9c46b",S3="#e66f51"))
pheatmap(c,
         annotation = annotation,
         scale = "row",
         cutree_rows = 4,
         show_colnames =F,
         color = colorRampPalette(c("navy", "white", "red"))(400),
         annotation_colors = ann_colors[1],
         cluster_cols =F,
         fontsize = 10,
         fontsize_row=7,
         fontsize_col=7,
         breaks=bk)
dev.off()
dev.off()


####SURVIVAL ANALYSIS####
setwd("E:/R/ESO3/SURVIVAL")
surv <- import("original data_TCGA.xlsx",sheet = 2)
surv3 <- import(file = 'original data_GSE53625.xlsx', sheet = 2)
surv<-surv[which(surv$SUBTYPE=="ESCA_ESCC"),]
rownames(surv) <- gsub("-",".",rownames(surv))

colnames(counts) <- substr(colnames(counts),1,12)
comsample<-intersect(rownames(surv),colnames(counts))
surv3<-surv3[comsample1,]
surv<-surv3
surv3$OS_Time<-surv3$OS_Time/12*365
surv$group<-group$group
surv$group <- factor(surv$group)


surv.expr$OS_STATUS<-substr(surv.expr$OS_STATUS,1,1)
surv.expr$OS_STATUS<-as.numeric(surv.expr$OS_STATUS)
write.table(surv,"surv.txt",sep = "\t",row.names = T,col.names = NA,quote = F)


library(survival)
library(survminer)
library(ggplot2)
surv3$group<-group$group
fitd <- survdiff(Surv(OS_MONTHS, OS_STATUS=="1") ~ group,
                 data      = surv,
                 na.action = na.exclude)
fitd <- survdiff(Surv(OS_Time, OS_Status=="1") ~ group,
                 data      = surv3,
                 na.action = na.exclude)
pValue <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit <- survfit(Surv(OS_MONTHS, OS_STATUS=="1")~ group, data =surv)
fit2 <- survfit(Surv(OS_Time, OS_Status=="1")~ group, data =surv3)
fit3<-survfit(Surv(PFS_MONTHS, PFS_STATUS)~ group, data =surv)
p.lab <- paste0("P", ifelse(pValue < 0.001, " < 0.001", paste0(" = ",round(pValue, 3))))
ggsurvplot(fit2,
           data = surv3,
           pval = p.lab,
           conf.int = TRUE,
           risk.table = TRUE, 
           risk.table.col = "strata",
           palette = "jco", 
           legend.labs = c("MPC1", "MPC2","MPC3"),
           size = 1,
           xlim = c(0,80),
           break.time.by = 20, 
           legend.title = "",
           surv.median.line = "hv", 
           ylab = "Survival probability (%)", 
           xlab = "Time (Months)",
           ncensor.plot = TRUE, 
           ncensor.plot.height = 0.25,
           risk.table.y.text = FALSE)
dev.off()



####CLINICAL ANNOTATION####
setwd("E:/R/ESO3/COMPLEXPHEATMAP")
library(ComplexHeatmap)
library(rio)
library(circlize)
identical(group$sample,rownames(surv))
surv$group<-group$group
surv<-surv%>%rownames_to_column("sample")
surv3$group<-group$group
surv3<-surv3%>%rownames_to_column("sample")
export(surv3,"geo.csv")
export(surv,"tcga.csv")
rank<-import("geo.csv")
ha1 = HeatmapAnnotation(df = rank[,c(8,17,23)],
                        gap = unit(2, "mm"),
                        show_legend = T,
                        simple_anno_size = unit(0.7, "cm"))

abc<-t(keggEs1)
Heatmap(keggEs2, top_annotation = ha1)
draw(ha1)
dev.off()

####COX####
setwd("E:/R/ESO3/COX")
library(survival)
library(forestplot)
library(tidyverse)
comcox<-intersect(rownames(surv3),colnames(counts3))
surv3<-surv3[comcox,]

surv$sample<- rownames(surv)
surv$sample <- gsub("-",".",surv$sample)
rownames(surv) <- surv$sample
surv <- surv[,-1]
surv <- surv[,-2]
fix(surv)
surv.expr<-surv
surv.expr<-surv3[,c(6,9,11,15,16,19:20,22,23)]
surv$AJCC_PATHOLOGIC_TUMOR_STAGE<-gsub("IIC","II",surv$AJCC_PATHOLOGIC_TUMOR_STAGE)
library(forestmodel)
library(survival) 
library(dplyr)
library(forestplot)
surv[surv == ""]<-NA
surv.expr[surv.expr == "yes"]<- 1
surv<-surv[complete.cases(surv$AJCC_PATHOLOGIC_TUMOR_STAGE), ] 
surv$group<-group$group
surv.expr<-surv.expr[,c(4:7,19:21,25,29:30,37)]
surv.expr[surv.expr == "no"]<- 0
surv.expr<-read.table(file = 'geo.txt', sep = '\t', header = TRUE,row.names = 1)
# install.packages("ezcox")
library(ezcox)
surv.expr$group<-group$group
surv.expr$group<-as.factor(surv.expr$group)
for(i in 9) {surv.expr[, i] = as.numeric(surv.expr[, i])}

colnames(surv.expr)<- c("age", "t","Sex","time","tobacco_use","tnm_stage","n_stage","alcohol_use","death","group")
cox <- ezcox(surv,time ="OS_MONTHS", status = "OS_STATUS",
             covariates =c("AGE","SEX","AJCC_PATHOLOGIC_TUMOR_STAGE","PATH_M_STAGE","PATH_N_STAGE","PATH_T_STAGE","RACE","group"))
cox <- ezcox(surv.expr,time ="time", status = "death",
             covariates =c("age","Sex","t","tobacco_use","alcohol_use","tnm_stage","n_stage","group"))

cox1<-cox[,-c(1:2)]
fix(cox1)
surv.expr$group<-paste0("S", group$group)
ins <- function(x) {c(x, rep(NA, ncol(cox1)-1))}
result<-rbind(c("a","b","c","d","e","f","g","h","i","j"),
              cox1[1, ],
              ins("SEX"),
              ins("Female"),
              cox1[2, ],
              ins("STAGE"),
              ins("I"),
              cox1[3:5, ],
              ins("M_STAGE"),
              ins("M0"),
              cox1[6,],
              ins("N_STAGE"),
              ins("N0"),
              cox1[7:9,],
              ins("T_stage"),
              ins("T1"),
              cox1[10:12,],
              ins("RACE"),
              ins("Asian"),
              cox1[13:14,],
              ins("group"),
              ins("S1"),
              cox1[15:16,],
              
              c(NA, NA, NA, NA, NA,NA,NA,NA,NA,NA)
)
result<-rbind(c("a","b","c","d","e","f","g","h","i","j"),
              cox1[1, ],
              ins("SEX"),
              ins("Female"),
              cox1[2, ],
              ins("tobacco_use"),
              ins("no"),
              cox1[6, ],
              ins("alcohol_use"),
              ins("no"),
              cox1[7,],
              ins("tnm_stage"),
              ins("I"),
              cox1[8:9,],
              ins("t_stage"),
              ins("T1"),
              cox1[3:5,],
              ins("n_stage"),
              ins("N0"),
              cox1[10:12,],
              ins("group"),
              ins("S1"),
              cox1[13:14,],
              
              c(NA, NA, NA, NA, NA,NA,NA,NA,NA,NA)
)
fix(result)

result[1,]<-colnames(result)
result<-result[,-c(2,4,5)]
forestplot(result,
           mean=result[,3],
           lower=result[,4],
           upper=result[,5],
           zero=1,
           boxsize=0.6,
           graph.pos= 4 ,
           hrzl_lines=list("1" = gpar(lty=1,lwd=2),
                           "2" = gpar(lty=2),
                           "30"= gpar(lwd=2,lty=1,columns=c(1:3,5:8)) ),
           graphwidth = unit(.25,"npc"),
           xlab="Hazard ratios",
           xticks=c(0.4,1,3,5,7,9) ,
           is.summary=c(T,T,T,F,F,T,F,F,T,F,F,T,F,F,F,T,F,F,F,F,T,F,F,F,F,T,F,F,F),
           txt_gp=fpTxtGp(
             label=gpar(cex=2.5),
             ticks=gpar(cex=2.5),
             xlab=gpar(cex=2.5),
             title=gpar(cex=2.5)),
           lwd.zero=1,
           lwd.ci=1.5,
           lwd.xaxis=2,
           lty.ci=1.5,
           ci.vertices =T,
           ci.vertices.height=0.2,
           clip=c(0.1,8),
           ineheight=unit(10, 'mm'),
           line.margin=unit(8, 'mm'),
           lineheight = unit(11.5,'mm'),
           colgap=unit(0.5, 'mm'),
           fn.ci_norm="fpDrawDiamondCI",
           col=fpColors(box="#00A896", lines="#02C39A", zero = "black"),)
dev.off()
write.table(result,"georesult.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
surv.expr$AJCC_PATHOLOGIC_TUMOR_STAGE<-gsub("IA","I",surv.expr$AJCC_PATHOLOGIC_TUMOR_STAGE)
surv.expr$group<-paste0("S", surv.expr$group)
surv$OS_STATUS<-substr(surv$OS_STATUS,1,1)
surv$OS_STATUS<-as.numeric(surv$OS_STATUS)
surv.expr<-surv
surv.expr<-subset(surv,surv$OS_MONTHS>=1)
surv.expr$death<-gsub("yes","1",surv.expr$death)
surv.expr$death<-as.numeric(surv.expr$death)
surv.expr$`death at fu`<-as.numeric(surv.expr$`death at fu`)
surv.expr$group<-paste0("S", group$group)
surv$AJCC_PATHOLOGIC_TUMOR_STAGE<-gsub("IVA","IV",surv$AJCC_PATHOLOGIC_TUMOR_STAGE)
coxphmodel1 <- coxph(Surv(OS_MONTHS, OS_STATUS) ~ AGE+SEX+STAGE+group, surv)
coxphmodel1 <- coxph(Surv(survival.time.months., death.at.fu) ~ age+Sex+tnm_stage+group, surv.expr)
forest_model(coxphmodel1,
             format_options = forest_model_format_options(colour = "#00A896",      
                                                          shape = 15,          
                                                          text_size = 5,      
                                                          point_size = 5,        
                                                          banded = TRUE),
             factor_separate_line = TRUE) 



####CIBERSORT####
setwd("E:/R/ESO3/CIBERSORT")   

library(e1071)
library(parallel)
library(preprocessCore)
library(dplyr)
library()
source("CIBERSORT.R")   
sig_matrix <- "LM22.txt"   
mixture_file = 'GSE121931_expression.txt'   
res_cibersort <- CIBERSORT(sig_matrix, mixture_file, perm=100, QN=TRUE)
save(res_cibersort,file = "res_cibersortgeo2.Rdata")   
group$group<-paste0("S", group$group)
load("res_cibersort.Rdata")
res_cibersort <- res_cibersort[,1:22]   
ciber.res <- res_cibersort[,colSums(res_cibersort) > 0] 

ciber.res<-ciber.res[order(group$group),]
ciber.res<-as.data.frame(ciber.res)
ciber.res<-arrange(ciber.res,group$group,`B cells memory`)


library(RColorBrewer)
mycol <- ggplot2::alpha(rainbow(ncol(ciber.res)), 0.7) 

par(bty="o", mgp = c(2.5,0.3,0), mar = c(2.1,4.1,2.1,10.1),tcl=-.25,las = 1,xpd = F)
barplot(as.matrix(t(ciber.res)),
        border = NA, 
        names.arg = rep("",nrow(ciber.res)), 
        yaxt = "n", 
        ylab = "Relative percentage",
        col = mycol)

axis(side = 2, at = c(0,0.2,0.4,0.6,0.8,1), 
     labels = c("0%","20%","40%","60%","80%","100%"))
legend(par("usr")[2]-5, 
       par("usr")[4], 
       legend = colnames(ciber.res), 
       xpd = T,
       fill = mycol ,
       cex = 0.6, 
       border = NA, 
       y.intersp = 1,
       x.intersp = 0.2,
       bty = "n")
dev.off()  

library(tidyverse)

res_cibersort1<-res_cibersort[,c(9,21)]
a<-res_cibersort1
a<-as.data.frame(a)

identical(rownames(a),rownames(group))
b <- group
class(b$group)
a$group <- b$group
a <- a %>% rownames_to_column("sample")
library(ggsci)
library(tidyr)
library(ggpubr)

b <- gather(a,key=CIBERSORT,value = Proportion,-c(group,sample))

ggboxplot(b, x = "CIBERSORT", y = "Proportion",
          fill = "group", palette = "lancet")+
  stat_compare_means(aes(group = group),
                     method = "wilcox.test",
                     label = "p.signif",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))+
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle=45, hjust=1)) 
dev.off()


####MUTATION####
setwd("E:/R/ESO3/TMB")
library(TCGAbiolinks)
library(maftools)
library(tidyverse)
library(readxl)
library(readr)

colnames(surv)[37]<-"Tumor_Sample_Barcode"
mut2 <- read.maf("TCGA.ESCA.varscan.9dab6855-ba5e-4afb-b3a0-f034a6f82eb2.DR-10.0.somatic.maf",isTCGA = T,clinicalData=surv )
sample <- subset(surv, SUBTYPE=="ESCA_ESCC")$Tumor_Sample_Barcode
mut1<- subsetMaf(maf=mut2, tsb=sample, isTCGA=TRUE)
surv$group<-group$group
surv$sample <- rownames(surv)
sample1<-subset(surv,group=="S1")$Tumor_Sample_Barcode
sample2<-subset(surv,group=="S2")$Tumor_Sample_Barcode
sample3<-subset(surv,group=="S3")$Tumor_Sample_Barcode
mut3<- subsetMaf(maf=mut1, tsb=sample1, isTCGA=TRUE)
mut4<-subsetMaf(maf=mut1, tsb=sample2, isTCGA=TRUE)
mut5<-subsetMaf(maf=mut1, tsb=sample3, isTCGA=TRUE)
pt.vs.rt3 <- mafCompare(m1 = mut4, m2 = mut5, m1Name = '2', m2Name = '3',minMut = 5)
write.table(pt.vs.rt2[[1]],"1v3.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
maf.gene<-mut1[[4]]
mut1@clinical.data$group<-group$group
oncoplot(maf =mut1, fontSize = 0.6 ,
         showTumorSampleBarcodes = F ,clinicalFeatures = 'group',
         sortByAnnotation = TRUE,)

a <- mut1@data %>% 
  .[,c("Hugo_Symbol","Variant_Classification","Tumor_Sample_Barcode")] %>% 
  as.data.frame() %>% 
  mutate(Tumor_Sample_Barcode = substring(.$Tumor_Sample_Barcode,1,12))
gene <- as.character(unique(a$Hugo_Symbol))
sample <- as.character(unique(a$Tumor_Sample_Barcode))

mat <- as.data.frame(matrix("",length(gene),length(sample),
                            dimnames = list(gene,sample)))
mat_0_1 <- as.data.frame(matrix(0,length(gene),length(sample),
                                dimnames = list(gene,sample)))




for (i in 1:nrow(a)){
  mat[as.character(a[i,1]),as.character(a[i,3])] <- as.character(a[i,2])
}
for (i in 1:nrow(a)){
  mat_0_1[as.character(a[i,1]),as.character(a[i,3])] <- 1
}

mat_0_1<-mat_0_1%>%t()%>%as.data.frame()
rownames(mat_0_1) <- gsub("-",".",rownames(mat_0_1))
mat_0_1<- mat_0_1 %>% rownames_to_column('sample')
group<- group%>%rownames_to_column('sample')
mat_0_12<-inner_join(mat_0_1,group,by=c("sample"="sample"))
mat_0_S3<-mat_0_12[which(mat_0_12$group=="S3"),]
mat_0_1<-mat_0_1%>%column_to_rownames("sample")
mat_0_1<-mat_0_1[,-5694]
gene<-as.character(rownames(reactomees))
my_comparisons <- list( c("S1", "S2"), c("S2", "S3"), c("S1", "S3") )
mat_0_1<-mat_0_1[-1,]
gene_countall<-inner_join(gene_countall,gene_count3,by=c("gene"="gene"))
gene_countall<-gene_countall%>%column_to_rownames("gene")
gene_countall<-gene_countall[gene,]
write.table(gene_countall,"gene_countall.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
Pvaluef<-c(rep(0,ncol(mat_0_1)))
design<-data.frame(rownames(mat_0_1),mat_0_12$group)
for(i in 1:ncol(mat_0_1)){
  ab<-as.numeric(mat_0_1[1:nrow(mat_0_1),i])
  b<-design$group.group
  aa<-data.frame(ab,b)
  y1=fisher.test(ab,b)
  Pvaluef[i]<-y1$p.value
}
mat_0_1<-mat_0_1%>%t()%>%as.data.frame()

mat_0_1$P<-Pvaluef
mat_0_11<-mat_0_1[which(mat_0_1$P<=0.05),]
pvalue<-as.data.frame(mat_0_11$P)
rownames(pvalue)<-rownames(mat_0_11)
write.table(pvalue,"oncology.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
gene_count <- data.frame(gene=rownames(mat_0_1),
                         count=as.numeric(apply(mat_0_1,1,sum))) %>%
  arrange(desc(count))



rownames(mat_0_S3)<-mat_0_S3$sample
mat_0_S3<-mat_0_S3[,-c(1,5695)]
gene_count3 <- data.frame(gene=rownames(mat_0_S3),
                          count=as.numeric(apply(mat_0_S3,1,sum)))%>%
  arrange(desc(count))
mat_0_S3<-mat_0_S3%>%t()%>%as.data.frame()
gene_count<-full_join(gene_count,gene_count3, by= "gene")
gene_count<-gene_count[rownames(pvalue),]


gene_count <- data.frame(gene=rownames(mat_0_1),
                         count=as.numeric(apply(mat_0_1,1,sum))) %>%
  arrange(desc(count))
gene_count<- gene_count%>%column_to_rownames("gene")
rownames(gene_count) <- gsub("-",".",rownames(gene_count))
gene_count<- gene_count%>%rownames_to_column("sample")
gene_count<-inner_join(gene_count,group,by=c("sample"="sample"))
aov<-aov(count~group,data=gene_count)
summary(aov)


gene_count[is.na(gene_count)]<- 0
gene_count$count.x1<- 26-gene_count$count.x
gene_count$count.y1<- 34-gene_count$count.y
gene_count$count2<- 30-gene_count$count
gene_count$col<-90
gene_count0<-gene_count[c(1:8),]
pval=c()
for (i in 1:nrow(gene_count)){
  compare<-chisq.test(gene_count[i,])$p.value
  pval=c(pval,compare)
}
gene_count$p<-pval
colnames(gene_count)<-c("S1","S2","S3")
compare<-chisq.test(  gene_count[8,])
gene_top <- gene_count$gene[1:20] # 修改数字，代表TOP多少
save(mat,mat_0_1,file = "TMBESCC.rda")


reactomees<-mat_0_1
reactomees<-as.data.frame(reactomees)
Pvaluekw<-c(rep(0,ncol(reactomees)))
design<-data.frame(rownames(reactomees),mat_0_12$group)
for(i in 1:ncol(reactomees)){
  ab<-as.numeric(reactomees[1:nrow(reactomees),i])
  b<-design$mat_0_12.group
  aa<-data.frame(ab,b)
  y1=fisher.test(ab,b)
  Pvaluekw[i]<-y1$p.value
}
reactomees<-t(reactomees)
reactomees<-as.data.frame(reactomees)
reactomees$P<-Pvaluekw
reactomees<-reactomees[which(reactomees$P < 0.05),]
mat_0_1<-reactomees%>%t()%>%as.data.frame()


oncoplot(maf = mut1,
         top = 30,   
         fontSize = 0.6,   
         showTumorSampleBarcodes = F,
         sortByAnnotation = T,
         clinicalFeatures = "group")   
dev.off()


####CLINICAL TABLE####
setwd("E:/R/ESO3/CLINICAL")
# remove.packages("gtsummary")
# install.packages("rvcheck") 

library(gt)
library(webshot)
webshot::install_phantomjs()
rvcheck::update_all(check_R=FALSE,which=c("CRAN","BioC","github"))
library(gtsummary)
library(tidyverse)
surv$AJCC_PATHOLOGIC_TUMOR_STAGE<-substr(surv$AJCC_PATHOLOGIC_TUMOR_STAGE,1,9)
surv$AJCC_PATHOLOGIC_TUMOR_STAGE[which(surv$AJCC_PATHOLOGIC_TUMOR_STAGE=="STAGE IVA")] <- "STAGE IV"
surv3$group<-group$group
clin <- surv %>% select(AGE,SEX,AJCC_PATHOLOGIC_TUMOR_STAGE,PATH_M_STAGE,PATH_N_STAGE,PATH_T_STAGE,RACE,group)
tbl_summary(
  clin,
  by = group, 
  missing = "no" 
) %>%
  add_n() %>% 
  modify_header(label = "**Variable**") %>% 
  bold_labels()%>%   
  as_gt() 
tbl_summary %>%
  as_flex_table() %>%
  flextable::save_as_docx(tb1)
library(flextable)
install.packages("flextable")
####XCELL####
setwd("E:/R/ESO3/XCELL")
library(xCell)
library(ggpubr)
library(tidyverse)


celltypeuse<-xCell.data$spill$K
rs<-xCellAnalysis(counts4,parallel.sz=10) 

rs <- as.data.frame(rs)
rs<-rs%>%t()%>%as.data.frame()
rs <- rs[,1:64]   
rs <- rs[,colSums(rs) > 0] 

rs<-rs[order(group$group),]
rs<-arrange(rs,group$group,`DC`)
library(RColorBrewer)
mycol <- ggplot2::alpha(rainbow(ncol(rs)), 0.7) 
par(bty="o", mgp = c(2.5,0.3,0), mar = c(2.1,4.1,2.1,10.1),tcl=-.25,las = 1,xpd = F)
barplot(as.matrix(t(rs)),
        border = NA, 
        names.arg = rep("",nrow(rs)),
        
        ylab = "value",
        col = mycol)

axis(side = 2, at = c(0,0.2,0.4,0.6,0.8,1), 
     labels = c("0%","20%","40%","60%","80%","100%"))
legend(par("usr")[2]-5, 
       par("usr")[4]+0.2, 
       legend = colnames(rs), 
       xpd = T,
       fill = mycol ,
       cex = 0.3, 
       border = NA, 
       y.intersp = 1,
       x.intersp = 0.1,
       bty = "n")
dev.off()   


reactomees<-rs
reactomees<-as.data.frame(reactomees)
Pvaluekw<-c(rep(0,ncol(reactomees)))
design<-data.frame(rownames(reactomees),group$group)
for(i in 1:ncol(reactomees)){
  ab<-as.numeric(reactomees[1:nrow(reactomees),i])
  b<-design$group.group
  aa<-data.frame(ab,b)
  y1=kruskal.test(ab~b,data=aa)
  Pvaluekw[i]<-y1$p.value
}
reactomees<-t(reactomees)
reactomees<-as.data.frame(reactomees)
reactomees$P<-Pvaluekw
reactomees<-reactomees[which(reactomees$P < 0.05),]
reactomees<-reactomees[,-ncol(reactomees)]



a<-t(reactomees)
a<-as.data.frame(a)

identical(rownames(a),rownames(group))
b <- group
class(b$group)
a$group <- b$group
a <- a %>% rownames_to_column("sample")


library(ggsci)
library(tidyr)
library(ggpubr)
b <- gather(a,key=xCell,value = Expression,-c(group,sample))

ggboxplot(b, x = "xCell", y = "Expression",
          fill = "group", palette = "lancet")+
  stat_compare_means(aes(group = group),
                     method = "wilcox.test",
                     label = "p.signif",
                     hide.ns=T,
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))+
  theme(text = element_text(size=10),
        axis.text.x = element_text(angle=45, hjust=1)) 

dev.off()

####IOBR####

setwd("E:/R/ESO3/IOBR")

library(IOBR)
library(EPIC)
library(estimate) 
library(tidyverse)
library(tidyHeatmap)
library(maftools)
library(ggpubr)
library(ggplot2)
sig_res<-calculate_sig_score(pdata           = NULL,
                             eset            = counts3,
                             signature       = signature_collection,
                             method          = "ssGSEA",
                             mini_gene_count = 2)
sig_res<-sig_res%>%column_to_rownames("ID")
sig_res<-sig_res[,-1]

reactomees<-sig_res
reactomees<-as.data.frame(reactomees)
Pvaluekw<-c(rep(0,ncol(reactomees)))
design<-data.frame(rownames(reactomees),group$group)
for(i in 1:ncol(reactomees)){
  ab<-as.numeric(reactomees[1:nrow(reactomees),i])
  b<-design$group.group
  aa<-data.frame(ab,b)
  y1=kruskal.test(ab~b,data=aa)
  Pvaluekw[i]<-y1$p.value
}
reactomees<-t(reactomees)
reactomees<-as.data.frame(reactomees)
reactomees$P<-Pvaluekw
reactomees<-reactomees[which(reactomees$P < 0.05),]
reactomees<-reactomees[-c(28:34,36:60),]
reactomees<-reactomees[-39,]
reactomees<-reactomees[,-ncol(reactomees)]
a<-t(reactomees)
a<-as.data.frame(a)
group<-group%>%column_to_rownames("sample")

identical(rownames(a),rownames(group))
b <- group
class(b$group)
a$group <- b$group
a <- a %>% rownames_to_column("sample")
library(ggsci)
library(tidyr)
library(ggpubr)

b <- gather(a,key=IOBR,value = Proportion,-c(group,sample))

ggboxplot(b, x = "IOBR", y = "Proportion",
          fill = "group", palette = "lancet")+
  stat_compare_means(aes(group = group),
                     method = "wilcox.test",
                     label = "p.signif",
                     hide.ns=T,
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))+
  theme(text = element_text(size=18),
        axis.text.x = element_text(angle=45, hjust=1)) 
dev.off()


####MUTATIONSIGNATURE####
setwd("E:/R/ESO3/MUTATIONSIGNATURE")
# BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
library(maftools)
library(tidyverse)
laml.titv = titv(maf = mut1, plot = FALSE, useSyn = TRUE)
plotTiTv(res = laml.titv)
library(BSgenome.Hsapiens.UCSC.hg38, quietly=TRUE)
luad.tnm <- trinucleotideMatrix(maf=mut1, ref_genome="BSgenome.Hsapiens.UCSC.hg38")
apobec_enrich <- plotApobecDiff(tnm=luad.tnm, maf=mut1)
maf.sign = estimateSignatures(mat = luad.tnm, nTry = 6)
plotCophenetic(res = maf.sign)
maf.sig = extractSignatures(mat = luad.tnm, n = 4)
plotSignatures(nmfRes = maf.sig, title_size = 1.3)
laml.og30.cosm = compareSignatures(nmfRes = maf.sig, sig_db = "legacy")
library('pheatmap')
pheatmap::pheatmap(mat=laml.og30.cosm$cosine_similarities, cluster_rows=FALSE, main="cosine similarity against validated signatures")
dev.off()
com.sig<-maf.sig[[2]]
com.sig<-t(com.sig)
rownames(com.sig)<-gsub("-",".",rownames(com.sig))
com.sig<-com.sig[order(group$group),]
com.sig<-as.data.frame(com.sig)
com.sig<-com.sig[order(group$group,com.sig$Signature_1),]


library(RColorBrewer)
mycol <- ggplot2::alpha(rainbow(ncol(com.sig)), 0.7)
par(bty="o", mgp = c(2.5,0.3,0), mar = c(2.1,4.1,2.1,10.1),tcl=-.25,las = 1,xpd = F)
barplot(as.matrix(t(com.sig)),
        border = NA, 
        names.arg = rep("",nrow(com.sig)), 
        yaxt = "n", 
        ylab = "Relative percentage",
        col = mycol)

axis(side = 2, at = c(0,0.2,0.4,0.6,0.8,1), 
     labels = c("0%","20%","40%","60%","80%","100%"))
legend(par("usr")[2]-5, 
       par("usr")[4], 
       legend = colnames(com.sig), 
       xpd = T,
       fill = mycol ,
       cex = 0.6, 
       border = NA, 
       y.intersp = 1,
       x.intersp = 0.2,
       bty = "n")
dev.off()   
reactomees<-com.sig
reactomees<-as.data.frame(reactomees)
Pvaluekw<-c(rep(0,ncol(reactomees)))
reactomees1<-reactomees%>%rownames_to_column("sample")
reactomees1<-inner_join(reactomees1,group,by=c("sample"="sample"))
design<-data.frame(rownames(reactomees),reactomees1$group)
for(i in 1:ncol(reactomees)){
  ab<-as.numeric(reactomees[1:nrow(reactomees),i])
  b<-design$reactomees1.group
  aa<-data.frame(ab,b)
  y1=aov(ab~b,data=aa)
  Pvaluekw[i]<-y1$p.value
}
y=aov(Signature_4~group,data=reactomees1)
summary(y)
reactomees<-t(reactomees)
reactomees<-as.data.frame(reactomees)
reactomees$P<-Pvaluekw
reactomees<-reactomees[which(reactomees$P < 0.05),]

####CCLE####
setwd("E:/R/ESO3/CCLE")
library(rio)        
library(data.table)
cclecount<-import("originnal data_Cell line.xlsx",sheet=1)
cclecount<-cclecount%>%column_to_rownames("CCLE_ID")
cclecount<-cclecount%>%t()%>%as.data.frame()
library(NMF)
library(doMPI)
library(tidyverse)
library(pheatmap)
library(dplyr)
library(ggplot2)
library(ggpubr)
keggES11<-gsva(as.matrix(cclecount),reactset,abs.ranking=TRUE,method="gsva")
ccle1 <- nmf(keggES11,2:7,"lee",nrun=10)
plot(ccle1)
V.ccle1 <- randomize(keggES11)
ccle.random1 <- nmf(V.ccle1, 2:7,"lee",nrun=10)
plot(ccle1, ccle.random1)
ccle2<-nmf(keggES11,3,"lee",nrun=10)
groupccle<- predict(ccle2)
groupccle<-as.data.frame(groupccle)
identical(group$group,groupccle$groupccle)
groupccle<- as.data.frame(groupccle)
groupccle$groupccle<-paste0("MPC",groupccle$groupccle)

reactomees<-ccle%>%t()%>%as.data.frame()
reactomees<-as.data.frame(reactomees)
Pvaluekw<-c(rep(0,ncol(reactomees)))
design<-data.frame(rownames(reactomees),groupccle$groupccle)
for(i in 1:ncol(reactomees)){
  ab<-as.numeric(reactomees[1:nrow(reactomees),i])
  b<-design$groupccle.groupccle
  aa<-data.frame(ab,b)
  y1=kruskal.test(ab~b,data=aa)
  Pvaluekw[i]<-y1$p.value
}
reactomees<-reactomees[complete.cases(reactomees),]
reactomees<-t(reactomees)
reactomees<-as.data.frame(reactomees)
reactomees$P<-Pvaluekw
reactomees<-reactomees[which(reactomees$P < 0.05),]
reactomees<-reactomees[,-ncol(reactomees)]

groupccle <- groupccle %>% 
  rownames_to_column("sample")
annotation <- groupccle %>% arrange(groupccle) %>% column_to_rownames("sample")
a <- groupccle %>% arrange(groupccle) 
b <- t(reactomees) %>% 
  as.data.frame() %>% 
  rownames_to_column("sample") 
c <- inner_join(a,b,"sample") %>% .[,-2] %>% column_to_rownames("sample") %>% t(.)
range(c)
bk <- c(seq(-2,-0.01,by=0.01),seq(0,2,by=0.01))

ann_colors = list(
  groupccle= c(MPC1="#2a9d8c",MPC2= "#e9c46b",MPC3="#e66f51"))
pheatmap(c,
         annotation = annotation,
         scale = "row",
         
         show_colnames =F,
         color = colorRampPalette(c("#addcca", "white","#F95A37"))(400),
         annotation_colors = ann_colors[1],
         cluster_cols =F,
         
         fontsize = 10,
         fontsize_row=7,
         fontsize_col=7,
         breaks=bk)
a <- groupccle %>% arrange(groupccle)%>% mutate(sample=substring(.$sample,1,12))
b <- t(reactomees) %>% 
  as.data.frame() %>% 
  rownames_to_column("sample")%>% mutate(sample=substring(.$sample,1,12))
c <- inner_join(a,b,"sample") %>% .[,-2] %>% column_to_rownames("sample") %>% t(.)
range(c)
bk <- c(seq(-2,-0.01,by=0.01),seq(0,2,by=0.01))
pheatmap(c,
         annotation = annotation,
         scale = "row",
         cutree_rows = 3,
         show_colnames =F,
         color = colorRampPalette(c("navy", "white", "red"))(400),
         cluster_cols =F,
         fontsize = 10,
         fontsize_row=7,
         fontsize_col=7)

reactomees<-reactomees%>%t()%>%as.data.frame()

a<-reactomees
a<-a%>%rownames_to_column("sample")
groupccle<-groupccle%>%rownames_to_column("sample")
b<-groupccle
c<-inner_join(a,b,by=c("sample"="sample"))
c$groupccle<-as.factor(c$groupccle)
my_comparisons <- list( c("MPC1", "MPC2"), c("MPC2", "MPC3"), c("MPC1", "MPC3") )
p<-c%>%
  ggplot(aes(groupccle, oleylcarnitine,fill=groupccle)) +
  geom_dotplot(binwidth = 0.16,binaxis='y', stackdir='center', 
               position=position_dodge(0.9)) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar", color="black", width = .07, 
               position = position_dodge(.9), size = .8)+
  stat_summary(fun = "median", geom = "crossbar", 
               mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.3)+
  stat_compare_means(comparisons=my_comparisons,
                     method = "wilcox.test",
                     label = "p.signif",
                     hide.ns=T,
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))
p+  theme_classic() 
write.table(reactomees,'result.txt',sep = "\t",row.names = T,col.names = NA,quote = F )
save(groupccle,file = 'groupcclenew.rda')


####IMMUNE COHORT####
setwd("E:/R/ESO3/IMMU")
library(GSVA)
library(limma)
library(GSEABase)
library(org.Hs.eg.db)
library(clusterProfiler)
keggEs7=gsva(as.matrix(eset),meset,abs.ranking=TRUE,method="gsva")


library(NMF)
library(doMPI)

res1 <- nmf(keggEs7,2:7,"lee",nrun=10)
plot(res1)
V.random7 <- randomize(keggEs7)
res.random7 <- nmf(V.random7, 2:7,"lee",nrun=10)
plot(res1, res.random7)
res8<-nmf(keggEs7,3,"lee",nrun=10)
plot(res8)
group1 <- predict(res8)
group1<-as.data.frame(group1)

table(group1$group1)
identical(group$group,group1$group1)
group<-group1
group<- as.data.frame(group)
group$group <- factor(group$group,levels=c(1,2,3))
save(group,file = 'groupme.rda')
group$group<-paste0("MPC",group$group)
imcli$group<-group$group
group$response<-imcli$treatment.response.for.gene.profiling
pdata$group<-group$group
pdata<-pdata%>%rownames_to_column("sample")
pdata<-inner_join(pdata,group,by=c("sample"="sample"))
grouptab<- xtabs(~ group+`Immune phenotype`,data =pdata)
grouptab<-as.data.frame(grouptab)
grouptab<-grouptab[-c(1:3),]
fix(grouptab)
library(RColorBrewer)
library(ggplot2)
p<-ggplot(data=grouptab,aes(group,Freq,fill=Immune.phenotype))+  
  geom_bar(stat="identity", position="fill",color="black", 
           width=0.6,size=0.25)+    
  xlab("Sample") + ylab("Expression values")+ 
  theme(    axis.title=element_text(size=15,face="plain",color="black"),    
            axis.text = element_text(size=12,face="plain",color="black"),    
            legend.title=element_text(size=14,face="plain",color="black"),    
            legend.position = "top"  )

p+  theme_classic() 











library(ggpubr)
library(rstatix)
p_bar <- ggbarplot(grouptab,x="group",y="Freq" ,fill = "Immune.phenotype")
p_bar <- ggbarplot(grouptab, x = 'group', y="Freq" ,fill = "Immune.phenotype", add = 'mean_sd’,
                  color = ‘gray30', position = position_dodge(0.6), width = 0.5, size = 1, legend = 'top') +
  scale_fill_manual(values = c('#E7B800', '#00AFBB','#8b2fbd'))+
  labs(x='group', y='frequency',fill='Immune.phenotype')+
  stat_compare_means(aes(group = Immune.phenotype),
                     method = "wilcox.test",
                     label = "p.signif" ,
                     hide.ns = T,
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))+
  theme(text = element_text(size=10))
fisher.test(pdata$TC_Level,pdata$group)
p_bar
dev.off()
library(IMvigor210CoreBiologies)
library(tidyverse)
data(cds)
expMatrix <- counts(cds)
eff_length2 <- fData(cds)[,c("entrez_id","length","symbol")]
rownames(eff_length2) <- eff_length2$entrez_id
head(eff_length2)
feature_ids <- rownames(expMatrix)
expMatrix <- expMatrix[feature_ids %in% rownames(eff_length2),]
mm <- match(rownames(expMatrix),rownames(eff_length2))
eff_length2 <- eff_length2[mm,]

x <- expMatrix/eff_length2$length
eset <- t(t(x)/colSums(x))*1e6
summary(duplicated(rownames(eset)))

eset <- IOBR::anno_eset(eset = eset,
                        annotation = eff_length2,
                        symbol = "symbol",
                        probe = "entrez_id",
                        method = "mean")
tumor_type <- "blca"
if(max(eset)>100) eset <- log2(eset+1)
pdata <- pData(cds)
colnames(pdata) <- gsub(colnames(pdata),pattern = " ",replacement = "_")
pdata <- rownames_to_column(pdata[,c("binaryResponse",
                                     "FMOne_mutation_burden_per_MB",
                                     "Neoantigen_burden_per_MB",
                                     "censOS","os")],var = "ID")
colnames(pdata)<-c("ID","BOR_binary","TMB","TNB","status","time")
pdata<-pdata[!is.na(pdata$BOR_binary),]
pdata$BOR_binary<-ifelse(pdata$BOR_binary=="CR/PR","R","NR")
save(expMatrix,pdata,file = "expcli_IMvigor210.Rdata")
comimmu<-intersect(pdata$ID,colnames(eset))
eset<-eset[,comimmu]
pdata<-pdata%>%rownames_to_column("ID")
group<-group%>%rownames_to_column("sample")
pdata<-inner_join(pdata,group,by=c("ID"="sample"))
pdata$group<-paste0("MP" , pdata$group)
interimmu<-intersect(pdata$ID,colnames(eset))
eset<-eset[,interimmu]
identical(pdata$ID,group$sample)
imcli$group<-group$group

pdata<-pdata[complete.cases(pdata),]
pdata$group<-as.numeric(pdata$group)
class<-aggregate(pdata$`Neoantigen burden per MB`,         
                 list(pdata$group),
                 median)
pdata$s3 <- ifelse(pdata$`Neoantigen burden per MB` > class[3,2],"H3","L3")
export(pdata,"pdatacombine.csv")
pdata1<-import("pdatacombine.csv")
export(pdata,"pdata.csv")
pdata<-import("pdata.csv")


library(survival)
library(survminer)
library(ggplot2)
res9<-nmf(keggEs7,3,"lee",nrun=10)
res8 <- nmf(keggEs7,2:7,"lee",nrun=10)
plot(res8)
group <- predict(res9)
group<-as.data.frame(group)
pdata$group<-group$group
fitd <- survdiff(Surv(os, censOS=="1") ~ s1,
                 data      = pdata1,
                 na.action = na.exclude)
pValue <- 1 - pchisq(fitd$chisq, length(fitd$n) - 1)
fit <- survfit(Surv(os, censOS=="1")~ s1, data =pdata1)
table(pdata1$s2)
p.lab <- paste0("P", ifelse(pValue < 0.001, " < 0.001", paste0(" = ",round(pValue, 3))))
ggsurvplot(fit,
           data = pdata,
           pval = p.lab,
           conf.int = TRUE, 
           risk.table = TRUE, 
           risk.table.col = "strata",
           palette = "jco", 
           legend.labs = c("H1S1","H2S2","H3S3","H4S4","H5S5","H6S6"), # 图例
           size = 1,
           xlim = c(0,20), 
           break.time.by = 5, 
           legend.title = "",
           surv.median.line = "hv", 
           ylab = "Survival probability (%)", 
           xlab = "Time (Months)",
           ncensor.plot = TRUE, 
           ncensor.plot.height = 0.25,
           risk.table.y.text = FALSE)
dev.off()
ggboxplot(pdata, 
          x="group", 
          y="TMB", 
         
          color = "group",
          palette = c('#2a9d8c','#e9c46b',"#e66f51"),
          legend.title = "type",
          order=c("MPC1","MPC2","MPC3"),
          add = "jitter")+ 
  stat_boxplot(geom = "errorbar",width=0.15,aes(color =group))+ 
  stat_compare_means(comparisons = my_comparisons,
                     method = "wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("***", "**", "*", "ns")),
                     label = "p.signif")

dev.off()  

p<-pdata%>%
  ggplot(aes(group, TMB,fill=group)) +
  geom_dotplot(binaxis='y', stackdir='center', 
               position=position_dodge(0.9)) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar", color="black", width = .07, 
               position = position_dodge(.9), size = .8)+
  stat_summary(fun = "median", geom = "crossbar", 
               mapping = aes(ymin = ..y.., ymax = ..y..), width = 0.3)+
  stat_compare_means(comparisons=my_comparisons,
                     method = "wilcox.test",
                     label = "p.signif",
                     hide.ns=T,
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))
p+  theme_classic() 



eset<-eset%>%t()%>%as.data.frame()
identical(pdata$ID,rownames(eset))
eset$group<-group$group
eset$group<-paste0("MP",eset$group)
ggplot(eset, aes(x = group, y = CD274, fill = group)) +
  geom_violin(trim = F) + scale_fill_brewer(palette = "Dark2") + 
  theme_classic() + geom_boxplot(width = 0.1, fill = "white") + 
  theme(text = element_text(size = 18), 
        axis.text.x = element_text(size = 18), 
        plot.title = element_text(hjust = 0.5), 
        legend.position = c(0.15,0.9)) + 
  xlab("group") +  ylab("PD-L1")+
  stat_compare_means(comparisons=my_comparisons,
                     method = "wilcox.test",
                     label = "p.signif",
                     hide.ns=T,
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))
dev.off() 
####DRUG SENSITIVITY####
setwd("E:/R/ESO3/DRUG")
install.packages("pRRophetic_0.5.tar.gz", repos = NULL, dependencies = TRUE)
install.packages("oncoPredict")
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
BiocManager::install(c('sva', 'car', 'genefilter'))
rm(list = ls()) 
options(stringsAsFactors = F)
library(oncoPredict)
library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)
library(data.frame)
library(rio)
library(tidyverse)
library(pheatmap)
th=theme(axis.text.x = element_text(angle = 45,vjust = 0.5))
group$group<-paste0("MPC",group$group)
GDSC2_Expr = readRDS('GDSC2_Expr (RMA Normalized and Log Transformed).rds')
GDSC2_Res = readRDS("GDSC2_Res.rds")
GDSC2_Res <- exp(GDSC2_Res) 
testExpr<- counts3

testExpr[1:4,1:4]  
testExpr<-as.matrix(testExpr)
dim(testExpr) 

calcPhenotype(trainingExprData = GDSC2_Expr,
              trainingPtype = GDSC2_Res,
              testExprData = testExpr,
              batchCorrect = 'eb', 
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'rawData' )
drugcounts<-import("DrugPredictions.csv")
drugcounts<-drugcounts%>%column_to_rownames("V1")
drugcounts<-drugcounts%>%t()%>%as.data.frame()

reactomees<-t(drugcounts)
reactomees<-as.data.frame(reactomees)
Pvaluekw<-c(rep(0,ncol(reactomees)))
design<-data.frame(rownames(reactomees),group$group)
for(i in 1:ncol(reactomees)){
  ab<-as.numeric(reactomees[1:nrow(reactomees),i])
  b<-design$group.group
  aa<-data.frame(ab,b)
  y1=kruskal.test(ab~b,data=aa)
  Pvaluekw[i]<-y1$p.value
}
reactomees<-t(reactomees)
reactomees<-as.data.frame(reactomees)
reactomees$P<-Pvaluekw
reactomees<-reactomees[which(reactomees$P < 0.05),]
reactomees<-reactomees[,-ncol(reactomees)]






fix(predictdata)

my_comparisons=list(c("MPC1","MPC2"),
                    c("MPC1","MPC3"),
                    c("MPC2","MPC3"))
export(reactomees,"filter.csv")
ggboxplot(predictdata, 
          x="cluster", 
          y="IC50", 
          color = "cluster",
          palette = c('#2a9d8c','#e9c46b',"#e66f51"),
          legend.title = "type",
          order=c("MPC1","MPC2","MPC3"),
          add = "jitter")+ 
  stat_boxplot(geom = "errorbar",width=0.15,aes(color =cluster))+ 
  stat_compare_means(comparisons = my_comparisons,
                     method = "wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("***", "**", "*", "ns")),
                     label = "p.signif")


dev.off()

####CORRPLOT####
install.packages("corrplot")
library(corrplot)
corr<-round(cor(t(keggEs2),method = "pearson"),3)
sampleDist <- dist(keggEs2) 
sampleDistMatrix <- as.matrix(sampleDist) 
colnames(sampleDistMatrix) <- NULL 
library(RColorBrewer)
colors <- colorRampPalette(rev(brewer.pal(9,"Blues")))(255) 
bk1 = c(seq(-0.5,-0.01,by=0.01),seq(0,0.5,by=0.01))
pheatmap(corr, 
         
         color = colorRampPalette(c("#15559a", "white","#F95A37" ))(200),
         breaks = bk1) 
corrplot(corr,order="hclust",addrect=5,
         col=colorRampPalette(c("#204051", "white","#e75a5f" ))(400),
         bg="khaki1",
         tl.cex=0.1)

