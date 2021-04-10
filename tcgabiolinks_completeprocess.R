###1.Download data from TCGA
library(TCGAbiolinks)
library(dplyr)
library(DT)
library(SummarizedExperiment)
data_type <- "Gene Expression Quantification"
data_category <- "Transcriptome Profiling"
workflow_type <- "HTSeq - Counts"
cancer_type = "TCGA-LUAD"
#legacy default is FALSE, downloading from hg38
query_TranscriptomeCounts <- GDCquery(project = cancer_type,
                                      data.category = data_category,
                                      data.type = data_type,
                                      workflow.type = workflow_type)
#getresults:to get the sample barcodes from "query_TranscriptomeCounts"
samplesDown <- getResults(query_TranscriptomeCounts,cols=c("cases"))#594 in total
#to select the solid tumor samples in "samplesDown"
#TP：Primary Solid Tumor
#NT：Solid Tissue Normal
dataSmTP <- TCGAquery_SampleTypes(barcode = samplesDown,
                                  typesample = "TP")#533
dataSmNT <- TCGAquery_SampleTypes(barcode = samplesDown,
                                  typesample = "NT")#59
#to filter the qualified barcodes according to groups
queryDown <- GDCquery(project = cancer_type, 
                      data.category = data_category,
                      data.type = data_type, 
                      workflow.type = workflow_type, 
                      barcode = c(dataSmTP, dataSmNT))
#to download the data that comply with the conditions
GDCdownload(queryDown,
            method = "api",
            directory = "GDCdata",
            files.per.chunk = 6)#592 in total
#method：api or client，api is faster but easy to be interrupted
#files.per.chunk：For example, when use api, if files.per.chunk = 6, it would download 6 files at a time to reduce error
###2.Data optitions
#to transfer the result of GDCquery() into R-accessible pattern
dataPrep1 <- GDCprepare(query = queryDown,save = TRUE, save.filename =
                          "1-dataprep1_LUAD_case.rda")

###3.Data pre-processing
#to remove the abnormal data
dataPrep2 <- TCGAanalyze_Preprocessing(object = dataPrep1,
                                       cor.cut = 0.6,
                                       datatype = "HTSeq - Counts")#592
write.csv(dataPrep2,file = "2-dataPrep2_LUAD.csv",quote = FALSE)
###4.To select the barcodes that tumor purity > 60%
#There are estimate, absolute, lump, ihc, cpe to decide
#cpe: The mean value of purity after normalizing the contents in all methods to make them have the same mean value and standard diviation
purityDATA <- TCGAtumor_purity(colnames(dataPrep1), 0, 0, 0, 0, 0.6)#381 in total
#filtered: the barcodes of normal tissue (control)
#pure_barcodes: the barcodes of tumor tissue
Purity.LUAD <- purityDATA$pure_barcodes#322
normal.LUAD <- purityDATA$filtered#59
###5.Merge the expression set of tumor and normal
puried_data <-dataPrep2[,c(Purity.LUAD,normal.LUAD)]#共381个
###6.Gene annotation
rownames(puried_data)<-rowData(dataPrep1)$external_gene_name
write.csv(puried_data,file = "3-puried.LUAD.csv",quote = FALSE)
###7.Expression set for differential analysis after normalization and filteration
BiocManager::install("EDASeq")
library(EDASeq)
#normalization
dataNorm <- TCGAanalyze_Normalization(tabDF = puried_data,
                                      geneInfo = geneInfo,
                                      method = "gcContent")
#remove low count genes
dataFilt_LUAD_final <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                             method = "quantile", 
                                             qnt.cut =  0.25)
a <- c(rep("tumor",322),rep("normal",59))
a <- t(as.data.frame(a))
colnames(a)<-colnames(dataFilt_LUAD_final)
dataFilt_LUAD_final <- rbind(dataFilt_LUAD_final,a)
write.csv(dataFilt_LUAD_final,file = "4-TCGA_LUAD_final.csv",quote = FALSE)
###8.Differential analysis
#tumor group
mat1 <- dataFilt_LUAD_final[,1:322]
mat1 <- log(mat1+1)
#normal group
mat2 <- dataFilt_LUAD_final[,323:381]
mat2 <- log(mat2+1)
#analysis 
BiocManager::install("edgeR")
library(edgeR)
Data_DEGs <- TCGAanalyze_DEA(mat1 = mat1,
                             mat2 = mat2,
                             Cond1type = "Tumor",
                             Cond2type = "Normal",
                             pipeline="limma",
                             batch.factors = c("TSS"),
                             voom = TRUE,
                             contrast.formula = "Mycontrast=Tumor-Normal")
write.csv(Data_DEGs,file = "5-LUAD_DEGs.csv",quote = FALSE)
###9.Enrichment analysis
#to set the threshold of logFC 
Data_DEGs_high_expr <- Data_DEGs[Data_DEGs$logFC >=1,]
Genelist <- rownames(Data_DEGs_high_expr)
ansEA <- TCGAanalyze_EAcomplete(TFname="DEA genes Tumor Vs Normal",
                                Genelist)
#visualization
TCGAvisualize_EAbarplot(tf  = rownames(ansEA$ResBP),
                        GOBPTab = ansEA$ResBP,
                        GOCCTab = ansEA$ResCC,
                        GOMFTab = ansEA$ResMF,
                        PathTab = ansEA$ResPat,
                        nRGTab = Genelist,
                        nBar = 20, #the number of barplot
                        filename = "TCGAvisualize_EAbarplot_Output.pdf")
###10.Clustering analysis
BiocManager::install("ConsensusClusterPlus")
library(ConsensusClusterPlus)
res.hc <- TCGAanalyze_Clustering(Data_DEGs, 
                                 method = "hclust", 
                                 methodHC = "ward.D2")
plot(res.hc)
###11.Enrichment using ClusterProfiler
library(clusterProfiler)
library(org.Hs.eg.db)
Data_DEGs$symbol <- row.names(Data_DEGs)
gene <- row.names(Data_DEGs)
class(gene)
gene.df<-bitr(gene, fromType = "SYMBOL", 
              toType = c("ENSEMBL","ENTREZID"),
              OrgDb = org.Hs.eg.db)
head(gene.df)
ego_cc<-enrichGO(gene       = Data_DEGs$symbol,
                 OrgDb      = org.Hs.eg.db,
                 keyType = "SYMBOL",
                 ont        = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05)
ego_bp<-enrichGO(gene       = gene.df$ENTREZID,
                 OrgDb      = org.Hs.eg.db,
                 keyType = "ENTREZID",
                 ont        = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05)
ego_mf<-enrichGO(gene       = gene.df$ENTREZID,
                 OrgDb      = org.Hs.eg.db,
                 keyType = "ENTREZID",
                 ont        = "MF",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05)

dotplot(ego_cc,showCategory = 20,title="The GO_CC enrichment analysis of all DEGs ")

library(stringr)
kk<-enrichKEGG(gene      =gene.df$ENTREZID,
               organism = 'hsa',
               pvalueCutoff = 0.05)
kk=DOSE::setReadable(kk, OrgDb='org.Hs.eg.db',keyType='ENTREZID')
barplot(kk,showCategory = 25, title="The KEGG enrichment analysis of all DEGs")
cnetplot(kk, foldChange=gene, circular = TRUE, colorEdge = TRUE)
emapplot(kk)

genelist <- Data_DEGs$logFC
names(genelist) <- Data_DEGs[,7]
genelist <- sort(genelist, decreasing = TRUE)
gsemf <- gseGO(genelist,
               OrgDb = org.Hs.eg.db,
               keyType = "SYMBOL",
               ont="BP",
               pvalueCutoff = 0.9)
gseaplot(gsemf, geneSetID="GO:0000280")

library(pathview)
library(gage)
library(gageData)
data("kegg.sets.hs")
data("sigmet.idx.hs")
kegg.sets.hs =  kegg.sets.hs[sigmet.idx.hs]
head(kegg.sets.hs,3)
gene.df<-bitr(gene, fromType = "SYMBOL", 
              toType = c("ENSEMBL","ENTREZID"),
              OrgDb = org.Hs.eg.db)
head(Data_DEGs)
foldchanges = Data_DEGs$logFC
names(foldchanges)= gene.df$ENTREZID
#There was an error when running as follow
#Error in names(foldchanges) = gene.df$ENTREZID : 
#'names'属性的长度[14370]必需和矢量的长度[12963]一样
#duplicated(gene.df$SYMBOL) to find the duplicated value
#And I found some duplicated rows
#To solve:
#gene.df <- gene.df %>% distinct(SYMBOL,.keep_all = TRUE)
# Then removed
head(foldchanges)
keggres = gage(foldchanges, gsets = kegg.sets.hs, same.dir = TRUE)
lapply(keggres, head)

keggrespathways = data.frame(id=rownames(keggres$greater), keggres$greater) %>% 
  tbl_df() %>% 
  filter(row_number()<=10) %>% 
  .$id %>% 
  as.character()
keggrespathways
keggresids = substr(keggrespathways, start=1, stop=8)
keggresids
plot_pathway = function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa", new.signature=FALSE)
tmp = sapply(keggresids, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, species="hsa"))
threshold<-as.factor(
  (Data_DEGs$logFC>1|Data_DEGs$logFC<(-1) & 
     Data_DEGs$P.Value<0.05 ))
library(ggplot2)
ggplot(Data_DEGs,aes(x= logFC,
                     y= -1*log10(P.Value),colour=threshold))+xlab("log2 fold-change")+ylab("-log10 p-value")+geom_point() 
###12.Survival analysis
library(survminer)
library(survival)
#to download clinic information
clin.LUAD <- GDCquery_clinic("TCGA-LUAD", "clinical")
#to check the influence of "gender"
TCGAanalyze_survival(clin.LUAD,
                     clusterCol="gender",
                     risk.table = FALSE,
                     xlim = c(100,1000),
                     ylim = c(0.4,1),
                     conf.int = FALSE,
                     pvalue = TRUE,
                     color = c("Dark2"))
##to investigate the influence of single gene
#for example,METTL14
dataFilt_LUAD_final <- read.csv("4-TCGA_LUAD_final.csv",row.names = 1,check.names=FALSE)
samplesTP <- TCGAquery_SampleTypes(colnames(dataFilt_LUAD_final), typesample = c("TP"))
METTL14 <- dataFilt_LUAD_final[c("METTL14"),samplesTP]
#to modify sample name
#to transfer like this: TCGA-91-6840-01A-11R-1949-07
names(METTL14) <- sapply(strsplit(names(METTL14),'-'),function(x) paste0(x[1:3],collapse="-"))
METTL14 <- t(METTL14)
#to merge the gene and clinic data
clin.LUAD$"METTL14" <- METTL14[match(clin.LUAD$submitter_id,rownames(METTL14)),]
#to select the needed information for survival analyze
df<-subset(clin.LUAD,select =c(submitter_id,vital_status,days_to_death,days_to_last_follow_up,METTL14))
#to remove NA
df <- df[!is.na(df$METTL14),]
#to classify into high expression or low expression
df$METTL14 <- as.numeric(df$METTL14)
df$exp <- ''
df[df$METTL14 >= mean(df$METTL14),]$exp <- "H"
df[df$METTL14 < mean(df$METTL14),]$exp <- "L"
#analyse by TCGAanalyze_survival function
TCGAanalyze_survival(df,
                     legend = "METTL14",
                     clusterCol="exp",
                     risk.table = FALSE,
                     conf.int = FALSE,
                     pvalue = TRUE,
                     color = c("Dark2"))
###13.Differential gene - volcano and heatmap 
##to load the pre-processing data
TCGA_LUAD_data <- read.csv(file = "4-TCGA_LUAD_final.csv",
                           header = T,
                           row.names = 1,
                           check.names = FALSE )#to assure the column names not to be corrected automatically
#to get the barcodes
samplesNT <- TCGAquery_SampleTypes(colnames(TCGA_LUAD_data), typesample = c("NT"))
samplesTP <- TCGAquery_SampleTypes(colnames(TCGA_LUAD_data), typesample = c("TP"))
#pair the tumor and normal barcodes
paired <- intersect(substr(samplesNT, 1, 12),substr(samplesTP, 1, 12))
length(paired)
#to merge two groups using the 1st~12th data
NT <- data.frame(NT1=substr(samplesNT,1,12),NT2=samplesNT)
TP <- data.frame(TP1 =substr(samplesTP,1,12),TP2=samplesTP)
TP_NT <- merge(TP,NT,by.x = "TP1",by.y = "NT1")
head(TP_NT,3)
#to get paired normal barcodes
TP <- TP_NT$TP2
#to get paired tumor barcodes
NT <- TP_NT$NT2
##downloading and pre-processing
#Some data has been saved so just load it into workspace.
#queryDown <- GDCquery(project = "TCGA-LUAD",
#data.category = "Transcriptome Profiling",
#data.type = "Gene Expression Quantification",
#workflow.type = "HTSeq - Counts")
#GDCdownload(queryDown,
#method = "api",
#directory = "GDCdata",
#files.per.chunk = 6)
#dataPrep1 <- GDCprepare(query = queryDown, save = F)
#load the data
load("I:/tcgabiolinks/TCGAbiolinks/1-dataprep1_LUAD_case.rda")
#dataPrep2 <- TCGAanalyze_Preprocessing(object = dataPrep1,
#                                       cor.cut = 0.6,
#                                       datatype = "HTSeq - Counts")
#dataPrep1 <- data
#purityDATA <- TCGAtumor_purity(colnames(dataPrep1), 0, 0, 0, 0, 0.6)
#Purity.LUAD <- purityDATA$pure_barcodes#322
#normal.LUAD <- purityDATA$filtered#59
#puried_data <-dataPrep2[,c(Purity.LUAD,normal.LUAD)]
#write the result
#rownames(puried_data)<-rowData(dataPrep1)$external_gene_name
#to normalize library size and GC abundance
#dataNorm <- TCGAanalyze_Normalization(tabDF = puried_data,
#                                     geneInfo = geneInfo,
#                                      method = "gcContent")
#过滤低count基因
#dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
 #                                 method = "quantile", 
  #                                qnt.cut =  0.25)
#write.csv(dataFilt,file = "paired_TCGA_LUAD_final.csv",quote = FALSE)
##differential analysis
#dataFilt_LUAD_final <- read.csv(file = "4-TCGA_LUAD_final.csv",header = T,row.names = 1,check.names = FALSE)
#tumor vs normal
mat1 <- dataFilt_LUAD_final[,1:322]#322
mat2 <- dataFilt_LUAD_final[,323:381]#59
DEG_LUAD_edgeR <- TCGAanalyze_DEA(mat1 = mat1,
                                  mat2 = mat2,
                                  Cond1type = "Tumor",
                                  Cond2type = "Normal",
                                  pipeline="edgeR",
                                  batch.factors = c("TSS"),
                                  voom = FALSE,
                                  contrast.formula = "Mycontrast=Tumor-Normal",
                                  fdr.cut = 0.01,
                                  logFC.cut = 1,
                                  method = "glmLRT")
write.csv(DEG_LUAD_edgeR,file = "5-paired_DEG_by_edgeR.csv")
##to increase average gene expression in different grouping conditions  
#to get the barcodes respectively
samplesNT <- TCGAquery_SampleTypes(colnames(dataFilt_LUAD_final), typesample = c("NT"))
samplesTP <- TCGAquery_SampleTypes(colnames(dataFilt_LUAD_final), typesample = c("TP"))
dataDEGsFilt <- DEG_LUAD_edgeR[abs(DEG_LUAD_edgeR$logFC) >= 1,]
str(dataDEGsFilt)
#to get the expression set respectively
dataTP <- dataFilt_LUAD_final[,samplesTP]
dataTN <- dataFilt_LUAD_final[,samplesNT]
#to import the parameter
dataDEGsFiltLevel <- TCGAanalyze_LevelTab(FC_FDR_table_mRNA = dataDEGsFilt,
                                          typeCond1 ="Tumor",
                                          typeCond2 = "Normal",
                                          TableCond1 = dataTP,
                                          TableCond2 = dataTN)
head(dataDEGsFiltLevel, 2)
##PCA analysis
pca <-TCGAvisualize_PCA(dataFilt = dataFilt_LUAD_final, 
                        dataDEGsFiltLevel = dataDEGsFiltLevel,
                        ntopgenes = 100,
                        group1 = samplesTP, 
                        group2 = samplesNT)
##Heatmap of deg
#to get the expression set
datDEGs <- dataFilt_LUAD_final[match(rownames(DEG_LUAD_edgeR),rownames(dataFilt_LUAD_final)),]
#to obtain clinic information
query <- GDCquery(project = "TCGA-LUAD",
                  data.category = "Clinical",
                  file.type = "xml",
                  barcode = substr(colnames(datDEGs),1,12))
GDCdownload(query)  
clinical <- GDCprepare_clinic(query,"patient")
#to match the expression set and clinic information
datDEGs_test_barcodes <- as.data.frame(substr(colnames(datDEGs),1,12),ncol=1)
colnames(datDEGs_test_barcodes) <-"LUAD_patient_barcode"
m <- clinical[match(datDEGs_test_barcodes[,1], clinical[ , 1]),]
str(m)
#check the duplicated data
table(duplicated(m))
#heatmap of raw data
pheatmap::pheatmap(datDEGs,scale = "row",show_rownames = F,show_colnames = F)
#add metadata
col.mdat <- data.frame(Sex=m$gender,
                       status=m$vital_status,
                       group=c(rep("tumor",322),rep("normal",59)))
#annotations(row) should match samples(column)
rownames(col.mdat) <- colnames(datDEGs) 
#to set the legend limit
bk <- c(seq(-1,6,by=0.01))
#plot by pheatmap
pheatmap::pheatmap(datDEGs,scale = "row",show_rownames = F,show_colnames = F,
                   annotation_col = col.mdat,
                   border_color=NA,
                   main = "Heatmap by pheatmap(edgeR)",
                   filename = "Heatmap_by_pheatmap.pdf",
                   color =c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),
                            colorRampPalette(colors = c("white","red"))(length(bk)/2)),#to set the legends color
                   legend_breaks=seq(-1,6,2),
                   breaks=bk )
##volcano
#highlight the genes that logfc>=8
DEG.LUAD.filt<-DEG_LUAD_edgeR[which(abs(DEG_LUAD_edgeR$logFC) >= 8), ]
str(DEG_LUAD_edgeR)
TCGAVisualize_volcano(DEG_LUAD_edgeR$logFC, 
                      DEG_LUAD_edgeR$FDR,
                      filename = "TumorvsNormal_FC8.edgeR.pdf", 
                      xlab = "logFC",
                      names = rownames(DEG_LUAD_edgeR), 
                      show.names = "highlighted",
                      x.cut = 1, 
                      y.cut = 0.01, 
                      highlight = rownames(DEG_LUAD_edgeR)[which(abs(DEG_LUAD_edgeR$logFC) >= 8)],
                      highlight.color = "orange",
                      title = "volcano plot by edgeR")
# title = "volcano plot by limma")
###14.Infiltration analysis by CIBERSORT
#after the puried_data
condition <- factor(c(rep("tumor",322),rep("normal",59)), levels = c("tumor","normal"))
condition
table(condition)
coldata <- data.frame(row.names=colnames(puried_data), condition)
coldata
table(coldata)
dds <- DESeqDataSetFromMatrix(puried_data, coldata, design= ~ condition)
dds <- dds[rowSums(counts(dds)) > 1, ] 
dds <- estimateSizeFactors(dds) 
normalized_counts <- counts(dds,normalized=T) 
dds$condition <- relevel(dds$condition,ref = "normal")
dds <- DESeq(dds)
dds
vsdata <- vst(dds, blind=FALSE)
exprSet=assay(vsdata)
write.table(exprSet,file = "exprset.txt")


library(preprocessCore)
source("cibersort.R")
LM22.file <- "LM22.txt"
exp.file <- "exprset.txt"
TME.results = CIBERSORT(LM22.file, exp.file, perm = 1000, QN = TRUE)
write.csv(TME.results,file = "TME_cibersort.csv")

###15.Infiltration analysis by MCPcounter
install_github("ebecht/MCPcounter",ref="master", subdir="Source",force = T)
exprset <- read.csv(file = "4-TCGA_LUAD_final.csv",header = T,row.names = 1)
exprset <- log2(exprset+1)

MCPcounter.estimate <- MCPcounter.estimate(
  exprset,
  featuresType=c('affy133P2_probesets','HUGO_symbols','ENTREZ_ID')[2],           
  probesets=read.table(curl('http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/probesets.txt'),sep='\t',stringsAsFactors=FALSE,colClasses='character'),
  genes=read.table(curl('http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt'),sep='\t',stringsAsFactors=FALSE,header=TRUE,colClasses='character',check.names=FALSE)
)

group <- factor(c(rep("tumor",322),rep("normal",59)))
group <- t(as.data.frame(group))
colnames(group) <- colnames(MCPcounter.estimate)
MCPcounter.estimate <- rbind(MCPcounter.estimate,group)
write.csv(MCPcounter.estimate,file = "mcp_result_BRC.csv")

data<-data.frame(Sample<-c(rep('control1',3),rep('control2',3),rep('control3',3),rep('treat1',3),rep('treat2',3),rep('treat3',3),rep('treat4',3)), contion<-rep(c('Cell','Tissue','Organ'),7), value<-c(503,264,148,299,268,98,363,289,208,108,424,353,1,495,168,152,367,146,48,596,143))

colnames(data)=c('sample',"contion","value")

ggplot(mcp,
       mapping = aes(Sample,value,fill=contion))
+geom_bar(stat='identity',position='fill') 
+labs(x = 'Sample',y = 'frequnency') 
+theme(axis.title =element_text(size = 16),axis.text =element_text(size = 14, color = 'black'))
+theme(axis.text.x = element_text(angle = 45, hjust = 1))


#to show interested genes high or low expression in one plot 
condition <- factor(c(rep("tumor",322),rep("control",59)), levels = c("tumor","control"))
condition
table(condition)

exprset <- read.csv("4-TCGA_LUAD_final.csv",header = T,row.names = 1)
exprset <- log2(exprset+1)

m6agenes <- exprset[c("METTL3","METTL14","WTAP","RBM15","ZC3H13",
                      "YTHDF1","YTHDF2","YTHDF3","YTHDC1","YTHDC2","IGF2BP1","IGF2BP2","HNRNPC",
                      "FTO","ALKBH5"),]
m6agenes <- as.data.frame(t(m6agenes))
m6agenes$condition <- c(rep("tumor",322),rep("control",59))
#there is a loop
#a <- as.data.frame(m6agenes$ALKBH5)
#b <- as.data.frame(m6agenes$condition)
#c <- as.data.frame(rep("ALKBH5",381))
#d <- cbind(a,c,b)
#colnames(d) <- c("Expression","gene","group")
#dat <- d
#dat <- rbind(dat,d)
#dat$gene = factor(dat$gene, levels=c("METTL3","METTL14","WTAP","RBM15","ZC3H13",
                                     #"YTHDF1","YTHDF2","YTHDF3","YTHDC1","YTHDC2","IGF2BP1","IGF2BP2","HNRNPC",
                                     #"FTO","ALKBH5"))
#write.csv(dat,file = "dat_luad.csv")
#I've done previous steps,so just load the data
dat <- read.csv("7-dat_luad.csv",row.names = 1)
compare_means( Expression ~ group, data = dat, group.by = "gene")
p <- ggboxplot(dat, x = "group", y = "Expression",
               color = "group", palette = "jco",
               add = "jitter",
               facet.by = "gene",
               short.panel.labs = FALSE)
p + stat_compare_means(label = "p.format")
pdf("LUAD.pdf",width=20,height=20)
dev.off
#when "dev.off()" doesn't work use "while (!is.null(dev.list()))  dev.off()"

#heatmap by complexheatmap
condition <- factor(c(rep("tumor",322),rep("control",59)), levels = c("control","tumor"),ordered = F)
group <- data.frame(group=condition)
rownames(group)=colnames(m6agenes)
samples <- condition
m6agenes <- as.matrix(m6agenes)
#ATTENTION: before run, class the type of input data set, it must be "num"
heat <- Heatmap(m6agenes, 
                col = colorRampPalette(c('navy', 'white', 'firebrick3'))(100), #set the high and low colors
                heatmap_legend_param = list(grid_height = unit(10,'mm')),  #set the height of legends
                show_row_names = TRUE,  #show the gene name
                top_annotation = HeatmapAnnotation(Group = samples, 
                                                   simple_anno_size = unit(2, 'mm'), 
                                                   col = list(Group = c('control' = '#00DAE0', 'tumor' = '#FF9289')),  #set the sample groups color
                                                   show_annotation_name = FALSE), 
                column_names_gp = gpar(fontsize = 10), 
                row_names_gp = gpar(fontsize = 10),
                cluster_rows = F,
                cluster_columns = F)
print(heat)
pdf("LUAD_m6aheatmap.pdf",width=20,height=15)
dev.off
##when "dev.off()" doesn't work use "while (!is.null(dev.list()))  dev.off()"


#16.relationship analysis
fr <- m6agenes
fr <- as.data.frame(t(fr))
ggscatterstats(data = fr, 
               y =METTL14, 
               x = METTL3,
               centrality.para = "mean",
               margins = "both",
               xfill = "red", 
               yfill = "blue", 
               marginal.type = "histogram",
               title = "Relationship")

corr.result<-cor(fr,method = 'pearson')
corr.p<-ggcorrplot::cor_pmat(fr)
ggcorrplot(
  corr = corr.result,
  type = 'full',
  p.mat = corr.p,#P-Value
  sig.level = 0.05 #show the p-value >0.05
)

ggcorrplot(corr.result,
           method = "circle",
           hc.order = T,
           hc.method = "ward.D",
           outline.color = "white",
           ggtheme = theme_bw(),
           type = "upper",
           colors = c("#6D9EC1","white","#E46726"),
           lab = T,
           lab_size = 2,
           p.mat = corr.p,
           insig = "blank")
