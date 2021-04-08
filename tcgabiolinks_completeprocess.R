###1.数据下载
library(TCGAbiolinks)
library(dplyr)
library(DT)
library(SummarizedExperiment)
data_type <- "Gene Expression Quantification"
data_category <- "Transcriptome Profiling"
workflow_type <- "HTSeq - Counts"
cancer_type = "TCGA-BRCA"
#legacy默认为FALSE，从hg38下载数据
query_TranscriptomeCounts <- GDCquery(project = cancer_type,
                                      data.category = data_category,
                                      data.type = data_type,
                                      workflow.type = workflow_type)
#getResults(query, rows, cols)根据指定行名或列名从query——transcriptome中获取结果,此处用来获得样本的barcode
samplesDown <- getResults(query_TranscriptomeCounts,cols=c("cases"))#共594个
#从samplesDown中筛选出TP（实体肿瘤）样本的barcodes
#TP：Primary Solid Tumor
#NT：Solid Tissue Normal
dataSmTP <- TCGAquery_SampleTypes(barcode = samplesDown,
                                  typesample = "TP")#533个
dataSmNT <- TCGAquery_SampleTypes(barcode = samplesDown,
                                  typesample = "NT")#59个
#按照样本分组重新筛选符合条件的barcodes
queryDown <- GDCquery(project = cancer_type, 
                      data.category = data_category,
                      data.type = data_type, 
                      workflow.type = workflow_type, 
                      barcode = c(dataSmTP, dataSmNT))
#barcode：根据输入的数据进行过滤
GDCdownload(queryDown,
            method = "api",
            directory = "GDCdata",
            files.per.chunk = 6)#共592个文件
#method：api或者client，api速度快但可能出现下载中断
#files.per.chunk：使用api方法时一次只下载n个文件，将一个大文件拆分为n个小文件进行下载，减少下载出错
###2.数据处理
#GDCprepare()将前面GDCquery()的结果准备成R语言可处理的SE文件
dataPrep1 <- GDCprepare(query = queryDown, save = TRUE, save.filename =
                          "dataprep1_BRCA_case.rda")

###3.数据预处理
#去除异常值
dataPrep2 <- TCGAanalyze_Preprocessing(object = dataPrep1,
                                       cor.cut = 0.6,
                                       datatype = "HTSeq - Counts")#592个
write.csv(dataPrep2,file = "dataPrep2_BRCA.csv",quote = FALSE)
###4.筛选肿瘤纯度大于60%的barcodes
#使用estimate, absolute, lump, ihc, cpe来衡量
#cpe是派生的共识度量，是将所有方法的标准含量归一化后的均值纯度水平，以使它们具有相等的均值和标准差
purityDATA <- TCGAtumor_purity(colnames(dataPrep1), 0, 0, 0, 0, 0.6)#381个
#purityDATA$pure_barcodes <- dataSmTP#purity533个
#purityDATA$filtered <- dataSmNT#normal59个
#filterde为被过滤的数据（正常组织的barcodes），pure_barcodes为我们要的肿瘤数据
Purity.BRCA <- purityDATA$pure_barcodes#322个
normal.BRCA <- purityDATA$filtered#59个
###5.合并肿瘤和正常组织的表达矩阵
puried_data <-dataPrep2[,c(Purity.BRCA,normal.BRCA)]#共381个
###6.进行基因注释
rownames(puried_data)<-rowData(dataPrep1)$external_gene_name
write.csv(puried_data,file = "puried.BRCA.csv",quote = FALSE)
###7.表达矩阵标准化和过滤，得到用于差异分析的矩阵
BiocManager::install("EDASeq")
library(EDASeq)
#标准化
dataNorm <- TCGAanalyze_Normalization(tabDF = puried_data,
                                      geneInfo = geneInfo,
                                      method = "gcContent")
#过滤低表达基因
dataFilt_BRCA_final <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                             method = "quantile", 
                                             qnt.cut =  0.25)
a <- c(rep("tumor",340),rep("normal",50))
a <- t(as.data.frame(a))
colnames(a)<-colnames(dataFilt_BRCA_final)
dataFilt_BRCA_final <- rbind(dataFilt_BRCA_final,a)
write.csv(dataFilt_BRCA_final,file = "TCGA_LIHC_final.csv",quote = FALSE)
###8.差异分析
#dataFilt_LIHC_final <- read.csv("TCGA_LUAD_final.csv", header = T,check.names = FALSE）
#rownames(dataFilt_LIHC_final) <- dataFilt_LIHC_final[,1]
#dataFilt_LIHC_final <- dataFilt_LIHC_final[,-1]
#肿瘤分组
mat1 <- dataFilt_LIHC_final[,1:931]
mat1 <- log(mat1+1)
#正常分组
mat2 <- dataFilt_LIHC_final[,932:1044]
mat2 <- log(mat2+1)
#差异分析
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
write.csv(Data_DEGs,file = "LUAD_DEGs.csv",quote = FALSE)
###9.富集分析
#筛选表达差异的基因，设置logfc
Data_DEGs_high_expr <- Data_DEGs[Data_DEGs$logFC >=1,]
Genelist <- rownames(Data_DEGs_high_expr)
ansEA <- TCGAanalyze_EAcomplete(TFname="DEA genes Tumor Vs Normal",
                                Genelist)
#可视化
TCGAvisualize_EAbarplot(tf  = rownames(ansEA$ResBP),
                        GOBPTab = ansEA$ResBP,
                        GOCCTab = ansEA$ResCC,
                        GOMFTab = ansEA$ResMF,
                        PathTab = ansEA$ResPat,
                        nRGTab = Genelist,
                        nBar = 20, #显示条形图的数量
                        filename = "TCGAvisualize_EAbarplot_Output.pdf")
###10.聚类分析
BiocManager::install("ConsensusClusterPlus")
library(ConsensusClusterPlus)
res.hc <- TCGAanalyze_Clustering(Data_DEGs, 
                                 method = "hclust", 
                                 methodHC = "ward.D2")
plot(res.hc)
###11.clusterprofiler分析
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
                 
                 ont        = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05)
ego_mf<-enrichGO(gene       = gene.df$ENTREZID,
                 OrgDb      = org.Hs.eg.db,
                 
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
#当时直接运行了结果报错，
#Error in names(foldchanges) = gene.df$ENTREZID : 
#'names'属性的长度[14370]必需和矢量的长度[12963]一样
#duplicated(gene.df$SYMBOL)然后查重一下
#果然有重复的
#gene.df <- gene.df %>% distinct(SYMBOL,.keep_all = TRUE)
#用于data.frame去重
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
###12.生存分析
library(survminer)
library(survival)
#下载临床数据
clin.LUAD <- GDCquery_clinic("TCGA-LUAD", "clinical")
#性别对生存曲线的影响
TCGAanalyze_survival(clin.LUAD,
                     clusterCol="gender",
                     risk.table = FALSE,
                     xlim = c(100,1000),
                     ylim = c(0.4,1),
                     conf.int = FALSE,
                     pvalue = TRUE,
                     color = c("Dark2"))
##单个基因表达对生存曲线的影响
#提取肿瘤样本中A4GALT的表达
#以A4GALT为例
dataFilt_LUAD_final <- read.csv("TCGA_LUAD_final.csv",row.names = 1,check.names=FALSE)
samplesTP <- TCGAquery_SampleTypes(colnames(dataFilt_LUAD_final), typesample = c("TP"))
METTL14 <- dataFilt_LUAD_final[c("METTL14"),samplesTP]
#修改样本名称
#样本名像这样TCGA-91-6840-01A-11R-1949-07
names(METTL14) <- sapply(strsplit(names(METTL14),'-'),function(x) paste0(x[1:3],collapse="-"))
METTL3 <- t(METTL3)
#合并基因和临床数据
clin.LUAD$"METTL3" <- METTL3[match(clin.LUAD$submitter_id,rownames(METTL3)),]
#选择生存分析需要的数据
df<-subset(clin.LUAD,select =c(submitter_id,vital_status,days_to_death,days_to_last_follow_up,METTL3))
#去除基因表达缺失的值
df <- df[!is.na(df$METTL3),]
#根据基因表达情况分组
df$METTL3 <- as.numeric(df$METTL3)
df$exp <- ''
df[df$METTL3 >= mean(df$METTL3),]$exp <- "H"
df[df$METTL3 < mean(df$METTL3),]$exp <- "L"
#使用TCGAanalyze_survival函数分析
TCGAanalyze_survival(df,
                     legend = "METTL3",
                     clusterCol="exp",
                     risk.table = FALSE,
                     conf.int = FALSE,
                     pvalue = TRUE,
                     color = c("Dark2"))
###13.差异基因的热图和火山图
##加载预处理文件
TCGA_LUAD_data <- read.csv(file = "TCGA_LUAD_final.csv",
                           header = T,
                           row.names = 1,
                           check.names = FALSE )#保证列名不发生自动更正
#获取正常和肿瘤的barcodes
samplesNT <- TCGAquery_SampleTypes(colnames(TCGA_LUAD_data), typesample = c("NT"))
samplesTP <- TCGAquery_SampleTypes(colnames(TCGA_LUAD_data), typesample = c("TP"))
#配对正常和肿瘤的barcodes
paired <- intersect(substr(samplesNT, 1, 12),substr(samplesTP, 1, 12))
length(paired)
#用前12位数据作桥梁，建立肿瘤和正常样本的数据框
NT <- data.frame(NT1=substr(samplesNT,1,12),NT2=samplesNT)
TP <- data.frame(TP1 =substr(samplesTP,1,12),TP2=samplesTP)
TP_NT <- merge(TP,NT,by.x = "TP1",by.y = "NT1")
head(TP_NT,3)
#获取配对正常组织的barcodes
TP <- TP_NT$TP2
#获取配对肿瘤组织的barcodes
NT <- TP_NT$NT2
##数据下载和预处理
#由于之前已经保存了一些数据所以处理部分可以忽略
#queryDown <- GDCquery(project = "TCGA-LUAD",
                      #data.category = "Transcriptome Profiling",
                      #data.type = "Gene Expression Quantification",
                      #workflow.type = "HTSeq - Counts")
#GDCdownload(queryDown,
            #method = "api",
            #directory = "GDCdata",
            #files.per.chunk = 6)
#dataPrep1 <- GDCprepare(query = queryDown, save = F)
#直接读入之前保存好的数据
load("H:/tcgabiolinks/LUAD_case.rda")
dataPrep2 <- TCGAanalyze_Preprocessing(object = dataPrep1,
                                       cor.cut = 0.6,
                                       datatype = "HTSeq - Counts")
dataPrep1 <- data
purityDATA <- TCGAtumor_purity(colnames(dataPrep1), 0, 0, 0, 0, 0.6)
Purity.LUAD <- purityDATA$pure_barcodes#322
normal.LUAD <- purityDATA$filtered#59
puried_data <-dataPrep2[,c(Purity.LUAD,normal.LUAD)]
#写入结果
rownames(puried_data)<-rowData(dataPrep1)$external_gene_name
#进行文库大小和GC丰度标准化
dataNorm <- TCGAanalyze_Normalization(tabDF = puried_data,
                                      geneInfo = geneInfo,
                                      method = "gcContent")
#过滤低count基因
dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                             method = "quantile", 
                                             qnt.cut =  0.25)
#write.csv(dataFilt,file = "paired_TCGA_LUAD_final.csv",quote = FALSE)
##差异表达分析
#dataFilt_LUAD_final <- read.csv(file = "paired_TCGA_LUAD_final.csv",header = T,row.names = 1,check.names = FALSE)
#肿瘤vs正常
mat1 <- dataFilt_LUAD_final[,1:322]#肿瘤322个
mat2 <- dataFilt_LUAD_final[,323:381]#正常59个
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
#write.csv(DEG.LUAD.edgeR,file = "paired_DEG_by_edgeR.csv")
##增加不同分组条件下的gene平均表达量
#获取各自的barcodes
samplesNT <- TCGAquery_SampleTypes(colnames(dataFilt_LUAD_final), typesample = c("NT"))
samplesTP <- TCGAquery_SampleTypes(colnames(dataFilt_LUAD_final), typesample = c("TP"))
dataDEGsFilt <- DEG_LUAD_edgeR[abs(DEG_LUAD_edgeR$logFC) >= 1,]
str(dataDEGsFilt)
#获取各自的表达矩阵
dataTP <- dataFilt_LUAD_final[,samplesTP]
dataTN <- dataFilt_LUAD_final[,samplesNT]
#传入参数
dataDEGsFiltLevel <- TCGAanalyze_LevelTab(FC_FDR_table_mRNA = dataDEGsFilt,
                                          typeCond1 ="Normal",
                                          typeCond2 = "Tumor",
                                          TableCond1 = dataTN,
                                          TableCond2 = dataTP)
head(dataDEGsFiltLevel, 2)
##主成分分析
pca <-TCGAvisualize_PCA(dataFilt = dataFilt_LUAD_final, 
                        dataDEGsFiltLevel = dataDEGsFiltLevel,
                        ntopgenes = 100,
                        group1 = samplesNT, 
                        group2 = samplesTP)
##差异基因热图
#获取表达矩阵
datDEGs <- dataFilt_LUAD_final[match(rownames(DEG_LUAD_edgeR),rownames(dataFilt_LUAD_final)),]
#获取临床信息
query <- GDCquery(project = "TCGA-LUAD",
                  data.category = "Clinical",
                  file.type = "xml",
                  barcode = substr(colnames(datDEGs),1,12))
GDCdownload(query)  
clinical <- GDCprepare_clinic(query,"patient")
#根据表达矩阵样本的barcodes匹配临床信息
datDEGs_test_barcodes <- as.data.frame(substr(colnames(datDEGs),1,12),ncol=1)
colnames(datDEGs_test_barcodes) <-"LUAD_patient_barcode"
m <- clinical[match(datDEGs_test_barcodes[,1], clinical[ , 1]),]
str(m)
#查看重复信息
table(duplicated(m))
#原始数据作图
pheatmap::pheatmap(datDEGs,scale = "row",show_rownames = F,show_colnames = F)
#增加metadata信息
col.mdat <- data.frame(Sex=m$gender,
                       status=m$vital_status,
                       group=c(rep("tumor",322),rep("normal",59)))
#保证列注释信息的行名与样本名（对应列）一致
rownames(col.mdat) <- colnames(datDEGs) 
#设置图例范围
bk <- c(seq(-1,6,by=0.01))
#作图
pheatmap::pheatmap(datDEGs,scale = "row",show_rownames = F,show_colnames = F,
                   annotation_col = col.mdat,
                   border_color=NA,
                   main = "Heatmap by pheatmap(edgeR)",
                   filename = "Heatmap_by_pheatmap.pdf",
                   color =c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),
                            colorRampPalette(colors = c("white","red"))(length(bk)/2)),#设置图例的颜色,
                   legend_breaks=seq(-1,6,2),
                   breaks=bk )
##差异分析火山图
#突出显示logfc>=8的基因名
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
###14.肿瘤浸润分析
#接puried_data
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


library(preprocessCore)
source("cibersort.R")
LM22.file <- "LM22.txt"
exp.file <- "exprset.txt"
TME.results = CIBERSORT(LM22.file, exp.file, perm = 1000, QN = TRUE)
write.csv(TME.results,file = "TME_cibersort.csv")


condition <- factor(c(rep("tumor",931),rep("control",113)), levels = c("tumor","control"))
condition
table(condition)

exprset <- read.csv("TCGA_BRCA_final.csv",header = T,row.names = 1)
exprset <- log2(exprset+1)

m6agenes <- exprset[c("METTL3","METTL14","WTAP","RBM15","ZC3H13",
                      "YTHDF1","YTHDF2","YTHDF3","YTHDC1","YTHDC2","IGF2BP1","IGF2BP2","HNRNPC",
                      "FTO","ALKBH5"),]
m6agenes <- as.data.frame(t(m6agenes))
m6agenes$condition <- c(rep("tumor",931),rep("control",113))

a <- as.data.frame(m6agenes$ALKBH5)
b <- as.data.frame(m6agenes$condition)
c <- as.data.frame(rep("ALKBH5",1044))
d <- cbind(a,c,b)
colnames(d) <- c("Expression","gene","group")
#dat <- d
dat <- rbind(dat,d)
dat$gene = factor(dat$gene, levels=c("METTL3","METTL14","WTAP","RBM15","ZC3H13",
                                     "YTHDF1","YTHDF2","YTHDF3","YTHDC1","YTHDC2","IGF2BP1","IGF2BP2","HNRNPC",
                                     "FTO","ALKBH5"))
write.csv(dat,file = "dat_brca.csv")
compare_means( Expression ~ group, data = dat, group.by = "gene")
p <- ggboxplot(dat, x = "group", y = "Expression",
               color = "group", palette = "jco",
               add = "jitter",
               facet.by = "gene",
               short.panel.labs = FALSE)
p + stat_compare_means(label = "p.format")
pdf("LIHC.pdf",width=20,height=20)
dev.off
#4.当3步无效时使用while (!is.null(dev.list()))  dev.off()

condition <- factor(c(rep("tumor",322),rep("control",59)), levels = c("control","tumor"),ordered = F)
group <- data.frame(group=condition)
rownames(group)=colnames(m6agenes)
samples <- condition
m6agenes <- as.matrix(m6agenes)
heat <- Heatmap(m6agenes, 
                col = colorRampPalette(c('navy', 'white', 'firebrick3'))(100), #定义热图由低值到高值的渐变颜色
                heatmap_legend_param = list(grid_height = unit(10,'mm')),  #图例高度设置
                show_row_names = TRUE,  #不展示基因名称
                top_annotation = HeatmapAnnotation(Group = samples, 
                                                   simple_anno_size = unit(2, 'mm'), 
                                                   col = list(Group = c('control' = '#00DAE0', 'tumor' = '#FF9289')),  #定义样本分组的颜色
                                                   show_annotation_name = FALSE), 
                column_names_gp = gpar(fontsize = 10), 
                row_names_gp = gpar(fontsize = 10),
                cluster_rows = F,
                cluster_columns = F)
print(heat)
pdf("cd10_m6aheatmap.pdf",width=20,height=15)
dev.off
#4.当3步无效时使用while (!is.null(dev.list()))  dev.off()


#相关性分析
fr <- m6agenes
fr <- log(fr+1)
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
  sig.level = 0.05 #P-Value大于0.05的在图中标记出来
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
