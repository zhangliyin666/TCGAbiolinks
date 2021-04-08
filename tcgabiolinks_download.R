library(BiocManager)
BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)
TCGAbiolinks::getGDCprojects()$project_id
cancer_type = "TCGA-LUAD"
clincal <- GDCquery_clinic(project = cancer_type,type = "clinical")
clinical <- clincal
View(clinical)

#下载RNA-seq的数据
library(dplyr)
library(DT)
library(SummarizedExperiment)
data_type <- "Gene Expression Quantification"
data_category <- "Transcriptome Profiling"
workflow_type <- "HTSeq - Counts"
query_TranscriptomeCounts <- GDCquery(project = cancer_type,
                                      data.category = data_category,
                                      data.type = data_type,
                                      workflow.type = workflow_type)
GDCdownload(query_TranscriptomeCounts,method = "api")
expdat <- GDCprepare(query = query_TranscriptomeCounts)
count_matix = assay(expdat)
View(count_matix)
write.csv(count_matix,file = "TCGAbiolinks_LUAD_counts.csv")

#下载miRNA数据
TCGAbiolinks:::getProjectSummary("TCGA-LUAD")
query_mi <- GDCquery(project = cancer_type,
                    data.category = "Transcriptome Profiling",
                    data.type = "miRNA Expression Quantification",
                    workflow.type = "BCGSC miRNA Profiling")
GDCdownload(query_mi,method = "api",files.per.chunk = 50)
expdat <- GDCprepare(query = query_mi)                    
#count_matrix = assay(expdat)
#会产生一个报错，估计是因为miRNA格式和mRNA不太一样
#Error in (function (classes, fdef, mtable)  : 
#unable to find an inherited method for function ‘assay’ for signature ‘"data.frame", "missing"’

###手动调整格式
GDCdownload(query_mi)
expdat <- GDCprepare(query = query_mi)
write.csv(expdat,file = paste(cancer_type,"miRNAs.csv",sep = "-"))
row.names(expdat) <- as.character(expdat[,1])
expdat <- expdat[,-1]
col_name<-unlist(lapply(colnames(expdat), FUN = function(x) {return(strsplit(x, split = "TCGA",fixed = T)[[1]][2])}))
col_name<-col_name[!duplicated(col_name)]
rpkm_names<-paste("reads_per_million_miRNA_mapped_TCGA",col_name,sep = "")
count_names<-paste("read_count_TCGA",col_name,sep = "")
write.csv(expdat[,rpkm_names],file = paste(cancer_type,"miRNAs_RPKM.csv",sep = "-"))
write.csv(expdat[,count_names],file = paste(cancer_type,"miRNAs_counts.csv",sep = "-"))
