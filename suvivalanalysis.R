TCGAbiolinks::getGDCprojects()$project_id
clin.LUAD <- GDCquery_clinic("TCGA-LUAD", "clinical",save.csv = TRUE)
#write.csv(clin.LUAD,file = "LUAD_clinical.csv")

library(survminer)
TCGAanalyze_survival(clin.LUAD,
                     clusterCol="gender",
                     risk.table = FALSE,
                     xlim = c(100,1000),
                     ylim = c(0.4,1),
                     conf.int = FALSE,
                     pvalue = TRUE,
                     color = c("Dark2"))

