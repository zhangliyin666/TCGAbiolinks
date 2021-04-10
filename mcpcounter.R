install_github("ebecht/MCPcounter",ref="master", subdir="Source",force = T)
exprset <- read.csv(file = "TCGA_LUAD_final.csv",header = T,row.names = 1)
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

