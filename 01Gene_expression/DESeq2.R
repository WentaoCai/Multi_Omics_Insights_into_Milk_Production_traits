library("DESeq2")
#input read count data
data <- read.table("103_mammary_gene.expression.count.txt",header=T)
countdata<-data
rownames(countdata)<-countdata$ID
countdata<-countdata[,-1]
#input covariate
coldata <- read.table("Covariate.txt",sep="\t",header=T)
head(coldata)
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design= ~Stage+Project)
colData(dds)
dds <- DESeq(dds)
res <- results(dds, contrast=c("Stage","Non-lactation","Lactation"))
write.csv(res,"DEGs.csv")