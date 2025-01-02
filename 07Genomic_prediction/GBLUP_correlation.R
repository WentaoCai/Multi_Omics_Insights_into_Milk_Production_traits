args <- commandArgs(trailingOnly = TRUE)
phe_file = args[1]
n1=paste0("V",args[2])
pred_file = args[3]
phe=read.csv(phe_file,header = F)
predict=read.table(pred_file,header=T)
mergefile=merge(phe,predict, by.x = "V1", by.y = "IID")
Result=cor(mergefile[[n1]],mergefile[,ncol(mergefile)])
cat(Result,"\n")