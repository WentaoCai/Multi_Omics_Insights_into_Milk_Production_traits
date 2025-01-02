args <- commandArgs(trailingOnly = TRUE)
phe_file = args[1]
n1=paste0("V",args[2])
pred_file = args[3]
phe=read.csv(phe_file,header = F)
predict=read.table(pred_file,header=F)
colnames(predict)=c("ID","value")
mergefile=merge(phe,predict, by.x = "V1", by.y = "ID")
Result=cor(mergefile[[n1]],mergefile[,"value"])
cat(Result,"\n")
