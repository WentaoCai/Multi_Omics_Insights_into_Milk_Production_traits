data=read.table("Gene_expression_6642samples_FPKM.txt",header = T)
row.names(data)=data$ID
data=data[,-1]
covariate=read.table("expression_covariate.txt",header=T,sep="\t")
covariate1=covariate[match(colnames(data), covariate$ID),]
data1=data[ ,colnames(data) %in% covariate1$ID]

###loged and scaled RPKM
Gene_Expression_T_scaled <- log(data1+0.25)
Gene_Expression_T_scaled <- apply(Gene_Expression_T_scaled,1,scale)
rownames(Gene_Expression_T_scaled) <- colnames(data1)
Gene_Expression_T_scaled=Gene_Expression_T_scaled[, colSums(is.na(Gene_Expression_T_scaled)) != nrow(Gene_Expression_T_scaled)]
head(Gene_Expression_T_scaled[,1:10])
################################################################################
######################detecting tissue-specific genes using t statistics########
#step 1. group tissues into visual classes, on any subsets
All_samples <- covariate1$Tissue_name
All_samples_e <- covariate1$Tissue_categories
Samples_sets <- unstack(All_samples,All_samples~as.factor(All_samples_e))
str(Samples_sets)

#step 2. compute the t-statisitcs while correcting for study, age ans sex covariables 
Category <- names(Samples_sets)
tissues_unique <- lapply(Samples_sets,unique)
str(tissues_unique)


##before run, you should change the output dir
for(i in 1:length(Samples_sets)){
  Tissue_class <- Category[i]
  No_tissues_category <- length(tissues_unique[[i]])
  ###only one tissue in the category 
  if(No_tissues_category==1){
    
    X <- Gene_Expression_T_scaled
    Sample_class <- ifelse(covariate1$Tissue_categories%in%Tissue_class,1,-1)
    
    Gene_number <- dim(X)[2]
    Myres <- array(NA,dim = c(Gene_number,2))
    
    ##for loop to compute the t-statistics for each gene in a tissue
    for(i in 1:Gene_number){
      Y <- as.numeric(X[,i])
      myfit <- lm(Y~Sample_class+factor(covariate1$Bioproject)+factor(covariate1$Breed_class)+factor(covariate$Sex))
      a <- summary(myfit)
      t <- coef(a)["Sample_class",c("t value","Pr(>|t|)")]
      Myres[i,c(1,2)] <- t
    }
    row.names(Myres) <- colnames(X)
    write.table(Myres,file = paste(Tissue_class,"_t.txt",sep = ""),append = F,quote = F,row.names = T,col.names = F)
  }
  ##when multiple tissues in the category
  if(No_tissues_category>1){
    Tissue_analyzed <- tissues_unique[[i]]
    for(j in 1:length(Tissue_analyzed)){
      target_tissue <- Tissue_analyzed[j]
      Similar_Tissues <- Tissue_analyzed[!Tissue_analyzed%in%target_tissue]
      X <-Gene_Expression_T_scaled
      index <- which(covariate1$Tissue_name%in%Similar_Tissues)
      X <- X[-index,]
      covariate_rm_similar <- covariate1[-index,]
      
      Sample_class <- ifelse(covariate_rm_similar$Tissue_name%in%target_tissue,1,-1)
      Gene_number <- dim(X)[2]
      Myres <- array(NA,dim = c(Gene_number,2))
      ##for loop to compute the t-statistics for each gene in a tissue
      for(i in 1:Gene_number){
        Y <- as.numeric(X[,i])
        myfit <- lm(Y~Sample_class+factor(covariate_rm_similar$Bioproject)+factor(covariate_rm_similar$Breed_class)+factor(covariate_rm_similar$Sex))
        a <- summary(myfit)
        t <- coef(a)["Sample_class",c("t value","Pr(>|t|)")]
        Myres[i,c(1,2)] <- t
      }
      row.names(Myres) <- colnames(X)
      write.table(Myres,file = paste(target_tissue,"_t.txt",sep = ""),append = F,quote = F,row.names = T,col.names = F)
    }
  }
}

