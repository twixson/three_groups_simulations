#############################
# competitors script
### run first: generate_data.R 
#############################
library(edgeR)
library(limma)
library(qvalue)
library(MASS)
library(dplyr)
library(DESeq2)
#############################

######
# GWAS 
######
logit_pval <- rep(NA,num_genes)
for(i in 1:num_genes){
  logit_pval[i] <- summary(glm(Y_GWAS ~ X_GWAS[,i],
                               family = binomial()))$coeff[8]}
post_probs_null_GWAS <- qvalue(logit_pval)$lfdr


######
# edgeR
######
groups_factor <- factor(PD_indicator)
x <- matrix(Y_RNA, nrow=num_genes)         
# store simulated count data as x
y <- DGEList(counts=x,group=groups_factor) 
# convert x, design into DGEList object for edgeR
y <- calcNormFactors(y)                    
# calculate normalization factors to scale raw library sizes
design <- model.matrix(~ PD_indicator)     
# create design matrix
y <- estimateDisp(y,design)                
# estimate dispersions for each 'tag' (gene)

fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef=2)
edger_pval <- lrt$table$PValue
post_probs_null_edgeR <- qvalue(lrt$table$PValue)$lfdr


######
# voom
######
v <- voom(x, design, plot=F)
vfit <- lmFit(v, design)
efit <- eBayes(vfit)
voom_pval <- efit$p.value[,2]
post_probs_null_voom <- qvalue(efit$p.value[,2])$lfdr  


###### 
# DEseq2
###### 
des_x <- x    # matrix of counts
# !!!!columns of the count matrix MUST match, 
# !!!!    in order, the rows of the column data !!!!
colnames_vector <- rep(NA, dim(des_x)[2])
j_control <- j_PD <- 0
for(i in 1:length(colnames_vector)){
  if(PD_indicator[i] == 0){
    j_control <- j_control + 1
    colnames_vector[i] <- paste0("control_", j_control)
  } else {
    j_PD <- j_PD + 1
    colnames_vector[i] <- paste0("PD_", j_PD)
  }
}
colnames(des_x) <- colnames_vector
des_coldata <- data.frame(condition=factor(PD_indicator))
rownames(des_coldata) <- colnames(des_x)
des_data <- DESeqDataSetFromMatrix(countData = des_x,
                                   colData = des_coldata,
                                   design = ~ condition)
des_fit <- DESeq(des_data)
#resultsNames(des_fit) # lists the coefficients
des_res <- results(des_fit, name="condition_1_vs_0")
des_res$pvalue[is.na(des_res$pvalue)] <- 1
post_probs_null_DEseq2 <- qvalue(des_res$pvalue)$lfdr


###### 
# combine results from GWAS and edgeR 
######
gwas_edger1 <- -2*(log(logit_pval)+log(edger_pval))
gwas_edger2 <- pchisq(gwas_edger1,lower.tail = F, df = 4)
post_probs_null_GWAS_edgeR <- qvalue(gwas_edger2)$lfdr 


###### 
# combine results from gwas and voom
###### 
gwas_voom1 <- -2*(log(logit_pval)+log(voom_pval))
gwas_voom2 <- pchisq(gwas_voom1,lower.tail = F, df = 4)
post_probs_null_GWAS_voom <- qvalue(gwas_voom2)$lfdr 


###### 
# combine results from gwas and deseq2
###### 
gwas_deseq1 <- -2*(log(logit_pval)+log(des_res$pvalue))
gwas_deseq2 <- pchisq(gwas_deseq1,lower.tail = F, df = 4)
post_probs_null_GWAS_DEseq2 <- qvalue(gwas_deseq2)$lfdr  

