library(DGCA)
library(dplyr)
CyDat <- read.table("CyDat.txt", head = TRUE,stringsAsFactors = TRUE)
adjusted_counts <- read.table("adjusted_counts.v3.txt")
gene_list <- colnames(adjusted_counts)
length(gene_list)
samp.F1 <- list()
samp.M1 <- list()
samp.F2 <- list()
samp.M2 <- list()
samp.F3 <- list()
samp.M3 <- list()
samp.M1case_cor <- list()
samp.M2case_cor <- list()
samp.M3case_cor <- list()
samp.F1case_cor <- list()
samp.F2case_cor <- list()
samp.F3case_cor <- list()

mtDNA <- c("BR","LB","FHL")
  for(i in 1:100){
    for(j in 1:3){
    s<- sample(gene_list,size = 500,replace = FALSE)
    CyExpr.s <- as.matrix(adjusted_counts[,colnames(adjusted_counts) %in% s])
    control <- matrix(as.numeric(as.character(CyExpr.s[row.names(CyExpr.s) %in% CyDat[CyDat$Mito == "SD" & CyDat$Sex == "Male",]$ID,])),ncol=dim(CyExpr.s)[[2]],dimnames = list(CyDat[CyDat$Mito == "SD" & CyDat$Sex == "Male",]$ID,colnames(CyExpr.s)))
    case <- matrix(as.numeric(as.character(CyExpr.s[row.names(CyExpr.s) %in% CyDat[CyDat$Mito == mtDNA[j] & CyDat$Sex == "Male",]$ID,])),ncol=dim(CyExpr.s)[[2]],dimnames = list(CyDat[CyDat$Mito == mtDNA[j] & CyDat$Sex == "Male",]$ID,colnames(CyExpr.s)))
    x <- c(rep("control",dim(control)[[1]]),rep("case",dim(case)[[1]]));designMat <- model.matrix(~x +0);colnames(designMat) = c("case", "control");M <- rbind(control,case);M <- t(M);M <- log(M+1)
    if(mtDNA[j] == "BR"){
      tmp0 <-  ddcorAll(inputMat = M, design = designMat,compare = c("control", "case"),corrType = "spearman",adjust = "fdr",nPerms = 0);tmp0$run <-i ; tmp0$mtDNA <- mtDNA[j];tmp0$Sex <- "Males";samp.M1[[i]] <- tmp0 %>% filter(pValDiff < 0.05);samp.M1case_cor[[i]] <- tmp0 %>% filter(case_pVal < 0.05)}
    if(mtDNA[j] == "LB"){
      tmp0 <-  ddcorAll(inputMat = M, design = designMat,compare = c("control", "case"),corrType = "spearman",adjust = "fdr",nPerms = 0);tmp0$run <-i ; tmp0$mtDNA <- mtDNA[j];tmp0$Sex <- "Males";samp.M2[[i]] <- tmp0 %>% filter(pValDiff < 0.05);samp.M2case_cor[[i]] <- tmp0 %>% filter(case_pVal < 0.05)}
    if(mtDNA[j] == "FHL"){
      tmp0 <-  ddcorAll(inputMat = M, design = designMat,compare = c("control", "case"),corrType = "spearman",adjust = "fdr",nPerms = 0);tmp0$run <- i ; tmp0$mtDNA <- mtDNA[j];tmp0$Sex <- "Males";samp.M3[[i]] <- tmp0 %>% filter(pValDiff < 0.05);samp.M3case_cor[[i]] <- tmp0 %>% filter(case_pVal < 0.05)}
    s<- sample(gene_list,size = 500,replace = FALSE)
    CyExpr.s <- as.matrix(adjusted_counts[,colnames(adjusted_counts) %in% s])
    control <- matrix(as.numeric(as.character(CyExpr.s[row.names(CyExpr.s) %in% CyDat[CyDat$Mito == "SD" & CyDat$Sex == "Female",]$ID,])),ncol=dim(CyExpr.s)[[2]],dimnames = list(CyDat[CyDat$Mito == "SD" & CyDat$Sex == "Female",]$ID,colnames(CyExpr.s)))
    case <- matrix(as.numeric(as.character(CyExpr.s[row.names(CyExpr.s) %in% CyDat[CyDat$Mito == mtDNA[j] & CyDat$Sex == "Female",]$ID,])),ncol=dim(CyExpr.s)[[2]],dimnames = list(CyDat[CyDat$Mito == mtDNA[j] & CyDat$Sex == "Female",]$ID,colnames(CyExpr.s)))
    x <- c(rep("control",dim(control)[[1]]),rep("case",dim(case)[[1]]));designMat <- model.matrix(~x +0);colnames(designMat) = c("case", "control");M <- rbind(control,case);M <- t(M);M <- log(M+1)
    if(mtDNA[j] == "BR"){
      tmp0 <-  ddcorAll(inputMat = M, design = designMat,compare = c("control", "case"),corrType = "spearman",adjust = "fdr",nPerms = 0);tmp0$run <- i; tmp0$mtDNA <- mtDNA[j];tmp0$Sex <- "Females";samp.F1[[i]] <- tmp0 %>% filter(pValDiff < 0.05);samp.F1case_cor[[i]] <- tmp0 %>% filter(case_pVal < 0.05)}
    if(mtDNA[j] == "LB"){
      tmp0 <-  ddcorAll(inputMat = M, design = designMat,compare = c("control", "case"),corrType = "spearman",adjust = "fdr",nPerms = 0);tmp0$run <- i; tmp0$mtDNA <- mtDNA[j];tmp0$Sex <- "Females";samp.F2[[i]] <- tmp0 %>% filter(pValDiff < 0.05);samp.F2case_cor[[i]] <- tmp0 %>% filter(case_pVal < 0.05)}
    if(mtDNA[j] == "FHL"){
      tmp0 <-  ddcorAll(inputMat = M, design = designMat,compare = c("control", "case"),corrType = "spearman",adjust = "fdr",nPerms = 0);tmp0$run <- i; tmp0$mtDNA <- mtDNA[j];tmp0$Sex <- "Females";samp.F3[[i]] <- tmp0 %>% filter(pValDiff < 0.05);samp.F3case_cor[[i]] <- tmp0 %>% filter(case_pVal < 0.05)}    
    }
  }      
