library(VennDiagram)
library(edgeR)

# Remove poorly sequenced samples identified in MDS plot
CyDat <- CyDat[!CyDat$ID %in% c("SE6480_SA72437","SE6480_SA72439","SE6484_SA72395"),] 
adjusted_counts <- adjusted_counts[!rownames(adjusted_counts) %in% c("SE6480_SA72437","SE6480_SA72439","SE6484_SA72395"),]


mtDNA <- c("SD","BR","LB","FHL")
sbgList <- list()
sbg.idList <- list()
Tsbg.idList<- list()
slg.idList<- list()
adjusted_counts <- as.matrix(adjusted_counts)
for(i in mtDNA){
  dat <- CyDat %>% filter(Mito == i)
  m1 <- matrix(as.numeric(as.character(adjusted_counts[row.names(adjusted_counts) %in% CyDat[CyDat$Mito == i,]$ID,])),ncol=dim(adjusted_counts)[[2]],dimnames = list(CyDat[CyDat$Mito == i,]$ID,colnames(adjusted_counts)));
  m1 <- t(m1)
  dat$Sex.num <- ifelse(dat$Sex.num == 1, 0, 1)
  designMat <- model.matrix(~dat$Sex.num) # 0= Female, 1= Male
  rownames(designMat) <- colnames(m1)
  keep.exprs <- filterByExpr(m1,design = designMat, min.count=10) # remove genes with low expression
  m1 <- m1[keep.exprs,]
  dgList <- DGEList(counts=m1, genes=rownames(m1))
  dgList <- calcNormFactors(dgList)
  dgList <- estimateGLMCommonDisp(dgList, design=designMat)
  dgList <- estimateGLMTrendedDisp(dgList, design=designMat)
  dgList <- estimateGLMTagwiseDisp(dgList, design=designMat)
    pdf(paste(i,"_MDS.pdf", sep = ""))
    plotMDS(dgList, main=i)
    dev.off()
  fit <- glmQLFit(dgList, designMat)
  lrt <- glmQLFTest(fit, coef=2)
  sbgList[[i]] <- lrt
  min_cutoff <- quantile(sbgList[i][[1]]$table %>% filter(PValue < 0.05) %>% pull(logFC), 0.99)
  tmp  <- sbgList[i][[1]]$table %>% filter(PValue < 0.05, abs(logFC) > min_cutoff) 
  slg.idList [[i]] <- unique(row.names(tmp))
  sbg.idList [[i]] <- decideTestsDGE(sbgList[i][[1]], p=0.01,adjust.method = "bonf")
  sbg.idList [[i]] <- rownames(sbgList[i][[1]])[as.logical(sbg.idList[[i]])]
  Tsbg.idList [[i]] <- decideTestsDGE(sbgList[i][[1]], p=1.0e-12,adjust.method = "bonf")
  Tsbg.idList [[i]] <- rownames(sbgList[i][[1]])[as.logical(Tsbg.idList[[i]])]}

# Sex-biased genes
SDsbg.id <- sbg.idList["SD"][[1]];BRsbg.id <- sbg.idList["BR"][[1]];LBsbg.id <- sbg.idList["LB"][[1]];FHLsbg.id <- sbg.idList["FHL"][[1]];
# Top sex-biased genes  
SDtsbg.id <- Tsbg.idList["SD"][[1]];BRtsbg.id <- Tsbg.idList["BR"][[1]];LBtsbg.id <- Tsbg.idList["LB"][[1]];FHLtsbg.id <- Tsbg.idList["FHL"][[1]];
# Sex-limited genes
SDslg.id <- slg.idList["SD"][[1]];BRslg.id <- slg.idList["BR"][[1]];LBslg.id <- slg.idList["LB"][[1]];FHLslg.id <- slg.idList["FHL"][[1]];

# Investigate gene groupings with Venn diagrams
sbg.overlap<- VennDiagram::calculate.overlap(list(SD=SDsbg.id, BR=BRsbg.id, LB=LBsbg.id, FHL= FHLsbg.id))
hsSBG <- c("a1","a2","a3","a7","a8","a13","a14");BR.sbg <- "a14";SD.sbg <- "a9";LB.sbg <- "a1";FHL.sbg <- "a3";cSBG <- "a6";
BR.lSBG.id <- as.vector(unlist(sbg.overlap[BR.sbg]));SD.lSBG.id <- as.vector(unlist(sbg.overlap[SD.sbg]));LB.lSBG.id <- as.vector(unlist(sbg.overlap[LB.sbg]));FHL.lSBG.id <- as.vector(unlist(sbg.overlap[FHL.sbg]));hsSBG.id <- as.vector(unlist(sbg.overlap[hsSBG]));cSBG.id <- as.vector(unlist(sbg.overlap[cSBG]));
tsbg.overlap<- VennDiagram::calculate.overlap(list(SD=SDtsbg.id, BR=BRtsbg.id, LB=LBtsbg.id, FHL= FHLtsbg.id))
slg.overlap<- VennDiagram::calculate.overlap(list(SD=SDslg.id, BR=BRslg.id, LB=LBslg.id, FHL= FHLslg.id))

