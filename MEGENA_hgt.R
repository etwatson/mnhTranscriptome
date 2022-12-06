
rm(list = ls()) # rm R working space

setwd("/media/mcclintock/etw/cybrid_rna/")
library(MEGENA)
library(dplyr)


# input parameters
n.cores <- 18; # number of cores/threads to call for PCP
doPar <-TRUE; # do we want to parallelize?
method = "spearman" # method for correlation. either pearson or spearman. 
FDR.cutoff = 0.05 # FDR threshold to define significant correlations upon shuffling samples. 
module.pval = 0.05 # module significance p-value. Recommended is 0.05. 
hub.pval = 0.05 # connectivity significance p-value based random tetrahedral networks
cor.perm = 10; # number of permutations for calculating FDRs for all correlation pairs. 
hub.perm = 1000; # number of permutations for calculating connectivity significance p-value. 

# annotation to be done on the downstream
annot.table=NULL
id.col = 1
symbol.col= 2
###########

#data(Sample_Expression)
# load("f1dc.Rdata")
# load("m1dc.Rdata")
# load("f2dc.Rdata")
# load("m2dc.Rdata")
# load("f3dc.Rdata")
# load("m3dc.Rdata")

rho_l <- list(M1dc,M2dc,M3dc,F1dc,F2dc,F3dc)
names(rho_l) <- c("Males.BR","Males.LB","Males.FHL","Females.BR","Females.LB","Females.FHL")

summary.out <- list()
MEGENA.out <- list()
net.out <- list()

for(i in 1:length(rho_l)){
  ijw <- rho_l[i][[1]]
  ijw <- ijw %>% filter(pValDiff_adj < 0.1, case_cor < 1,control_cor <1, case_pVal >0, control_pVal>0)
  ijw <- ijw %>% select(Gene1, Gene2, zScoreDiff);
  names(ijw) <- c("row","col","weight")
  ijw$weight <- abs(ijw$weight)
  run.par = doPar & (getDoParWorkers() == 1) 
  if(run.par){
    cl <- parallel::makeCluster(n.cores)
    registerDoParallel(cl)
    cat(paste("number of cores to use:",getDoParWorkers(),"\n",sep = ""))}
  pfn_res = MEGENA::calculate.PFN(ijw,num.cores = 6,doPar = TRUE)
  pfn_res$weight = (pfn_res$weight/max(pfn_res$weight)) * 0.999999999
  pfn_res<- pfn_res[!duplicated(pfn_res),]
  g = igraph::graph.data.frame(pfn_res, directed = FALSE)
  MEGENA.output <- do.MEGENA(g,
                             mod.pval = module.pval,hub.pval = hub.pval,remove.unsig = FALSE,
                             min.size = 10,max.size = vcount(g)/2,
                             doPar = doPar,num.cores = n.cores,n.perm = hub.perm,
                             save.output = FALSE)
  if(getDoParWorkers() > 1){
    env <- foreach:::.foreachGlobals
    rm(list=ls(name=env), pos=env)}
  summary.output <- MEGENA.ModuleSummary(MEGENA.output,
                                         mod.pvalue = module.pval,hub.pvalue = hub.pval,
                                         min.size = 10,max.size = vcount(g)/2,
                                         annot.table = annot.table,id.col = id.col,symbol.col = symbol.col,
                                         output.sig = TRUE)
  
  summary.out[[i]] <- summary.output;
  MEGENA.out[[i]] <- MEGENA.output;
  net.out[[i]] <- g;
}

names(summary.out) <- c("Males.BR","Males.LB","Males.FHL","Females.BR","Females.LB","Females.FHL")
names(MEGENA.out) <- c("Males.BR","Males.LB","Males.FHL","Females.BR","Females.LB","Females.FHL")
names(net.out) <- c("Males.BR","Males.LB","Males.FHL","Females.BR","Females.LB","Females.FHL")


#### Hypergeometric test

dc.cls <- c("0/+","+/0","+/-","-/+") #define differential correlation classes of interest
mat.l <- list()
hgt.l <- list()
listNames <- vector()
for (i in 1:6){ # loop through experimental groupings
  n_mod <- nrow(summary.out[i][[1]]$module.table) # number of modules
  N = length(unique(c(rho_l[i][[1]]$Gene1,rho_l[i][[1]]$Gene2))) # total number of genes
  sex = strsplit(names(rho_l[i]),split = "\\.")[[1]][1]
  pop = strsplit(names(rho_l[i]),split = "\\.")[[1]][2]

    for (j in 1:n_mod){                             # loop through modules
    mod = names(summary.out[i][[1]]$modules[j]) # get module ID
    modGen <- summary.out[i][[1]]$modules[j][[1]] # Get module genes IDs
    n = length(summary.out[i][[1]]$modules[j][[1]]) # Get number of genes in module

      for(v in 1:4){ #loop through classes
      #get number of DC genes in class
      K = length(unique(c(rho_l[i][[1]] %>% filter(pValDiff_adj < 0.1, case_cor < 1,control_cor <1, case_pVal >0, control_pVal>0, Classes == dc.cls[v]) %>% pull(Gene1),
                          rho_l[i][[1]] %>% filter(pValDiff_adj < 0.1, case_cor < 1,control_cor <1, case_pVal >0, control_pVal>0, Classes == dc.cls[v]) %>% pull(Gene2))))
      # Get DC gene IDs for given Mito x Sex groupings and class of DC
      K.id <- unique(c(rho_l[i][[1]] %>% filter(pValDiff_adj < 0.1, case_cor < 1,control_cor <1, case_pVal >0, control_pVal>0, Classes == dc.cls[v]) %>% pull(Gene1),
                       rho_l[i][[1]] %>% filter(pValDiff_adj < 0.1, case_cor < 1,control_cor <1, case_pVal >0, control_pVal>0, Classes == dc.cls[v]) %>% pull(Gene2)))
      k = length(K.id[K.id %in% modGen]) # Get number of DC genes in class
      print(paste("Now running hypergeometric test for ",pop," population ",sex,", class = ",dc.cls[v],", Module = ",mod,sep = ""))
      mat.l[[length(mat.l) + 1]] <- matrix(c(k,n-k,K-k,N+k-n-K), ncol =2, nrow =2,dimnames = list(c("Genes in class","Genes not in class"),c("In module","Not in module")))
      hgt.l[[length(hgt.l) + 1]] <-fisher.test(matrix(c(k,n-k,K-k,N+k-n-K), ncol =2, nrow =2), alternative = "greater")
      w <- paste(sex,pop,mod,dc.cls[v], sep = ".")
      listNames <- c(listNames,w)}
  }
}  
names(mat.l) <- listNames; names(hgt.l) <- listNames;    


hgt.tab <- data.frame()
for(i in 1:length(listNames)){
  mod.n <- sum(mat.l[i][[1]][,1])
  p <- hgt.l[i][[1]]$p.value
  or <- hgt.l[i][[1]]$estimate[[1]]
  s <- strsplit(listNames[i],"\\.")[[1]][1]
  mt <- strsplit(listNames[i],"\\.")[[1]][2]
  m <- strsplit(listNames[i],"\\.")[[1]][3]
  c <- strsplit(listNames[i],"\\.")[[1]][4]
  hgt.tab <- rbind(hgt.tab,data.frame(mito=mt, sex=s, module=m, module.size=mod.n, class=c, odds.ratio=or, p.value=p))
}
save(hgt.tab,hgt.l,mat.l,summary.out,MEGENA.out,net.out, file = "DC_enrichment.RData")
