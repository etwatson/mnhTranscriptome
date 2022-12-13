library(WGCNA)
library(dplyr)

# Estimate connectivity, scaled connectivity, cluster coefficient, 
# maximum adjacency ratio, network density, network centralization,
# and network heterogeneity

# identify modules of interest
moi <- hgt.tab %>% filter(module.size > 50, class == '0/+', p.value < 0.05) %>% select(mito, sex, module)
moi$ID <- paste(moi$sex, moi$mito, moi$module, sep = ".")

subnetwork.stats <- list()
for(i in 1:nrow(moi)){
	graph = net.out[paste(moi[i,2],moi[i,1], sep = ".")][[1]]  #subgraph
	vids = summary.out[paste(moi[i,2],moi[i,1], sep = ".")][[1]]$modules[moi[i,3]][[1]] #module
	sub.PFN <- induced.subgraph(graph, vids)
	sub.adj <- get.adjacency(sub.PFN,sparse=FALSE, attr = 'weight')
	sub.fnc <- fundamentalNetworkConcepts(sub.adj)
  	subnetwork.stats[[i]] <- sub.fnc
  	}
names(subnetwork.stats) <- moi$ID
