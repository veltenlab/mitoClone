## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----example,fig.width = 6, fig.height=4--------------------------------------
library(mitoClone)
library(pheatmap)

LudwigFig5.Counts <- readRDS(url("http://steinmetzlab.embl.de/mutaseq/fig5_mc_out.RDS"))

LudwigFig5 <- mutationCallsFromBlacklist(LudwigFig5.Counts,min.af=0.05, min.num.samples=5, universal.var.cells = 0.5 * length(LudwigFig5.Counts), binarize = 0.1)

LudwigFig5.Meta <- data.frame(row.names = rownames(LudwigFig5@N), Clone = gsub("_.*","",gsub("Donor1_","",rownames(LudwigFig5@N))))

clustered <- quick_cluster(LudwigFig5, binarize = T, drop_empty = T, clustering.method = "ward.D2", annotation_col = LudwigFig5.Meta,show_colnames=F,fontsize_row = 7)


## ----fig7orig,fig.width = 6, fig.height=4-------------------------------------
## this is the Fig7 A/B data it is only the single cells from the Tumor samples (see bottom arrow of Fig7A) - It looks very messay and the healthy comparison is bulk sequencing (let me know if you think I should follow up)
LudwigFig7.Counts <- readRDS(url("http://steinmetzlab.embl.de/mutaseq/fig7_nucleotide_counts_per_position.RDS"))


supervised.variants <- c("9000 T>C","4819 G>A","11127 G>A","7229 A>G","5921 G>A","15735 G>A","15044 G>A","14866 G>A","6735 G>A","12651 T>C","14805 G>A")
shifted.pos <- as.integer(gsub(" .+","",supervised.variants)) - 1
basechange <- gsub("\\d+ ","",supervised.variants)
supervised.variants.shifted <- paste(shifted.pos, basechange)
hand.selected <- pullcounts.vars(LudwigFig7.Counts, supervised.variants, shift=-1)
hand.selected <- mutationCallsFromMatrix(t(hand.selected$M), t(hand.selected$N), cluster = rep(T, length(supervised.variants)))

clustered.hand <- quick_cluster(hand.selected, binarize = T, drop_empty = T, clustering.method = "ward.D2",fontsize_row = 7, show_colnames=F, cutree_cols = 12, cluster_distance_col = "manhattan")
 
orig.anno <- data.frame(row.names = names(cutree(clustered.hand$tree_col, k=12)), Clone = as.factor(cutree(clustered.hand$tree_col, k=12)))


## ----fig7new,fig.width = 6, fig.height=4--------------------------------------
LudwigFig7 <- mutationCallsFromBlacklist(LudwigFig7.Counts,min.af=0.1, min.num.samples=3, universal.var.cells = 0.25 * length(LudwigFig7.Counts), binarize = 0.1)
LudwigFig7 <- muta_cluster(LudwigFig7, cores=10,tempfolder = paste0(getwd(),"/LudwigFig7"))

LudwigFig7@metadata$original <- factor(NA, levels = unique(orig.anno$Clone))
LudwigFig7@metadata[rownames(orig.anno), "original"] <- orig.anno$Clone
 
LudwigFig7 <- clusterMetaclones(LudwigFig7, plot = F)
 
plotClones(LudwigFig7)
 

