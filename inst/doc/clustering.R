## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(mitoClone)
P1 <- mutationCallsFromMatrix(as.matrix(M_P1), as.matrix(N_P1))
P2 <- mutationCallsFromMatrix(as.matrix(M_P2), as.matrix(N_P2))

## ----select-------------------------------------------------------------------
P2@cluster["X9010GC"] <- F
P2@cluster["X2392TC"] <- F

## ----runPhiscs----------------------------------------------------------------
P1 <- muta_cluster(P1, cores=4, tempfolder = paste0(getwd(),"/P342debug"))
P2 <- muta_cluster(P2, cores=4, tempfolder = paste0(getwd(),"/P101debug"))

## ----plotTree-----------------------------------------------------------------
plotTree(P1, file = "P1.ps")
plotTree(P2, file = "P2.ps")

## ----clusterClonesP1, fig.width=8,fig.height=6--------------------------------
P1 <- clusterMetaclones(P1)

## ----clusterClonesP2, fig.width=8,fig.height=6--------------------------------
P2 <- clusterMetaclones(P2)

## ----plotClonesP1, fig.width=8,fig.height=6-----------------------------------
plotClones(P1)

## ----plotClonesP2, fig.width=8,fig.height=6-----------------------------------
plotClones(P2)

