## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----positionList, eval=F-----------------------------------------------------
#  positions <- read.csv(url("http://steinmetzlab.embl.de/mutaseq/position_list.txt"), col.names = c("chr","pos"),header=F, stringsAsFactors = F)

## ----countTables, eval=F------------------------------------------------------
#  library(deepSNV)
#  require(VGAM)
#  
#  
#  bam <- list.files("bam", full.names = T)
#  
#  getCounts <- function(chr,pos, bamfilelist) {
#    if (chr == "chr23") chr <- "chrX"
#    if (chr == "chr24") chr <- "chrY"
#    out <- lapply(bamfilelist, bam2R,chr =chr,start=pos, stop=pos)
#    out <- do.call(rbind,out)[,c("A","C","G","T","N","-","INS","DEL")]
#    out
#  }
#  
#  countTables <- mapply(getCounts, positions$chr, positions$pos, MoreArgs = list(bamfilelist= bam), SIMPLIFY = F)

## ----simpleModel, eval=F------------------------------------------------------
#  likeli.simple <- function(params, data) {
#    params <- 1/(1+exp(-params))
#    -sum(mapply(dbetabinom, data[,1], apply(data,1,sum), MoreArgs = list(prob = params[1], rho = params[2], log=T)))
#  
#  }

## ----complexModel, eval=F-----------------------------------------------------
#  
#  likeli.clas <- function(params, data, clones) {
#    params <- 1/(1+exp(-params))
#    out <- 0
#    for ( i in 1:max(clones)) {
#    add <- -sum(mapply(dbetabinom, data[clones==i,1], apply(data[clones==i,],1,sum), MoreArgs = list(prob = params[i], rho = params[max(clones)+1], log=T)))
#    if(!is.na(add)) out <- out +add
#    }
#    out
#  
#  }

## ----getclone, eval=F---------------------------------------------------------
#  clone <- apply(P1@mainClone,1,which.max)

## ----fitModels, eval=F--------------------------------------------------------
#  
#  getAIC <- function(counts, clones) {
#  
#    test <- counts[,c("A","C","G","T","INS","DEL")]
#    colsums <- apply(test,2,sum)
#    if (max(colsums) > sum(colsums) - 100) return(NA) else {
#      test <- test[,rank(colsums) >= length(colsums)-1]
#      if(!is.matrix(test)) return(NA) else {
#        totalreads <- apply(test,1,sum)
#        ar <- apply(test,2,function(x) mean(x/totalreads,na.rm=T))[1]
#        ar.byclone <- sapply(unique(clones), function(cl) {
#          apply(test[clones == cl,],2,function(x) mean(x/totalreads[clones==cl],na.rm=T))
#        })[1,]
#        logit <- function(p) log(p/(1-p))
#  
#        #start with the max likelihood estimate of a binomal distribution to speed up the fit
#        ar <- logit(ar)
#        ar.byclone <- logit(ar.byclone)
#        ar.byclone[is.infinite(ar.byclone)] <- sign(ar.byclone)[is.infinite(ar.byclone)] * 100
#        ar.byclone[is.na(ar.byclone)] <- 0
#        optim.simple <- optim(c(ar,0), likeli.simple, data=test, control = list(reltol = 0.00001))
#        optim.complex <- optim(c(ar.byclone,0), likeli.clas, data=test, control = list(reltol = 0.00001))
#        return(optim.simple$value - optim.complex$value - 2 * max(clones))
#      }
#  
#  
#    }
#  
#  }
#  
#  
#  deltaAIC <- lapply(countTables, getAIC)

