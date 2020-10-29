#'Inference of mutational trees by of single cell mutational status
#'
#'From data on the observed mutational status of single cells at a number of genomic sites, computes a likely phylogenetic tree using PhISCS (https://github.com/sfu-compbio/PhISCS) and associates single cells with leaves of the tree.
#'The function \code{\link{clusterMetaclones}} should be called on the output in order to group mutations into clones using a likelihood-based approach.
#'@param mutcalls object of class \code{\link{mutationCalls}}.
#'@param fn false negative rate, i.e. the probability of only observing the reference allele if there is a mutation. #add gene-wise
#'@param fp false positive, i.e. the probability of observing the mutant allele if there is no mutation.
#'@param cores number of cores to use for PhISCS (defaults to 1)
#'@param time maximum time to be used for PhISCS optimization, in seconds (defaults to 10000)
#'@param tempfolder temporary folder to use for PhISCS output
#'@param python_env Any shell commands to execute in order to make the gurobi python package available. The easiest solution is running R from an environment where the gurobi python package is avaiable. In some settings (e.g. RStudio Server), this parameter can be used instead. \code{muta_clone} executes PhISCS using a \code{system} call to python. The value of this parameter is prepended to the call. If you have a conda environment \code{myenv} that contains gurobipy, \code{source activate myenv} can work. Occassionally RStudio Server modifies your PATH so that that the conda and source commands are not available. In that case you can for example use \code{export PATH=/path/to/conda/:$PATH; source activate myenv}. easybuild users can \code{module load anaconda/v3; source activate myenv} 
#'@param force_recalc Rerun PhISCS even if the \code{tempfolder} contains valid PhISCS output
#'@return an object of class \code{\link{mutationCalls}}, with an inferred tree structure and cell to clone assignment added.
#'@export
muta_cluster <- function(mutcalls, fn = 0.1, fp = 0.02, cores = 1, time =10000, tempfolder = tempdir(), python_env = '', force_recalc = F) {

  #prepare data and run PhISCS
  suppressWarnings(dir.create(tempfolder))

  usedata <- mutcalls@ternary[,mutcalls@cluster]
  ##dir.create(file.path(tempfolder,'out'))
  write.table(usedata, file = file.path(tempfolder,"in.txt"), quote = F, sep="\t",
              row.names = gsub("[><_]","",rownames(usedata)),col.names = gsub("[><_]","",colnames(usedata)))


  base <- system.file("extdata/python/PhISCS-I",package = "mitoClone")
  command <- sprintf("%spython %s -SCFile %s -fn %.2f -fp %.2f -o %s -threads %d -time %d --drawTree",
                     ifelse(python_env =="", "", paste0(python_env,"; ")),
                     base, file.path(tempfolder,"in.txt"),
                     fn, fp,
                     file.path(tempfolder,"out"),
                     cores, time)
  if (!file.exists(file.path(tempfolder, "out", "in.CFMatrix")) | force_recalc) {
    message("Now running the following command:", command)
    tryCatch(system(command), error = function(e) stop("PhISCS error: ",e,"Make sure that the gurobi python package is available and consider specifying python_env."))
  } else {
    message("Results found, skipping PhISCS run")
  }



  #read in the result and create tree data structure
  physics <- read.table(file.path(tempfolder, "out","in.CFMatrix"),header=T, row.names = 1,sep="\t")
  txtCon <- textConnection(unique(apply(physics, 1, paste, collapse = ",")))
  clones <- read.csv(txtCon, header=F , col.names = colnames(physics))
  mutcalls@tree <- clones2tree(clones)

  nodes.order <- getNodes(mutcalls@tree)[-1]

  #retrieve the assignment for every cell
  #compute likelihood of assignments
  clone.names <- apply(clones[,nodes.order],1,function(x) max(which(x==1)))
  clone.names <- nodes.order[clone.names]
  clone.names[is.na(clone.names)] <- "root"

  clones <- as.matrix(clones)
  clone.names -> rownames(clones)

  cell2clone.prob <- apply(usedata[,nodes.order],1, function(cell) {
    apply(clones[,nodes.order], 1, function(clone){
        #likelihood
        sum(log10(ifelse(cell == "0" & clone == 1, fn,
                         ifelse(cell == "1" & clone == 0, fp,
                                ifelse(cell == "?", 1, (1-fn) * (1-fp))))))
      })
    })
    mutcalls@cell2clone <- t(apply(cell2clone.prob, 2, function(x) {
      10^x / sum(10^x)
    }))


    #finally, determine which mutations to group:
    evaluate_likelihood <- function(data, idealized) {
      sum(log10(ifelse(data == "0" & idealized == 1, fn,
                       ifelse(data == "1" & idealized == 0, fp,
                              ifelse(data == "?", 1, (1-fn) * (1-fp))))))
    }

    ref <- evaluate_likelihood(usedata[,colnames(physics)], physics)
    mutcalls@treeLikelihoods <- sapply(colnames(physics), function(node1) {
      sapply(c(colnames(physics),"root"), function(node2) {
        newmodel <- physics
        if (node2 == "root") newmodel[,node1] <- rep(0, nrow(newmodel)) else newmodel[,node1] <- newmodel[,node2]
        evaluate_likelihood(usedata[,colnames(physics)], newmodel) - ref
      })
    })



    return(mutcalls)


}

#'Quick clustering of mutations
#'
#'Performs a quick hierarchical clustering on a object of class \code{\link{mutationCalls}}. See \code{\link{muta_cluster}} for an alternative that infers mutational trees and uses sound models of dropout.
#'@param mutcalls object of class \code{\link{mutationCalls}}.
#'@param binarize If \code{FALSE}, will use raw allele frequencies for the clustering. If \code{TRUE}, will use binarized mutation/reference/dropout calls.
#'@param ... Parameters passed to \code{\link{pheatmap::pheatmap}}
#'@return The result of running \code{\link{pheatmap::pheatmap}}
#'@export
quick_cluster <- function(mutcalls, binarize = F,drop_empty =T,  ...) {

    if (drop_empty) mutcalls@ternary <- mutcalls@ternary[apply(mutcalls@ternary,1,function(x) any(x=="1")),]
    if (binarize ) converted <- t(apply(mutcalls@ternary, 2, function(x) ifelse(x == "1", 1, ifelse(x=="0",-1,0))))
    if (!binarize) converted <- t(mutcalls@M / (mutcalls@M + mutcalls@N))

  clustered <- pheatmap::pheatmap(converted[mutcalls@cluster,],...)

}
# clusterMetaclones <- function(mutcalls, nclust = 3) {
#   #determine clustering
#
#   grouped <- pheatmap::pheatmap(mutcalls@treeLikelihoods,
#                      clustering_distance_cols = as.dist(1-cor(mutcalls@treeLikelihoods)),
#                      clustering_distance_rows = as.dist(1-cor(t(mutcalls@treeLikelihoods))),
#                      cutree_rows = nclust)
#   mutcalls@mut2clone <- cutree(grouped$tree_row, k=nclust)
#   use <- mutcalls@mut2clone[colnames(mutcalls@cell2clone)]
#   mutcalls@mainClone <- sapply(unique(use), function(mainClone) {
#     apply(mutcalls@cell2clone[,use == mainClone], 1, sum)
#   }) #there is a problem here because cell2clone does not have the exact same length a mut2clone (identical muts were already excluded)
#   mutcalls
# }

#'Cluster mutations into clones - following the tree structure
#'
#'PhISCS orders all mutations into a hierarchical mutational tree; in many cases, the exact order of the acquisition of individual mutations in not unanimously determined from the data. This function computes the change in likelihood of the infered clonal assignment if two mutations are merged into a clone. Hierarchical clustering is then used to determine the clonal structure. The result is visualized and should be fine-tuned using the \code{min.lik} parameter.
#'@param mutcalls mutcalls object of class \code{\link{mutationCalls}} for which \code{\link{muta_cluster}} has been run
#'@param min.lik specifies the minimum difference in likelihood required
#'@param plot whether dendrograms should be plotted.
#'
#'@export
clusterMetaclones <- function(mutcalls, min.lik = 1, plot = T) {
  #split the tree into branches with no further splits
  branches <- getBranches(mutcalls@tree)
  mutcalls@mut2clone <- as.integer(rep(0, nrow(mutcalls@treeLikelihoods)))
  names(mutcalls@mut2clone) <- rownames(mutcalls@treeLikelihoods)

  par(mfrow= c(ceiling(sqrt(length(branches))), ceiling(sqrt(length(branches)))))

  for (i in 1:length(branches)) {
    if (sum(branches[[i]]!="root" ) <= 1) {
      mutcalls@mut2clone[branches[[i]]] <- as.integer(max(mutcalls@mut2clone) + 1)
    } else {
      #d <- as.dist(1-cor(t(mutcalls@treeLikelihoods[branches[[i]],])))
      ub <- branches[[i]][branches[[i]]!="root"]
      d <- dist( t(mutcalls@treeLikelihoods[ub,ub]) )
      #pheatmap::pheatmap(mutcalls@treeLikelihoods[branches[[i]],branches[[i]]],
      #                   clustering_distance_cols = as.dist(1-cor(mutcalls@treeLikelihoods[branches[[i]],branches[[i]]])),
      #                   clustering_distance_rows = as.dist(1-cor(t(mutcalls@treeLikelihoods[branches[[i]],branches[[i]]]))))
      cl<- hclust(d)
      tryCatch(plot(cl), error=function(e) "Failed to plot")
      mutcalls@mut2clone[branches[[i]]] <- as.integer(max(mutcalls@mut2clone) + cutree(cl, h = nrow(mutcalls@M) * min.lik))
    }

  }
  use <- mutcalls@mut2clone[colnames(mutcalls@cell2clone)]
  mutcalls@mainClone <- sapply(unique(use), function(mainClone) {
    if (sum(use == mainClone) == 1) mutcalls@cell2clone[,use == mainClone] else apply(mutcalls@cell2clone[,use == mainClone], 1, sum)
  }) #there is a problem here because cell2clone does not have the exact same length a mut2clone (identical muts were already excluded)

  return(mutcalls)
}


getBranches <- function(tree) {
  current <- c()
  while (length(tree$children) == 1) {
    current <- c(current, tree$mutation)
    tree <- tree$children[[1]]
  }
  if (length(tree$children) == 0) return(c(current, tree$mutation))
  result <- removeDepth(lapply(tree$children, getBranches))
  result[[length(result)+1]] <- c(current, tree$mutation)
  return(result)
}

removeDepth <- function(list) {
  out <- list()
  for (i in 1:length(list)) {
    if (is.list(list[[i]])) {
      for (j in 1:length(list[[i]])) {
        out[[length(out)+1]] <- list[[i]][[j]]
      }
    } else out[[length(out)+1]] <- list[[i]]
  }
  out
}

