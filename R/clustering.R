#'Clustering of single cells by mutation status
#'
#'From data on the observed mutational status of single cells at a number of genomic sites, computes a likely phylogenetic tree using PhISCS (https://github.com/sfu-compbio/PhISCS) and associates single cells with leaves of the tree.
#'The function \code{\link{clusterMetaclones}} should be called on the output in order to group mutations into clones using a likelihood-based approach.
#'@param mutcalls object of class \code{\link{mutationCalls}}.
#'@param fn false negative rate, i.e. the probability of only observing the reference allele if there is a mutation. #add gene-wise
#'@param fp false positive, i.e. the probability of observing the mutant allele if there is no mutation.
#'@param png file name to plot the tree into, or NULL if no plotting is required
#'@param cores number of cores to use for PhISCS (defaults to 1)
#'@param time maximum time to be used for PhISCS optimization, in seconds (defaults to 10000)
#'@param tempfolder temporary folder to use for PhISCS output
#'@param python_env Any shell commands to execute in order to make the gurobi python package availabler, such as \code{source activate myenv}
#'@return an object of class \code{\link{mutationCalls}}, with an inferred tree structure and cell to clone assignment added.
#'@export
muta_cluster <- function(mutcalls, fn = 0.1, fp = 0.02, png = NULL, cores = 1, time =10000, tempfolder = tempdir(), python_env = "") {

  #prepare data and run PhISCS
  usedata <- mutcalls@ternary[,mutcalls@cluster]
  write.table(usedata, file = file.path(tempfolder,"in.txt"), quote = F, sep="\t",
              row.names = gsub("[><_]","",rownames(usedata)),col.names = gsub("[><_]","",colnames(usedata)))


  base <- system.file("extdata/python/PhISCS-I",package = "mitoseq")
  command <- sprintf("%spython %s -SCFile %s -fn %.2f -fp %.2f -o %s -threads %d -time %d --drawTree",
                     ifelse(python_env =="", "", paste0(python_env,"; ")),
                     base, file.path(tempfolder,"in.txt"),
                     fn, fp,
                     file.path(tempfolder,"out"),
                     cores, time)
  message("Now running the following command:", command)
  tryCatch(system(command), error = function(e) stop("PhISCS error: ",e,"Make sure that the gurobi python package is available and consider specifying python_env."))


  #read in the result and create tree data structure
  physics <- read.table(file.path(tempfolder, "out","in.CFMatrix"),header=T, row.names = 1)
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

    ref <- evaluate_likelihood(usedata, physics)
    mutcalls@treeLikelihoods <- sapply(colnames(physics), function(node1) {
      sapply(c(colnames(physics),"root"), function(node2) {
        newmodel <- physics
        if (node2 == "root") newmodel[,node1] <- rep(0, nrow(newmodel)) else newmodel[,node1] <- newmodel[,node2]
        evaluate_likelihood(usedata, newmodel) - ref
      })
    })



    return(mutcalls)


}



#'Cluster mutations into clones
#'
#'PhISCS orders all mutations into a hierarchical mutational tree; in many cases, the exact order of the acquisition of individual mutations in not unanimously determined from the data. This function computes the change in likelihood of the infered clonal assignment if two mutations are merged into a clone. Hierarchical clustering is then used to determine the clonal structure. The result is visualized and should be fine-tuned using the \code{nclust} parameter.
#'@param mutcalls mutcalls object of class \code{\link{mutationCalls}} for which \code{\link{muta_cluster}} has been run
#'@param nclust number of metaclones expected
#'
#'@export
clusterMetaclones <- function(mutcalls, nclust = 3) {
  #determine clustering

  grouped <- pheatmap::pheatmap(mutcalls@treeLikelihoods,
                     clustering_distance_cols = as.dist(1-cor(mutcalls@treeLikelihoods)),
                     clustering_distance_rows = as.dist(1-cor(t(mutcalls@treeLikelihoods))),
                     cutree_rows = nclust)
  mutcalls@mut2clone <- cutree(grouped$tree_row, k=nclust)
  use <- mutcalls@mut2clone[colnames(mutcalls@cell2clone)]
  mutcalls@mainClone <- sapply(unique(use), function(mainClone) {
    apply(mutcalls@cell2clone[,use == mainClone], 1, sum)
  }) #there is a problem here because cell2clone does not have the exact same length a mut2clone (identical muts were already excluded)
  mutcalls
}
