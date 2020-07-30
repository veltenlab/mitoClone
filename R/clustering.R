#'Clustering of single cells by mutation status
#'
#'From data on the observed mutational status of single cells at a number of genomic sites, performs the following analysis:
#'a) Computes a likely phylogenetic tree using PhISCS (https://github.com/sfu-compbio/PhISCS)
#'b) Groups mutations in clusters using a likelihood based approach
#'c) Associates single cells with the clusters and computes a likelihood score
#'@param mutcalls object of class \code{\link{mutationCalls}}.
#'@param fn false negative rate, i.e. the probability of only observing the reference allele if there is a mutation. #add gene-wise
#'@param fp false positive, i.e. the probability of observing the mutant allele if there is no mutation.
#'@param png file name to plot the tree into, or NULL if no plotting is required
#'@param cores number of cores to use for PhISCS (defaults to 1)
#'@param time maximum time to be used for PhISCS optimization, in seconds (defaults to 10000)
#'@param tempfolder temporary folder to use for PhISCS output
#'@param python_env Any shell commands to execute in order to make the gurobi python package availabler, such as \code{source activate myenv}
#'@return Tbd
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

    return(mutcalls)


}
