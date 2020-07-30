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
#'@param python_env Any shell commands to execute in order to be able to run PhISCS-I, such as \code{source activate myenv}
#'@return Tbd
#'@export
muta_cluster <- function(mutcalls, fn = 0.1, fp = 0.02, png = NULL, cores = 1, time =10000, tempfolder = tempdir(), python_env = "") {
  write.table(mutcalls@ternary, file = file.path(tempfolder,"in.txt"), quote = F, sep="\t",
              row.names = gsub("[><_]","",rownames(forphiscs)),col.names = gsub("[><_]","",colnames(forphiscs)))


  base <- system.file("extdata/python/PhISCS-I",package = "mitoseq")
  command <- sprintf("%s%s -SCFile %s -fn %.2f -fp %.2f -o %s -threads %d -time %d --drawTree",
                     ifelse(pythonenv =="", "", paste0(pythonenv,"; ")),
                     base, file.path(tempfolder,"in.txt"),
                     fn, fp,
                     file.path(tempfolder,"out"),
                     cores, time)
  print(command)
}
