#'mutationCalls class
#'
#'To create this class from a list of bam files (where each bam file corresponds to a single cell), use \code{\link{mutationCallsFromBamList}}
#'To create this class from a bam file containing cell identity as a tag, use \code{\link{mutationCallsFromSingleBam}}
#'To create this class if you already have the matrices of mutation counts, use its contstructor, i.e. \code{mutationCallsFromMatrix(M = data1, N = data2)}
#'
#'@slot M A matrix of read counts mapping to the \emph{mutant} allele. Columns are genomic sites and rows and single cells.
#'@slot N A matrix of read counts mapping to the \emph{mutant} allele. Columns are genomic sites and rows and single cells.
#'@slot ternary Discretized version describing the mutational status of each gene in each cell, where 1 signfiies mutant, 0 signifies reference, and ? signifies dropout
#'@slot cluster Boolean vector of length \code{ncol(M)} specifying if the given mutation should be included for clustering (\code{TRUE}) or only used for annotation.
#'@slot metadata Metadata frame for annotation of single cells (used for plotting). Row names should be the same as in \code{M}
#'@slot tree Inferred mutation tree
#'@slot cell2clone Probability matrix of single cells and their assignment to clones.
#'@slot mut2clone Maps mutations to main clones
#'@slot mainClone Probability matrix of single cells and their assignment to main clones
#'@slot treeLikelihoods Likelihood matrix underlying the inference of main clones, see \code{\link{clusterMetaclones}}
#'@export
mutationCalls <- setClass(
  "mutationCalls",
  slots = c(
    M = "matrix",
    N = "matrix",
    metadata = "data.frame",
    ternary = "matrix",
    cluster = "logical",
    tree= "list",
    cell2clone = "matrix",
    mut2clone = "integer",
    mainClone = "matrix",
    treeLikelihoods = "matrix"

  ),
  validity = function(object) {
    if (!identical(dim(object@M), dim(object@N))) {
      return("Matrices M and N must have identical dimensions")
    }
    return(TRUE)
  }

)
#'mutationCalls constructor
#'
#'To be used when allele-specific count matrices are available.
#'@param M A matrix of read counts mapping to the \emph{mutant} allele. Columns are genomic sites and rows and single cells.
#'@param N A matrix of read counts mapping to the \emph{mutant} allele. Columns are genomic sites and rows and single cells.
#'@param cluster If \code{NULL}, only mutations with coverage in 20% of the cells or more will be used for the clustering, and all other mutations will be used for cluster annotation only. Alternatively, a boolean vector of length \code{ncol(M)} that specifies the desired behavior for each genomic site.
#'@param binarize A function that turns M and N into a ternary matrix of mutation calls in single cells, where where -1 signfiies reference, 0 signifies dropout, and 1 signifies wild-type. The default was found to work well on mitochondrial data.
#'@return An object of class \code{\link{mutationCalls}}.
#'@export
mutationCallsFromMatrix <- function(M, N, cluster=NULL, metadata = data.frame(row.names = rownames(M)), binarize = function(M,N) {alleleRatio <- M/(M+N); apply(alleleRatio, 2, function(x) ifelse(is.na(x),"?", ifelse(x>0.05,"1","0")))}) {
  out <- new("mutationCalls", M=M, N=N, metadata = metadata, ternary=binarize(M,N))

  if (!is.null(cluster)) out@cluster <- cluster else {
    out@cluster <- apply(out@ternary!="?", 2, mean) > 0.2 #& apply(out@ternary=="1", 2, mean) > 0.04 #the last filter was not used when I made the figure, there was a filter on the allele freq. in RNA. Should maybe include this in the other routines? But this works as well
  }
  out
}

#'Plot clonal assignment of single cells
#'
#'Creates a heatmap of single cell mutation calls, clustered using PhISCS.
#'@param mutcalls object of class \code{\link{mutationCalls}}.
#'@param what One of the following: \emph{alleleRatio}: The fraction of reads mapping to the mutant allele or \emph{ternary}: Ternarized mutation status
#'@param show boolean vector specifying for each mutation if it should be plotted on top of the heatmap as metadata; defaults to mutations not used for the clustering \code{!mutcalls@cluster}
#'@param ... any arguments passed to \code{\link{pheatmap::pheatmap}}
#'@export
plotClones <- function(mutcalls, what = "alleleRatio", show = !mutcalls@cluster, ...) {
  if (what == "alleleRatio") plotData <- mutcalls@M / (mutcalls@M + mutcalls@N)
  if (what == "ternary") plotData <- apply(mutcalls@ternary, 2, function(x) ifelse(x == "1", 1, ifelse(x=="?", 0, -1)))
  plotData <- t(plotData[,getNodes(mutcalls@tree)[-1]]) #how to order rows?
  annos <- data.frame(row.names = rownames(mutcalls@M), mutcalls@ternary[,show],
                      mutcalls@metadata)
  if (length(mutcalls@mut2clone) > 0) {
    annos$mainClone <- as.factor(apply(mutcalls@mainClone, 1, which.max))
    annos$confidence <- apply(mutcalls@mainClone, 1, max)
    plotData <- plotData[,order(annos$mainClone)]
  }

  pheatmap::pheatmap(plotData, cluster_cols = F, cluster_rows = F, show_colnames = F,
           color = colorRampPalette(rev(c("#9B0000","#FFD72E","#FFD72E","#00009B")))(100),
           annotation_col = annos, ...)

}


#'Plot clonal tree
#'
#'Uses graphviz to plot the clonal tree
#'@param mutcalls object of class \code{\link{mutationCalls}}.
#'@param file file name of postscript output file
#'@export
plotTree <- function(mutcalls, file = "mytree.ps") {
  gv <- toGraphviz(mutcalls@tree)
  tmp <- tempfile()
  writeLines(gv, con = tmp)
  system(sprintf("dot -Tps %s > %s", tmp, file),wait = T)
  file.remove(tmp)
}


#'mutationCalls accessors
#'
#'Retrieves the full matrix of likelihoods associating single cells with clones
#'@param  mutcall object of class \code{\link{mutationCalls}}.
#'@param mainClones Retrieve likelihoods associated with the main Clones. Defaults to \code{TRUE} if \code{\link{clusterMetaclones}} has been run.
#'@export

getCloneLikelihood <- function(mutcall, mainClones =length(mutcall@mut2clone) > 0)

#' @describeIn getCloneLikelihood Retrieve the most likely clone associate with each cell.
getMainClone <- function(mutcall, mainClones = length(mutcall@mut2clone) > 0) as.factor(apply(getCloneLikelihood(mutcall, mainClones = mainClones), 1, which.max))

#' @describeIn getCloneLikelihood Retrieve the likelihood of the most likely clone for each cell.
getConfidence <- function(mutcalls, mainClones = length(mutcall@mut2clone) > 0) as.factor(apply(getCloneLikelihood(mutcall, mainClones = mainClones), 1, max))


