#'mutationCalls class
#'
#'To create this class from a list of bam files (where each bam file corresponds to a single cell), use \code{\link{mutationCallsFromBamList}}
#'To create this class from a bam file containing cell identity as a tag, use \code{\link{mutationCallsFromSingleBam}}
#'To create this class if you already have the matrices of mutation counts, use its contstructor, i.e. \code{mutationCallsFromMatrix(M = data1, N = data2)}
#'
#'@slot M A matrix of read counts mapping to the \emph{mutant} allele. Columns are genomic sites and rows and single cells.
#'@slot N A matrix of read counts mapping to the \emph{mutant} allele. Columns are genomic sites and rows and single cells.
#'@ternary Discretized version describing the mutational status of each gene in each cell, where 1 signfiies mutant, 0 signifies reference, and ? signifies dropout
#'@cluster Boolean vector of length \code{ncol(M)} specifying if the given mutation should be included for clustering (\code{TRUE}) or only used for annotation.
mutationCalls <- setClass(
  "mutationCalls",
  slots = c(
    M = "matrix",
    N = "matrix",
    ternary = "matrix",
    cluster = "logical"
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
mutationCallsFromMatrix <- function(M, N, cluster=NULL, binarize = function(M,N) {alleleRatio <- M/(M+N); apply(alleleRatio, 2, function(x) ifelse(is.na(x),"?", ifelse(x>0.05,"1","0")))}) {
  out <- new("mutationCalls", M=M, N=N, ternary=binarize(M,N))

  if (!is.null(cluster)) out@cluster <- cluster else {
    out@cluster <- apply(out@ternary!=0, 2, mean) > 0.2
  }
  out
}

