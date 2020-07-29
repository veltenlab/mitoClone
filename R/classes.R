#'mutationCalls class
#'
#'To create this class from a list of bam files (where each bam file corresponds to a single cell), use \code{\link{mutationCallsFromBamList}}
#'To create this class from a bam file containing cell identity as a tag, use \code{\link{mutationCallsFromSingleBam}}
#'To create this class if you already have the matrices of mutation counts, use the contstructor, i.e. \code{new("mutationCalls", M = data1, N = data2)}
#'
#'@slot M A matrix of read counts mapping to the \emph{mutant} allele. Columns are genomic sites and rows and single cells.
#'@slot N A matrix of read counts mapping to the \emph{mutant} allele. Columns are genomic sites and rows and single cells.
#'@slot any other information that we want to keep?
mutationCalls <- setClass(
  "mutationCalls",
  slots = c(
    M = "matrix",
    N = "matrix"
  ),
  validity = function(object) {
    if (!identical(dim(object@M), dim(object@N))) {
      return("Matrices M and N must have identical dimensions")
    }
    return(TRUE)
  }

)

#we can specify some plotting methods for this object if we want? so as to plot some QC stat from the mutation calling, etc.

