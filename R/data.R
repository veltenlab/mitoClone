#' Allele Count tables
#'
#'M: Mutant allele counts; N: Reference allele counts. P1: Patient 1; P2: Patient 2
#' @format a data frame of variable sites (columns) across single cells (rows)
#' @name data
NULL

#'@rdname data
"M_P1"

#'@rdname data
"N_P1"

#'@rdname data
"M_P2"

#'@rdname data
"N_P2"

#' Mitochondrial blacklist
#'
#' Blacklist of variants that are likely not true somatic mutations
#' @format A list with three entries:
#' #' \itemize{
#'   \item \emph{three}: Regions of the mitochondrial genome that are within 1 nt of a 3-mer homopolymer (e.g. AAA)
#'   \item \emph{mutaseq}: Mutations in the mitochondrial genome that were reoccuring across patients (present in more than one individual in the MutaSeq dataset)
#'   \item \emph{masked}: Regions of the mitochondrial genome that are soft-masked in the UCSC or Ensembl annotations
#' }
"blacklists"
