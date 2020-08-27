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
#'   \item \emph{three} **explanations ben**
#'   \item \emph{mutaseq} **explanations ben**
#'   \item \emph{masked} **explanations ben**
#' }
"blacklists"
