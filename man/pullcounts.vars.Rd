% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mutationCalling.R
\name{pullcounts.vars}
\alias{pullcounts.vars}
\title{Pull variant counts}
\usage{
pullcounts.vars(mc.out, vars, cells = NULL, shift = 0)
}
\arguments{
\item{mc.out}{A list of base call matrices (one matrix per cell) as produced by \code{\link{baseCountsFromBamList}}}

\item{vars}{Character vector of variants to pull, in format 5643G>T}

\item{cells}{Character vector for cells to select, or NULL if all cells from the input are to be used}
}
\value{
A list with two entries, M (count table on the variant allele) and N (count table on the reference allele)
}
\description{
Pull variant counts
}
