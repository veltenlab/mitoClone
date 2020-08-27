#'Create a mutationCalls object from a list of single-cell BAM files
#'
#'Uses the \code{deepSNV} package to count nucleotide frequencies at every position in the mitochondrial genome for every cell and passes the result to the \code{\link{mutationCallsFromDeepSNV}} function for filtering variants.
#'@param bamfiles A character vector specifyign the bam file paths
#'@return A list of base count matrices which can serve as an input to \code{\link{mutationCallsFromBaseCounts}}
#'@export
baseCountsFromBamList <- function(bamfiles) {
NULL
}


#'Create a mutationCalls object from a list of single-cell BAM files
#'
#'More explanations of what happens
#'@param bam Path to the bam file
#'@param tag Name of the bam file tag
#'@return A list of base count matrices which can serve as an input to \code{\link{mutationCallsFromBaseCounts}}
#'@export
baseCountsFromSingleBam <- function(bam, tag = "XQ:C") {
NULL
}

#'Create a mutationCalls object from nucleotide base calls
#'
#'Identifies relevant mitochondrial somatic variants from raw counts of nucleotide frequencies. Applies two sets of filters: In the first step, filters on coverage and minimum allele frequency to include potentially noisy variants; in the second step, filters against a blacklist of variants that were observed in several individuals and that therefore are unlikely to represent true somatic variants (e.g. RNA editing events).
#'@param BaseCounts A list of base call matrices (one matrix per cell) as produced by \code{\link{baseCountsFromSingleBam}} or \code{\link{baseCountsFromBamList}}
#'@param lim.cov Minimal coverage required per cell for a cell to be classified as covered
#'@param min.af Minimal allele frequency for a cell to be classified as mutant
#'@param min.num.samples Minimal number of cells required to be classified as covered and mutant according to the thresholds set in \code{lim.cov} and \code{min.af}. Usually specified as a fraction of the total number of cells.
#'@param universal.var.cells Maximum number of cells required to be classified as mutant according to the threshold set in \code{min.af}.  Usually specified as a fraction of the total number of cells; serves to avoid e.g. germline variants.
#'@param blacklists Blacklists to use **explanations ben**
#'@param max.var.na Final filtering step: Remove all mutations with no coverage in more than this fraction of cells
#'@param max.cell.na Final filtering step: Remove all cells with no coverage in more than this fraction of mutations
#'@return An object of class \code{\link{mutationCalls}}
#'@export
mutationCallsFromBaseCounts <- function(BaseCounts,lim.cov=20, min.af=0.2, min.num.samples=0.01*length(BaseCounts), universal.var.cells=0.95*length(BaseCounts), blacklists.use = c("mutaseq","masked","three"), max.var.na = 0.5, max.cell.na = 0.95) {
  varaf <- parallel::mclapply(BaseCounts,function(x){
    ## focus on A,G,C,T
    x <- x[,1:4]
    ## find cell that have less than 100 cov over agct at a given pos
    zeroes <- rowSums(x) < lim.cov
    ## af calc
    #x.af <- x/rowSums(x)
    x.af <- x / (x+apply(x,1,max))
    x.af <- reshape2::melt(x.af)
    colnames(x.af) <- c('pos','nt','af')
    ## remove reference af's
    x.af <- x.af[!(mito.dna[x.af$pos] == x.af$nt),]
    ## remove N site
    x.af <- x.af[!(mito.dna[x.af$pos] == 'N'),]
    x.af$name <- paste0(x.af$pos,' ',mito.dna[x.af$pos],'>',x.af$nt)
    ## find dominant NT
    x.af$af[x.af$pos %in% which(zeroes)] <- NA
    x <- x.af$af
    names(x) <- x.af$name
    return(x)
  }, mc.cores=10) #remove parallelism here
  varaf <- do.call(cbind, varaf)
  ## you could allow for only sites with coverage! currently you filter at a rate of 10% cells dropping out max
  ##varaf <- varaf[rowSums(is.na(varaf))/length(mc.out) < max.fraction.na,]
  varaf <- varaf[rowSums(varaf > min.af,na.rm=TRUE) >= min.num.samples,]

  is.names <- sapply(blacklists.use, function(x) typeof(blacklists[[x]]) == "character")
  #part 2 - filter based on the blacklist
  if(sum(is.names) > 0){
    removal.names.list <- unique(unlist(blacklists[blacklists.use[is.names]]))
    varaf <- varaf[!row.names(varaf) %in% removal.names.list,]
  }
  if(sum(!is.names) > 0){
    removal.ranges.list <- unique(unlist(GenomicRanges::GRangesList(blacklists[blacklists.use[!is.names]])))
    varaf <- varaf[-c(S4Vectors::queryHits(GenomicRanges::findOverlaps(mut2gr(row.names(varaf)),removal.ranges.list))),]
  }
  #if(drop.empty){
  varaf <- varaf[rowSums(varaf,na.rm=T) > 0,] #colSums(varaf,na.rm=T) > 0
  #}

  varaf <- varaf[!rowSums(varaf >= min.af,na.rm=TRUE) >= universal.var.cells,]
  ## vars must have less than X % NA's
  varaf <- varaf[rowSums(is.na(varaf)) < max.var.na*NCOL(varaf),]
  ## cells must have less than X % NA's
  varaf <- varaf[,colSums(is.na(varaf)) < max.cell.na*NROW(varaf)]
  return(varaf)
}


### convert mutation to GRanges
mut2gr <- function(mut) {
  gr <- GenomicRanges::GRanges(paste0('chrM:',as.numeric(gsub(" [A-Z].*","",mut))))
  gr$ref <- sapply(mut,function(x) { unlist(strsplit(gsub("\\d+ ","",x),">"))[[1]] })
  gr$alt <- sapply(mut,function(x) { unlist(strsplit(gsub("\\d+ ","",x),">"))[[2]] })
  return(gr)
}


#'Make this private
#'
#'@export
pullcounts.vars <- function(mc.out,vars){
  var.gr <- mut2gr(vars)
  varcounts <- sapply(mc.out,function(x){
    ## focus on A,G,C,T
    x.ref <- sapply(seq(var.gr),function(y){
      x[BiocGenerics::start(var.gr)[[y]],var.gr$ref[[y]]]
    })
    x.alt <- sapply(seq(var.gr),function(y){
      x[BiocGenerics::start(var.gr)[[y]],var.gr$alt[[y]]]
    })
    return(list('ref'=x.ref,'alt'=x.alt))
  })
  var.ref <- do.call(cbind,varcounts['ref',])
  var.alt <- do.call(cbind,varcounts['alt',])
  row.names(var.alt) <- row.names(var.ref) <- vars
  return(list('M'=var.alt,'N'=var.ref))
}
### convert mutation to GRanges
