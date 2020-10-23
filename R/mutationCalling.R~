#'Create a mutationCalls object from a list of single-cell BAM files
#'
#'Uses the \code{deepSNV} package to count nucleotide frequencies at every position in the mitochondrial genome for every cell and passes the result to the \code{\link{mutationCallsFromDeepSNV}} function for filtering variants.
#'@param bamfiles A character vector specifyign the bam file paths
#'@param sites Vector specifying genomic regions, defaults to the entire mitochondrial genome. Can be modified to query nuclear mutations
#'@return A list of base count matrices which can serve as an input to \code{\link{mutationCallsFromBlacklist}} or \code{\link{mutationCallsFromCohort}}
#'@export
baseCountsFromBamList <- function(bamfiles, sites = "chrM:1-16659") {
NULL
}


#'Create a mutationCalls object from a list of single-cell BAM files
#'
#'More explanations of what happens
#'@param bam Path to the bam file
#'@param sites Vector specifying genomic regions, defaults to the entire mitochondrial genome. Can be modified to query nuclear mutations
#'@param tag Name of the bam file tag
#'@return A list of base count matrices which can serve as an input to \code{\link{mutationCallsFromBlacklist}} or \code{\link{mutationCallsFromCohort}}
#'@export
baseCountsFromSingleBam <- function(bam, sites = "chrM:1-16659", tag = "XQ:C") {
NULL
}

#'Create a mutationCalls objects from nucleotide base calls and defines a blacklist (cohort)
#'
#'Identifies relevant mitochondrial somatic variants from raw counts of nucleotide frequencies measured in single cells from several individuals. Applies two sets of filters: In the first step, filters on coverage to include potentially noisy variants; in the second step, compares allele frequencies between patients to remove variants that were observed in several individuals and that therefore are unlikely to represent true somatic variants (e.g. RNA editing events). The blacklist derived from the MutaSeq dataaset is available in \code{\link{blacklist}} and can be used on single individuals using \code{\link{mutationCallsFromBlacklist}}
#'@param BaseCounts A list of base call matrices (one matrix per cell) as produced by \code{\link{baseCountsFromSingleBam}} or \code{\link{baseCountsFromBamList}}
#'@param patient A character vector associating each cell / entry in the \code{BaseCount} list with a patient
#'@param MINREADS Minimum number of reads on a site in a single cell to qualify the site as covered
#'@param MINCELL Minimum number of cells across the whole data set to cover a site
#'@param MINFRAC Fraction of reads on the mutant allele to provisionally classify a cell as mutant
#'@param MINCELLS.PATIENT Minimum number of mutant cells per patient to classify the mutation as relevant in that patient, AND
#'@param MINRELATIVE.PATIENT Minimum fraction of mutant cells per patient to classify the mutation as relevant in that patient
#'@param MINRELATIVE.OTHER Minimum fraction of mutant cells identified in a second patient for the mutation to be excluded
#'@return A list of \code{\link{mutationCalls}} objects (one for each \code{patient}) and an entry \code{blacklist} containing a blacklist of sites with variants in several individuals
#'@export
mutationCallsFromCohort <- function(BaseCounts, patient, MINREADS = 5, MINCELL = 20, MINFRAC = 0.1, MINCELLS.PATIENT = 10, MINRELATIVE.PATIENT = 0.01, MINRELATIVE.OTHER = 0.1) {
  nuc.count.per.position.array <- array(data = 0,
                                        dim = c(length(BaseCounts),
                                                nrow(BaseCounts[[1]]),
                                                ncol(BaseCounts[[1]])),
                                        dimnames = list(names(BaseCounts),
                                                        paste0("X", 1:nrow(BaseCounts[[1]])),
                                                        colnames(BaseCounts[[1]])
                                        ))

  for (i in 1:length(BaseCounts)) nuc.count.per.position.array[i,,] <- BaseCounts[[i]]

  #determine the overall reference
  sum.overall <- apply(nuc.count.per.position.array, c(2,3), sum)
  reference <- colnames(sum.overall)[apply(sum.overall, 1, which.max)]

  mt.reads.per.cell <- apply(nuc.count.per.position.array, 1, sum)


  #turn the array into a binary array of VARIANTS
  variant_calls <- lapply(1:length(reference), function(pos) {
    #which variants exist at this site?
    support <- apply(nuc.count.per.position.array[,pos,] > MINREADS,2,sum )
    support <- support[!names(support) %in% c(reference[pos], "N")]
    use <- names(support)[support > MINCELL]

    if (length(use) == 0) NULL else {
      out <- matrix(data =NA, ncol = length(use), nrow = nrow(nuc.count.per.position.array[,pos,]),
                    dimnames = list(rownames(nuc.count.per.position.array[,pos,]), paste0(pos,reference[pos],">",use)))
      for (i in 1:length(use)) {
        pos_sum <- apply(nuc.count.per.position.array[,pos,],1,sum)
        condition_mut <- nuc.count.per.position.array[,pos,use[i]] > MINREADS &  nuc.count.per.position.array[,pos,use[i]] > MINFRAC * pos_sum
        condition_ref <- nuc.count.per.position.array[,pos,reference[pos]] > MINREADS &  nuc.count.per.position.array[,pos,reference[pos]] > MINFRAC * pos_sum

        out[, i] <- ifelse(condition_mut,
                           ifelse(condition_ref, "BOTH", "MUT"),
                           ifelse(condition_ref, "WT", "DROP")
        )
      }
      out
    }

  })

  variant_calls <- do.call(cbind, variant_calls)

  #now check:
  #how often does each variant exist per patient?
  varcount.bypatient <- sapply(unique(patient), function(pa) {
    apply(variant_calls[patient == pa, ] ,2, function(x) sum(x %in% c("BOTH","MUT")))
  })


  patient.count <-  as.vector(table(patient)[colnames(varcount.bypatient)])
  names(patient.count) <- colnames(varcount.bypatient)
  varcount.bypatient.fraction <- t(t(varcount.bypatient) / patient.count)

  #throw out anything with less than
  #a) 10 cells and 1% support in any patient
  filter <- apply(varcount.bypatient, 1, max) > MINCELLS.PATIENT & apply(varcount.bypatient, 1, function(x) max(x) / patient.count[which.max(x)] ) > MINRELATIVE.PATIENT

  #b) support of more than 10 cells or 10% the level in a second patient
  patientfilter <- filter & apply(varcount.bypatient.fraction, 1, function(x) sum(x > MINRELATIVE.OTHER*max(x)) ) == 1 &
    apply(varcount.bypatient, 1, function(x) sum(x >= MINCELLS.PATIENT) ) == 1

  #what are the variants that occur in multiple patients abundantly?
  multipatient <- filter & apply(varcount.bypatient.fraction, 1, function(x) sum(x > MINRELATIVE.OTHER*max(x)) ) > 1 &
    apply(varcount.bypatient, 1, function(x) sum(x >= MINCELLS.PATIENT) ) > 1 & !grepl(">-",rownames(varcount.bypatient))

  mutation.bypatient <- colnames(varcount.bypatient)[apply(varcount.bypatient[patientfilter,],1,which.max)]

  variant_calls_selected <- variant_calls[,patientfilter]

  #now, prepare return values.
  mutation.bypatient <- mutation.bypatient[!grepl("->",colnames(variant_calls_selected))]
  variant_calls_selected <- variant_calls_selected[,!grepl("->",colnames(variant_calls_selected))]

  out <- lapply(unique(patient), function(pa) {
    #a, retrieve matrices of allele counts for patient specific variants
    if (sum(mutation.bypatient == pa) == 0) return(NULL)
    MN <- pullcounts.vars(BaseCounts[patient == pa], colnames(variant_calls_selected)[mutation.bypatient == pa])
    #b, create mutationCalls object
    o <- mutationCallsFromMatrix(t(MN$M), t(MN$N))

  })

  names(out) <- unique(patient)

  #discuss this part with ben - not consistent with his blacklist.
  out$blacklist <- rownames(varcount.bypatient[multipatient,])
  out$blacklist <- gsub("(\\d+)(.+)","\\1 \\2", out$blacklist)
  return(out)
}


#'Create a mutationCalls object from nucleotide base calls using a blacklist (single individual)
#'
#'Identifies relevant mitochondrial somatic variants from raw counts of nucleotide frequencies. Applies two sets of filters: In the first step, filters on coverage and minimum allele frequency to include potentially noisy variants; in the second step, filters against a blacklist of variants that were observed in several individuals and that therefore are unlikely to represent true somatic variants (e.g. RNA editing events). These blacklists are created using \code{\link{mutationCallsFromCohort}}
#'@param BaseCounts A list of base call matrices (one matrix per cell) as produced by \code{\link{baseCountsFromSingleBam}} or \code{\link{baseCountsFromBamList}}
#'@param lim.cov Minimal coverage required per cell for a cell to be classified as covered
#'@param min.af Minimal allele frequency for a cell to be classified as mutant
#'@param min.num.samples Minimal number of cells required to be classified as covered and mutant according to the thresholds set in \code{lim.cov} and \code{min.af}. Usually specified as a fraction of the total number of cells.
#'@param universal.var.cells Maximum number of cells required to be classified as mutant according to the threshold set in \code{min.af.universal}.  Usually specified as a fraction of the total number of cells; serves to avoid e.g. germline variants.
#'@param min.af.universal Minimal allele frequency for a cell to be classified as mutant, in the context of removing universal variants. Defaults to \code{min.af}, but can be set to lower values.
#'@param blacklists Blacklists to use **explanations ben**
#'@param max.var.na Final filtering step: Remove all mutations with no coverage in more than this fraction of cells
#'@param max.cell.na Final filtering step: Remove all cells with no coverage in more than this fraction of mutations
#'@param ... Parameters passed to \code{\link{mutationCallsFromMatrix}}
#'@return An object of class \code{\link{mutationCalls}}
#'@export
mutationCallsFromBlacklist <- function(BaseCounts,lim.cov=20, min.af=0.2, min.num.samples=0.01*length(BaseCounts), min.af.universal =min.af, universal.var.cells=0.95*length(BaseCounts), blacklists.use = blacklists, max.var.na = 0.5, max.cell.na = 0.95, ...) {
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

  is.names <- sapply(blacklists.use, function(x) typeof(x) == "character")
  #part 2 - filter based on the blacklist
  if(sum(is.names) > 0){
    removal.names.list <- unique(unlist(blacklists.use[is.names]))
    varaf <- varaf[!row.names(varaf) %in% removal.names.list,]
  }
  if(sum(!is.names) > 0){
    removal.ranges.list <- unique(unlist(GenomicRanges::GRangesList(blacklists.use[!is.names])))
    varaf <- varaf[-c(S4Vectors::queryHits(GenomicRanges::findOverlaps(mut2gr(row.names(varaf)),removal.ranges.list))),]
  }
  #if(drop.empty){
  varaf <- varaf[rowSums(varaf,na.rm=T) > 0,] #colSums(varaf,na.rm=T) > 0
  #}

  varaf <- varaf[!rowSums(varaf >= min.af.universal,na.rm=TRUE) >= universal.var.cells,]
  ## vars must have less than X % NA's
  varaf <- varaf[rowSums(is.na(varaf)) < max.var.na*NCOL(varaf),]
  ## cells must have less than X % NA's
  varaf <- varaf[,colSums(is.na(varaf)) < max.cell.na*NROW(varaf)]

  MN <- pullcounts.vars(BaseCounts, rownames(varaf), colnames(varaf))
  mutationCallsFromMatrix(t(MN$M), t(MN$N), ...)
}


### convert mutation to GRanges
mut2gr <- function(mut) {
  gr <- GenomicRanges::GRanges(paste0('chrM:',as.numeric(gsub(" [A-Z].*","",mut))))
  gr$ref <- sapply(mut,function(x) { unlist(strsplit(gsub("\\d+ ","",x),">"))[[1]] })
  gr$alt <- sapply(mut,function(x) { unlist(strsplit(gsub("\\d+ ","",x),">"))[[2]] })
  return(gr)
}


#'Pull variant counts
#'
#'@param BaseCounts A list of base call matrices (one matrix per cell) as produced by \code{\link{baseCountsFromSingleBam}} or \code{\link{baseCountsFromBamList}}
#'@param vars Character vector of variants to pull, in format 5643G>T
#'@param cells Character vector for cells to select, or NULL if all cells from the input are to be used
#'@return A list with two entries, M (count table on the variant allele) and N (count table on the reference allele)
#'@export
pullcounts.vars <- function(mc.out,vars, cells=NULL, shift=0){
  pos <- as.integer(gsub(" *[ACGT].+","",vars)) + shift
  ref <- gsub("\\d+ *([ACGT])>(.+)","\\1",vars)
  alt <- gsub("\\d+ *([ACGT])>(.+)","\\2",vars)
  N <- sapply(mc.out, function(cell) {
    mapply(function(p,x) cell[p,x], pos, ref)
  })
  M <- sapply(mc.out, function(cell) {
    mapply(function(p,x) cell[p,x], pos, alt)
  })
  if(!is.matrix(M)) {
    M <- matrix(M, ncol = length(M),dimnames = list(vars, names(M)))
    N <- matrix(N, ncol = length(N),dimnames = list(vars,names(N)))
  }
  rownames(M) <- vars -> rownames(N)
  if (is.null(cells)) return(list(M = M, N = N)) else return(list(M = M[,cells], N = N[,cells]))
  # var.gr <- mut2gr(vars)
  # varcounts <- sapply(mc.out,function(x){
  #   ## focus on A,G,C,T
  #   x.ref <- sapply(seq(var.gr),function(y){
  #     x[BiocGenerics::start(var.gr)[[y]],var.gr$ref[[y]]]
  #   })
  #   x.alt <- sapply(seq(var.gr),function(y){
  #     x[BiocGenerics::start(var.gr)[[y]],var.gr$alt[[y]]]
  #   })
  #   return(list('ref'=x.ref,'alt'=x.alt))
  # })
  # var.ref <- do.call(cbind,varcounts['ref',])
  # var.alt <- do.call(cbind,varcounts['alt',])
  # row.names(var.alt) <- row.names(var.ref) <- vars
  # return(list('M'=var.alt[,cells],'N'=var.ref[,cells]))
}
### convert mutation to GRanges
