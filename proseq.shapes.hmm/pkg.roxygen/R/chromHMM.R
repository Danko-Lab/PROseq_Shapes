#' identify regions contained within target set
#'
#' @param container container regions
#' @param regions regions to test for inclusion in container regions
#' @return logical vector indicating inclusion of each region
#' @export
contained.mask <- function(container, regions) {
  N = dim(regions)[1]
  contained = vector(mode="logical", length=N)

  # split containers by chromosome
  conts.chrom = lapply(levels(container[,1]), function(chrom)
    container[container[,1] == chrom, 2:3])
  names(conts.chrom) <- levels(container[,1])

  #
  foreach.bed(regions, function(i, chrom, start, end, strand) {
    cset = conts.chrom[[chrom]]

    contained[i] <<- any(cset[,1] <= start & cset[,2] >= end)
  })

  return(contained)
}

#' Identify regions overlapping within target set
#'
#' @param setSrc regions to process
#' @param setSel regions used to test overlap
#' @return logical vector indicating overlap of each source region with selection set
#' @export
overlaped.mask <- function(setSrc, setSel) {
  N = dim(setSrc)[1]
  selected = vector(mode="logical", length=N)

  # split containers by chromosome
  conts.chrom = lapply(levels(setSel[,1]), function(chrom)
    setSel[setSel[,1] == chrom, 2:3])
  names(conts.chrom) <- levels(setSel[,1])

  #
  foreach.bed(setSrc, function(i, chrom, start, end, strand) {
    cset = conts.chrom[[chrom]]

    selected[i] <<- any(start < cset[,2] & end > cset[,1])
  })

  return(selected)
}

#' compute fraction of regions contained within target set
#'
#' @param container container regions
#' @param regions regions to test for inclusion in container regions
#' @return fraction of source regions contained inside target set
#' @export
contained.fraction <- function(container, regions) {
  contained = contained.mask(container, regions)
  
  # result
  sum(contained)/length(contained)
}

#' identify pairs contained within target set
#'
#' @param pair.plus bed data.frame with plus side of pairs
#' @param pair.minus bed data.frame with minus side of pairs
#' @param container container regions
#' @return logical vector indicating inclusion of each pair
#' @export pairs.in.regions
pairs.in.regions <- function(pair.plus, pair.minus, container) {
  mask.plus = contained.mask(container, pair.plus)
  mask.minus = contained.mask(container, pair.minus)

  return(mask.plus & mask.minus)
}


#' Compute number of covered bases
#'
#' @param bed BED data.frame
#' @return total number of covered bases
#' @export
bitmask.total <- function(bed) {
  bed.chroms = lapply(levels(bed[,1]), function(chrom) bed[bed[,1] == chrom,])
  total = 0
  names(bed.chroms) <- levels(bed[,1])
  
  for (chrom in levels(bed[,1])) {
    subset = bed.chroms[[chrom]]
    mask = vector(mode="logical", length=max(subset[,3]))

    foreach.bed(subset, function(i, chrom, start, end, strand) {
      start = start + 1
      mask[start:end] <<- T
    })

    total = total + sum(mask)
  }

  return(total)
}

#' Load chromHMM regions for given cellLine
#'
#' @param cellLine cell line to use (gm12878 or k562)
#' @param prom.states character string with state numbers to use as promoters
#' @param enh.states character string with state numbers to use as enhancers
#' @param basePath Root directory for chromHMM bed files
#' @return list with promoter(proms) and enhancer(enh) sets
#' @export
chromHMM.regions <- function(cellLine, prom.states = "12", enh.states = "4567", basePath = "/usr/data/GROseq.parser/hg19/") {
  chrom.hmm = NULL

  if (cellLine == "gm12878")
    chrom.hmm = read.table(paste(basePath, "gm12878/chromhmm/wgEncodeBroadHmmGm12878HMM.bed.gz", sep='/'), colClasses = c("factor", "integer", "integer", "factor", "integer", "factor", "integer", "integer", "factor"))
  else
    chrom.hmm = read.table(paste(basePath, "k562/chromhmm/wgEncodeBroadHmmK562HMM.bed.gz", sep='/'), colClasses = c("factor", "integer", "integer", "factor", "integer", "factor", "integer", "integer", "factor"))

  # split into promoters and enhancers
  prom.mask = paste("^[", prom.states, "]_", sep='')
  enh.mask = paste("^[", enh.states, "]_", sep='')

  proms = grepl(prom.mask, chrom.hmm[,4])
  enh = grepl(enh.mask, chrom.hmm[,4])

  #
  return(list(proms = chrom.hmm[proms,], enh = chrom.hmm[enh,]))
}
