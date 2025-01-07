#
# Normalization
#

#' @export
total.reads <- function(bw) {
  total = 0

  if (is.null(bw))
    return(NA)
  
  for (chrom in bw$chroms) {
    tmp = chromStepSum.bigWig(bw, chrom, 1, 0)
    total = total + sum(tmp)
  }

  return(total)
}

#' Compute the total read count for each bigWig file in bwSet's list
#'
#' @param bwSet A list of bigWig objects; for example, one produced by \code{mainBwSet.grocap}
#' @return Single step sum of values across all chromosomes, for each bigWig file
#' @export
compute.normalization <- function(bwSet) {
  lapply(bwSet, total.reads)
}

# norm = compute.normalization(bwSet)
# save(norm, file="totals.norm.Rdata")
