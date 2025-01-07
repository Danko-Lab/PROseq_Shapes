#' Collect GROcap and CAGE (or other bigWig dataset) in a set of HMM TSS predictions
#'
#' If predictions are larger than 2*hwin bp, then their GRO-cap weighted center
#' is computed and a 2*hwin window around that center is used.
#'
#' @param preds BED data-frame with HMM identified TSS regions
#' @param bwGROcap bigWig object with GRO-cap data (in the strand matching 'preds')
#' @param bwCAGE bigWig object with CAGE (or other) data (in the strand matching 'preds')
#' @param hwin half-window in bp (see details)
#' @param snd.win extension applied to the collection window for the second dataset (bwCAGE). Window is extended in the strand sense.
#' @return Nx2 matrix with data for each prediction
#' @export
values.preds <- function(preds, bwGROcap, bwCAGE, hwin = 50, snd.win = 0) {
  N = dim(preds)[1]

  res = matrix(data=0, nrow=N, ncol=2)
  
  foreach.bed(preds, function(i, chrom, start, end, strand) {
    reads = abs(queryByStep.bigWig(bwGROcap, chrom, start, end, step=1))
    L = end - start

    grocap = 0
    cage = 0
    if (L < 2*hwin) {
      grocap = sum(reads)

      if (strand == '+')
        end = end + snd.win
      else
        start = start - snd.win
      
      cage = sum(abs(queryByStep.bigWig(bwCAGE, chrom, start, end, step=1, default.null = F)))
    } else {
      offset = sum(reads*(1:L))/sum(reads) - 1
      ilow = as.integer(start + offset - hwin)
      ihigh = as.integer(start + offset + hwin)

      grocap = sum(abs(queryByStep.bigWig(bwGROcap, chrom, ilow, ihigh, step=1)))

      if (strand == '+')
        ihigh = ihigh + snd.win
      else
        ilow = ilow - snd.win

      cage = sum(abs(queryByStep.bigWig(bwCAGE, chrom, ilow, ihigh, step=1, default.null = F)))
    }

    res[i,] <<- c(grocap, cage)
  })

  return(res)
}
