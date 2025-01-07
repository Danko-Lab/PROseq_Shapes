#
# Plots and statistics over post-processed results
#
#

#' step type coverage
#' @export
step.coverage <- function(preds, chroms, bwSet, bwBck, scale.factor, step = 10, n.cores = 1) {
  make.mask <- function(preds, L) {
    res = vector(mode="logical", length=L)
    
    foreach.bed(preds, function(i, chrom, start, end, strand) {
      is = (start + 1) %/% step + 1
      ie = end %/% step

      res[is:ie] <<- TRUE
    })

    return(res)
  }

  strand.totals <- function(preds.chrom, bwTAP, bwNoTAP, chrom, strand) {
    tap = abs(chromStepSum.bigWig(bwTAP, chrom, step, 0))
    notap = abs(chromStepSum.bigWig(bwNoTAP, chrom, step, 0))

    mask = make.mask(preds.chrom[preds.chrom[,6] == strand,],
      length(tap))

    # get coverate
    #

    # TAP+ > 0
    s1 = tap > 0
    t1p = sum(s1)
    t1pc = sum(s1[mask])

    # TAP+ > 0 && TAP- > 0
    s2 = tap > notap & notap > 0
    t2p = sum(s2)
    t2pc = sum(s2[mask])

    # TAP+ > 0 && TAP- == 0
    s3 = tap > notap & notap == 0
    t3p = sum(s3)
    t3pc = sum(s3[mask])

    # TAP+ > 0 && TAP+ < TAP-
    s4 = tap > 0 & tap < notap
    t4p = sum(s4)
    t4pc = sum(s4[mask])

    return(list(total = c(t1p, t2p, t3p, t4p),
                cover = c(t1pc, t2pc, t3pc, t4pc)))
  }
  
    
  # types:
  # TAP+ > 0
  # TAP+ > TAP- & TAP- > 0
  # TAP+ > TAP- & TAP- == 0
  # TAP+ > 0 & TAP- < 0
  totals = c(0, 0, 0, 0)
  totals.covered = c(0, 0, 0, 0)
  
  chrom.totals = mclapply(chroms, function(chrom) {
    preds.chrom = preds[preds[,1] == chrom,]

    # collect data
    # - plus strand
    cat("  *", chrom, "(+)\n")
    tp = strand.totals(preds.chrom, bwSet$GROcap.plus, bwBck$GROcap.plus,
      chrom, '+')

    totals = totals + tp$total
    totals.covered = totals.covered + tp$cover

    # - minus strand
    cat("  *", chrom, "(-)\n")
    tm = strand.totals(preds.chrom, bwSet$GROcap.minus, bwBck$GROcap.minus,
      chrom, '-')

    return(list(totals = tp$total + tm$total,
                covered = tp$cover + tm$cover))
  }, mc.cores = n.cores)

  # aggregate
  return(Reduce(function(r1, r2) { list(totals = r1$totals + r2$totals, covered = r1$covered + r2$covered) }, chrom.totals, list(totals = 0, covered = 0)))
}

#' Library coverage
#' @export
read.coverage <- function(preds, chroms, bwSet, scale.factor, step = 10, n.cores = 1) {
  make.mask <- function(preds, L) {
    res = vector(mode="logical", length=L)
    
    foreach.bed(preds, function(i, chrom, start, end, strand) {
      is = (start + 1) %/% step + 1
      ie = end %/% step

      res[is:ie] <<- TRUE
    })

    return(res)
  }

  strand.totals <- function(preds.chrom, bwTAP, chrom, strand) {
    tap = abs(chromStepSum.bigWig(bwTAP, chrom, step, 0))

    mask = make.mask(preds.chrom[preds.chrom[,6] == strand,],
      length(tap))

    total = sum(tap)
    total.covered = sum(tap[mask])

    return(list(total = total, cover = total.covered))
  }

  # chrom data
  chrom.totals = mclapply(chroms, function(chrom) {
    preds.chrom = preds[preds[,1] == chrom,]

    # collect data
    # - plus strand
    cat("  *", chrom, "(+)\n")
    tp = strand.totals(preds.chrom, bwSet$GROcap.plus, chrom, '+')

    # - minus strand
    cat("  *", chrom, "(-)\n")
    tm = strand.totals(preds.chrom, bwSet$GROcap.minus, chrom, '-')

    return(list(totals = tp$total + tm$total,
                covered = tp$cover + tm$cover))
  }, mc.cores = n.cores)

  # merge
  return(Reduce(function(r1, r2) { list(totals = r1$totals + r2$totals, covered = r1$covered + r2$covered) }, chrom.totals, list(totals = 0, covered = 0)))
}
