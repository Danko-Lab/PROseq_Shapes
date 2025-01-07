
#' Classify pairs (in order, [-, +]) as SS, SU, US, UU
#'
#' @export
stable.unstable.classify <- function(pair.plus, pair.minus, bwSet, cage.thresh = 8) {
  pair.plus.values = values.preds(pair.plus, bwSet$GROcap.plus, bwSet$CAGENuc.bp.plus)
  pair.minus.values = values.preds(pair.minus, bwSet$GROcap.minus, bwSet$CAGENuc.bp.minus)

  # normaly high.values == median(log10(...))
  #
  # here, we'll use the top 80% ...

  tmp.plus = log10(pair.plus.values[,1])
  tmp.minus = log10(pair.minus.values[,1])

  thresh.plus = quantile(tmp.plus, prob = 0.2)
  thresh.minus = quantile(tmp.minus, prob = 0.2)

  #high.values = log10(pair.plus.values[,1]) > median(log10(pair.plus.values[,1])) & log10(pair.minus.values[,1]) > median(log10(pair.minus.values[,1]))

  high.values = tmp.plus > thresh.plus & tmp.minus > thresh.minus

  # stability thresholds
  #

  # redef stable/unstable
  w.cage.plus = pair.plus.values[high.values, 2] >= cage.thresh
  w.cage.minus = pair.minus.values[high.values, 2] >= cage.thresh
  # leave buffer
  wo.cage.plus = pair.plus.values[high.values, 2] == 0
  wo.cage.minus = pair.minus.values[high.values, 2] == 0

  # define classes
  #
  # order of pairs is (-, +)
  #
  stable.stable = w.cage.plus & w.cage.minus
  stable.unstable = wo.cage.plus & w.cage.minus
  unstable.stable = w.cage.plus & wo.cage.minus
  unstable.unstable = wo.cage.plus & wo.cage.minus

  
  # get actual coordinates
  ss.plus = pair.plus[high.values,][stable.stable,]
  ss.minus = pair.minus[high.values,][stable.stable,]

  su.plus = pair.plus[high.values,][stable.unstable,]
  su.minus = pair.minus[high.values,][stable.unstable,]
  
  us.plus = pair.plus[high.values,][unstable.stable,]
  us.minus = pair.minus[high.values,][unstable.stable,]
  
  uu.plus = pair.plus[high.values,][unstable.unstable,]
  uu.minus = pair.minus[high.values,][unstable.unstable,]

  return(list(ss.plus = ss.plus, ss.minus = ss.minus,
              su.plus = su.plus, su.minus = su.minus,
              us.plus = us.plus, us.minus = us.minus,
              uu.plus = uu.plus, uu.minus = uu.minus,
              thresh.plus = thresh.plus, thresh.minus = thresh.minus))
}

#' Write stable-unstable BED files
#'
#' @export
stable.unstable.write <- function(prefix, classified.lst) {
  write.bed <- function(bed, filename) {
    write.table(bed, file=filename, sep='\t', quote=F, col.names=F, row.names=F)
  }
  make.name <- function(type) {
    paste(prefix, type, "bed", sep='.')
  }

  write.bed(classified.lst$ss.plus, make.name("SS_plus"))
  write.bed(classified.lst$ss.minus, make.name("SS_minus"))
  write.bed(classified.lst$su.plus, make.name("SU_plus"))
  write.bed(classified.lst$su.minus, make.name("SU_minus"))
  write.bed(classified.lst$us.plus, make.name("US_plus"))
  write.bed(classified.lst$us.minus, make.name("US_minus"))
  write.bed(classified.lst$uu.plus, make.name("UU_plus"))
  write.bed(classified.lst$uu.minus, make.name("UU_minus"))
}

#' Filter singletons
#'
#' @export
single.filter <- function(singles, bwSet, thresh.plus, thresh.minus, cage.thresh = 8) {
  single.plus = singles[singles[,6] == '+',]
  single.minus = singles[singles[,6] == '-',]

  # get scores
  single.plus.values = values.preds(single.plus, bwSet$GROcap.plus, bwSet$CAGENuc.bp.plus)
  single.minus.values = values.preds(single.minus, bwSet$GROcap.minus, bwSet$CAGENuc.bp.minus)

  # filter by thresholds
  tmp.plus = log10(single.plus.values[,1])
  tmp.minus = log10(single.minus.values[,1])

  #
  singletons = rbind(
      single.plus[tmp.plus > thresh.plus,],
      single.minus[tmp.minus > thresh.minus,])

  # classify
  # redef stable/unstable
  w.cage.plus = single.plus.values[, 2] >= cage.thresh
  w.cage.minus = single.minus.values[, 2] >= cage.thresh
  # leave buffer
  wo.cage.plus = single.plus.values[, 2] == 0
  wo.cage.minus = single.minus.values[, 2] == 0

  stable.plus = tmp.plus > thresh.plus & w.cage.plus
  unstable.plus = tmp.plus > thresh.plus & wo.cage.plus

  stable.minus = tmp.minus > thresh.minus & w.cage.minus
  unstable.minus = tmp.minus > thresh.minus & wo.cage.minus

  stable.singletons = rbind(
      single.plus[stable.plus,],
      single.minus[stable.minus,])
  unstable.singletons = rbind(
      single.plus[unstable.plus,],
      single.minus[unstable.minus,])

  return(list(all = singletons, stable = stable.singletons, unstable = unstable.singletons))
}
