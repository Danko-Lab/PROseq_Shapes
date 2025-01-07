
#' Take a vector of read values (in 10bp steps) and produce a thresholded vector of zeros and ones
#'
#' @param reads list of read vectors
#' @param scale.factor scale factor
#' @param thresh a step is peaked if value is above this value
#' @param w number of steps to look for background values (pos +/- w)
#' @param k minimum number of background steps to run log2 ratio test
#' @return data suitable for the Naive HMMs
#' @export
processReads <- function(reads, scale.factor, thresh, testMode = FALSE, w = 10, k = 2) {
  if (testMode) {
    reads.scale = abs(reads[0:100000]) * scale.factor
  } else {
    reads.scale = abs(reads) * scale.factor
  }
  meanReadVal = mean(reads.scale)
  readLength = length(reads.scale)
  X = rep(0, readLength) # default no data
  
  useful = reads.scale > 0
  if (any(useful)) {
    X[useful] = sapply(which(useful), function(idx) {
      # val = log2(reads.scale[idx] / meanReadVal)
      
      # if (val > thresh)
      #  return(1)
      
      # return(0)
      
      i = max(idx - w, 1)
      j = min(idx + w, readLength)
      bck.idxs = which(reads.scale[i:j] > 0)
      if (length(bck.idxs) >= k) {
        bck = mean(reads.scale[i:j][bck.idxs])
        val = log2(reads.scale[idx] / bck)
        
        if (val > thresh) {
          return(1)
        }
      }
      return(0)
    })
  }
  
  # aboveThreshold = reads.scale > thresh
  # aboveThreshold = log2(reads.scale / meanReadVal) > thresh
  # X[aboveThreshold[1:readLength]] = 1
  
  return(X)
}

#' Take a vector of read values (in 10bp steps) and produce a thresholded vector of zeros and ones
#'
#' @param reads list of read vectors
#' @param scale.factor scale factor
#' @param thresh a step is peaked if value is above this value
#' @return data suitable for the Naive HMMs
#' @export
processReads.globalMean <- function(reads, scale.factor, thresh, testMode = FALSE) {
  if (testMode) {
    reads.scale = abs(reads[0:100000]) * scale.factor
  } else {
    reads.scale = abs(reads) * scale.factor
  }
  meanReadVal = mean(reads.scale)
  readLength = length(reads.scale)
  X = rep(0, readLength) # default no data
  
  useful = reads.scale > 0
  if (any(useful)) {
    X[useful] = sapply(which(useful), function(idx) {
      val = log2(reads.scale[idx] / meanReadVal)
      
      if (val > thresh)
        return(1)
      
      return(0)
    })
  }
  return(X)
}

#' Take a list of vectors of read counts (in 10bp steps) and produce a combined dataset of all processed vectors
#' needed by the HMM
#'
#' @param tss list of read vectors
#' @param plusGenestart list of read vectors
#' @param plusGenebody list of read vectors
#' @param plusGeneend list of read vectors
#' @param plusAftergene list of read vectors
#' @param scale.factor scale factor
#' @param log2.thresh a step is peaked if value is above this value
#' @return data suitable for the Naive HMMs
#' @export
combine5.sequence.to.data <- function(tss, plusGenestart, plusGenebody, plusGeneend, plusAftergene, 
                                      scale.factor, w = 10, k = 2, thresh = 1, testMode = FALSE) {
  
  # browser()
  E1 = processReads(tss, scale.factor, thresh, testMode)
  E2 = processReads(plusGenestart, scale.factor, thresh, testMode)
  E3 = processReads(plusGenebody, scale.factor, thresh, testMode)
  E4 = processReads(plusGeneend, scale.factor, thresh, testMode)
  E5 = processReads(plusAftergene, scale.factor, thresh, testMode)

  cat("tss mean:", mean(E1), "\n")
  cat("plusGenestart mean:", mean(E2), "\n")
  cat("plusGenebody mean:", mean(E3), "\n")
  cat("plusGeneend mean:", mean(E4), "\n")
  cat("plusAftergene mean:", mean(E5), "\n")

  return(rbind(E1, E2, E3, E4, E5))
}

#' Take a list of vectors of read counts (in 10bp steps) and produce a combined dataset of all processed vectors
#' needed by the HMM
#'
#' @param tss list of read vectors
#' @param plusGenestart list of read vectors
#' @param plusGenebody list of read vectors
#' @param plusGeneend list of read vectors
#' @param plusAftergene list of read vectors
#' @param plusNeg list of read vectors
#' @param scale.factor scale factor
#' @param log2.thresh a step is peaked if value is above this value
#' @return data suitable for the Naive HMMs
#' @export
combine6.sequence.to.data <- function(tss, plusGenestart, plusGenebody, plusGeneend, plusAftergene, plusNeg, 
                                      scale.factor, w = 10, k = 2, thresh = 1, testMode = FALSE) {
  
  browser()
  E1 = processReads(tss, scale.factor, thresh, testMode)
  E2 = processReads(plusGenestart, scale.factor, thresh, testMode)
  E3 = processReads(plusGenebody, scale.factor, thresh, testMode)
  E4 = processReads(plusGeneend, scale.factor, thresh, testMode)
  E5 = processReads(plusAftergene, scale.factor, thresh, testMode)
  E6 = processReads(plusNeg, scale.factor, thresh, testMode)
  
  cat("tss mean:", mean(E1), "\n")
  cat("plusGenestart mean:", mean(E2), "\n")
  cat("plusGenebody mean:", mean(E3), "\n")
  cat("plusGeneend mean:", mean(E4), "\n")
  cat("plusAftergene mean:", mean(E5), "\n")
  cat("plusNeg mean:", mean(E6), "\n")
  
  return(rbind(E1, E2, E3, E4, E5, E6))
}

#' Take a list of vectors of read counts (in 10bp steps) and produce a combined dataset of all processed vectors
#' needed by the HMM
#'
#' @param plusGenestart list of read vectors
#' @param plusGenebody list of read vectors
#' @param plusGeneend list of read vectors
#' @param plusNeg list of read vectors
#' @param scale.factor scale factor
#' @param log2.thresh a step is peaked if value is above this value
#' @return data suitable for the Naive HMMs
#' @export
combine4.sequence.to.data <- function(plusGenestart, plusGenebody, plusGeneend, plusNeg, 
                                      scale.factor, w = 10, k = 2, thresh = 1, testMode = FALSE) {
  
  # browser()
  E1 = processReads(plusGenestart, scale.factor, thresh, testMode)
  E2 = processReads(plusGenebody, scale.factor, thresh, testMode)
  E3 = processReads(plusGeneend, scale.factor, thresh, testMode)
  E4 = processReads(plusNeg, scale.factor, thresh, testMode)
  
  cat("plusGenestart mean:", mean(E1), "\n")
  cat("plusGenebody mean:", mean(E2), "\n")
  cat("plusGeneend mean:", mean(E3), "\n")
  cat("plusNeg mean:", mean(E4), "\n")
  
  return(rbind(E1, E2, E3, E4))
}

#' Take a list of vectors of read counts (in 10bp steps) and produce a combined dataset of all processed vectors
#' needed by the HMM
#'
#' @param plusGenestart list of read vectors
#' @param plusGenebody list of read vectors
#' @param plusGeneend list of read vectors
#' @param scale.factor scale factor
#' @param log2.thresh a step is peaked if value is above this value
#' @return data suitable for the Naive HMMs
#' @export
combine3.sequence.to.data <- function(plusGenestart, plusGenebody, plusGeneend, 
                                      scale.factor, w = 10, k = 2, thresh = 1, testMode = FALSE) {
  
  # browser()
  E1 = processReads(plusGenestart, scale.factor, thresh, testMode)
  E2 = processReads(plusGenebody, scale.factor, thresh, testMode)
  E3 = processReads(plusGeneend, scale.factor, thresh, testMode)
  
  return(rbind(E1, E2, E3))
}

#' Take a list of vectors of read counts (in 10bp steps) and produce a combined dataset of all processed vectors
#' needed by the HMM
#'
#' @param plusGenestart list of read vectors
#' @param plusGenebody list of read vectors
#' @param scale.factor scale factor
#' @param log2.thresh a step is peaked if value is above this value
#' @return data suitable for the Naive HMMs
#' @export
combine2.sequence.to.data <- function(plusGenestart, plusGenebody, 
                                      scale.factor, w = 10, k = 2, thresh = 1, testMode = FALSE) {
  
  # browser()
  E1 = processReads(plusGenestart, scale.factor, thresh, testMode)
  E2 = processReads(plusGenebody, scale.factor, thresh, testMode)

  return(rbind(E1, E2))
}

#' Take a vector of read counts (in 10bp steps) and produce the (Y,X) pair
#' needed by the HMM
#'
#' Peaks identification makes use of steps in the local region marked as depleated. If a sufficient number is present, a step is called a peak if the TAP+ read count is above the mean TAP+ reads in neighboring background (depleated) steps.
#'
#' @param reads.in.steps.tap reads in TAP+ condition
#' @param reads.in.steps.notap reads in TAP- condition
#' @param scale.factor TAP- to TAP+ scale factor
#' @param w number of steps to look for background values (pos +/- w)
#' @param k minimum number of background steps to run log2 ratio test
#' @param log2.thresh a step is peaked if ratio to background windows is above this value
#' @return data suitable for the Naive HMMs
#' @export
sequence.to.data <- function(reads.in.steps.tap, reads.in.steps.notap, scale.factor, w = 10, k = 2, log2.thresh = 1) {
  reads.tap = abs(reads.in.steps.tap)
  reads.notap = abs(reads.in.steps.notap) * scale.factor

  ratio = log2(reads.tap / reads.notap)
  
  useful = reads.tap > 0 & reads.notap > 0
  increased = reads.tap > reads.notap

  L = length(reads.tap)

  # create X
  #
  X = rep(1, L) # default no data
  X[increased] = 2 # enriched
  X[useful & ratio < 0] = 3 # depleated (but > 0)

  # create Y
  #
  Y = rep(0, L)
  if (any(increased & useful)) {
    Y[increased & useful] = sapply(which(increased & useful), function(idx) {
      i = max(idx - w, 1)
      j = min(idx + w, L)
      bck.idxs = which(X[i:j] == 3)
      if (length(bck.idxs) >= k) {
        bck = mean(reads.tap[i:j][bck.idxs])
        val = log2(reads.tap[idx] / bck)
        
        if (val > log2.thresh)
          return(1)
      }
      return(0)
    })
  }

  return(rbind(Y, X))
}


#' Write a BED data frame as a prediction track (with colors)
#'
#' @param bed BED data frame
#' @param filename output BED filename
#' @param track track/description string
#' @param color.plus RGB color (comma separated) for plus strand
#' @param color.minus RGB color (comma separated) for minus strand
#' @export
write.track.bed <- function(bed, filename, track="hmm.preds", color.plus = "197,0,11", color.minus = "0,132,209") {
  colors = rep(color.plus, dim(bed)[1])
  colors[bed[,6] == '-'] = color.minus

  bed = cbind(bed[,1], as.integer(bed[,2]), as.integer(bed[,3]), bed[,4:6],
    as.integer(bed[,2]), as.integer(bed[,3]), colors)

  fout = file(filename, "w")
  cat("track name=", track, " description=\"", track, "\"  itemRgb=On\n", sep='', file=fout)
  write.table(bed, file=fout, col.names=F, row.names=F, quote=F, sep='\t')
  close(fout)
}
