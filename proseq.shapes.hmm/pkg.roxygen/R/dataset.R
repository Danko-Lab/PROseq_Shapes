#
# BigWig Datasets
#

#' Load the main datasets
#'
#' @param basePath Root directory for bigWig files
#' @param fileName Directory name for cell specific files
#' @return List of bigWig objects
#' @export
BwSet <- function(basePath = "/usr/data/", fileName) {

  # load bigWigs
  path.merge <- function(basePath, fileName) {
    paste(basePath, "/", fileName, sep='')
  }
  
  list(  
       tss = load.bigWig(path.merge(basePath, "tss.bw")),
       plusGenestart = load.bigWig(path.merge(basePath, "plus-genestart.bw")),
       minusGenestart = load.bigWig(path.merge(basePath, "minus-genestart.bw")),
       plusGenebody = load.bigWig(path.merge(basePath, "plus-genebody.bw")),
       minusGenebody = load.bigWig(path.merge(basePath, "minus-genebody.bw")),
       plusGeneend = load.bigWig(path.merge(basePath, "plus-geneend.bw")),
       minusGeneend = load.bigWig(path.merge(basePath, "minus-geneend.bw")),
       plusAftergene = load.bigWig(path.merge(basePath, "plus-aftergene.bw")),
       minusAftergene = load.bigWig(path.merge(basePath, "minus-aftergene.bw")),
       plusNeg = load.bigWig(path.merge(basePath, "plus-neg.bw")),
       minusNeg = load.bigWig(path.merge(basePath, "minus-neg.bw")),
       genestart = load.bigWig(path.merge(basePath, "genestart.bw")),
       genebody = load.bigWig(path.merge(basePath, "genebody.bw")),
       geneend = load.bigWig(path.merge(basePath, "geneend.bw")),
       aftergene = load.bigWig(path.merge(basePath, "aftergene.bw")),
       testplusGenebody = load.bigWig(path.merge(basePath, "test-plus-genebody.bw")),
       testminusGenebody = load.bigWig(path.merge(basePath, "test-minus-genebody.bw"))
  )
}
