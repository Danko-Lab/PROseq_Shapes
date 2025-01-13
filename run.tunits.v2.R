
library(tunits)
source("hmm.prototypes.R")

step = 50

args <- commandArgs(trailingOnly=TRUE)
chrom = args[1]
cell_type = args[2]
dataPath = args[3]

seedVal = 3000
set.seed(seedVal)

# chrom = "chr7"
# chrom=c("chr1",  "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17","chr18", "chr19", "chr2",  "chr3",  "chr4" , "chr5" , "chr6" , "chr7" , "chr8", "chr9", "chrX")
# chrom=c("chr20", "chr21")

bigwigPath = paste(dataPath, "/bigwigs/LAB_V2_CNN_V4/bigwigs_all_positions_50bp_", cell_type, "_", chrom, sep='')
bwBodyPlus.path = paste(bigwigPath, "/", "plus-genebody.bw", sep='')
bwBodyMinus.path = paste(bigwigPath, "/", "minus-genebody.bw", sep='')
ref.params.path = NA

#
# prepare data
#
dregBED.covar <- function(dataset, dreg.bed) {
  N = length(dataset)
  covars = vector(mode="list", length=N)
  
  for (i in 1:N) {
    chrom = dataset[[i]]$chrom
    positions = dataset[[i]]$positions
    step = dataset[[i]]$step
    
    bed.i = dreg.bed[dreg.bed[,1] == chrom,]
    
    covars.i = 1 - sapply(positions, function(pos) {
      pos = (pos - 1) * step
      idx = which(bed.i[,2] <= pos & bed.i[,3] > pos)
      
      if (length(idx) >= 1)
        return(max(bed.i[idx, 4])) # use 4 for bedgraph. 5 for 6 bed
      return(0)
    })
    
    covars[[i]] = covars.i
  }
  
  return(covars)
}

dreg.clamp <- function(lst) lapply(lst, function(values) {
  1 - pmax(0, pmin(1 - values, 1))
})

all.dset = load.dataset(chrom, list(bw.plus = bwBodyPlus.path, bw.minus = bwBodyMinus.path), step, transform=TRUE)
# all.dset.start = load.dataset(chrom, list(bw.plus = bwStartPlus.path, bw.minus = bwStartMinus.path), step, transform=TRUE)
all.dset.body = load.dataset(chrom, list(bw.plus = bwBodyPlus.path, bw.minus = bwBodyMinus.path), step, transform=TRUE)
# all.dset.aftergene = load.dataset(chrom, list(bw.plus = bwAftergenePlus.path, bw.minus = bwAftergeneMinus.path), step, transform=TRUE)

# One or two sets of covars with strand info
covars = dregBED.covar(all.dset, read.table(paste(bigwigPath, "/", "plus-genestart-thresh01-merge.wig", sep='')))
saveRDS(covars, file = paste("covars.", chrom, ".plus-genestart-thresh01-merge.rds", sep=''))
# covars = readRDS(file = paste("covars.", chrom, ".plus-genestart-thresh01-merge.rds", sep=''))
covars.clamp.plus.start = dreg.clamp(covars)
# covars.clamp = dreg.clamp(covars)
covars = dregBED.covar(all.dset, read.table(paste(bigwigPath, "/", "minus-genestart-thresh01-merge.wig", sep='')))
saveRDS(covars, file = paste("covars.", chrom, ".minus-genestart-thresh01-merge.rds", sep=''))
# covars = readRDS(file = paste("covars.", chrom, ".minus-genestart-thresh01-merge.rds", sep=''))
covars.clamp.minus.start = dreg.clamp(covars)

covars.clamp = dreg.clamp(covars)
covars.clamp[[1]] = rbind(covars.clamp.plus.start[[1]])
covars.clamp[[2]] = rbind(covars.clamp.minus.start[[2]])

all.dset[[1]]$data = all.dset.body[[1]]$data[2,]
all.dset[[2]]$data = all.dset.body[[2]]$data[2,]

train.dset = train.dataset(all.dset)

# extended
all.dset.c <- all.dset
for (i in 1:length(all.dset)) {
  all.dset.c[[i]]$covar = covars.clamp[[i]]
}

hmmName = "hmm3states1slot1covarContinuous.startpriors.LAB_V2_CNN_V4.gamma"
hmm.dreg = hmm3states1slot1covarContinuous()
nStates = 3

# if present, pre-load starting parameters
if (!is.na(ref.params.path)) {
  load(ref.params.path) # defines hmm.params
  restore.params.qhmm(hmm.dreg, hmm.params)
}

# run
trace.dreg = em.qhmm(hmm.dreg, train.dset, covar.lst = covars.clamp, n_threads = 2)
if (length(trace.dreg$loglik) == 2) {
  # give it another change in case some error occured (I should fix this on the QHMM side ...
  trace.dreg = em.qhmm(hmm.dreg, train.dset, covar.lst = covars.clamp, n_threads = 2)
}

# save parameters
hmm.params = collect.params.qhmm(hmm.dreg)
save(hmm.params, file=paste(chrom, ".params.Rdata", sep=''))

# save predictions

outputFilePrefix = paste(hmmName, ".", chrom, ".", seedVal, ".preds", sep='')
outputFile = paste(outputFilePrefix, ".bed", sep='')
outputFilePlus = paste(outputFilePrefix, ".plus.bed", sep='')
outputFileMinus = paste(outputFilePrefix, ".minus.bed", sep='')
outputAftergeneFile = paste(outputFilePrefix, ".aftergene.bed", sep='')
outputAftergeneFilePlus = paste(outputFilePrefix, ".aftergene.plus.bed", sep='')
outputAftergeneFileMinus = paste(outputFilePrefix, ".aftergene.minus.bed", sep='')

# 1. just body
# preds.dreg = decode.dataset(hmm.dreg, all.dset.c, 2:nStates, covar.lst = covars.clamp)
preds.dreg = decode.dataset(hmm.dreg, all.dset.c, 2, covar.lst = covars.clamp)
write.track(preds.dreg, outputFilePrefix, outputFilePrefix, outputFile)

commandLine = paste("cat", outputFile, "| awk 'BEGIN { OFS = \"\t\" } ; ($6 == \"+\") {print $0}' >", outputFilePlus)
text = system(commandLine, intern=TRUE)

commandLine = paste("cat", outputFile, "| awk 'BEGIN { OFS = \"\t\" } ; ($6 == \"-\") {print $0}' >", outputFileMinus)
text = system(commandLine, intern=TRUE)

# 2. after gene
preds.aftergene = decode.dataset(hmm.dreg, all.dset.c, 3, covar.lst = covars.clamp)
write.track(preds.aftergene, outputFilePrefix, outputFilePrefix, outputAftergeneFile)

commandLine = paste("cat", outputAftergeneFile, "| awk 'BEGIN { OFS = \"\t\" } ; ($6 == \"+\") {print $0}' >", outputAftergeneFilePlus)
text = system(commandLine, intern=TRUE)

commandLine = paste("cat", outputAftergeneFile, "| awk 'BEGIN { OFS = \"\t\" } ; ($6 == \"-\") {print $0}' >", outputAftergeneFileMinus)
text = system(commandLine, intern=TRUE)

# 2. full version
# preds.full.dreg = decode.dataset(hmm.dreg, all.dset.c, 2:3, covar.lst = covars.clamp)
# write.track(preds.full.dreg, chrom, chrom, paste(chrom, ".", seedVal, ".preds.full.bed", sep=''))

# 3. extended version
# write.extended.track(preds.full.dreg, preds.dreg, chrom, chrom, paste(chrom, ".", seedVal, ".preds.ext.bed", sep=''))
