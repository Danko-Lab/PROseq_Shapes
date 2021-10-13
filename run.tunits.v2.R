#
# New HMM from work on dense k562 work
#
#

library(proseq.shapes.hmm)
library(tunits)
source("common.R")
source("hmm.prototypes.R")
# dataPath = "/home/paul/WorkingMemory/Box Sync/a4_PROseq_shapes/data"
# dataPath = "/local/workdir/prm88/a4_PROseq_shapes/data"

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


# bigwigPath = paste(dataPath, "/bigwigs/LAB_V1_CNN_V3/bigwigs_all_positions_50bp_chr7_ChROseq_merged", sep='')
# bigwigPath = paste(dataPath, "/bigwigs/LAB_V1_CNN_V3/bigwigs_all_positions_50bp_", chrom, sep='')
# bigwigPath = paste(dataPath, "/bigwigs/LAB_V1_CNN_V3/bigwigs_all_positions_50bp_", chrom, "_preds_on_simulation_filtered_002", sep='')
bigwigPath = paste(dataPath, "/bigwigs/LAB_V2_CNN_V4/bigwigs_all_positions_50bp_", cell_type, "_", chrom, sep='')
# bed.path = paste(bigwigPath, "/", "genestart-thresh001-merge.wig", sep='')
# bed.path = "G1.dREG.peak.full.bed.gz"
# bwPlus.path = paste(dataPath, "/", "GSM1480321_K562_GROcap_wTAP_plus.bigWig", sep='')
# bwMinus.path = paste(dataPath, "/", "GSM1480321_K562_GROcap_wTAP_minus.bigWig", sep='')
# bwPlus.path = paste(dataPath, "/", "G1_plus.bw", sep='')
# bwMinus.path = paste(dataPath, "/", "G1_minus.bw", sep='')

# bwStartPlus.path = paste(bigwigPath, "/", "plus-genestart.bw", sep='')
# bwStartMinus.path = paste(bigwigPath, "/", "minus-genestart.bw", sep='')
bwBodyPlus.path = paste(bigwigPath, "/", "plus-genebody.bw", sep='')
bwBodyMinus.path = paste(bigwigPath, "/", "minus-genebody.bw", sep='')
# bwEndPlus.path = paste(bigwigPath, "/", "plus-geneend.bw", sep='')
# bwEndMinus.path = paste(bigwigPath, "/", "minus-geneend.bw", sep='')
# bwAftergenePlus.path = paste(bigwigPath, "/", "plus-aftergene.bw", sep='')
# bwAftergeneMinus.path = paste(bigwigPath, "/", "minus-aftergene.bw", sep='')
# ref.params.path = "chr7.params.poisson.3000.Rdata"
# ref.params.path = "chr7.params.gamma.3000.Rdata"
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

# all.dset = all.dset.body
# train.dset = train.dataset(all.dset)

# Single set of covars
# covars = dregBED.covar(all.dset, read.table(bed.path))
# saveRDS(covars, file = paste("covars.", chrom, ".G1.dREG.peak.full.rds", sep=''))
# saveRDS(covars, file = paste("covars.", chrom, ".simulated_002.dREG.peak.full.rds", sep=''))
# saveRDS(covars, file = paste("covars.", chrom, ".genestart-thresh001-merge.rds", sep=''))
# saveRDS(covars, file = paste("covars.", chrom, ".genestart-thresh001-merge.simulated_002.rds", sep=''))

# covars = readRDS(file = paste("covars.", chrom, ".G1.dREG.peak.full.rds", sep=''))
# covars = readRDS(file = paste("covars.", chrom, ".simulated_002.dREG.peak.full.rds", sep=''))
# covars = readRDS(file = paste("covars.", chrom, ".genestart-thresh001-merge.rds", sep=''))
# covars = readRDS(file = paste("covars.", chrom, ".genestart-thresh001-merge.simulated_002.rds", sep=''))
# covars.clamp = dreg.clamp(covars)
# covars.clamp.tss = dreg.clamp(covars)

# TSSs from predictions
# covars = dregBED.covar(all.dset, read.table(paste(bigwigPath, "/", "tss-thresh08-merge.wig", sep='')))
# saveRDS(covars, file = paste("covars.", chrom, ".tss-thresh08-merge.rds", sep=''))
# covars = readRDS(file = paste("covars.", chrom, ".tss-thresh08-merge.rds", sep=''))
# covars.clamp.start = dreg.clamp(covars)
# covars.clamp = dreg.clamp(covars)

# Two sets of covars
# covars = dregBED.covar(all.dset, read.table(paste(bigwigPath, "/", "genestart-thresh001-merge.wig", sep='')))
# saveRDS(covars, file = paste("covars.", chrom, ".genestart-thresh001-merge.rds", sep=''))
# covars = readRDS(file = paste("covars.", chrom, ".genestart-thresh001-merge.rds", sep=''))
# covars.clamp.start = dreg.clamp(covars)
# covars.clamp = dreg.clamp(covars)
# covars = dregBED.covar(all.dset, read.table(paste(bigwigPath, "/", "geneend-thresh01-merge.wig", sep='')))
# saveRDS(covars, file = paste("covars.", chrom, "geneend-thresh01-merge.rds", sep=''))
# covars = readRDS(file = paste("covars.", chrom, "aftergene-thresh07-merge.rds", sep=''))
# covars.clamp.end = dreg.clamp(covars)
# covars.clamp[[1]] = rbind(covars.clamp.start[[1]], covars.clamp.end[[1]])
# covars.clamp[[2]] = rbind(covars.clamp.start[[2]], covars.clamp.end[[2]])

# Three sets of covars
# covars = dregBED.covar(all.dset, read.table(paste(bigwigPath, "/", "aftergene-thresh07-merge.wig", sep='')))
# saveRDS(covars, file = paste("covars.", chrom, "aftergene-thresh07-merge.rds", sep=''))
# covars = readRDS(file = paste("covars.", chrom, "aftergene-thresh07-merge.rds", sep=''))
# covars.clamp.after = dreg.clamp(covars)
# covars.clamp[[1]] = rbind(covars.clamp.start[[1]], covars.clamp.end[[1]], covars.clamp.after[[1]])
# covars.clamp[[2]] = rbind(covars.clamp.start[[2]], covars.clamp.end[[2]], covars.clamp.after[[2]])

# Using TSS priors
# covars = dregBED.covar(all.dset, read.table(bed.path))
# saveRDS(covars, file = paste("covars.", chrom, ".G1.dREG.peak.full.rds", sep=''))
# covars = readRDS(file = paste("covars.", chrom, ".G1.dREG.peak.full.rds", sep=''))

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

# covars = dregBED.covar(all.dset, read.table(paste(bigwigPath, "/", "plus-geneend-thresh01-merge.wig", sep='')))
# saveRDS(covars, file = paste("covars.", chrom, "plus-geneend-thresh01-merge.rds", sep=''))
# covars = readRDS(file = paste("covars.", chrom, "plus-geneend-thresh01-merge.rds", sep=''))
# covars.clamp.plus.end = dreg.clamp(covars)
# covars = dregBED.covar(all.dset, read.table(paste(bigwigPath, "/", "minus-geneend-thresh01-merge.wig", sep='')))
# saveRDS(covars, file = paste("covars.", chrom, "minus-geneend-thresh01-merge.rds", sep=''))
# covars = readRDS(file = paste("covars.", chrom, "minus-geneend-thresh01-merge.rds", sep=''))
# covars.clamp.minus.end = dreg.clamp(covars)

# covars.clamp = dreg.clamp(covars)
# covars.clamp[[1]] = rbind(covars.clamp.plus.start[[1]], covars.clamp.plus.end[[1]])
# covars.clamp[[2]] = rbind(covars.clamp.minus.start[[2]], covars.clamp.minus.end[[2]])
# covars.clamp[[1]] = rbind(covars.clamp.tss[[1]], covars.clamp.plus.end[[1]])
# covars.clamp[[2]] = rbind(covars.clamp.tss[[2]], covars.clamp.minus.end[[2]])

# Three sets of covars with strand info
# covars = dregBED.covar(all.dset, read.table(paste(bigwigPath, "/", "plus-aftergene-thresh02-merge.wig", sep='')))
# saveRDS(covars, file = paste("covars.", chrom, "plus-aftergene-thresh02-merge.rds", sep=''))
# covars = readRDS(file = paste("covars.", chrom, "plus-aftergene-thresh02-merge.rds", sep=''))
# covars.clamp.plus.after = dreg.clamp(covars)
# covars = dregBED.covar(all.dset, read.table(paste(bigwigPath, "/", "minus-aftergene-thresh02-merge.wig", sep='')))
# saveRDS(covars, file = paste("covars.", chrom, "minus-aftergene-thresh02-merge.rds", sep=''))
# covars = readRDS(file = paste("covars.", chrom, "minus-aftergene-thresh02-merge.rds", sep=''))
# covars.clamp.minus.after = dreg.clamp(covars)
# covars = dregBED.covar(all.dset, read.table(paste(bigwigPath, "/", "plus-stable-unstable-thresh07-merge.wig", sep='')))
# saveRDS(covars, file = paste("covars.", chrom, "plus-stable-unstable-thresh07-merge.rds", sep=''))
# covars = readRDS(file = paste("covars.", chrom, "plus-stable-unstable-thresh07-merge.rds", sep=''))
# covars.clamp.plus.stable.unstable = dreg.clamp(covars)
# covars = dregBED.covar(all.dset, read.table(paste(bigwigPath, "/", "minus-stable-unstable-thresh07-merge.wig", sep='')))
# saveRDS(covars, file = paste("covars.", chrom, "minus-stable-unstable-thresh07-merge.rds", sep=''))
# covars = readRDS(file = paste("covars.", chrom, "minus-stable-unstable-thresh07-merge.rds", sep=''))
# covars.clamp.minus.stable.unstable = dreg.clamp(covars)

# covars = dregBED.covar(all.dset, read.table(paste(bigwigPath, "/", "plus-unstable-thresh01-merge.wig", sep='')))
# saveRDS(covars, file = paste("covars.", chrom, "plus-unstable-thresh01-merge.rds", sep=''))
# covars = readRDS(file = paste("covars.", chrom, "plus-unstable-thresh01-merge.rds", sep=''))
# covars.clamp.plus.unstable = dreg.clamp(covars)
# covars = dregBED.covar(all.dset, read.table(paste(bigwigPath, "/", "minus-unstable-thresh01-merge.wig", sep='')))
# saveRDS(covars, file = paste("covars.", chrom, "minus-unstable-thresh01-merge.rds", sep=''))
# covars = readRDS(file = paste("covars.", chrom, "minus-unstable-thresh01-merge.rds", sep=''))
# covars.clamp.minus.unstable = dreg.clamp(covars)

# covars.clamp = dreg.clamp(covars)
# covars.clamp[[1]] = rbind(covars.clamp.plus.start[[1]])
# covars.clamp[[2]] = rbind(covars.clamp.minus.start[[2]])

# covars.clamp[[1]] = rbind(covars.clamp.plus.start[[1]], covars.clamp.plus.after[[1]])
# covars.clamp[[2]] = rbind(covars.clamp.minus.start[[2]], covars.clamp.minus.after[[2]])

# covars.clamp[[1]] = rbind(covars.clamp.plus.start[[1]], covars.clamp.plus.end[[1]], covars.clamp.plus.after[[1]])
# covars.clamp[[2]] = rbind(covars.clamp.minus.start[[2]], covars.clamp.minus.end[[2]], covars.clamp.minus.after[[2]])

# covars.clamp[[1]] = rbind(covars.clamp.plus.stable[[1]], covars.clamp.plus.unstable[[1]], covars.clamp.plus.end[[1]])
# covars.clamp[[2]] = rbind(covars.clamp.minus.stable[[2]], covars.clamp.minus.unstable[[2]], covars.clamp.minus.end[[2]])

# Get list of BigWig files
# bwSet = BwSet(basePath = commonDatasetPath)
# plusGenestart.slot = step.bpQuery.bigWig(bwSet$plusGenestart, chrom, step = step, start = NULL, end = NULL, op = "avg")
# plusGenebody.slot = step.bpQuery.bigWig(bwSet$plusGenebody, chrom, step = step, start = NULL, end = NULL, op = "avg")
# plusGeneend.slot = step.bpQuery.bigWig(bwSet$plusGeneend, chrom, step = step, start = NULL, end = NULL, op = "avg")
# minusGenestart.slot = step.bpQuery.bigWig(bwSet$minusGenestart, chrom, step = step, start = NULL, end = NULL, op = "avg")
# minusGenebody.slot = step.bpQuery.bigWig(bwSet$minusGenebody, chrom, step = step, start = NULL, end = NULL, op = "avg")
# minusGeneend.slot = step.bpQuery.bigWig(bwSet$minusGeneend, chrom, step = step, start = NULL, end = NULL, op = "avg")

# all.dset[[1]]$data = plusGenebody.slot[1:3182772]
# all.dset[[2]]$data = minusGenebody.slot[1:3182772]
# all.dset[[1]]$data = rbind(plusGenebody.slot[1:3182772], plusGeneend.slot[1:3182772])
# all.dset[[2]]$data = rbind(minusGenebody.slot[1:3182772], minusGeneend.slot[1:3182772])
# all.dset[[1]]$data = rbind(plusGenestart.slot[1:3182772], plusGenebody.slot[1:3182772])
# all.dset[[2]]$data = rbind(minusGenestart.slot[1:3182772], minusGenebody.slot[1:3182772])
# all.dset[[1]]$data = rbind(plusGenestart.slot[1:3182772], plusGenebody.slot[1:3182772], plusGeneend.slot[1:3182772])
# all.dset[[2]]$data = rbind(minusGenestart.slot[1:3182772], minusGenebody.slot[1:3182772], minusGeneend.slot[1:3182772])


# all.dset[[1]]$data = rbind(all.dset.start[[1]]$data[2,], all.dset.body[[1]]$data[2,])
# all.dset[[2]]$data = rbind(all.dset.start[[2]]$data[2,], all.dset.body[[2]]$data[2,])

# all.dset[[1]]$data = rbind(all.dset.start[[1]]$data[2,], all.dset.body[[1]]$data[2,])
# all.dset[[2]]$data = rbind(all.dset.start[[2]]$data[2,], all.dset.body[[2]]$data[2,])

# all.dset[[1]]$data = rbind(all.dset.start[[1]]$data[2,], all.dset.body[[1]]$data[2,], all.dset.end[[1]]$data[2,])
# all.dset[[2]]$data = rbind(all.dset.start[[2]]$data[2,], all.dset.body[[2]]$data[2,], all.dset.end[[2]]$data[2,])

all.dset[[1]]$data = all.dset.body[[1]]$data[2,]
all.dset[[2]]$data = all.dset.body[[2]]$data[2,]

# all.dset[[1]]$data = rbind(all.dset.body[[1]]$data[2,], all.dset.aftergene[[1]]$data[2,])
# all.dset[[2]]$data = rbind(all.dset.body[[2]]$data[2,], all.dset.aftergene[[2]]$data[2,])

train.dset = train.dataset(all.dset)

# extended
all.dset.c <- all.dset
for (i in 1:length(all.dset)) {
  all.dset.c[[i]]$covar = covars.clamp[[i]]
}

# create HMM instance
# hmm.dreg = splithmm3.hmm(0.153, with.shortcut = TRUE, no.egrps = TRUE, poisson.decay = TRUE, use.negbinom = TRUE)
# hmm.dreg = splithmm4.hmm(with.shortcut = TRUE, no.egrps = FALSE, use.negbinom = TRUE)

# hmmName = "hmm2states1slot1covarContinuous.tsspriors.LAB_V1_CNN_V3"
# hmm.dreg = hmm2states1slot1covarContinuous()
# nStates = 2

# hmmName = "hmm2states2slots1covarContinuous.tsspriors"
# hmm.dreg = hmm2states2slots1covarContinuous()
# nStates = 2

# hmmName = "hmm3states1slot2covarsContinuous.startsAndEndsWStrands.LAB_V1_CNN_V3.gamma"
# hmmName = "hmm3states1slot2covarsContinuous.startsAndAfterWStrands.LAB_V1_CNN_V3.gamma"
# hmmName = "hmm3states1slot2covarsContinuous.startsAndAfterWStrandsAndShortcut.LAB_V1_CNN_V3.gamma"
# hmmName = "hmm3states1slot2covarsContinuous.startsAndEnds.LAB_V2_CNN_V4.poisson"
# hmmName = "hmm3states1slot2covarsContinuous.startsAndAfter.LAB_V2_CNN_V4.gamma"
# hmm.dreg = hmm3states1slot2covarsContinuous()
# nStates = 3

# hmmName = "hmm3states1slot3covarsContinuous.startsEndsAfter.LAB_V1_CNN_V3.gamma"
# hmmName = "hmm3states1slot3covarsContinuous.startsEndsAfterWStrands.LAB_V1_CNN_V3.gamma"
# hmm.dreg = hmm3states1slot3covarsContinuous()
# nStates = 3

# hmmName = "hmm3states1slot1covarContinuous.tsspriors.LAB_V1_CNN_V3.poisson.simulated_002"
# hmmName = "hmm3states1slot1covarContinuous.startpriors.LAB_V1_CNN_V3.poisson"
# hmmName = "hmm3states1slot1covarContinuous.startpriors.LAB_V1_CNN_V3.poisson.ChROseq_merge"
# hmmName = "hmm3states1slot1covarContinuous.startpriors.LAB_V1_CNN_V3.gamma.neg"
# hmmName = "hmm3states1slot1covarContinuous.startpriors.LAB_V1_CNN_V3.gamma"
# hmmName = "hmm3states1slot1covarContinuous.startpriorsWStrands.LAB_V1_CNN_V3.gamma"
# hmmName = "hmm3states1slot1covarContinuous.startpriorsWStrandsAndShortcut.LAB_V1_CNN_V3.gamma"
# hmmName = "hmm3states1slot1covarContinuous.startpriors.LAB_V1_CNN_V3.gamma.simulated_002"
# hmmName = "hmm3states1slot1covarContinuous.tsspriors.LAB_V1_CNN_V3.gamma.simulated_002"
# hmmName = "hmm3states1slot1covarContinuous.startpriors.LAB_V2_CNN_V4.gamma.GAN_input"
# hmmName = "hmm3states1slot1covarContinuous.stable.LAB_V2_CNN_V4.gamma"
# hmmName = "hmm3states1slot1covarContinuous.stable-unstable-thresh07.LAB_V2_CNN_V4.gamma"
hmmName = "hmm3states1slot1covarContinuous.startpriors.LAB_V2_CNN_V4.gamma"
hmm.dreg = hmm3states1slot1covarContinuous()
nStates = 3

# hmmName = "hmm3states2slots1covarContinuous.tsspriors.LAB_V1_CNN_V3"
# hmm.dreg = hmm3states2slots1covarContinuous()
# nStates = 3

# hmmName = "hmm4states1slot3covarsContinuous.stableunstableends.LAB_V2_CNN_V4.gamma"
# hmm.dreg = hmm4states1slot3covarsContinuous()
# nStates = 4

# hmmName = "hmm4states1slot2covarsContinuous.tssends.LAB_V2_CNN_V4.gamma"
# hmm.dreg = hmm4states1slot2covarsContinuous()
# nStates = 4


# hmmName = "hmm4states1slot1covarContinuous.tsspriors.LAB_V2_CNN_V4.gamma"
# hmm.dreg = hmm4states1slot1covarContinuous()
# nStates = 4

#
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

# Evaluate
library(groHMM)
library(rtracklayer)
library(GenomicRanges)

# activeTranscriptsFile = "../../data/labels/LAB_V2/gencode_high_confidence_chr7.bed"
activeTranscriptsFile = paste("../../data/labels/LAB_V2/gencode_high_confidence_", chrom, ".bed", sep='')
# activeTranscriptsFile = "gencode_high_confidence_chr7.bed"
# activeTranscriptsFile = "refGene_high_confidence_10Kb.bed"
# activeTranscriptsFile = "refGene_high_confidence_10Kb_chr7.bed"
# activeTranscriptsFile ="../../data/labels/LAB_V1/refGene_and_GROcap_chr7.bed"
inputFile = outputFile

# Only keep parts of input file that intersect with active transcripts file
commandLine = paste("bedtools intersect -s -a", inputFile, "-b", activeTranscriptsFile, "-wa > inputFile_intersection.bed")
text = system(commandLine, intern=TRUE)

# Try readBed instead and see if I get different results
# activeTranscripts = readBed("K562_active_transcripts_chr7_sorted.bed")
activeTranscripts <- import(activeTranscriptsFile, format = "BED", genome = "hg19")
activeTranscripts$gene_id = activeTranscripts$name
ca_activeTranscripts <- makeConsensusAnnotations(activeTranscripts)

# tunits eval
# splithmm3Transcripts <- import("splithmm3.chr7.preds.bed", format = "BED", genome = "hg19")
# splithmm3Transcripts$gene_id = splithmm3Transcripts$name
# ca_splithmm3Transcripts <- makeConsensusAnnotations(splithmm3Transcripts)
# e <- evaluateHMMInAnnotations(ca_splithmm3Transcripts, ca_activeTranscripts)
# e$eval
predicted.transcripts <- import("inputFile_intersection.bed", format = "BED", genome = "hg19")
predicted.transcripts$gene_id = predicted.transcripts$name
ca_predicted.transcripts <- makeConsensusAnnotations(predicted.transcripts)

e <- evaluateHMMInAnnotations(ca_predicted.transcripts, ca_activeTranscripts)
cat("\nResults for:", inputFile, "\n\n")
# cat("merged disassociated\ttotal\terrorRate\ttxSize\n", as.character(e$eval), "\n")

splitEval = strsplit(as.character(e$eval)," ")
cat("Merged:", splitEval[[1]][1], "\n")
cat("Disassociated:", splitEval[[2]][1], "\n")
cat("Total:", splitEval[[3]][1], "\n")
cat("Error rate:", splitEval[[4]][1], "\n")
cat("txSize:", splitEval[[5]][1], "\n")

# Obtain similarity measure
# commandLine = paste("sort -k1,1 -k2,2n inputFile_intersection.bed > outfile_sorted.bed")
commandLine = paste("sort -k1,1 -k2,2n ", inputFile, " > outfile_sorted.bed")
text = system(commandLine, intern=TRUE)

commandLine = paste("bedtools merge -s -c 6 -o distinct -i ", activeTranscriptsFile, " | awk 'BEGIN { OFS = \"\\t\" } ; {print $1, $2, $3, \".\", \".\", $4}' > activeTranscriptsFileMerged.bed")
text = system(commandLine, intern=TRUE)

commandLine = paste("bedtools jaccard -a outfile_sorted.bed -b activeTranscriptsFileMerged.bed")
jaccardResult = read.table(text = system(commandLine, intern=TRUE))
cat("Jaccard similarity:", as.character(jaccardResult$V3[2]), "\n")

commandLine = paste("echo $(bedtools intersect -s -c -b activeTranscriptsFileMerged.bed -a outfile_sorted.bed | awk 'BEGIN { OFS = \"\\t\" } ; ($10 > 1) {print $0}' | wc -l) | bc -l")
runTogetherResult = read.table(text = system(commandLine, intern=TRUE))
cat("Genes run together:", runTogetherResult[["V1"]][1], "\n")

commandLine = paste("echo $(bedtools intersect -s -c -a activeTranscriptsFileMerged.bed -b outfile_sorted.bed | awk 'BEGIN { OFS = \"\\t\" } ; ($7 > 1) {print $0}' | wc -l) | bc -l")
fragmentationResult = read.table(text = system(commandLine, intern=TRUE))
cat("Genes broken up:", fragmentationResult[["V1"]][1], "\n")

# cat("Jaccard similarity:", as.character(jaccardResult$V3[2]), "\n", "Genes broken up:", fragmentationResult[["V1"]][1], "\n", "Genes run together:", runTogetherResult[["V1"]][1], "\n")
#cat("Genes broken up:", fragmentationResult[["V1"]][1], "\nGenes run together:", runTogetherResult[["V1"]][1], "\n")

density = getTxDensity(ca_predicted.transcripts, ca_activeTranscripts)

cat("density$FivePrimeFP", density$FivePrimeFP, "\n")
cat("density$TP", density$TP, "\n")
cat("density$PostTTS", density$PostTTS, "\n")
cat("density$TUA", density$TUA, "\n")

columnNames = c("input_file", "merged", "disassociated", "total", "error_rate", "tx_size", "jaccard", "run_together", "broken_up", "FivePrimeFP", "TP", "PostTTS", "TUA")
line = c(inputFile, splitEval[[1]][1], splitEval[[2]][1], splitEval[[3]][1], splitEval[[4]][1], splitEval[[5]][1], as.character(jaccardResult$V3[2]), runTogetherResult[["V1"]][1], fragmentationResult[["V1"]][1], density$FivePrimeFP, density$TP, density$PostTTS, density$TUA)
mline = matrix(line, 1, length(line))
write.table(mline, file="evaluation_metrics.csv", row.names=FALSE, col.names=FALSE, sep="\t", eol="\n", quote=FALSE, dec=".", append=TRUE)

cat("\n")

