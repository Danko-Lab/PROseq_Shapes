#
# PRO-seq shapes Transcription Unit HMM
#

#' Create HMM with 11 states and 11 observation slots (this one goes up to 11)
#'
#' @param tm matrix of valid transitions
#' @param free.params.vec vector of free parameters
#' @export
make.qhmm2Shapes11states11slots <- function(tm = NULL, free.params.vec = c((1 - 0.9)/2, (1 - 0.9)/2, 0.9, 0.01, 0.45, 0.45)) {
  full.params.vec <- function(a, b) {
    c(1 - (a + b), a, b)
  }
  
  # Data: (E1, E2, E3, E4, E5, E6)
  #
  # E1: TSS predictions
  # E2: Plus gene start predictions
  # E3: Plus gene body predictions
  # E4: Plus gene end predictions
  # E5: Plus after gene predictions
  # E6: Plus neg predictions
  
  # states: PNEG, TSS, PG.1, PG.2, PG.3, PAG :: total 6

  n.states = 6
  TSS = 1
  PG.S = 2
  PG.B = 3
  PG.E = 4
  PAG = 5
  PNEG = 6
  
  # valid transitions
  tm = matrix(data=0, nrow=n.states, ncol=n.states)
  
  tm[TSS, PNEG] = 1
  # tm[TSS, TSS] = 2
  tm[TSS, PG.S] = 2
  
  # tm[PG.S, PG.S] = 1
  tm[PG.S, PG.B] = 1
  
  tm[PG.B, PG.B] = 1
  tm[PG.B, PG.E] = 2
  
  tm[PG.E, PNEG] = 1
  tm[PG.E, PG.E] = 2
  tm[PG.E, PAG] = 3
  
  tm[PAG, PNEG] = 1
  tm[PAG, PAG] = 2

  tm[PNEG, PNEG] = 1
  tm[PNEG, TSS] = 2
  tm[PNEG, PG.S] = 3
  
  # Emmision groups
  nSlots = 6
  egrp = new.emission.groups(n.states, nSlots)
  egrp = add.emission.groups(egrp, group=c(nSlots, 2:4))
  # egrp = add.emission.groups(egrp, group=c(nSlots, 5:7))
  
  #
  hmm <- new.qhmm(list(rep(1, nSlots), NULL),
                  tm, rep("discrete", n.states),
                  rep(list(rep("discrete", nSlots)), n.states),
                  emission.groups = egrp)
  cat("Transition matrix:", hmm$valid.transitions, "\n")

  # set initial probabilities
  set.initial.probs.qhmm(hmm, c(rep(0, n.states - 1), 1))
  
  # Transition and emission probabilities should be set randomly, 
  # so that we can have multiple start points for the em algoithm
  # set.seed(1)
  
  # set transition parameters
  # randomly generating n numbers that sum to 1 
  # (where n is the number of possible state transitions)
  # https://stackoverflow.com/questions/11003967/generate-3-random-number-that-sum-to-1-in-r
  # rowSums(tm>0)[i] gives the number of transitions for state i
  # diff(c(0, sort(runif(n-1)), 1)) gives n random numbers that sum to 1
  # (this gives an error if n < 1, so makes sure eash state always has at least one transition)
  # diff(c(0, sort(runif((rowSums(tm>0)[i])-1)), 1)) gives us a random number for each state, all adding to 1
  for (i in 1:n.states) {
    set.transition.params.qhmm(hmm, i, diff(c(0, sort(runif((rowSums(tm>0)[i])-1)), 1)))
  }
  # set.transition.params.qhmm(hmm, PNEG, diff(c(0, sort(runif((rowSums(tm>0)[PNEG])-1)), 1)))
  # set.transition.params.qhmm(hmm, TSS, diff(c(0, sort(runif((rowSums(tm>0)[TSS])-1)), 1)))
  # set.transition.params.qhmm(hmm, PG.S, diff(c(0, sort(runif((rowSums(tm>0)[PG.S])-1)), 1)))
  # set.transition.params.qhmm(hmm, PG.B, diff(c(0, sort(runif((rowSums(tm>0)[PG.B])-1)), 1)))
  # set.transition.params.qhmm(hmm, PG.E, diff(c(0, sort(runif((rowSums(tm>0)[PG.E])-1)), 1)))
  # set.transition.params.qhmm(hmm, PAG, diff(c(0, sort(runif((rowSums(tm>0)[PAG])-1)), 1)))
  
  # Set emission parameters & options
  # (initial probabilities based on distribution of predictions in each of the observation slots)
  tssParams = c(0.16, 0.84)
  tssParams.reverse = c(0.84, 0.16)
  set.emission.params.qhmm(hmm, 1:n.states, tssParams.reverse, slot = 1)
  set.emission.params.qhmm(hmm, TSS, tssParams, slot = 1)
  
  plusGenestartParams = c(0.04, 0.96)
  plusGenestartParams.reverse = c(0.96, 0.04)
  set.emission.params.qhmm(hmm, 1:n.states, plusGenestartParams.reverse, slot = 2)
  set.emission.params.qhmm(hmm, PG.S, plusGenestartParams, slot = 2)
  
  plusGenebodyParams = c(0.17, 0.83)
  plusGenebodyParams.reverse = c(0.83, 0.17)
  set.emission.params.qhmm(hmm, 1:n.states, plusGenebodyParams.reverse, slot = 3)
  set.emission.params.qhmm(hmm, PG.B, plusGenebodyParams, slot = 3)
  
  plusGeneendParams = c(0.11, 0.89)
  plusGeneendParams.reverse = c(0.89, 0.11)
  set.emission.params.qhmm(hmm, 1:n.states, plusGeneendParams.reverse, slot = 4)
  set.emission.params.qhmm(hmm, PG.E, plusGeneendParams, slot = 4)
  
  plusAftergene = c(0.12, 0.88)
  plusAftergene.reverse = c(0.88, 0.12)
  set.emission.params.qhmm(hmm, 1:n.states, plusAftergene.reverse, slot = 5)
  set.emission.params.qhmm(hmm, PAG, plusAftergene, slot = 5)
  
  plusNeg = c(0.29, 0.71)
  plusNeg.reverse = c(0.71, 0.29)
  set.emission.params.qhmm(hmm, 1:n.states, plusNeg.reverse, slot = 6)
  set.emission.params.qhmm(hmm, PNEG, plusNeg, slot = 6)
  
  # Alternatively, set emission params randomly
  # (see notes above for what diff(c(0, sort(runif(n - 1)), 1)) is doing)
  nEmissionVals = 2  # number of possible emission values (0 and 1 in this case)
  for (i in 1:n.states) {
    for (slotNum in 1:nSlots) {
      # set.emission.params.qhmm(hmm, i, diff(c(0, sort(runif(nEmissionVals - 1)), 1)), slot = slotNum)
      set.emission.option.qhmm(hmm, i, "offset", 0, slot = slotNum)
    }
  }
  
  return(hmm)
}

#' Create HMM with 4 states and 4 observation slots
#'
#' @param free.params.vec vector of free parameters
#' @export
make.qhmm2Shapes4states5slots <- function(free.params.vec = c((1 - 0.9)/2, (1 - 0.9)/2, 0.9, 0.01, 0.45, 0.45)) {
  full.params.vec <- function(a, b) {
    c(1 - (a + b), a, b)
  }
  
  # Data: (E1, E2, E3, E4, E5)
  #
  # E1: TSS predictions
  # E2: Plus gene start predictions
  # E3: Plus gene body predictions
  # E4: Plus gene end predictions
  # E5: Plus after gene predictions
  
  # states: PG.S, PG.B, PG.E, PNEG :: total 4
  
  n.states = 4
  PG.S = 1
  PG.B = 2
  PG.E = 3
  PNEG = 4
  
  # valid transitions
  # if (is.null(tm)) {
  tm = matrix(data=0, nrow=n.states, ncol=n.states)
  
  tm[PG.S, PG.S] = 1
  tm[PG.S, PG.B] = 2
  
  tm[PG.B, PNEG] = 1
  tm[PG.B, PG.B] = 2
  tm[PG.B, PG.E] = 3
  
  tm[PG.E, PNEG] = 1
  tm[PG.E, PG.E] = 2
  
  tm[PNEG, PNEG] = 1
  tm[PNEG, PG.S] = 2
  
  
  # Emmision groups
  nSlots = 5
  egrp = new.emission.groups(n.states, nSlots)
  egrp = add.emission.groups(egrp, group=c(nSlots, 2:4))
  # egrp = add.emission.groups(egrp, group=c(nSlots, 5:7))
  
  #
  hmm <- new.qhmm(list(rep(1, nSlots), NULL),
                  tm, rep("discrete", n.states),
                  rep(list(rep("discrete", nSlots)), n.states),
                  emission.groups = egrp)
  
  # set initial probabilities
  set.initial.probs.qhmm(hmm, c(rep(0, n.states - 1), 1))
  
  # Transition and emission probabilities should be set randomly, 
  # so that we can have multiple start points for the em algoithm
  # set.seed(1)
  
  # set transition parameters
  # Give a very high probability to self transitions
  set.transition.params.qhmm(hmm, PG.S, c(0.98, 0.02))
  set.transition.params.qhmm(hmm, PG.B, c(0.01, 0.95, 0.04))
  set.transition.params.qhmm(hmm, PG.E, c(0.01, 0.99))
  set.transition.params.qhmm(hmm, PNEG, c(0.97, 0.03))
  
  # randomly generating n numbers that sum to 1 
  # (where n is the number of possible state transitions)
  # https://stackoverflow.com/questions/11003967/generate-3-random-number-that-sum-to-1-in-r
  # rowSums(tm>0)[i] gives the number of transitions for state i
  # diff(c(0, sort(runif(n-1)), 1)) gives n random numbers that sum to 1
  # (this gives an error if n < 1, so makes sure eash state always has at least one transition)
  # diff(c(0, sort(runif((rowSums(tm>0)[i])-1)), 1)) gives us a random number for each state, all adding to 1
  # for (i in 1:n.states) {
  #   set.transition.params.qhmm(hmm, i, diff(c(0, sort(runif((rowSums(tm>0)[i])-1)), 1)))
  # }
  
  # Set emission parameters & options
  # (initial probabilities based on distribution of predictions in each of the observation slots)
  # plusGenestartParams = c(0.04, 0.96)
  # plusGenestartParams.reverse = c(0.96, 0.04)
  # set.emission.params.qhmm(hmm, 1:n.states, plusGenestartParams.reverse, slot = 1)
  # set.emission.params.qhmm(hmm, PG.S, plusGenestartParams, slot = 1)
  
  # plusGenebodyParams = c(0.17, 0.83)
  # plusGenebodyParams.reverse = c(0.83, 0.17)
  # set.emission.params.qhmm(hmm, 1:n.states, plusGenebodyParams.reverse, slot = 2)
  # set.emission.params.qhmm(hmm, PG.B, plusGenebodyParams, slot = 2)
  
  # plusGeneendParams = c(0.11, 0.89)
  # plusGeneendParams.reverse = c(0.89, 0.11)
  # set.emission.params.qhmm(hmm, 1:n.states, plusGeneendParams.reverse, slot = 3)
  # set.emission.params.qhmm(hmm, PG.E, plusGeneendParams, slot = 3)
  
  # plusNeg = c(0.29, 0.71)
  # plusNeg.reverse = c(0.71, 0.29)
  # set.emission.params.qhmm(hmm, 1:n.states, plusNeg.reverse, slot = 4)
  # set.emission.params.qhmm(hmm, PNEG, plusNeg, slot = 4)
  
  # Alternatively, set emission params randomly
  # (see notes above for what diff(c(0, sort(runif(n - 1)), 1)) is doing)
  nEmissionVals = 2  # number of possible emission values (0 and 1 in this case)
  for (i in 1:n.states) {
    for (slotNum in 1:nSlots) {
      set.emission.params.qhmm(hmm, i, diff(c(0, sort(runif(nEmissionVals - 1)), 1)), slot = slotNum)
      set.emission.option.qhmm(hmm, i, "offset", 0, slot = slotNum)
    }
  }
  
  return(hmm)
}

#' Create HMM with 2 states and 5 observation slots
#'
#' @param free.params.vec vector of free parameters
#' @export
make.qhmm2Shapes2states5slots <- function(free.params.vec = c((1 - 0.9)/2, (1 - 0.9)/2, 0.9, 0.01, 0.45, 0.45)) {
  full.params.vec <- function(a, b) {
    c(1 - (a + b), a, b)
  }
  
  # Data: (E1, E2, E3, E4, E5)
  #
  # E1: TSS predictions
  # E2: Plus gene start predictions
  # E3: Plus gene body predictions
  # E4: Plus gene end predictions
  # E5: Plus after gene predictions
  
  # states: PG, PNEG :: total 2
  
  n.states = 2
  PG = 1
  PNEG = 2
  
  # valid transitions
  # if (is.null(tm)) {
  tm = matrix(data=0, nrow=n.states, ncol=n.states)
  
  tm[PG, PG] = 1
  tm[PG, PNEG] = 2
  
  tm[PNEG, PG] = 1
  tm[PNEG, PNEG] = 2
  
  
  # Emmision groups
  nSlots = 5
  # egrp = new.emission.groups(n.states, nSlots)
  # egrp = add.emission.groups(egrp, group=c(nSlots, 1))

  #
  # tnames = c("autocorr_covar", "autocorr")
  # tdist="dgamma"
  # hmm <- new.qhmm(list(rep(1, nSlots), 1),
  #                tm, tnames, 
  #                list(c("geometric", tdist, tdist, tdist, tdist), 
  #                     c("geometric", "poisson", "poisson", "poisson", "poisson")),
  #                emission.groups = egrp)
  
  hmm <- new.qhmm(list(rep(1, nSlots), NULL),
                 tm, 
                 rep("discrete", n.states),
                 rep(list(rep("discrete", nSlots)), n.states))
  
  # set initial probabilities
  set.initial.probs.qhmm(hmm, c(rep(0, n.states - 1), 1))
  
  # Transition and emission probabilities should be set randomly, 
  # so that we can have multiple start points for the em algoithm
  # set.seed(1)
  
  # set transition parameters
  # Give a very high probability to self transitions
  set.transition.params.qhmm(hmm, PG, c(0.98, 0.02))
  set.transition.params.qhmm(hmm, PNEG, c(0.02, 0.98))
  
  # randomly generating n numbers that sum to 1 
  # (where n is the number of possible state transitions)
  # https://stackoverflow.com/questions/11003967/generate-3-random-number-that-sum-to-1-in-r
  # rowSums(tm>0)[i] gives the number of transitions for state i
  # diff(c(0, sort(runif(n-1)), 1)) gives n random numbers that sum to 1
  # (this gives an error if n < 1, so makes sure eash state always has at least one transition)
  # diff(c(0, sort(runif((rowSums(tm>0)[i])-1)), 1)) gives us a random number for each state, all adding to 1
  # for (i in 1:n.states) {
  #   set.transition.params.qhmm(hmm, i, diff(c(0, sort(runif((rowSums(tm>0)[i])-1)), 1)))
  # }
  
  # Set emission parameters & options
  # (initial probabilities based on distribution of predictions in each of the observation slots)
  # plusGenestartParams = c(0.04, 0.96)
  # plusGenestartParams.reverse = c(0.96, 0.04)
  # set.emission.params.qhmm(hmm, 1:n.states, plusGenestartParams.reverse, slot = 1)
  # set.emission.params.qhmm(hmm, PG.S, plusGenestartParams, slot = 1)
  
  # plusGenebodyParams = c(0.17, 0.83)
  # plusGenebodyParams.reverse = c(0.83, 0.17)
  # set.emission.params.qhmm(hmm, 1:n.states, plusGenebodyParams.reverse, slot = 2)
  # set.emission.params.qhmm(hmm, PG.B, plusGenebodyParams, slot = 2)
  
  # plusGeneendParams = c(0.11, 0.89)
  # plusGeneendParams.reverse = c(0.89, 0.11)
  # set.emission.params.qhmm(hmm, 1:n.states, plusGeneendParams.reverse, slot = 3)
  # set.emission.params.qhmm(hmm, PG.E, plusGeneendParams, slot = 3)
  
  # plusNeg = c(0.29, 0.71)
  # plusNeg.reverse = c(0.71, 0.29)
  # set.emission.params.qhmm(hmm, 1:n.states, plusNeg.reverse, slot = 4)
  # set.emission.params.qhmm(hmm, PNEG, plusNeg, slot = 4)
  
  # Alternatively, set emission params randomly
  # (see notes above for what diff(c(0, sort(runif(n - 1)), 1)) is doing)
  nEmissionVals = 2  # number of possible emission values (0 and 1 in this case)
  for (i in 1:n.states) {
    for (slotNum in 1:nSlots) {
      set.emission.params.qhmm(hmm, i, diff(c(0, sort(runif(nEmissionVals - 1)), 1)), slot = slotNum)
      set.emission.option.qhmm(hmm, i, "offset", 0, slot = slotNum)
    }
  }
  
  return(hmm)
}

#' Create HMM with 4 states and 3 observation slots
#'
#' @export
make.qhmm2Shapes4states3slots <- function(free.params.vec = c(
  (1 - 0.9)/2, (1 - 0.9)/2,
  0.9, 0.01,
  0.45, 0.45)) {
  full.params.vec <- function(a, b) {
    c(1 - (a + b), a, b)
  }
  
  # Data: (E1, E2, E3)
  #
  # E1: Plus gene start predictions
  # E2: Plus gene body predictions
  # E3: Plus gene end predictions
  
  # states: B, PG.1, PG.2, PG.3 :: total 4
  #
  
  B = 1
  PG.1 = 2
  PG.2 = 3
  PG.3 = 4
  n.states = 4
  
  # valid transitions
  tm = matrix(data=0, nrow=n.states, ncol=n.states)
  
  tm[B, B] = 1
  tm[B, PG.1] = 2
  tm[B, PG.2] = 3
  # tm[B, MG.1] = 3
  
  tm[PG.1, PG.2] = 1
  
  tm[PG.2, PG.2] = 1
  tm[PG.2, PG.3] = 2
  
  tm[PG.3, B] = 1
  
  # tm[M2.1, M2.1] = 1
  # tm[M2.1, M2.2] = 2
  
  # tm[M2.2, M2.1] = 1
  # tm[M2.2, M2.2] = 2
  # tm[M2.2, M2.3] = 3
  
  # tm[M2.3, M2.3] = 1
  # tm[MG.3, B] = 2
  
  # Emmision groups
  nSlots = 3
  egrp = new.emission.groups(n.states, nSlots)
  egrp = add.emission.groups(egrp, group=c(nSlots, 2:4))
  # egrp = add.emission.groups(egrp, group=c(nSlots, 5:7))
  
  #
  hmm <- new.qhmm(list(rep(1, nSlots), NULL),
                  tm, rep("discrete", n.states),
                  rep(list(rep("discrete", nSlots)), n.states),
                  emission.groups = egrp)
  
  # set initial probabilities
  set.initial.probs.qhmm(hmm, c(1, rep(0, n.states - 1)))
  
  # set transition parameters
  # set transition parameters
  set.transition.params.qhmm(hmm, B, c(0.7, 0.2, 0.1))
  set.transition.params.qhmm(hmm, PG.1, 1)
  set.transition.params.qhmm(hmm, PG.2, c(0.6, 0.4))
  set.transition.params.qhmm(hmm, PG.3, 1)
  # set.transition.params.qhmm(hmm, MG.1, c(0.5, 0.5))
  # set.transition.params.qhmm(hmm, MG.2, c(0.45, 0.1, 0.45))
  # set.transition.params.qhmm(hmm, MG.3, c(0.5, 0.5))
  
  # Set emission parameters & options
  params = c(1, 0)
  params.1 = c(0, 1)
  
  # set.emission.params.qhmm(hmm, 1:n.states, params, slot = 1)
  set.emission.params.qhmm(hmm, B, c(0.8, 0.2), slot = 1)
  set.emission.params.qhmm(hmm, PG.1, c(0.2, 0.8), slot = 1)
  set.emission.params.qhmm(hmm, PG.2, c(0.5, 0.5), slot = 1)
  set.emission.params.qhmm(hmm, PG.3, c(0.6, 0.4), slot = 1)
  
  for (i in 1:n.states)
    set.emission.option.qhmm(hmm, i, "offset", 0, slot = 1)
  
  # set.emission.params.qhmm(hmm, 1:n.states, params, slot = 2)
  # set.emission.params.qhmm(hmm, PG.2, params.1, slot = 2)
  set.emission.params.qhmm(hmm, B, c(0.75, 0.25), slot = 2)
  set.emission.params.qhmm(hmm, PG.1, c(0.3, 0.7), slot = 2)
  set.emission.params.qhmm(hmm, PG.2, c(0.55, 0.45), slot = 2)
  set.emission.params.qhmm(hmm, PG.3, c(0.5, 0.5), slot = 2)
  
  for (i in 1:n.states)
    set.emission.option.qhmm(hmm, i, "offset", 0, slot = 2)
  
  # set.emission.params.qhmm(hmm, 1:n.states, params, slot = 3)
  # set.emission.params.qhmm(hmm, PG.2, params.1, slot = 3)
  set.emission.params.qhmm(hmm, B, c(0.9, 0.1), slot = 3)
  set.emission.params.qhmm(hmm, PG.1, c(0.2, 0.8), slot = 3)
  set.emission.params.qhmm(hmm, PG.2, c(0.35, 0.65), slot = 3)
  set.emission.params.qhmm(hmm, PG.3, c(0.8, 0.2), slot = 3)
  
  for (i in 1:n.states)
    set.emission.option.qhmm(hmm, i, "offset", 0, slot = 3)
  
  return(hmm)
}

#' Create HMM
#'
#' @export
make.qhmm2Old <- function(free.params.vec = c(
  (1 - 0.9)/2, (1 - 0.9)/2,
  0.9, 0.01,
  0.45, 0.45)) {
  full.params.vec <- function(a, b) {
    c(1 - (a + b), a, b)
  }
  
  # Data: (Y, X)
  #
  # Y: 0 :: no peak
  #    1 :: peak
  #
  # X: 1 :: [no signap TAP+ = 0],
  #    2 :: [enriched TAP+ > TAP-],
  #    3 :: [depleated TAP- > TAP+ > 0]
  #    
  
  # states: B, M1.[1..3], M2.1, M2.2, M2.3 :: total 7
  #
  
  B = 1
  M1.1 = 2
  M1.2 = 3
  M1.3 = 4
  M2.1 = 5
  M2.2 = 6
  M2.3 = 7
  n.states = 7
  
  # valid transitions
  tm = matrix(data=0, nrow=n.states, ncol=n.states)
  
  tm[B,B] = 1
  tm[B, M1.1] = 2
  tm[B, M2.1] = 3
  
  tm[M1.1, M1.2] = 1
  
  tm[M1.2, M1.2] = 1
  tm[M1.2, M1.3] = 2
  
  tm[M1.3, B] = 1
  
  tm[M2.1, M2.1] = 1
  tm[M2.1, M2.2] = 2
  
  tm[M2.2, M2.1] = 1
  tm[M2.2, M2.2] = 2
  tm[M2.2, M2.3] = 3
  
  tm[M2.3, M2.3] = 1
  tm[M2.3, B] = 2
  
  #
  egrp = new.emission.groups(n.states, 2)
  egrp = add.emission.groups(egrp, group=c(2, 2:4))
  egrp = add.emission.groups(egrp, group=c(2, 5:7))
  
  #
  hmm <- new.qhmm(list(c(1,1), NULL),
                  tm, rep("discrete", n.states),
                  rep(list(c("discrete", "discrete")), n.states),
                  emission.groups = egrp)
  
  # set initial probabilities
  set.initial.probs.qhmm(hmm, c(1, rep(0, n.states - 1)))
  
  # set transition parameters
  set.transition.params.qhmm(hmm, B, c(0.99, 0.005, 0.005))
  set.transition.params.qhmm(hmm, M1.1, 1)
  set.transition.params.qhmm(hmm, M1.2, c(0.5, 0.5))
  set.transition.params.qhmm(hmm, M1.3, 1)
  set.transition.params.qhmm(hmm, M2.1, c(0.5, 0.5))
  set.transition.params.qhmm(hmm, M2.2, c(0.45, 0.1, 0.45))
  set.transition.params.qhmm(hmm, M2.3, c(0.5, 0.5))
  
  # set emission parameters & options
  
  # . Y
  Y.params = c(1, 0)
  Y.params.1 = c(0, 1)
  set.emission.params.qhmm(hmm, 1:n.states, Y.params, slot = 1, fixed = c(T, T))
  set.emission.params.qhmm(hmm, M2.2, Y.params.1, slot = 1, fixed = c(T, T))
  
  for (i in 1:n.states)
    set.emission.option.qhmm(hmm, i, "offset", 0, slot = 1)
  
  # . X
  X.B.params = full.params.vec(free.params.vec[1], free.params.vec[2])
  X.M1.params = full.params.vec(free.params.vec[3], free.params.vec[4])
  X.M2.params = full.params.vec(free.params.vec[5], free.params.vec[6])
  
  set.emission.params.qhmm(hmm, B, X.B.params, slot = 2)
  set.emission.params.qhmm(hmm, M1.1:M1.3, X.M1.params, slot = 2)
  set.emission.params.qhmm(hmm, M2.1:M2.3, X.M2.params, slot = 2)
  
  return(hmm)
}

#' Decode HMM predictions
#'
decode.data.qhmm2Shapes <- function (hmm, data, start = 0, step = 50, 
                                     TSS.range = NA, PG.S.range = NA, PG.B.range = NA, PG.E.range = NA, 
                                     PAG.range = NA, PNEG.range = NA, PG.range = NA) {
  
  # browser()
  path = viterbi.qhmm(hmm, data)
  
  PG.blocks = NULL
  TSS.blocks = NULL
  PG.S.blocks = NULL
  PG.B.blocks = NULL
  PG.E.blocks = NULL
  PAG.blocks = NULL
  PNEG.blocks = NULL
  
  if (is.na(PG.range)) {
    if (!is.na(PG.S.range)) {
      PG.S.blocks = path.blocks.qhmm(path, PG.S.range)
    }
    if (!is.na(PG.B.range)) {
      PG.B.blocks = path.blocks.qhmm(path, PG.B.range)
    }
    if (!is.na(PG.E.range)) {
      PG.E.blocks = path.blocks.qhmm(path, PG.E.range)
    }
  } else {
    PG.blocks = path.blocks.qhmm(path, PG.range)
  }
  
  if (!is.na(TSS.range)) {
    TSS.blocks = path.blocks.qhmm(path, TSS.range)
  }
  if (!is.na(PAG.range)) {
    PAG.blocks = path.blocks.qhmm(path, PAG.range)
  }
  PNEG.blocks = path.blocks.qhmm(path, PNEG.range)

  # starts = as.integer((c(PG.blocks[1, ]) - 1) * step + start)
  # ends = as.integer(c(PG.blocks[2, ]) * step + start)
  starts = as.integer((c(PG.blocks[1, ], TSS.blocks[1, ], PG.S.blocks[1, ], PG.B.blocks[1, ], PG.E.blocks[1, ], PAG.blocks[1, ], PNEG.blocks[1, ]) - 1) * step + start)
  ends = as.integer(c(PG.blocks[2, ], TSS.blocks[2, ], PG.S.blocks[2, ], PG.B.blocks[2, ], PG.E.blocks[2, ], PAG.blocks[2, ], PNEG.blocks[2, ]) * step + start)
  
  PG.dim = dim(PG.blocks)[2]
  if (is.null(PG.blocks))
    PG.dim = 0
  
  TSS.dim = dim(TSS.blocks)[2]
  if (is.null(TSS.blocks))
    TSS.dim = 0
  
  PG.S.dim = dim(PG.S.blocks)[2]
  if (is.null(PG.S.blocks))
    PG.S.dim = 0
  
  PG.B.dim = dim(PG.B.blocks)[2]
  if (is.null(PG.B.blocks))
    PG.B.dim = 0
  
  PG.E.dim = dim(PG.E.blocks)[2]
  if (is.null(PG.E.blocks))
    PG.E.dim = 0
  
  PAG.dim = dim(PAG.blocks)[2]
  if (is.null(PAG.blocks))
    PAG.dim = 0
  
  PNEG.dim = dim(PNEG.blocks)[2]
  if (is.null(PNEG.blocks))
    PNEG.dim = 0
  
  types = c(rep("PG", PG.dim), rep("TSS", TSS.dim), rep("PG.S", PG.S.dim), rep("PG.B", PG.B.dim), rep("PG.E", PG.E.dim), rep("PAG", PAG.dim), rep("PNEG", PNEG.dim))
  data.frame(starts, ends, types)
}

#' Decode HMM predictions
#'
decode.data.qhmmOld <- function (hmm, data, start = 0, step = 10,
                              m1.range = 2:4, m2.range = NA,
                              m2.start = 5, m2.mid = c(5, 6),
                              m2.end = 7) {
  
  path = viterbi.qhmm(hmm, data)
  m1.blocks = path.blocks.qhmm(path, m1.range)
  
  m2.blocks = NULL
  if (is.na(m2.range))
    m2.blocks = path.blocks2.qhmm(path, m2.start, m2.mid, m2.end)
  else
    m2.blocks = path.blocks.qhmm(path, m2.range)
  
  starts = as.integer((c(m1.blocks[1, ], m2.blocks[1, ]) - 1) * step + start)
  ends = as.integer(c(m1.blocks[2, ], m2.blocks[2, ]) * step + start)
  
  m1.dim = dim(m1.blocks)[2]
  if (is.null(m1.blocks))
    m1.dim = 0
  m2.dim = dim(m2.blocks)[2]
  if (is.null(m2.blocks))
    m2.dim = 0
  types = c(rep("M1", m1.dim), rep("M2", m2.dim))
  data.frame(starts, ends, types)
}

#' Process chromosome with HMM
#'
process.chromosome.qhmm <- function (bwSet, chrom, scale.factor, step = 50, hmm = make.qhmm2Shapes2states5slots(), 
                                     covar.lst = NULL, enable.EM = FALSE, thresh = 1, free.param.lst = NULL, testMode = FALSE) {
  decode <- function(hmm, chrom, start, strand, data) {
    tmp = decode.data.qhmm2Shapes(hmm, data, start = start, step = step, TSS.range = NA, PG.S.range = NA, 
                                  PG.B.range = 2, PG.E.range = NA, PAG.range = NA, PNEG.range = 1, PG.range = 2)
    if (dim(tmp)[1] == 0) 
      return(NULL)
    else {
      bed = cbind(chrom, tmp, 0, strand)
      colnames(bed) <- c("chrom", "start", "end", "type", "score", "strand")
      return(bed[bed[,4] == 'PG',])
    }
  }
  decode.peaks <- function(hmm, chrom, start, strand, data) {
    tmp = decode.data.qhmm2Shapes(hmm, data, start = start, step = step, TSS.range = NA, PG.S.range = NA, 
                                  PG.B.range = 2, PG.E.range = NA, PAG.range = NA, PNEG.range = 1, PG.range = NA)
    if (dim(tmp)[1] == 0) 
      return(NULL)
    else {
      bed = cbind(chrom, tmp, 0, strand)
      colnames(bed) <- c("chrom", "start", "end", "type", "score", "strand")
      return(bed)
    }
  }
  combine.strands <- function(plus.lst, minus.lst) {
    if (length(plus.lst) == 1)
      return(rbind(plus.lst[[1]], minus.lst[[1]]))
    return(lapply(1:length(plus.lst), function(idx) {
      rbind(plus.lst[[idx]], minus.lst[[idx]])
    }))
  }

  em.trace.plus = NULL
  em.trace.minus = NULL
  em.params.plus = NULL
  em.params.minus = NULL
  start.params = collect.params.qhmm(hmm)

  N = max(1, length(free.param.lst))
  bed.plus.lst = vector(mode="list", length=N)
  bed.minus.lst = vector(mode="list", length=N)
  bed.peaks.plus.lst = vector(mode="list", length=N)
  bed.peaks.minus.lst = vector(mode="list", length=N)

  cat(chrom, "\n")
  cat(" * data (+)\n")
  # tap.plus = chromStepSum.bigWig(bwSet$GROcap.plus, chrom, step = 10, defaultValue = 0)
  # notap.plus = chromStepSum.bigWig(bwBck$GROcap.plus, chrom, step = 10, defaultValue = 0)
  # res.plus = sequence.to.data(tap.plus, notap.plus, scale.factor, log2.thresh = log2.thresh)
  # tap.plus = NULL
  # notap.plus = NULL
  # gc()
  
  # tss.slot = chromStepSum.bigWig(bwSet$tss, chrom, step = 10, defaultValue = 0)
  # plusGenestart.slot = chromStepSum.bigWig(bwSet$plusGenestart, chrom, step = 10, defaultValue = 0)
  # plusGenebody.slot = chromStepSum.bigWig(bwSet$plusGenebody, chrom, step = 10, defaultValue = 0)
  # plusGeneend.slot = chromStepSum.bigWig(bwSet$plusGeneend, chrom, step = 10, defaultValue = 0)
  # plusAftergene.slot = chromStepSum.bigWig(bwSet$plusAftergene, chrom, step = 10, defaultValue = 0)
  # plusNeg.slot = chromStepSum.bigWig(bwSet$plusNeg, chrom, step = 10, defaultValue = 0)
  
  # plusGenebody.slot = step.bpQuery.bigWig(bwSet$plusGenebody, chrom, step = 50, start = 10125, end = 159115875, op = "max")
  # plusGenebody.slot = step.bpQuery.bigWig(bwSet$plusGenebody, chrom, step = step, start = NULL, end = NULL, op = "max")
  tss.slot = step.bpQuery.bigWig(bwSet$tss, chrom, step = step, start = NULL, end = NULL, op = "avg")
  plusGenestart.slot = step.bpQuery.bigWig(bwSet$plusGenestart, chrom, step = step, start = NULL, end = NULL, op = "avg")
  plusGenebody.slot = step.bpQuery.bigWig(bwSet$plusGenebody, chrom, step = step, start = NULL, end = NULL, op = "avg")
  testplusGenebody.slot = step.bpQuery.bigWig(bwSet$testplusGenebody, chrom, step = step, start = NULL, end = NULL, op = "avg")
  plusGeneend.slot = step.bpQuery.bigWig(bwSet$plusGeneend, chrom, step = step, start = NULL, end = NULL, op = "avg")
  # plusAftergene.slot = step.bpQuery.bigWig(bwSet$plusAftergene, chrom, step = step, start = NULL, end = NULL, op = "avg")
  
  combinedData = plusGenebody.slot
  # combinedData = rbind(plusGenestart.slot, plusGenebody.slot)
  # combinedData = rbind(plusGenestart.slot, plusGenebody.slot, plusGeneend.slot)
  
  # combinedData = plusGenestart.slot
  # combinedData = testplusGenebody.slot
  # combinedData = processReads(plusGenebody.slot, scale.factor, thresh, testMode)
  
  # processed.plusGenestart.slot = processReads(plusGenestart.slot, scale.factor, thresh, testMode)
  # processed.plusGenebody.slot = processReads(plusGenebody.slot, scale.factor, thresh, testMode)
  # combinedData = rbind(processed.plusGenestart.slot, processed.plusGenebody.slot)
  
  # data.lst = list(plusGenestart.slot, plusGenebody.slot, plusGeneend.slot)
  # combinedData = multi.sequence.to.data(data.lst, scale.factor, thresh = thresh)
  
  # combinedData = combine2.sequence.to.data(plusGenestart.slot, plusGenebody.slot, 
  #                                         scale.factor, thresh = thresh)
  
  # combinedData = combine3.sequence.to.data(plusGenestart.slot, plusGenebody.slot, 
  #                                         plusGeneend.slot, scale.factor, thresh = thresh)
  
  # combinedData = combine4.sequence.to.data(plusGenestart.slot, plusGenebody.slot, 
  #                                         plusGeneend.slot, plusNeg.slot, 
  #                                         scale.factor, thresh = thresh, testMode = testMode)
  
  # combinedData = combine5.sequence.to.data(tss.slot, plusGenestart.slot, plusGenebody.slot, 
  #                                        plusGeneend.slot, plusAftergene.slot, 
  #                                        scale.factor, thresh = thresh, testMode = testMode)
  
  # combinedData = combine6.sequence.to.data(tss.slot, plusGenestart.slot, plusGenebody.slot, 
  #                                          plusGeneend.slot, plusAftergene.slot, plusNeg.slot, 
  #                                          scale.factor, thresh = thresh, testMode = testMode)
  
  # browser()
  plusGenestart.slot = NULL
  plusGenebody.slot = NULL
  testplusGenebody.slot = NULL
  plusGeneend.slot = NULL
  plusAftergene.slot = NULL
  plusNeg.slot = NULL
  gc()
  
  cat(" * process (+)\n")
  if (is.null(free.param.lst)) {
    if (enable.EM) {
      if (is.null(covar.lst)) {
        em.trace.plus = em.qhmm(hmm, list(combinedData))
      } else {
        em.trace.plus = em.qhmm(hmm, list(combinedData), covar.lst = covar.lst)
      }
      em.params.plus = collect.params.qhmm(hmm)
    }
    bed.plus = decode(hmm, chrom, 0, "+", combinedData)
    bed.peaks.plus = decode.peaks(hmm, chrom, 0, "+", combinedData)
    combinedData = NULL
    if (enable.EM)
      restore.params.qhmm(hmm, start.params)
    gc()

    # store in list
    bed.plus.lst[[1]] = bed.plus
    bed.peaks.plus.lst[[1]] = bed.peaks.plus
  } else {
    for (i in 1:N) {
      # make new HMM
      free.params.i = free.param.lst[[i]]
      hmm = make.qhmm2Shapes2states5slots(free.params.i)

      # parse data
      bed.plus = decode(hmm, chrom, 0, "+", combinedData)
      bed.peaks.plus = decode.peaks(hmm, chrom, 0, "+", combinedData)

      # store results
      bed.plus.lst[[i]] = bed.plus
      bed.peaks.plus.lst[[i]] = bed.peaks.plus
    }
    # clean up
    combinedData = NULL
    gc()
  }
  
  cat(" * data (-)\n")
  # minusGenestart.slot = chromStepSum.bigWig(bwSet$minusGenestart, chrom, step = 10, defaultValue = 0)
  # minusGenebody.slot = chromStepSum.bigWig(bwSet$minusGenebody, chrom, step = 10, defaultValue = 0)
  # minusGeneend.slot = chromStepSum.bigWig(bwSet$minusGeneend, chrom, step = 10, defaultValue = 0)
  # minusAftergene.slot = chromStepSum.bigWig(bwSet$minusAftergene, chrom, step = 10, defaultValue = 0)
  # minusNeg.slot = chromStepSum.bigWig(bwSet$minusNeg, chrom, step = 10, defaultValue = 0)

  # minusGenebody.slot = step.bpQuery.bigWig(bwSet$minusGenebody, chrom, step = 50, start = 10125, end = 159115875, op = "max")
  # minusGenebody.slot = step.bpQuery.bigWig(bwSet$minusGenebody, chrom, step = step, start = NULL, end = NULL, op = "max")
  minusGenestart.slot = step.bpQuery.bigWig(bwSet$minusGenestart, chrom, step = step, start = NULL, end = NULL, op = "avg")
  minusGenebody.slot = step.bpQuery.bigWig(bwSet$minusGenebody, chrom, step = step, start = NULL, end = NULL, op = "avg")
  testminusGenebody.slot = step.bpQuery.bigWig(bwSet$testminusGenebody, chrom, step = step, start = NULL, end = NULL, op = "avg")
  minusGeneend.slot = step.bpQuery.bigWig(bwSet$minusGeneend, chrom, step = step, start = NULL, end = NULL, op = "avg")
  # minusAftergene.slot = step.bpQuery.bigWig(bwSet$minusAftergene, chrom, step = step, start = NULL, end = NULL, op = "avg")
  
  combinedData = minusGenebody.slot
  # combinedData = rbind(minusGenestart.slot, minusGenebody.slot)
  # combinedData = rbind(minusGenestart.slot, minusGenebody.slot, minusGeneend.slot)
  
  # combinedData = minusGenestart.slot
  # combinedData = testminusGenebody.slot
  # combinedData = processReads(minusGenebody.slot, scale.factor, thresh, testMode)
  
  # processed.minusGenestart.slot = processReads(minusGenestart.slot, scale.factor, thresh, testMode)
  # processed.minusGenebody.slot = processReads(minusGenebody.slot, scale.factor, thresh, testMode)
  # combinedData = rbind(processed.minusGenestart.slot, processed.minusGenebody.slot)
  
  # combinedData = combine2.sequence.to.data(minusGenestart.slot, minusGenebody.slot, 
  #                                         scale.factor, thresh = thresh)
  
  # combinedData = combine3.sequence.to.data(minusGenestart.slot, minusGenebody.slot, 
  #                                         minusGeneend.slot, scale.factor, thresh = thresh)
  
  # combinedData = combine5.sequence.to.data(tss.slot, minusGenestart.slot, minusGenebody.slot, 
  #                                         minusGeneend.slot, minusAftergene.slot, 
  #                                         scale.factor, thresh = thresh, testMode = testMode)
  
  tss.slot = NULL
  minusGenestart.slot = NULL
  minusGenebody.slot = NULL
  testminusGenebody.slot = NULL
  minusGeneend.slot = NULL
  minusAftergene.slot = NULL
  minusNeg.slot = NULL
  gc()
  
  cat(" * process (-)\n")
  if (is.null(free.param.lst)) {
    if (enable.EM) {
      if (is.null(covar.lst)) {
        em.trace.minus = em.qhmm(hmm, list(combinedData))
      } else {
        em.trace.minus = em.qhmm(hmm, list(combinedData), covar.lst = covar.lst)
      }
      em.params.minus = collect.params.qhmm(hmm)
    }
    bed.minus = decode(hmm, chrom, 0, "-", combinedData)
    bed.peaks.minus = decode.peaks(hmm, chrom, 0, "-", combinedData)
    combinedData = NULL
    if (enable.EM)
      restore.params.qhmm(hmm, start.params)
    gc()

    # store in list
    bed.minus.lst[[1]] = bed.minus
    bed.peaks.minus.lst[[1]] = bed.peaks.minus
  } else {
    for (i in 1:N) {
      # make new HMM
      free.params.i = free.param.lst[[i]]
      hmm = make.qhmm2Shapes2states5slots(free.params.i)
      
      # parse data
      bed.minus = decode(hmm, chrom, 0, "-", combinedData)
      bed.peaks.minus = decode.peaks(hmm, chrom, 0, "-", combinedData)

      # store results
      bed.minus.lst[[i]] = bed.minus
      bed.peaks.minus.lst[[i]] = bed.peaks.minus
    }
    # clean up
    combinedData = NULL
    gc()
  }

  #
  # prepare final result
  #
  preds = combine.strands(bed.plus.lst, bed.minus.lst)
  peaks = combine.strands(bed.peaks.plus.lst, bed.peaks.minus.lst)
  
  if (!enable.EM)
    return(list(preds = preds, peaks = peaks))
  else
    return(list(preds = preds, peaks = peaks,
  	       em = list(trace.plus = em.trace.plus,
                     trace.minus = em.trace.minus,
   	  	   params.plus = em.params.plus,
    		   params.minus = em.params.minus)))
  
  # browser()
  # if (length(bed.plus.lst) == 1)
  #   preds = bed.plus.lst[[1]]
  # else {
  #   preds = lapply(1:length(bed.plus.lst), function(idx) {
  #   rbind(bed.plus.lst[[idx]])
  #   })
  # }
  
  # if (length(bed.peaks.plus.lst) == 1)
  #   peaks = bed.peaks.plus.lst[[1]]
  # else {
  #   peaks = lapply(1:length(bed.peaks.plus.lst), function(idx) {
  #   rbind(bed.peaks.plus.lst[[idx]])
  #   })
  # }
  
  # if (!enable.EM)
  #   return(list(preds = preds, peaks = peaks))
  # else
  #   return(list(preds = preds, peaks = peaks,
  # 	        em = list(trace.plus = em.trace.plus,
  #             	  	  params.plus = em.params.plus)))
}

#' Process entire genome with HMM
#'
#' @export
process.genome.qhmm <- function (bwSet, chroms, scale.factor, step = 50, hmm = make.qhmm2Shapes2states5slots(), 
                                 covar.lst = NULL, enable.EM = FALSE, thresh = 1, free.param.lst = NULL, testMode = FALSE) {
  combine.chroms <- function(curSet, new.lst) {
    if (is.null(curSet))
      return(new.lst)
    else {
      N = length(curSet)
      return(lapply(1:N, function(idx) {
        rbind(curSet[[idx]], new.lst[[idx]])
      }))
    }
  }
  
  bed = NULL
  peaks = NULL
  
  em.data = NULL

  for (chrom in chroms) {
    res.chrom = process.chromosome.qhmm(bwSet, chrom, scale.factor, step, hmm = hmm, covar.lst = covar.lst, 
      enable.EM = enable.EM, thresh = thresh, free.param.lst = free.param.lst, testMode = testMode)

    if (is.null(free.param.lst)) {
      bed.chrom = res.chrom$preds
      bed.peaks.chrom = res.chrom$peaks
    
      bed = rbind(bed, bed.chrom)
      peaks = rbind(peaks, bed.peaks.chrom)
      if (enable.EM) {
        em.data = c(em.data, list(res.chrom$em))
      }
    } else {
      bed.chrom.lst = res.chrom$preds
      bed.peaks.chrom.lst = res.chrom$peaks

      bed = combine.chroms(bed, bed.chrom.lst)
      peaks = combine.chroms(peaks, bed.peaks.chrom.lst)
    }
  }
  if (!enable.EM)
    return(list(preds = bed, peaks = peaks))
  else {
    names(em.data) <- chroms
    return(list(preds = bed, peaks = peaks, em = em.data))
  }
}
