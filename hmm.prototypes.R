#
# HMM model
#
library(rqhmm)

# Create HMM with 3 states and 2 observation slots
#
hmm3states2slotsDiscrete <- function() {
  
  # Data:
  #
  # E1: Gene start predictions
  # E2: Gene body predictions
  # states: PNEG, PG.S, PG.B :: total 3
  
  n.states = 3
  PNEG = 1
  PG.S = 2
  PG.B = 3
  
  # valid transitions
  vtbl = NULL
  
  # valid transitions
  vtbl = matrix(data=0, nrow=n.states, ncol=n.states)
  
  vtbl[PNEG, PNEG] = 1
  vtbl[PNEG, PG.S] = 2
  
  vtbl[PG.S, PG.S] = 1
  vtbl[PG.S, PG.B] = 2
  
  vtbl[PG.B, PG.B] = 1
  vtbl[PG.B, PNEG] = 2
  
  # transition groups
  # tgrps = list(2:n.states)  # share size distribution; it's a poor fit for all, but want to avoid getting a tiny 'T' state
  # tgrps = list(c(2,4)) # must exclude PG.B=3 state as it has a different number of transitions
  tgrps = NULL
  
  # emission groups
  nSlots = 2
  # egrps = NULL
  egrps = new.emission.groups(n.states, nSlots)
  egrps = add.emission.groups(egrp, group=c(nSlots, 2:3))
  # egrps = add.emission.groups(egrps, states = c(2, 3, 4), slots = c(2, 2, 2)) # share GROseq over T, P and D [scaled]
  
  tdist = "discrete"

  hmm <- new.qhmm(list(rep(1, nSlots), NULL),
                  vtbl, 
                  rep("discrete", n.states),
                  rep(list(rep(tdist, nSlots)), n.states),
                  transition.groups = tgrps,
                  emission.groups = egrps)
  
  # set initial parameters
  set.initial.probs.qhmm(hmm, c(1, rep(0, n.states - 1))) # start with background
  
  # Transition and emission probabilities should be set randomly, 
  # so that we can have multiple start points for the em algoithm
  # set.seed(1)
  
  # set transition parameters
  # Give a very high probability to self transitions
  # set.transition.params.qhmm(hmm, PNEG, c(0.98, 0.02))
  # set.transition.params.qhmm(hmm, PG.B, c(0.98, 0.02))
  
  # randomly generating n numbers that sum to 1 
  # (where n is the number of possible state transitions)
  # https://stackoverflow.com/questions/11003967/generate-3-random-number-that-sum-to-1-in-r
  # rowSums(tm>0)[i] gives the number of transitions for state i
  # diff(c(0, sort(runif(n-1)), 1)) gives n random numbers that sum to 1
  # (this gives an error if n < 1, so makes sure eash state always has at least one transition)
  # diff(c(0, sort(runif((rowSums(tm>0)[i])-1)), 1)) gives us a random number for each state, all adding to 1
  for (i in 1:n.states) {
    set.transition.params.qhmm(hmm, i, diff(c(0, sort(runif((rowSums(vtbl>0)[i])-1)), 1)))
  }
  
  # Set emission parameters & options
  # (initial probabilities based on distribution of predictions in each of the observation slots)
  # plusGenestartParams = c(0.04, 0.96)
  # plusGenestartParams.reverse = c(0.96, 0.04)
  # set.emission.params.qhmm(hmm, 1:n.states, plusGenestartParams.reverse, slot = 1)
  # set.emission.params.qhmm(hmm, PG.S, plusGenestartParams, slot = 1)
  
  # plusGenebodyParams = c(0.17, 0.83)
  # plusGenebodyParams.reverse = c(0.83, 0.17)
  # set.emission.params.qhmm(hmm, PNEG, plusGenebodyParams.reverse, slot = 1)
  # set.emission.params.qhmm(hmm, PG.B, plusGenebodyParams, slot = 1)
  
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
      # set.emission.params.qhmm(hmm, i, c(0.1, 0.1), slot = slotNum) # dgamma params
      set.emission.params.qhmm(hmm, i, diff(c(0, sort(runif(nEmissionVals - 1)), 1)), slot = slotNum)
      set.emission.option.qhmm(hmm, i, "offset", 0, slot = slotNum)
    }
  }
  
  return(hmm)
}

# Create HMM with 2 states and 1 observation slot
#
hmm2states1slotDiscrete <- function() {
  
  # Data:
  #
  # E1: Gene body predictions
  # states: PG.B, PNEG :: total 2
  
  n.states = 2
  PNEG = 1
  PG.B = 2
  
  # valid transitions
  vtbl = NULL
  
  # valid transitions
  vtbl = matrix(data=0, nrow=n.states, ncol=n.states)
  
  vtbl[PNEG, PNEG] = 1
  vtbl[PNEG, PG.B] = 2
  
  vtbl[PG.B, PG.B] = 1
  vtbl[PG.B, PNEG] = 2
  
  # transition groups
  # tgrps = list(2:n.states)  # share size distribution; it's a poor fit for all, but want to avoid getting a tiny 'T' state
  # tgrps = list(c(2,4)) # must exclude PG.B=3 state as it has a different number of transitions
  tgrps = NULL
  
  # emission groups
  nSlots = 1
  # egrps = new.emission.groups(n.states, nSlots)
  # egrps = add.emission.groups(egrps, states = c(2, 3, 4), slots = c(2, 2, 2)) # share GROseq over T, P and D [scaled]
  egrps = NULL
  
  tdist = "discrete"
  
  hmm <- new.qhmm(list(rep(1, nSlots), NULL),
                  vtbl, 
                  rep("discrete", n.states),
                  rep(list(rep(tdist, nSlots)), n.states),
                  transition.groups = tgrps,
                  emission.groups = egrps)
  
  # set initial parameters
  set.initial.probs.qhmm(hmm, c(1, rep(0, n.states - 1))) # start with background
  
  # Transition and emission probabilities should be set randomly, 
  # so that we can have multiple start points for the em algoithm
  # set.seed(1)
  
  # set transition parameters
  # Give a very high probability to self transitions
  # set.transition.params.qhmm(hmm, PNEG, c(0.98, 0.02))
  # set.transition.params.qhmm(hmm, PG.B, c(0.98, 0.02))
  
  # randomly generating n numbers that sum to 1 
  # (where n is the number of possible state transitions)
  # https://stackoverflow.com/questions/11003967/generate-3-random-number-that-sum-to-1-in-r
  # rowSums(tm>0)[i] gives the number of transitions for state i
  # diff(c(0, sort(runif(n-1)), 1)) gives n random numbers that sum to 1
  # (this gives an error if n < 1, so makes sure eash state always has at least one transition)
  # diff(c(0, sort(runif((rowSums(tm>0)[i])-1)), 1)) gives us a random number for each state, all adding to 1
  for (i in 1:n.states) {
    set.transition.params.qhmm(hmm, i, diff(c(0, sort(runif((rowSums(vtbl>0)[i])-1)), 1)))
  }
  
  # Set emission parameters & options
  # (initial probabilities based on distribution of predictions in each of the observation slots)
  # plusGenestartParams = c(0.04, 0.96)
  # plusGenestartParams.reverse = c(0.96, 0.04)
  # set.emission.params.qhmm(hmm, 1:n.states, plusGenestartParams.reverse, slot = 1)
  # set.emission.params.qhmm(hmm, PG.S, plusGenestartParams, slot = 1)
  
  # plusGenebodyParams = c(0.17, 0.83)
  # plusGenebodyParams.reverse = c(0.83, 0.17)
  # set.emission.params.qhmm(hmm, PNEG, plusGenebodyParams.reverse, slot = 1)
  # set.emission.params.qhmm(hmm, PG.B, plusGenebodyParams, slot = 1)
  
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

# NOTE: I would like to do the HMM with a 'scaled' neg. binomial
#       but I have not implemented that yet, so first run will 
#       use a DGamma instead
#


# 1. no RNA info (use TSS prior)
#
# state paths:
#    B -> T -> PP -> D -> B
#           -> B [optional link]
#
#  state  #pro  dist
#    B     [1]   [1]
#    T     [2]   [2]
#   PP     [3]   [2]
#    D     [4]   [2]
#
# #pro: [1] background level (Poisson)
#       [2] body level (NegBinom | DGamma)
#       [3] fixed scale up of body level (NegBinom | DGamma)
#       [4] fixed scale down of body level (NegBinom | DGamma)
#
# dist: [1] background distances (Geom0)
#       [2] body distances (Geom0)
#
splithmm1.hmm <- function(scaleTP, scalePD, with.shortcut = FALSE, use.negbinom = FALSE) {
  N = 4
  
  # B, T, P, D
  
  # valid transitions
  vtbl = NULL
  tnames = c("autocorr_covar", "autocorr", "autocorr", "autocorr")
  
  if (!with.shortcut) {
    #
    # no link from T to B directly
    #
    vtbl = rbind(
      c(1, 2, 0, 0), # B to: B or T
      c(0, 1, 2, 0), # T to: T or P
      c(0, 0, 1, 2), # P to: P or D
      c(2, 0, 0, 1)) # D to: D or B
  } else {
    #
    # no link from T to B directly
    #
    vtbl = rbind(
      c(1, 2, 0, 0), # B to: B or T
      c(3, 1, 2, 0), # T to: T or B or P
      c(0, 0, 1, 2), # P to: P or D
      c(2, 0, 0, 1)) # D to: D or B
  }

  tdist = "dgamma"
  if (use.negbinom)
    tdist = "neg_binomial"

  # transition groups
  tgrps = list(2:N)  # share size distribution; it's a poor fit for all, but want to avoid getting a tiny 'T' state
  if (with.shortcut)
   tgrps = list(3:N) # must exclude T=2 state as it has a different number of transitions

  # emission groups
  egrps = new.emission.groups(N, 2)
  egrps = add.emission.groups(egrps, states = c(2, 3, 4), slots = c(2, 2, 2)) # share GROseq over T, P and D [scaled]

  hmm = new.qhmm(list(c(1, 1), 1), # covar: TSS signal
    vtbl,
    tnames,
    list(c("geometric", "poisson"),  # B
         c("geometric", tdist),      # T
         c("geometric", tdist),      # P
         c("geometric", tdist)),     # D
    transition.groups = tgrps,
    emission.groups = egrps)
  
  # set initial parameters
  set.initial.probs.qhmm(hmm, c(1, 0, 0, 0)) # start with background
  
  # set transitions
  set.transition.params.qhmm(hmm, 2:N, 0.99)
  
  # set emissions
  geom.base = 0
  # background
  set.emission.params.qhmm(hmm, 1, 1/3000, slot = 1)
  set.emission.option.qhmm(hmm, 1, "base", geom.base, slot = 1)
  set.emission.params.qhmm(hmm, 1, 0.1, slot = 2) # lambda
  
  # transcribed
  for (i in 2:N) {
    set.emission.params.qhmm(hmm, i, 1/10, slot = 1)
    set.emission.option.qhmm(hmm, i, "base", geom.base, slot = 1)
  }
  
  if (use.negbinom) {
    set.emission.params.qhmm(hmm, 2:N, c(10, 0.1), slot = 2)
    # set scale factors
    stop("not implemented!")
  } else {
    set.emission.params.qhmm(hmm, 2:N, c(1, 1), slot = 2) # gamma params
    set.emission.option.qhmm(hmm, 3, "scale_private", scaleTP, slot = 2)
    set.emission.option.qhmm(hmm, 4, "scale_private", scalePD, slot = 2)
  }

  return(hmm)
}


# 2. w/ RNA info (use TSS prior)
#    B -> U
#      -> FEXON -> INTRON -> EXON -> INTRON
#                                 -> PPAUSE -> PDECAY -> B
#
#   state  #pro  dist  #rna
#     B     [1]   [1]   [1]
#     U     [2]   [2]   [1]
#   FEXON   [2]   [2]   [2]
#  INTRON   [3]   [2]   [1]
#   EXON    [3]   [2]   [2]
#  PPAUSE   [4]   [2]   [1]
#  PDECAY   [5]   [2]   [1]
#
# #pro: [1] background level (Poisson)
#       [2] start level (NegBinom | DGamma)
#       [3] body level (NegBinom | DGamma)
#       [4] fixed scale up of body level (NegBinom | DGamma)
#       [5] fixed scale down of body level (NegBinom | DGamma)
#
# dist: [1] background distances (Geom0)
#       [2] body distances (Geom0)
#
# #rna: [1] background level (Poisson)
#       [2] exon level (NegBinom | DGamma)
#
splithmm2.hmm <- function(scaleTP, scalePD, use.negbinom = FALSE) {
  N = 7
  
  # B, U, F, I, E, P, D
  
  # valid transitions
  vtbl = NULL
  tnames = c("autocorr_covar", rep("autocorr", N - 1))

  vtbl = rbind(
    c(1, 2, 3, 0, 0, 0, 0), # B to: B or U or F
    c(2, 1, 0, 0, 0, 0, 0), # U to: U or B
    c(0, 0, 1, 2, 0, 0, 0), # F to: F or I
    c(0, 0, 0, 1, 2, 0, 0), # I to: I or E
    c(0, 0, 0, 2, 1, 3, 0), # E to: E or I or P
    c(0, 0, 0, 0, 0, 1, 2), # P to: P or D
    c(2, 0, 0, 0, 0, 0, 1)) # D to: D or B

  tdist = "dgamma"
  if (use.negbinom)
    tdist = "neg_binomial"

  # emission groups
  egrps = new.emission.groups(N, 3)
  # distances
  egrps = add.emission.groups(egrps, states = 2:7, slots = rep(1, 6)) # share distance distribution between all non-background states
  # pro/gro
  egrps = add.emission.groups(egrps, states = 4:7, slots = rep(2, 4)) # share GROseq over E, I, P, D (last two scaled versions)
  egrps = add.emission.groups(egrps, states = 2:3, slots = c(2, 2)) # share GROseq over U and F
  # rna
  egrps = add.emission.groups(egrps, states = c(1, 2, 4, 6, 7), slots = c(3, 3, 3, 3, 3)) # share background RNAseq (B, U, I, P, D)
  egrps = add.emission.groups(egrps, states = c(3, 5), slots = c(3, 3)) # share exon RNAseq (F, E)

  hmm = new.qhmm(list(c(1, 1, 1), 1), # covar: TSS signal
    vtbl,
    tnames,
    list(c("geometric", "poisson", "poisson"),  # B
         c("geometric", tdist, "poisson"),      # U
         c("geometric", tdist, tdist),          # F
         c("geometric", tdist, "poisson"),      # I
         c("geometric", tdist, tdist),          # E
         c("geometric", tdist, "poisson"),      # P
         c("geometric", tdist, "poisson")),     # D
    transition.groups = list(c(3, 5)), # share size distribution between F and E (both exons)
    emission.groups = egrps)
  
  # set initial parameters
  set.initial.probs.qhmm(hmm, c(1, rep(0, N - 1))) # start with background
  
  # set transitions
  set.transition.params.qhmm(hmm, 2:N, 0.99)

  # set emissions
  geom.base = 0
  # background
  set.emission.params.qhmm(hmm, 1, 1/3000, slot = 1)
  set.emission.option.qhmm(hmm, 1, "base", geom.base, slot = 1)
  set.emission.params.qhmm(hmm, 1, 0.1, slot = 2) # lambda
  set.emission.params.qhmm(hmm, 1, 0.1, slot = 3) # lambda
  
  # transcribed
  
  # dist
  for (i in 2:N) {
    set.emission.params.qhmm(hmm, i, 1/10, slot = 1)
    set.emission.option.qhmm(hmm, i, "base", geom.base, slot = 1)
  }
  
  # pro/gro
  if (use.negbinom) {
    set.emission.params.qhmm(hmm, 2:N, c(10, 0.1), slot = 2)
    # set scale factors
    stop("not implemented!")
  } else {
    set.emission.params.qhmm(hmm, 2:N, c(1, 1), slot = 2) # gamma params
    set.emission.option.qhmm(hmm, 6, "scale_private", scaleTP, slot = 2)
    set.emission.option.qhmm(hmm, 7, "scale_private", scalePD, slot = 2)
  }
  
  # rna
  set.emission.params.qhmm(hmm, c(2,4,6,7) , 0.1, slot = 3) # lambda
  
  if (use.negbinom) {
    set.emission.params.qhmm(hmm, c(3, 5), c(10, 0.1), slot = 3)
  } else {
    set.emission.params.qhmm(hmm, c(3, 5), c(1, 1), slot = 3) # gamma params
  }

  return(hmm)
}

#
# 3. three state HMM
#
# B -> T -> D
#
# state paths:
#    B -> T -> D -> B
#           -> B (optional link)
#
#  state  #pro  dist
#    B     [1]   [1]
#    T     [2]   [2]
#    D     [3]   [2]
#
# #pro: [1] background level (Poisson)
#       [2] body level (NegBinom | DGamma)
#       [3] fixed scale down of body level (NegBinom | DGamma)
#
# dist: [1] background distances (Geom0)
#       [2] body distances (Geom0)
#
splithmm3.hmm <- function(scaleTD, with.shortcut = FALSE, no.egrps = FALSE, poisson.decay = FALSE, use.negbinom = FALSE) {
  N = 3
  
  if (poisson.decay & !no.egrps)
    stop("poisson.decay requires no.egrps")
  
  # B, T, D
  
  # valid transitions
  vtbl = NULL
  tnames = c("autocorr_covar", "autocorr", "autocorr")
  
  if (!with.shortcut) {
    #
    # no link from T to B directly
    #
    vtbl = rbind(
      c(1, 2, 0), # B to: B or T
      c(0, 1, 2), # T to: T or D
      c(2, 0, 1)) # D to: D or B
  } else {
    #
    # no link from T to B directly
    #
    vtbl = rbind(
      c(1, 2, 0), # B to: B or T
      c(3, 1, 2), # T to: T or D or B
      c(2, 0, 1)) # D to: D or B
  }

  tdist = "dgamma"
  if (use.negbinom)
    tdist = "neg_binomial"
  ddist = tdist
  if (poisson.decay)
    ddist = "poisson"

  # transition groups
  tgrps = list(2:N)  # share size distribution; it's a poor fit for all, but want to avoid getting a tiny 'T' state
  if (with.shortcut)
   tgrps = NULL # must exclude T=2 state as it has a different number of transitions

  # emission groups
  egrps = NULL
  if (!no.egrps) {
    egrps = new.emission.groups(N, 2)
    egrps = add.emission.groups(egrps, states = c(2, 3), slots = c(2, 2)) # share GROseq over T and D [scaled]
    egrps = add.emission.groups(egrps, states = c(2, 3), slots = c(1, 1)) # share distance over T and D
  }

  hmm = new.qhmm(list(c(1, 1), 1), # covar: TSS signal
    vtbl,
    tnames,
    list(c("geometric", "poisson"),  # B
         c("geometric", tdist),      # T
         c("geometric", ddist)),     # D
    transition.groups = tgrps,
    emission.groups = egrps)
  
  # set initial parameters
  set.initial.probs.qhmm(hmm, c(1, 0, 0)) # start with background
  
  # set transitions
  set.transition.params.qhmm(hmm, 2:N, 0.99)
  
  # set emissions
  geom.base = 0
  # background
  set.emission.params.qhmm(hmm, 1, 1/3000, slot = 1)
  set.emission.option.qhmm(hmm, 1, "base", geom.base, slot = 1)
  set.emission.params.qhmm(hmm, 1, 0.1, slot = 2) # lambda
  
  # transcribed
  for (i in 2:N) {
    set.emission.params.qhmm(hmm, i, 1/10, slot = 1)
    set.emission.option.qhmm(hmm, i, "base", geom.base, slot = 1)
  }
  
  if (use.negbinom) {
    if (!no.egrps) {
      set.emission.params.qhmm(hmm, 2:N, c(10, 0.1), slot = 2)
      # set scale factors
      stop("not implemented!")
    } else {
      if (poisson.decay) {
        set.emission.params.qhmm(hmm, 2, c(10, 0.1), slot = 2)
        set.emission.params.qhmm(hmm, 3, 0.1, slot = 2)
      } else {
        set.emission.params.qhmm(hmm, 2, c(10, 0.1), slot = 2)
        set.emission.params.qhmm(hmm, 2, c(0.1, 0.1), slot = 2)
      }
    }
  } else {
    if (!no.egrps) {
      set.emission.params.qhmm(hmm, 2:N, c(1, 1), slot = 2) # gamma params
      set.emission.option.qhmm(hmm, 3, "scale_private", scaleTD, slot = 2)
    } else {
      if (poisson.decay) {
        set.emission.params.qhmm(hmm, 2, c(1, 1), slot = 2) # gamma params
        set.emission.params.qhmm(hmm, 3, 0.1, slot = 2) # gamma params
      } else {
        set.emission.params.qhmm(hmm, 2, c(1, 1), slot = 2) # gamma params
        set.emission.params.qhmm(hmm, 3, c(1, 0.1), slot = 2) # gamma params 
      }
    }
  }

  return(hmm)
}

#
# 4. five state HMM
#
# B -> I -> P -> T -> D
#
# state paths:
#    B -> I -> P -> T -> D -> B
#           -....-> T -> D -> B [skippable pause]
#                     -> B (optional link)
#
#  state  #pro  dist
#    B     [1]   [1]
#    I     [2]   [2]
#    P     [3]   [2]
#    T     [4]   [3]
#    D     [5]   [4]
#
# #pro: [1] background level (Poisson)
#       [2] initiation level (Poisson)
#       [4] pause level (Poisson)
#       [4] body level (NegBinom | DGamma)
#       [5] decay level (Poisson)
#
# dist: [1] background distances (Geom0)
#       [2] initiation/pause distances (Geom0)
#       [3] body distances (Geom0)
#       [4] decay distances (Geom0)
#
splithmm4.hmm <- function(with.shortcut = TRUE, no.egrps = FALSE, use.negbinom = FALSE) {
  N = 5

  # B, I, P, T, D
  
  # valid transitions
  vtbl = NULL
  tnames = c("autocorr_covar", rep("autocorr", N - 1))
  
  if (!with.shortcut) {
    #
    # no link from T to B directly
    #
    vtbl = rbind(
      c(1, 2, 0, 0, 0), # B to: B or I
      c(0, 1, 2, 3, 0), # I to: I or P or T
      c(0, 0, 1, 2, 0), # P to: P or T
      c(0, 0, 0, 1, 2), # T to: T or D
      c(2, 0, 0, 0, 1)) # D to: D or B
  } else {
    #
    # no link from T to B directly
    #
    vtbl = rbind(
      c(1, 2, 0, 0, 0), # B to: B or I
      c(0, 1, 2, 3, 0), # I to: I or P or T
      c(0, 0, 1, 2, 0), # P to: P or T
      c(3, 0, 0, 1, 2), # T to: T or D or B
      c(2, 0, 0, 0, 1)) # D to: D or B
  }
  
  tdist = "dgamma"
  if (use.negbinom)
    tdist = "neg_binomial"
  ddist = "poisson"

  # transition groups
  
# NOT valid!!
#  tgrps = list(2:3, 4:5)  # share size distribution; it's a poor fit for all, but want to avoid getting a tiny 'T' state
#  if (with.shortcut)
#   tgrps = list(2:3) # must exclude T=4 state as it has a different number of transitions
  tgrps = list(4:5)
  if (with.shortcut)
    tgrps = NULL # must exclude T=4 state as it has a different number of transitions

  # emission groups
  egrps = NULL
  if (!no.egrps) {
    egrps = new.emission.groups(N, 2)
    egrps = add.emission.groups(egrps, states = c(2, 3), slots = c(1, 1)) # share distance over I and P
    egrps = add.emission.groups(egrps, states = c(2, 5), slots = c(2, 2)) # share reads over I and D
  }

  hmm = new.qhmm(list(c(1, 1), 1), # covar: TSS signal
    vtbl,
    tnames,
    list(c("geometric", "poisson"),  # B
         c("geometric", "poisson"),  # I
         c("geometric", "poisson"),  # P
         c("geometric", tdist),      # T
         c("geometric", ddist)),     # D
    transition.groups = tgrps,
    emission.groups = egrps)
  
  # set initial parameters
  set.initial.probs.qhmm(hmm, c(1, rep(0, N - 1))) # start with background

  # set transitions
  set.transition.params.qhmm(hmm, 2:3, 0.50)
  set.transition.params.qhmm(hmm, 4:N, 0.99)
  
  # set emissions
  geom.base = 0
  # background
  set.emission.params.qhmm(hmm, 1, 1/3000, slot = 1)
  set.emission.option.qhmm(hmm, 1, "base", geom.base, slot = 1)
  set.emission.params.qhmm(hmm, 1, 0.1, slot = 2) # lambda
  
  # initiation
  set.emission.params.qhmm(hmm, 2, 1/10, slot = 1)
  set.emission.option.qhmm(hmm, 2, "base", geom.base, slot = 1)
  set.emission.params.qhmm(hmm, 2, 0.2, slot = 2) # lambda
  
  # pause
  set.emission.params.qhmm(hmm, 3, 1/10, slot = 1)
  set.emission.option.qhmm(hmm, 3, "base", geom.base, slot = 1)
  set.emission.params.qhmm(hmm, 3, 10, slot = 2) # lambda
  
  # transcribed
  for (i in 4:N) {
    set.emission.params.qhmm(hmm, i, 1/10, slot = 1)
    set.emission.option.qhmm(hmm, i, "base", geom.base, slot = 1)
  }
  
  if (use.negbinom) {
    set.emission.params.qhmm(hmm, 4, c(10, 0.1), slot = 2)
    set.emission.params.qhmm(hmm, 5, 0.1, slot = 2)
  } else {
    set.emission.params.qhmm(hmm, 4, c(1, 1), slot = 2) # gamma params
    set.emission.params.qhmm(hmm, 5, 0.1, slot = 2) # gamma params
  }

  return(hmm)
}

# 5. three state HMM
#
# B -> T -> D
#
# state paths:
#    B -> T -> D -> B
#           -> B (optional link)
#
#  state  #pro  dist
#    B     [1]   [1]
#    T     [2]   [2]
#    D     [3]   [2]
#
# #pro: [1] background level (Poisson)
#       [2] body level (NegBinom | DGamma)
#       [3] fixed scale down of body level (NegBinom | DGamma)
#
# dist: [1] background distances (Geom0)
#       [2] body distances (Geom0)
#
splithmm5.hmm <- function(scaleTD = 1, with.shortcut = FALSE, no.egrps = FALSE, poisson.decay = FALSE, use.negbinom = FALSE) {
  N = 3
  
  if (poisson.decay & !no.egrps)
    stop("poisson.decay requires no.egrps")
  
  # B, T, D
  
  # valid transitions
  vtbl = NULL
  tnames = c("autocorr_covar", "autocorr", "autocorr")
  
  if (!with.shortcut) {
    #
    # no link from T to B directly
    #
    vtbl = rbind(
      c(1, 2, 0), # B to: B or T
      c(0, 1, 2), # T to: T or D
      c(2, 0, 1)) # D to: D or B
  } else {
    #
    # no link from T to B directly
    #
    vtbl = rbind(
      c(1, 2, 0), # B to: B or T
      c(3, 1, 2), # T to: T or D or B
      c(2, 0, 1)) # D to: D or B
  }
  
  tdist = "dgamma"
  if (use.negbinom)
    tdist = "neg_binomial"
  ddist = tdist
  if (poisson.decay)
    ddist = "poisson"
  
  # transition groups
  tgrps = list(2:N)  # share size distribution; it's a poor fit for all, but want to avoid getting a tiny 'T' state
  if (with.shortcut)
    tgrps = NULL # must exclude T=2 state as it has a different number of transitions
  
  # emission groups
  egrps = NULL
  if (!no.egrps) {
    egrps = new.emission.groups(N, 2)
    egrps = add.emission.groups(egrps, states = c(2, 3), slots = c(2, 2)) # share GROseq over T and D [scaled]
    egrps = add.emission.groups(egrps, states = c(2, 3), slots = c(1, 1)) # share distance over T and D
  }
  
  hmm = new.qhmm(list(c(1, 1), 1), # covar: TSS signal
                 vtbl,
                 tnames,
                 list(c("geometric", "poisson"),  # B
                      c("geometric", tdist),      # T
                      c("geometric", ddist)),     # D
                 transition.groups = tgrps,
                 emission.groups = egrps)
  
  # set initial parameters
  set.initial.probs.qhmm(hmm, c(1, 0, 0)) # start with background
  
  # set transitions
  set.transition.params.qhmm(hmm, 2:N, 0.99)
  
  # set emissions
  geom.base = 0
  # background
  set.emission.params.qhmm(hmm, 1, 1/3000, slot = 1)
  set.emission.option.qhmm(hmm, 1, "base", geom.base, slot = 1)
  set.emission.params.qhmm(hmm, 1, 0.1, slot = 2) # lambda
  
  # transcribed
  for (i in 2:N) {
    set.emission.params.qhmm(hmm, i, 1/10, slot = 1)
    set.emission.option.qhmm(hmm, i, "base", geom.base, slot = 1)
  }
  
  if (use.negbinom) {
    if (!no.egrps) {
      set.emission.params.qhmm(hmm, 2:N, c(10, 0.1), slot = 2)
      # set scale factors
      stop("not implemented!")
    } else {
      if (poisson.decay) {
        set.emission.params.qhmm(hmm, 2, c(10, 0.1), slot = 2)
        set.emission.params.qhmm(hmm, 3, 0.1, slot = 2)
      } else {
        set.emission.params.qhmm(hmm, 2, c(10, 0.1), slot = 2)
        set.emission.params.qhmm(hmm, 2, c(0.1, 0.1), slot = 2)
      }
    }
  } else {
    if (!no.egrps) {
      set.emission.params.qhmm(hmm, 2:N, c(1, 1), slot = 2) # gamma params
      set.emission.option.qhmm(hmm, 3, "scale_private", scaleTD, slot = 2)
    } else {
      if (poisson.decay) {
        set.emission.params.qhmm(hmm, 2, c(1, 1), slot = 2) # gamma params
        set.emission.params.qhmm(hmm, 3, 0.1, slot = 2) # gamma params
      } else {
        set.emission.params.qhmm(hmm, 2, c(1, 1), slot = 2) # gamma params
        set.emission.params.qhmm(hmm, 3, c(1, 0.1), slot = 2) # gamma params 
      }
    }
  }
  
  return(hmm)
}

# 6. 4 state HMM, 2 slots
#
# state paths:
#    B -> T -> PP -> D -> B
#           -> B [optional link]
#
#  state  #pro  dist
#   PNEG     [1]   [1]
#   PG.S     [2]   [2]
#   PG.B     [3]   [2]
#   PG.E     [4]   [2]
#
# #pro: [1] background level (Poisson)
#       [2] gene start level (NegBinom | DGamma)
#       [3] body level (NegBinom | DGamma)
#       [4] gene end level (NegBinom | DGamma)
#
# dist: [1] start predictions (Geom0)
#       [2] body predictions (Geom0)
#
splithmm6.hmm <- function(scaleTP = 0.153, scalePD = 0.153, use.negbinom = FALSE) {
  
  # PNEG, PG.S, PG.B, PG.E
  n.states = 4
  PNEG = 1
  PG.S = 2
  PG.B = 3
  PG.E = 4
  
  # valid transitions
  vtbl = NULL
  tnames = c(rep("autocorr", n.states))
  # tnames = c("autocorr_covar", rep("autocorr", n.states - 1))
  
  # valid transitions
  vtbl = matrix(data=0, nrow=n.states, ncol=n.states)
  
  vtbl[PNEG, PNEG] = 1
  vtbl[PNEG, PG.S] = 2
  
  vtbl[PG.S, PG.S] = 1
  vtbl[PG.S, PG.B] = 2
  
  vtbl[PG.B, PG.B] = 1
  vtbl[PG.B, PG.E] = 2
  vtbl[PG.B, PNEG] = 3
  
  vtbl[PG.E, PG.E] = 1
  vtbl[PG.E, PNEG] = 2
  
  # if (!with.shortcut) {
  #   #
  #   # no link from T to B directly
  #   #
  #   vtbl = rbind(
  #     c(1, 2, 0, 0), # B to: B or T
  #     c(0, 1, 2, 0), # T to: T or P
  #     c(0, 0, 1, 2), # P to: P or D
  #     c(2, 0, 0, 1)) # D to: D or B
  # } else {
  #   #
  #   # no link from T to B directly
  #   #
  #   vtbl = rbind(
  #     c(1, 2, 0, 0), # B to: B or T
  #     c(3, 1, 2, 0), # T to: T or B or P
  #     c(0, 0, 1, 2), # P to: P or D
  #     c(2, 0, 0, 1)) # D to: D or B
  # }
  
  # transition groups
  # tgrps = list(2:n.states)  # share size distribution; it's a poor fit for all, but want to avoid getting a tiny 'T' state
  # tgrps = list(c(2,4)) # must exclude PG.B=3 state as it has a different number of transitions
  tgrps = NULL
  
  # emission groups
  nSlots = 2
  # egrps = new.emission.groups(n.states, nSlots)
  # egrps = add.emission.groups(egrps, states = c(2, 3, 4), slots = c(2, 2, 2)) # share GROseq over T, P and D [scaled]
  egrps = NULL
  
  # tdist = "gamma"
  tdist = "geometric"
  
  hmm = new.qhmm(list(c(1, 1), NULL), 
                 vtbl,
                 tnames,
                 list(c(tdist, tdist),      # PNEG
                      c(tdist, tdist),      # PG.S
                      c(tdist, tdist),      # PG.B
                      c(tdist, tdist)),     # PG.E
                 transition.groups = tgrps,
                 emission.groups = egrps)
  
  # set initial parameters
  set.initial.probs.qhmm(hmm, c(1, rep(0, n.states - 1))) # start with background
  
  # set transitions
  set.transition.params.qhmm(hmm, 2:n.states, 0.99)
  
  # set emissions
  # Geometric
  geom.base = 0
  # for (i in 1:N) {
  #  set.emission.params.qhmm(hmm, i, 1/3000, slot = 1)
  #  set.emission.option.qhmm(hmm, i, "base", geom.base, slot = 1)
  #}
  
  # background11
  # set.emission.params.qhmm(hmm, 1, 1/3000, slot = 1)
  # set.emission.option.qhmm(hmm, 1, "base", geom.base, slot = 1)
  
  # Gamma
  # set.emission.params.qhmm(hmm, 1, c(1.23, 5.68), slot = 1) # lambda
  # set.emission.params.qhmm(hmm, 1, c(1.23, 5.68), slot = 2) # lambda
  
  # transcribed
  # for (i in 2:n.states) {
  #  set.emission.params.qhmm(hmm, i, 1/10, slot = 1)
  #  set.emission.option.qhmm(hmm, i, "base", geom.base, slot = 1)
  # }
  
  # set.emission.params.qhmm(hmm, 2:n.states, c(0.1, 0.1), slot = 2) # gamma params
  # set.emission.option.qhmm(hmm, 3, "scale_private", scaleTP, slot = 2)
  # set.emission.option.qhmm(hmm, 4, "scale_private", scalePD, slot = 2)
  
  # set.emission.params.qhmm(hmm, 2, c(1.6, 0.3), slot = 1) # gamma params
  # set.emission.params.qhmm(hmm, 3, c(1.6, 0.3), slot = 1) # gamma params 
  # set.emission.params.qhmm(hmm, 4, c(1.6, 0.3), slot = 1) # gamma params 
  
  # set.emission.params.qhmm(hmm, 2, c(1.6, 0.3), slot = 2) # gamma params
  # set.emission.params.qhmm(hmm, 3, c(1.6, 0.3), slot = 2) # gamma params 
  # set.emission.params.qhmm(hmm, 4, c(1.6, 0.3), slot = 2) # gamma params 
  
  # set.emission.params.qhmm(hmm, 2, 0.1, slot = 2) # gamma params
  # set.emission.params.qhmm(hmm, 3, 0.1, slot = 2) # gamma params
  # set.emission.params.qhmm(hmm, 4, 0.1, slot = 2) # gamma params
  
  # Set shape and scale params for gamma distribution
  for (i in 2:n.states) {
    for (slotNum in 1:nSlots) {
      set.emission.params.qhmm(hmm, i, c(1.6, 14.7), slot = slotNum)
    }
  }
  
  # Alternatively, set emission params randomly
  # (see notes in hmm.parse.R for what diff(c(0, sort(runif(n - 1)), 1)) is doing)
  # nEmissionVals = 2  # number of possible emission values (0 and 1 in this case)
  # for (i in 1:n.states) {
  #  for (slotNum in 1:nSlots) {
  #    set.emission.params.qhmm(hmm, i, diff(c(0, sort(runif(nEmissionVals - 1)), 1)), slot = slotNum)
  #    # set.emission.option.qhmm(hmm, i, "offset", 0, slot = slotNum)
  #  }
  # }
  
  return(hmm)
}

# 7. 2 states HMM, 1 slot
#
# state paths:
#    B -> T -> B
#
#  state  #pro  dist
#   PNEG     [1]   [1]
#   PG.B     [2]   [2]
#
# #pro: [1] background level (Gamm)
#       [2] body level (Gamma)
#
# dist: [1] body predictions (Geom0)
#
hmm2states1slotContinuous <- function() {
  
  # PNEG, PG.B
  n.states = 2
  PNEG = 1
  PG.B = 2
  
  # valid transitions
  vtbl = NULL
  tnames = c(rep("autocorr", n.states))
  # tnames = c("autocorr_covar", rep("autocorr", n.states - 1))
  
  # valid transitions
  vtbl = matrix(data=0, nrow=n.states, ncol=n.states)
  
  vtbl[PNEG, PNEG] = 1
  vtbl[PNEG, PG.B] = 2
  
  vtbl[PG.B, PG.B] = 1
  vtbl[PG.B, PNEG] = 2
  
  # transition groups
  # tgrps = list(2:n.states)  # share size distribution; it's a poor fit for all, but want to avoid getting a tiny 'T' state
  # tgrps = list(c(2,4)) # must exclude PG.B=3 state as it has a different number of transitions
  tgrps = NULL
  
  # emission groups
  nSlots = 1
  # egrps = new.emission.groups(n.states, nSlots)
  # egrps = add.emission.groups(egrps, states = c(2, 3, 4), slots = c(2, 2, 2)) # share GROseq over T, P and D [scaled]
  egrps = NULL
  
  tdist = "gamma"

  
  hmm = new.qhmm(list(1, NULL), 
                 vtbl,
                 tnames,
                 list(tdist, tdist),  
                 transition.groups = tgrps,
                 emission.groups = egrps)
  
  # set initial parameters
  set.initial.probs.qhmm(hmm, c(1, 0)) # start with background
  
  # set transitions
  # set.transition.params.qhmm(hmm, 1, 0.99)
  # set.transition.params.qhmm(hmm, 2, 0.99)
  
  # Give a very high probability to self transitions
  # set.transition.params.qhmm(hmm, PNEG, c(0.99, 0.01))
  # set.transition.params.qhmm(hmm, PG.B, c(0.99, 0.01))
  
  # Randomly generating n numbers that sum to 1 
  # (where n is the number of possible state transitions)
  # https://stackoverflow.com/questions/11003967/generate-3-random-number-that-sum-to-1-in-r
  # rowSums(tm>0)[i] gives the number of transitions for state i
  # diff(c(0, sort(runif(n-1)), 1)) gives n random numbers that sum to 1
  # (this gives an error if n < 1, so makes sure eash state always has at least one transition)
  # diff(c(0, sort(runif((rowSums(tm>0)[i])-1)), 1)) gives us a random number for each state, all adding to 1
  for (i in 1:n.states) {
    set.transition.params.qhmm(hmm, i, diff(c(0, sort(runif((rowSums(vtbl>0)[i])-1)), 1)))
  }
  
  # set emissions

  # Alternatively, set emission params randomly
  # (see notes above for what diff(c(0, sort(runif(n - 1)), 1)) is doing)
  nEmissionVals = 2  # number of possible emission values (0 and 1 in this case)
  for (i in 1:n.states) {
    for (slotNum in 1:nSlots) {
      # set.emission.params.qhmm(hmm, i, c(0.1, 0.1), slot = slotNum) # dgamma params
      set.emission.params.qhmm(hmm, i, diff(c(0, sort(runif(nEmissionVals - 1)), 1)), slot = slotNum)
      set.emission.option.qhmm(hmm, i, "offset", 0, slot = slotNum)
    }
  }
  
  return(hmm)
}

# 7. 2 states HMM, 1 slot with 1 set of priors
#
# state paths:
#    B -> T -> B
#
#  state  #pro  dist
#   PNEG     [1]   [1]
#   PG.B     [2]   [2]
#
# #pro: [1] background level (Gamm)
#       [2] body level (Gamma)
#
# dist: [1] body predictions (Geom0)
#
hmm2states1slot1covarContinuous <- function() {
  
  # PNEG, PG.B
  n.states = 2
  PNEG = 1
  PG.B = 2
  
  # valid transitions
  vtbl = NULL
  # tnames = c(rep("autocorr", n.states))
  tnames = c("autocorr_covar", rep("autocorr", n.states - 1))
  # tnames = c(rep("autocorr_covar", n.states))
  
  # valid transitions
  vtbl = matrix(data=0, nrow=n.states, ncol=n.states)
  
  vtbl[PNEG, PNEG] = 1
  vtbl[PNEG, PG.B] = 2
  
  vtbl[PG.B, PG.B] = 1
  vtbl[PG.B, PNEG] = 2
  
  # transition groups
  # tgrps = list(2:n.states)  # share size distribution; it's a poor fit for all, but want to avoid getting a tiny 'T' state
  # tgrps = list(c(2,4)) # must exclude PG.B=3 state as it has a different number of transitions
  tgrps = NULL
  
  # emission groups
  nSlots = 1
  # egrps = new.emission.groups(n.states, nSlots)
  # egrps = add.emission.groups(egrps, states = c(2, 3, 4), slots = c(2, 2, 2)) # share GROseq over T, P and D [scaled]
  egrps = NULL
  
  tdist = "gamma"
  
  hmm = new.qhmm(list(1, 1), 
                 vtbl,
                 tnames,
                 list(tdist, tdist),  
                 transition.groups = tgrps,
                 emission.groups = egrps)
  
  # set initial parameters
  set.initial.probs.qhmm(hmm, c(1, 0)) # start with background
  
  # set transitions
  set.transition.params.qhmm(hmm, 2:n.states, 0.99)
  
  # set emissions
  
  # Alternatively, set emission params randomly
  # (see notes above for what diff(c(0, sort(runif(n - 1)), 1)) is doing)
  nEmissionVals = 2  # number of possible emission values (0 and 1 in this case)
  for (i in 1:n.states) {
    for (slotNum in 1:nSlots) {
      # set.emission.params.qhmm(hmm, i, c(0.1, 0.1), slot = slotNum) # dgamma params
      set.emission.params.qhmm(hmm, i, diff(c(0, sort(runif(nEmissionVals - 1)), 1)), slot = slotNum)
      set.emission.option.qhmm(hmm, i, "offset", 0, slot = slotNum)
    }
  }
  
  return(hmm)
}

# 7. 2 states HMM, 1 slot with 2 sets of priors
#
# state paths:
#    B -> T -> B
#
#  state  #pro  dist
#   PNEG     [1]   [1]
#   PG.B     [2]   [2]
#
# #pro: [1] background level (Gamm)
#       [2] body level (Gamma)
#
# dist: [1] body predictions (Geom0)
#
hmm2states1slot2covarsContinuous <- function() {
  
  # PNEG, PG.B
  n.states = 2
  PNEG = 1
  PG.B = 2
  
  # valid transitions
  vtbl = NULL
  # tnames = c(rep("autocorr", n.states))
  # tnames = c("autocorr_covar", rep("autocorr", n.states - 1))
  tnames = c(rep("autocorr_covar", n.states))
  
  # valid transitions
  vtbl = matrix(data=0, nrow=n.states, ncol=n.states)
  
  vtbl[PNEG, PNEG] = 1
  vtbl[PNEG, PG.B] = 2
  
  vtbl[PG.B, PG.B] = 1
  vtbl[PG.B, PNEG] = 2
  
  # transition groups
  # tgrps = list(2:n.states)  # share size distribution; it's a poor fit for all, but want to avoid getting a tiny 'T' state
  # tgrps = list(c(2,4)) # must exclude PG.B=3 state as it has a different number of transitions
  tgrps = NULL
  
  # emission groups
  nSlots = 1
  # egrps = new.emission.groups(n.states, nSlots)
  # egrps = add.emission.groups(egrps, states = c(2, 3, 4), slots = c(2, 2, 2)) # share GROseq over T, P and D [scaled]
  egrps = NULL
  
  tdist = "gamma"
  
  hmm = new.qhmm(list(1, c(1,1)), 
                 vtbl,
                 tnames,
                 list(tdist, tdist),  
                 transition.groups = tgrps,
                 emission.groups = egrps)
  
  # set initial parameters
  set.initial.probs.qhmm(hmm, c(1, 0)) # start with background
  
  # set transitions
  # set.transition.params.qhmm(hmm, 2:n.states, 0.99)
  
  # set emissions
  
  # Alternatively, set emission params randomly
  # (see notes above for what diff(c(0, sort(runif(n - 1)), 1)) is doing)
  nEmissionVals = 2  # number of possible emission values (0 and 1 in this case)
  for (i in 1:n.states) {
    for (slotNum in 1:nSlots) {
      # set.emission.params.qhmm(hmm, i, c(0.1, 0.1), slot = slotNum) # dgamma params
      set.emission.params.qhmm(hmm, i, diff(c(0, sort(runif(nEmissionVals - 1)), 1)), slot = slotNum)
      set.emission.option.qhmm(hmm, i, "offset", 0, slot = slotNum)
    }
  }
  
  return(hmm)
}

# 7. 2 states HMM, 2 slots, Continuous emissions distributions w/ covars
#
# state paths:
#    B -> T -> B
#
#  state  #pro  dist
#   PNEG     [1]   [1]
#   PG.B     [2]   [2]
#
# #pro: [1] background level (Gamm)
#       [2] body level (Gamma)
#
# dist: [1] body predictions (Gamma)
#
hmm2states2slots1covarContinuous <- function() {
  
  # PNEG, PG.B
  n.states = 2
  PNEG = 1
  PG.B = 2
  
  # valid transitions
  vtbl = NULL
  # tnames = c(rep("autocorr", n.states))
  tnames = c("autocorr_covar", rep("autocorr", n.states - 1))
  # tnames = c("autocorr", "autocorr_covar")
  
  # valid transitions
  vtbl = matrix(data=0, nrow=n.states, ncol=n.states)
  
  vtbl[PNEG, PNEG] = 1
  vtbl[PNEG, PG.B] = 2
  
  vtbl[PG.B, PG.B] = 1
  vtbl[PG.B, PNEG] = 2
  
  # transition groups
  # tgrps = list(2:n.states)  # share size distribution; it's a poor fit for all, but want to avoid getting a tiny 'T' state
  # tgrps = list(c(2,4)) # must exclude PG.B=3 state as it has a different number of transitions
  tgrps = NULL
  
  # emission groups
  nSlots = 2
  # egrps = new.emission.groups(n.states, nSlots)
  # egrps = add.emission.groups(egrps, states = c(2, 3, 4), slots = c(2, 2, 2)) # share GROseq over T, P and D [scaled]
  egrps = NULL
  
  tdist = "gamma"

  hmm = new.qhmm(list(c(1, 1), 1), # covar: TSS signal
                 vtbl,
                 tnames,
                 list(c("geometric", tdist),  # PNEG
                      c("geometric", tdist)),     # PG.B
                 transition.groups = tgrps,
                 emission.groups = egrps)
  
  # set initial parameters
  set.initial.probs.qhmm(hmm, c(1, 0)) # start with background
  
  # set transitions
  set.transition.params.qhmm(hmm, 2:n.states, 0.99)
  # set.transition.params.qhmm(hmm, 1, 0.99)

  # set emissions
  geom.base = 0
  # background
  # set.emission.params.qhmm(hmm, 1, 1/3000, slot = 1)
  # set.emission.option.qhmm(hmm, 1, "base", geom.base, slot = 1)

  # transcribed
  for (i in 1:n.states) {
    set.emission.params.qhmm(hmm, i, 1/10, slot = 1)
    set.emission.option.qhmm(hmm, i, "base", geom.base, slot = 1)
  }
  
  
  # Alternatively, set emission params randomly
  # (see notes above for what diff(c(0, sort(runif(n - 1)), 1)) is doing)
  nEmissionVals = 2  # number of possible emission values (0 and 1 in this case)
  for (i in 1:n.states) {
    for (slotNum in 2:nSlots) {
      # set.emission.params.qhmm(hmm, i, c(0.1, 0.1), slot = slotNum) # dgamma params
      set.emission.params.qhmm(hmm, i, diff(c(0, sort(runif(nEmissionVals - 1)), 1)), slot = slotNum)
      set.emission.option.qhmm(hmm, i, "offset", 0, slot = slotNum)
    }
  }
  
  return(hmm)
}

# 7. 2 states HMM, 3 slots, Continuous emissions distributions w/ covars
#
# state paths:
#    B -> T -> B
#
#  state  #pro  dist
#   PNEG     [1]   [1]
#   PG.B     [2]   [2]
#
# #pro: [1] background level (Gamm)
#       [2] body level (Gamma)
#
# dist: [1] body predictions (Gamma)
#
hmm2states3slotsContinuousWCovars <- function() {
  
  # PNEG, PG.B
  n.states = 2
  PNEG = 1
  PG.B = 2
  
  # valid transitions
  vtbl = NULL
  # tnames = c(rep("autocorr", n.states))
  tnames = c("autocorr_covar", rep("autocorr", n.states - 1))
  # tnames = c("autocorr", "autocorr_covar")
  
  # valid transitions
  vtbl = matrix(data=0, nrow=n.states, ncol=n.states)
  
  vtbl[PNEG, PNEG] = 1
  vtbl[PNEG, PG.B] = 2
  
  vtbl[PG.B, PG.B] = 1
  vtbl[PG.B, PNEG] = 2
  
  # transition groups
  # tgrps = list(2:n.states)  # share size distribution; it's a poor fit for all, but want to avoid getting a tiny 'T' state
  # tgrps = list(c(2,4)) # must exclude PG.B=3 state as it has a different number of transitions
  tgrps = NULL
  
  # emission groups
  nSlots = 3
  # egrps = new.emission.groups(n.states, nSlots)
  # egrps = add.emission.groups(egrps, states = c(2, 3, 4), slots = c(2, 2, 2)) # share GROseq over T, P and D [scaled]
  egrps = NULL
  
  tdist = "gamma"

  hmm <- new.qhmm(list(rep(1, nSlots), 1),
                  vtbl,
                  tnames,
                  rep(list(rep(tdist, nSlots)), n.states),
                  transition.groups = tgrps,
                  emission.groups = egrps)
  
  # set initial parameters
  set.initial.probs.qhmm(hmm, c(1, 0)) # start with background
  
  # set transitions
  set.transition.params.qhmm(hmm, 2:n.states, 0.99)

  # set emissions
  
  # Alternatively, set emission params randomly
  # (see notes above for what diff(c(0, sort(runif(n - 1)), 1)) is doing)
  nEmissionVals = 2  # number of possible emission values (0 and 1 in this case)
  for (i in 1:n.states) {
    for (slotNum in 1:nSlots) {
      # set.emission.params.qhmm(hmm, i, c(0.1, 0.1), slot = slotNum) # dgamma params
      set.emission.params.qhmm(hmm, i, diff(c(0, sort(runif(nEmissionVals - 1)), 1)), slot = slotNum)
      set.emission.option.qhmm(hmm, i, "offset", 0, slot = slotNum)
    }
  }
  
  return(hmm)
}

# 7. 2 states HMM, 4 slots, Continuous emissions distributions w/ covars
#
# state paths:
#    B -> T -> B
#
#  state  #pro  dist
#   PNEG     [1]   [1]
#   PG.B     [2]   [2]
#
# #pro: [1] background level (Gamm)
#       [2] body level (Gamma)
#
# dist: [1] body predictions (Gamma)
#
hmm2states4slotsContinuousWCovars <- function() {
  
  # PNEG, PG.B
  n.states = 2
  PNEG = 1
  PG.B = 2
  
  # valid transitions
  vtbl = NULL
  # tnames = c(rep("autocorr", n.states))
  tnames = c("autocorr_covar", rep("autocorr", n.states - 1))
  
  # valid transitions
  vtbl = matrix(data=0, nrow=n.states, ncol=n.states)
  
  vtbl[PNEG, PNEG] = 1
  vtbl[PNEG, PG.B] = 2
  
  vtbl[PG.B, PG.B] = 1
  vtbl[PG.B, PNEG] = 2
  
  # transition groups
  # tgrps = list(2:n.states)  # share size distribution; it's a poor fit for all, but want to avoid getting a tiny 'T' state
  # tgrps = list(c(2,4)) # must exclude PG.B=3 state as it has a different number of transitions
  tgrps = NULL
  
  # emission groups
  nSlots = 4
  # egrps = new.emission.groups(n.states, nSlots)
  # egrps = add.emission.groups(egrps, states = c(2, 3, 4), slots = c(2, 2, 2)) # share GROseq over T, P and D [scaled]
  egrps = NULL
  
  tdist = "gamma"

  hmm <- new.qhmm(list(rep(1, nSlots), 1),
                  vtbl,
                  tnames,
                  rep(list(rep(tdist, nSlots)), n.states),
                  transition.groups = tgrps,
                  emission.groups = egrps)
  
  # set initial parameters
  set.initial.probs.qhmm(hmm, c(1, 0)) # start with background
  
  # set transitions
  set.transition.params.qhmm(hmm, 2:n.states, 0.99)

  # set emissions

  # Alternatively, set emission params randomly
  # (see notes above for what diff(c(0, sort(runif(n - 1)), 1)) is doing)
  nEmissionVals = 2  # number of possible emission values (0 and 1 in this case)
  for (i in 1:n.states) {
    for (slotNum in 1:nSlots) {
      # set.emission.params.qhmm(hmm, i, c(0.1, 0.1), slot = slotNum) # dgamma params
      set.emission.params.qhmm(hmm, i, diff(c(0, sort(runif(nEmissionVals - 1)), 1)), slot = slotNum)
      set.emission.option.qhmm(hmm, i, "offset", 0, slot = slotNum)
    }
  }
  
  return(hmm)
}

# 8. 2 states HMM, 3 slots
#
# state paths:
#    B -> T -> B
#
#  state  #pro  dist
#   PNEG     [1]   [1]
#   PG.B     [2]   [2]
#
# #pro: [1] background level (gamma)
#       [2] gene start level (gamma)
#       [3] body level (gamma)
#       [4] gene end level (gamma)
#
splithmm8.hmm <- function() {
  
  # PNEG, PG.B
  n.states = 2
  PNEG = 1
  PG.B = 2
  
  # valid transitions
  vtbl = NULL
  tnames = c(rep("autocorr", n.states))
  # tnames = c("autocorr_covar", rep("autocorr", n.states - 1))
  
  # valid transitions
  vtbl = matrix(data=0, nrow=n.states, ncol=n.states)
  
  vtbl[PNEG, PNEG] = 1
  vtbl[PNEG, PG.B] = 2
  
  vtbl[PG.B, PG.B] = 1
  vtbl[PG.B, PNEG] = 2
  
  # transition groups
  # tgrps = list(2:n.states)  # share size distribution; it's a poor fit for all, but want to avoid getting a tiny 'T' state
  # tgrps = list(c(2,4)) # must exclude PG.B=3 state as it has a different number of transitions
  tgrps = NULL
  
  # emission groups
  nSlots = 3
  # egrps = new.emission.groups(n.states, nSlots)
  # egrps = add.emission.groups(egrps, states = c(2, 3, 4), slots = c(2, 2, 2)) # share GROseq over T, P and D [scaled]
  egrps = NULL
  
  backgrounddist = "normal"
  tdist = "gamma"
  
  hmm <- new.qhmm(list(rep(1, nSlots), NULL),
                  vtbl,
                  tnames,
                  rep(list(rep(tdist, nSlots)), n.states),
                  transition.groups = tgrps,
                  emission.groups = egrps)
  
  # set initial parameters
  set.initial.probs.qhmm(hmm, c(1, rep(0, n.states - 1))) # start with background
  
  # set transitions
  set.transition.params.qhmm(hmm, 1:n.states, 0.99)
  
  # set emissions
  
  # Geometric
  geom.base = 0
  # for (i in 1:n.states) {
  #  set.emission.params.qhmm(hmm, i, 1/3000, slot = 1)
  #  set.emission.option.qhmm(hmm, i, "base", geom.base, slot = 1)
  # }
  
  # background11
  # set.emission.params.qhmm(hmm, 1, 1/3000, slot = 1)
  # set.emission.option.qhmm(hmm, 1, "base", geom.base, slot = 1)
  
  # transcribed
  # for (i in 2:n.states) {
  #  set.emission.params.qhmm(hmm, i, 1/10, slot = 1)
  #  set.emission.option.qhmm(hmm, i, "base", geom.base, slot = 1)
  # }
  
  # Set shape and scale params for gamma distribution
  set.emission.params.qhmm(hmm, 1, c(0.581222, 0.023408), slot = 1) # gamma
  set.emission.params.qhmm(hmm, 2, c(0.280569, 0.042092), slot = 1) # gamma
  
  set.emission.params.qhmm(hmm, 1, c(1.685649, 0.049756), slot = 2) # gamma
  set.emission.params.qhmm(hmm, 2, c(0.862261, 0.367080), slot = 2) # gamma
  
  set.emission.params.qhmm(hmm, 1, c(0.986111, 0.024571), slot = 3) # gamma
  set.emission.params.qhmm(hmm, 2, c(0.563369, 0.064000), slot = 3) # gamma
  
  # Normal
  # set.emission.params.qhmm(hmm, 1, c(0.01360517, 0.001050558), slot = 1)
  # set.emission.params.qhmm(hmm, 1, c(0.08387034, 0.006401999), slot = 2)
  # set.emission.params.qhmm(hmm, 1, c(0.02422957, 0.001753884), slot = 3)
  
  return(hmm)
}

# 4 states HMM, 1 slot, 1 covar
#
# state paths:
#    PNEG -> S -> B -> E -> PNEG
#
#  state  #pro  dist
#   PNEG     [1]   [1]
#   PG.s     [2]   [2]
#   PG.B     [3]   [2]
#   PG.E     [4]   [2]
#
# #pro: [1] background level (Gamma)
#
hmm4states1slot1covarContinuous <- function() {
  
  # PNEG, PG.S, PG.B, PG.E
  n.states = 4
  PNEG = 1
  PG.S = 2
  PG.B = 3
  PG.E = 4
  
  # valid transitions
  nCovars = 1
  # tnames = c(rep("autocorr", n.states))
  tnames = c("autocorr_covar", rep("autocorr", n.states - 1))
  
  # valid transitions
  vtbl = NULL
  vtbl = matrix(data=0, nrow=n.states, ncol=n.states)
  
  vtbl[PNEG, PNEG] = 1
  vtbl[PNEG, PG.S] = 2
  
  vtbl[PG.S, PG.S] = 1
  vtbl[PG.S, PG.B] = 2
  
  vtbl[PG.B, PG.B] = 1
  vtbl[PG.B, PG.E] = 2
  
  vtbl[PG.E, PG.E] = 1
  vtbl[PG.E, PNEG] = 2
  
  # transition groups
  # tgrps = list(2:n.states)  # share size distribution; it's a poor fit for all, but want to avoid getting a tiny 'T' state
  # tgrps = list(c(2,4)) # must exclude PG.B=3 state as it has a different number of transitions
  tgrps = NULL
  
  # emission groups
  nSlots = 1
  # egrps = new.emission.groups(n.states, nSlots)
  # egrps = add.emission.groups(egrps, states = c(2, 3, 4), slots = c(2, 2, 2)) # share GROseq over T, P and D [scaled]
  egrps = NULL
  
  # backgrounddist = "normal"
  tdist = "gamma"
  
  hmm <- new.qhmm(list(rep(1, nSlots), rep(1, nCovars)),
                  vtbl,
                  tnames,
                  rep(list(rep(tdist, nSlots)), n.states),
                  transition.groups = tgrps,
                  emission.groups = egrps)
  
  # set initial parameters
  set.initial.probs.qhmm(hmm, c(1, rep(0, n.states - 1))) # start with background
  
  # set transitions
  set.transition.params.qhmm(hmm, 2:n.states, 0.99)
  
  # Give a very high probability to self transitions
  # set.transition.params.qhmm(hmm, PNEG, c(0.99, 0.01))
  # set.transition.params.qhmm(hmm, PG.S, c(0.99, 0.01))
  # set.transition.params.qhmm(hmm, PG.B, c(0.99, 0.01))
  # set.transition.params.qhmm(hmm, PG.E, c(0.99, 0.01))
  
  # set emissions
  
  # Alternatively, set emission params randomly
  # (see notes above for what diff(c(0, sort(runif(n - 1)), 1)) is doing)
  nEmissionVals = 2  # number of possible emission values (0 and 1 in this case)
  for (i in 1:n.states) {
    for (slotNum in 1:nSlots) {
      # set.emission.params.qhmm(hmm, i, c(0.1, 0.1), slot = slotNum) # dgamma params
      set.emission.params.qhmm(hmm, i, diff(c(0, sort(runif(nEmissionVals - 1)), 1)), slot = slotNum)
      set.emission.option.qhmm(hmm, i, "offset", 0, slot = slotNum)
    }
  }
  
  return(hmm)
}

# 9. 4 states HMM, 3 slots
#
# state paths:
#    B -> T -> B
#
#  state  #pro  dist
#   PNEG     [1]   [1]
#   PG.s     [2]   [2]
#   PG.B     [3]   [2]
#   PG.E     [4]   [2]
#
# #pro: [1] background level (Poisson)
#       [2] gene start level (NegBinom | DGamma)
#       [3] body level (NegBinom | DGamma)
#       [4] gene end level (NegBinom | DGamma)
#
hmm4states3slotsContinuous <- function() {
  
  # PNEG, PG.S, PG.B, PG.E
  n.states = 4
  PNEG = 1
  PG.S = 2
  PG.B = 3
  PG.E = 4
  
  # valid transitions
  vtbl = NULL
  tnames = c(rep("autocorr", n.states))
  # tnames = c("autocorr_covar", rep("autocorr", n.states - 1))
  
  # valid transitions
  vtbl = matrix(data=0, nrow=n.states, ncol=n.states)
  
  vtbl[PNEG, PNEG] = 1
  vtbl[PNEG, PG.S] = 2
  
  vtbl[PG.S, PG.S] = 1
  vtbl[PG.S, PG.B] = 2
  
  vtbl[PG.B, PG.B] = 1
  vtbl[PG.B, PG.E] = 2
  
  vtbl[PG.E, PG.E] = 1
  vtbl[PG.E, PNEG] = 2
  
  # transition groups
  # tgrps = list(2:n.states)  # share size distribution; it's a poor fit for all, but want to avoid getting a tiny 'T' state
  # tgrps = list(c(2,4)) # must exclude PG.B=3 state as it has a different number of transitions
  tgrps = NULL
  
  # emission groups
  nSlots = 3
  # egrps = new.emission.groups(n.states, nSlots)
  # egrps = add.emission.groups(egrps, states = c(2, 3, 4), slots = c(2, 2, 2)) # share GROseq over T, P and D [scaled]
  egrps = NULL
  
  # backgrounddist = "normal"
  tdist = "gamma"
  
  hmm <- new.qhmm(list(rep(1, nSlots), NULL),
                  vtbl,
                  tnames,
                  rep(list(rep(tdist, nSlots)), n.states),
                  transition.groups = tgrps,
                  emission.groups = egrps)
  
  # set initial parameters
  set.initial.probs.qhmm(hmm, c(1, rep(0, n.states - 1))) # start with background
  
  # set transitions
  # set.transition.params.qhmm(hmm, 1:n.states, 0.99)
  
  # Give a very high probability to self transitions
  set.transition.params.qhmm(hmm, PNEG, c(0.99, 0.01))
  set.transition.params.qhmm(hmm, PG.S, c(0.99, 0.01))
  set.transition.params.qhmm(hmm, PG.B, c(0.99, 0.01))
  set.transition.params.qhmm(hmm, PG.E, c(0.99, 0.01))
  
  # Randomly generating n numbers that sum to 1 
  # (where n is the number of possible state transitions)
  # https://stackoverflow.com/questions/11003967/generate-3-random-number-that-sum-to-1-in-r
  # rowSums(tm>0)[i] gives the number of transitions for state i
  # diff(c(0, sort(runif(n-1)), 1)) gives n random numbers that sum to 1
  # (this gives an error if n < 1, so makes sure eash state always has at least one transition)
  # diff(c(0, sort(runif((rowSums(tm>0)[i])-1)), 1)) gives us a random number for each state, all adding to 1
  # for (i in 1:n.states) {
  #   set.transition.params.qhmm(hmm, i, diff(c(0, sort(runif((rowSums(vtbl>0)[i])-1)), 1)))
  # }
  
  # set emissions
  
  # Geometric
  geom.base = 0
  # for (i in 1:n.states) {
  #  set.emission.params.qhmm(hmm, i, 1/3000, slot = 1)
  #  set.emission.option.qhmm(hmm, i, "base", geom.base, slot = 1)
  # }
  
  # background11
  # set.emission.params.qhmm(hmm, 1, 1/3000, slot = 1)
  # set.emission.option.qhmm(hmm, 1, "base", geom.base, slot = 1)
  
  # transcribed
  # for (i in 2:n.states) {
  #  set.emission.params.qhmm(hmm, i, 1/10, slot = 1)
  #  set.emission.option.qhmm(hmm, i, "base", geom.base, slot = 1)
  # }
  
  # Set shape and scale params for gamma distribution
  set.emission.params.qhmm(hmm, 1, c(0.581222, 0.023408), slot = 1) # gamma
  set.emission.params.qhmm(hmm, 2, c(0.280569, 0.042092), slot = 1) # gamma
  set.emission.params.qhmm(hmm, 3, c(0.862261, 0.367080), slot = 1) # gamma
  set.emission.params.qhmm(hmm, 4, c(0.563369, 0.064000), slot = 1) # gamma
  
  set.emission.params.qhmm(hmm, 1, c(1.685649, 0.049756), slot = 2) # gamma
  set.emission.params.qhmm(hmm, 2, c(0.280569, 0.042092), slot = 2) # gamma
  set.emission.params.qhmm(hmm, 3, c(0.862261, 0.367080), slot = 2) # gamma
  set.emission.params.qhmm(hmm, 4, c(0.563369, 0.064000), slot = 2) # gamma
  
  set.emission.params.qhmm(hmm, 1, c(0.986111, 0.024571), slot = 3) # gamma
  set.emission.params.qhmm(hmm, 2, c(0.280569, 0.042092), slot = 3) # gamma
  set.emission.params.qhmm(hmm, 3, c(0.862261, 0.367080), slot = 3) # gamma
  set.emission.params.qhmm(hmm, 4, c(0.563369, 0.064000), slot = 3) # gamma
  
  # Normal
  # set.emission.params.qhmm(hmm, 1, c(0.01360517, 0.001050558), slot = 1)
  # set.emission.params.qhmm(hmm, 1, c(0.08387034, 0.006401999), slot = 2)
  # set.emission.params.qhmm(hmm, 1, c(0.02422957, 0.001753884), slot = 3)
  
  # Alternatively, set emission params randomly
  # (see notes above for what diff(c(0, sort(runif(n - 1)), 1)) is doing)
  nEmissionVals = 2  # number of possible emission values (0 and 1 in this case)
  # for (i in 1:n.states) {
  #  for (slotNum in 1:nSlots) {
  #    # set.emission.params.qhmm(hmm, i, c(0.1, 0.1), slot = slotNum) # dgamma params
  #    set.emission.params.qhmm(hmm, i, diff(c(0, sort(runif(nEmissionVals - 1)), 1)), slot = slotNum)
  #    set.emission.option.qhmm(hmm, i, "offset", 0, slot = slotNum)
  #  }
  # }
  
  return(hmm)
}

# 4 states HMM, 3 slots, 1 covar
#
# state paths:
#    PNEG -> S -> B -> E -> PNEG
#
#  state  #pro  dist
#   PNEG     [1]   [1]
#   PG.s     [2]   [2]
#   PG.B     [3]   [2]
#   PG.E     [4]   [2]
#
# #pro: [1] background level (Gamma)
#       [2] gene start level (Gamma)
#       [3] body level (Gamma)
#       [4] gene end level (Gamma)
#
hmm4states3slots1covarContinuous <- function() {
  
  # PNEG, PG.S, PG.B, PG.E
  n.states = 4
  PNEG = 1
  PG.S = 2
  PG.B = 3
  PG.E = 4
  
  # valid transitions
  nCovars = 1
  # tnames = c(rep("autocorr", n.states))
  tnames = c("autocorr_covar", rep("autocorr", n.states - 1))
  
  # valid transitions
  vtbl = NULL
  vtbl = matrix(data=0, nrow=n.states, ncol=n.states)
  
  vtbl[PNEG, PNEG] = 1
  vtbl[PNEG, PG.S] = 2
  
  vtbl[PG.S, PG.S] = 1
  vtbl[PG.S, PG.B] = 2
  
  vtbl[PG.B, PG.B] = 1
  vtbl[PG.B, PG.E] = 2
  
  vtbl[PG.E, PG.E] = 1
  vtbl[PG.E, PNEG] = 2
  
  # transition groups
  # tgrps = list(2:n.states)  # share size distribution; it's a poor fit for all, but want to avoid getting a tiny 'T' state
  # tgrps = list(c(2,4)) # must exclude PG.B=3 state as it has a different number of transitions
  tgrps = NULL
  
  # emission groups
  nSlots = 3
  # egrps = new.emission.groups(n.states, nSlots)
  # egrps = add.emission.groups(egrps, states = c(2, 3, 4), slots = c(2, 2, 2)) # share GROseq over T, P and D [scaled]
  egrps = NULL
  
  # backgrounddist = "normal"
  tdist = "gamma"
  
  hmm <- new.qhmm(list(rep(1, nSlots), rep(1, nCovars)),
                  vtbl,
                  tnames,
                  rep(list(rep(tdist, nSlots)), n.states),
                  transition.groups = tgrps,
                  emission.groups = egrps)
  
  # set initial parameters
  set.initial.probs.qhmm(hmm, c(1, rep(0, n.states - 1))) # start with background
  
  # set transitions
  set.transition.params.qhmm(hmm, 2:n.states, 0.99)
  
  # Give a very high probability to self transitions
  # set.transition.params.qhmm(hmm, PNEG, c(0.99, 0.01))
  # set.transition.params.qhmm(hmm, PG.S, c(0.99, 0.01))
  # set.transition.params.qhmm(hmm, PG.B, c(0.99, 0.01))
  # set.transition.params.qhmm(hmm, PG.E, c(0.99, 0.01))
  
  # set emissions
  
  # Alternatively, set emission params randomly
  # (see notes above for what diff(c(0, sort(runif(n - 1)), 1)) is doing)
  nEmissionVals = 2  # number of possible emission values (0 and 1 in this case)
  for (i in 1:n.states) {
    for (slotNum in 1:nSlots) {
      # set.emission.params.qhmm(hmm, i, c(0.1, 0.1), slot = slotNum) # dgamma params
      set.emission.params.qhmm(hmm, i, diff(c(0, sort(runif(nEmissionVals - 1)), 1)), slot = slotNum)
      set.emission.option.qhmm(hmm, i, "offset", 0, slot = slotNum)
    }
  }
  
  return(hmm)
}

# 10. 6 states HMM, 5 slots
#
# state paths:
#    B -> T -> B
#
#  state  #pro  dist
#   PNEG     [1]   [1]
#   TSS      [2]   [2]
#   PG.S     [2]   [2]
#   PG.B     [3]   [2]
#   PG.E     [4]   [2]
#   PAG      [4]   [2]
#
# #pro: [1] background level (Normal)
#       [2] TSS (Gamma)
#       [3] gene start level (Gamma)
#       [4] body level (Gamma)
#       [5] gene end level (Gamma)
#       [6] after gene (Gamma)
#
splithmm10.hmm <- function() {
  
  # PNEG, PG.S, PG.B, PG.E
  n.states = 6
  PNEG = 1
  TSS = 2
  PG.S = 3
  PG.B = 4
  PG.E = 5
  PAG = 6
  
  # valid transitions
  vtbl = NULL
  tnames = c(rep("autocorr", n.states))
  # tnames = c("autocorr_covar", rep("autocorr", n.states - 1))
  
  # valid transitions
  vtbl = matrix(data=0, nrow=n.states, ncol=n.states)
  
  vtbl[PNEG, PNEG] = 1
  vtbl[PNEG, TSS] = 2
  vtbl[PNEG, PG.S] = 3
  
  vtbl[TSS, TSS] = 1
  vtbl[TSS, PNEG] = 2
  vtbl[TSS, PG.S] = 3
  
  vtbl[PG.S, PG.S] = 1
  vtbl[PG.S, PG.B] = 2
  
  vtbl[PG.B, PG.B] = 1
  vtbl[PG.B, PG.E] = 2
  vtbl[PG.B, PAG] = 3
  vtbl[PG.B, PNEG] = 4
  
  vtbl[PG.E, PG.E] = 1
  vtbl[PG.E, PAG] = 2
  vtbl[PG.E, PNEG] = 3
  
  vtbl[PAG, PAG] = 1
  vtbl[PAG, PNEG] = 2
  
  # transition groups
  # tgrps = list(2:n.states)  # share size distribution; it's a poor fit for all, but want to avoid getting a tiny 'T' state
  # tgrps = list(c(2,4)) # must exclude PG.B=3 state as it has a different number of transitions
  tgrps = NULL
  
  # emission groups
  nSlots = 5
  # egrps = new.emission.groups(n.states, nSlots)
  # egrps = add.emission.groups(egrps, states = c(2, 3, 4), slots = c(2, 2, 2)) # share GROseq over T, P and D [scaled]
  egrps = NULL
  
  backgrounddist = "normal"
  tdist = "gamma"
  
  hmm <- new.qhmm(list(rep(1, nSlots), NULL),
                  vtbl,
                  tnames,
                  list(c(tdist, tdist, tdist, tdist, tdist),      # PNEG
                       c(tdist, tdist, tdist, tdist, tdist),      # TSS
                       c(tdist, tdist, tdist, tdist, tdist),      # PG.S
                       c(backgrounddist, tdist, tdist, tdist, tdist),      # PG.B
                       c(tdist, tdist, tdist, tdist, tdist),      # PG.E
                       c(tdist, tdist, tdist, tdist, tdist)),     # PAG
                  transition.groups = tgrps,
                  emission.groups = egrps)
  
  # set initial parameters
  set.initial.probs.qhmm(hmm, c(1, rep(0, n.states - 1))) # start with background
  
  # set transitions
  # set.transition.params.qhmm(hmm, 1:n.states, 0.99)
  
  # Give a very high probability to self transitions
  set.transition.params.qhmm(hmm, PNEG, c(0.98, 0.01, 0.01))
  set.transition.params.qhmm(hmm, TSS, c(0.98, 0.01, 0.01))
  set.transition.params.qhmm(hmm, PG.S, c(0.99, 0.01))
  set.transition.params.qhmm(hmm, PG.B, c(0.97, 0.01, 0.01, 0.01))
  set.transition.params.qhmm(hmm, PG.E, c(0.98, 0.01, 0.01))
  set.transition.params.qhmm(hmm, PAG, c(0.99, 0.01))
  
  # set emissions
  
  # Geometric
  geom.base = 0
  # for (i in 1:n.states) {
  #  set.emission.params.qhmm(hmm, i, 1/3000, slot = 1)
  #  set.emission.option.qhmm(hmm, i, "base", geom.base, slot = 1)
  # }
  
  # background11
  # set.emission.params.qhmm(hmm, 1, 1/3000, slot = 1)
  # set.emission.option.qhmm(hmm, 1, "base", geom.base, slot = 1)
  
  # transcribed
  # for (i in 2:n.states) {
  #  set.emission.params.qhmm(hmm, i, 1/10, slot = 1)
  #  set.emission.option.qhmm(hmm, i, "base", geom.base, slot = 1)
  # }
  
  # Set shape and scale params for gamma distribution
  set.emission.params.qhmm(hmm, 1, c(0.581222, 0.023408), slot = 1) # gamma
  set.emission.params.qhmm(hmm, 2, c(0.280569, 0.042092), slot = 1) # gamma
  set.emission.params.qhmm(hmm, 3, c(0.280569, 0.042092), slot = 1) # gamma
  set.emission.params.qhmm(hmm, 4, c(0.862261, 0.367080), slot = 1) # gamma
  set.emission.params.qhmm(hmm, 5, c(0.563369, 0.064000), slot = 1) # gamma
  set.emission.params.qhmm(hmm, 6, c(0.563369, 0.064000), slot = 1) # gamma
  
  set.emission.params.qhmm(hmm, 1, c(0.581222, 0.023408), slot = 2) # gamma
  set.emission.params.qhmm(hmm, 2, c(0.280569, 0.042092), slot = 2) # gamma
  set.emission.params.qhmm(hmm, 3, c(0.280569, 0.042092), slot = 2) # gamma
  set.emission.params.qhmm(hmm, 4, c(0.862261, 0.367080), slot = 2) # gamma
  set.emission.params.qhmm(hmm, 5, c(0.563369, 0.064000), slot = 2) # gamma
  set.emission.params.qhmm(hmm, 6, c(0.563369, 0.064000), slot = 2) # gamma
  
  # set.emission.params.qhmm(hmm, 1, c(1.685649, 0.049756), slot = 3) # gamma
  set.emission.params.qhmm(hmm, 2, c(0.280569, 0.042092), slot = 3) # gamma
  set.emission.params.qhmm(hmm, 3, c(0.280569, 0.042092), slot = 3) # gamma
  set.emission.params.qhmm(hmm, 4, c(0.862261, 0.367080), slot = 3) # gamma
  set.emission.params.qhmm(hmm, 5, c(0.563369, 0.064000), slot = 3) # gamma
  set.emission.params.qhmm(hmm, 6, c(0.563369, 0.064000), slot = 3) # gamma
  
  set.emission.params.qhmm(hmm, 1, c(0.986111, 0.024571), slot = 4) # gamma
  set.emission.params.qhmm(hmm, 2, c(0.280569, 0.042092), slot = 4) # gamma
  set.emission.params.qhmm(hmm, 3, c(0.280569, 0.042092), slot = 4) # gamma
  set.emission.params.qhmm(hmm, 4, c(0.862261, 0.367080), slot = 4) # gamma
  set.emission.params.qhmm(hmm, 5, c(0.563369, 0.064000), slot = 4) # gamma
  set.emission.params.qhmm(hmm, 6, c(0.563369, 0.064000), slot = 4) # gamma
  
  set.emission.params.qhmm(hmm, 1, c(0.986111, 0.024571), slot = 5) # gamma
  set.emission.params.qhmm(hmm, 2, c(0.280569, 0.042092), slot = 5) # gamma
  set.emission.params.qhmm(hmm, 3, c(0.280569, 0.042092), slot = 5) # gamma
  set.emission.params.qhmm(hmm, 4, c(0.862261, 0.367080), slot = 5) # gamma
  set.emission.params.qhmm(hmm, 5, c(0.563369, 0.064000), slot = 5) # gamma
  set.emission.params.qhmm(hmm, 6, c(0.563369, 0.064000), slot = 5) # gamma
  
  # Normal
  # set.emission.params.qhmm(hmm, 1, c(0.01360517, 0.001050558), slot = 1)
  # set.emission.params.qhmm(hmm, 1, c(0.08387034, 0.006401999), slot = 2)
  set.emission.params.qhmm(hmm, 1, c(0.02422957, 0.001753884), slot = 3)
  
  return(hmm)
}

# 11. 3 state HMM, 2 slots
#
# state paths:
#    B -> T -> PP -> D -> B
#           -> B [optional link]
#
#  state  #pro  dist
#   PNEG     [1]   [1]
#   PG.S     [2]   [2]
#   PG.B     [3]   [2]
#
# #pro: [1] background level (gamma)
#       [2] gene start level (gamma)
#       [3] body level (gamma)
#
# dist: [1] start predictions (Geom0)
#       [2] body predictions (Geom0)
#
hmm3states2slotsContinuous <- function() {
  
  # PNEG, PG.S, PG.B
  n.states = 3
  PNEG = 1
  PG.S = 2
  PG.B = 3
  
  # valid transitions
  vtbl = NULL
  tnames = c(rep("autocorr", n.states))
  # tnames = c("autocorr_covar", rep("autocorr", n.states - 1))
  
  # valid transitions
  vtbl = matrix(data=0, nrow=n.states, ncol=n.states)
  
  vtbl[PNEG, PNEG] = 1
  vtbl[PNEG, PG.S] = 2
  
  vtbl[PG.S, PG.S] = 1
  vtbl[PG.S, PG.B] = 2
  
  vtbl[PG.B, PG.B] = 1
  vtbl[PG.B, PNEG] = 2
  
  # transition groups
  # tgrps = list(2:n.states)  # share size distribution; it's a poor fit for all, but want to avoid getting a tiny 'T' state
  # tgrps = list(c(2,4)) # must exclude PG.B=3 state as it has a different number of transitions
  tgrps = NULL
  
  # emission groups
  nSlots = 2
  # egrps = new.emission.groups(n.states, nSlots)
  # egrps = add.emission.groups(egrps, states = c(2, 3, 4), slots = c(2, 2, 2)) # share GROseq over T, P and D [scaled]
  egrps = NULL
  
  # backgrounddist = "normal"
  tdist = "gamma"
  
  hmm = new.qhmm(list(c(1, 1), NULL), 
                 vtbl,
                 tnames,
                 list(c(tdist, tdist),      # PNEG
                      c(tdist, tdist),      # PG.S
                      c(tdist, tdist)),     # PG.B
                 transition.groups = tgrps,
                 emission.groups = egrps)
  
  # set initial parameters
  set.initial.probs.qhmm(hmm, c(1, rep(0, n.states - 1))) # start with background
  
  # set transitions
  # set.transition.params.qhmm(hmm, 1:n.states, 0.99)
  
  # Give a very high probability to self transitions
  # set.transition.params.qhmm(hmm, PNEG, c(0.99, 0.01))
  # set.transition.params.qhmm(hmm, PG.S, c(0.99, 0.01))
  # set.transition.params.qhmm(hmm, PG.B, c(0.99, 0.01))
  
  # Randomly generating n numbers that sum to 1 
  # (where n is the number of possible state transitions)
  # https://stackoverflow.com/questions/11003967/generate-3-random-number-that-sum-to-1-in-r
  # rowSums(tm>0)[i] gives the number of transitions for state i
  # diff(c(0, sort(runif(n-1)), 1)) gives n random numbers that sum to 1
  # (this gives an error if n < 1, so makes sure eash state always has at least one transition)
  # diff(c(0, sort(runif((rowSums(tm>0)[i])-1)), 1)) gives us a random number for each state, all adding to 1
  for (i in 1:n.states) {
    set.transition.params.qhmm(hmm, i, diff(c(0, sort(runif((rowSums(vtbl>0)[i])-1)), 1)))
  }
  
  # Set emissions
  # Set shape and scale params for gamma distribution
  # set.emission.params.qhmm(hmm, 1, c(0.38, 0.03), slot = 1) # gamma
  # set.emission.params.qhmm(hmm, 2, c(0.25, 1.08), slot = 1) # gamma
  # set.emission.params.qhmm(hmm, 3, c(0.84, 0.38), slot = 1) # gamma
  
  # set.emission.params.qhmm(hmm, 1, c(1.71, 0.05), slot = 2) # gamma
  # set.emission.params.qhmm(hmm, 2, c(0.25, 1.08), slot = 2) # gamma
  # set.emission.params.qhmm(hmm, 3, c(0.84, 0.38), slot = 2) # gamma
  
  # Normal
  # set.emission.params.qhmm(hmm, 1, c(0.01360517, 0.001050558), slot = 1)
  # set.emission.params.qhmm(hmm, 1, c(0.0869, 0.0064), slot = 2)
  
  # Alternatively, set emission params randomly
  # (see notes above for what diff(c(0, sort(runif(n - 1)), 1)) is doing)
  nEmissionVals = 2  # number of possible emission values (0 and 1 in this case)
  for (i in 1:n.states) {
    for (slotNum in 1:nSlots) {
      # set.emission.params.qhmm(hmm, i, c(0.1, 0.1), slot = slotNum) # dgamma params
      set.emission.params.qhmm(hmm, i, diff(c(0, sort(runif(nEmissionVals - 1)), 1)), slot = slotNum)
      set.emission.option.qhmm(hmm, i, "offset", 0, slot = slotNum)
    }
  }
  
  return(hmm)
}

# 11. 3 state HMM, 2 slots, 1 covar
#
# state paths:
#    PNEG -> PG.B -> PG.D -> PNEG
#
#  state    #pro  
#   PNEG     [1]   
#   PG.B     [2]   
#   PG.D     [3]   
#
# #pro: [1] background level (gamma)
#       [2] gene body level (gamma)
#       [3] after gene level (gamma)
#
hmm3states2slots1covarContinuous <- function() {
  
  # PNEG, PG.B, PG.D
  n.states = 3
  PNEG = 1
  PG.B = 2
  PG.D = 3
  
  # valid transitions
  vtbl = NULL
  # tnames = c(rep("autocorr", n.states))
  tnames = c("autocorr_covar", rep("autocorr", n.states - 1))
  
  # valid transitions
  vtbl = matrix(data=0, nrow=n.states, ncol=n.states)
  
  vtbl[PNEG, PNEG] = 1
  vtbl[PNEG, PG.B] = 2
  
  vtbl[PG.B, PG.B] = 1
  vtbl[PG.B, PG.D] = 2
  
  vtbl[PG.D, PG.D] = 1
  vtbl[PG.D, PNEG] = 2
  
  # transition groups
  tgrps = NULL
  tgrps = list(2:n.states)  # share distribution; it's a poor fit for all, but want to avoid getting a tiny 'T' state

  # emission groups
  nSlots = 2
  # egrps = new.emission.groups(n.states, nSlots)
  # egrps = add.emission.groups(egrps, states = c(2, 3, 4), slots = c(2, 2, 2)) # share GROseq over T, P and D [scaled]
  egrps = NULL
  
  tdist = "gamma"
  
  hmm = new.qhmm(list(c(1, 1), 1), 
                 vtbl,
                 tnames,
                 list(c(tdist, tdist),      # PNEG
                      c(tdist, tdist),      # PG.B
                      c(tdist, tdist)),     # PG.D
                 transition.groups = tgrps,
                 emission.groups = egrps)
  
  # set initial parameters
  set.initial.probs.qhmm(hmm, c(1, rep(0, n.states - 1))) # start with background
  
  # set transitions
  set.transition.params.qhmm(hmm, 2:n.states, 0.99)
  
  # Randomly generating n numbers that sum to 1 
  # (where n is the number of possible state transitions)
  # https://stackoverflow.com/questions/11003967/generate-3-random-number-that-sum-to-1-in-r
  # rowSums(tm>0)[i] gives the number of transitions for state i
  # diff(c(0, sort(runif(n-1)), 1)) gives n random numbers that sum to 1
  # (this gives an error if n < 1, so makes sure eash state always has at least one transition)
  # diff(c(0, sort(runif((rowSums(tm>0)[i])-1)), 1)) gives us a random number for each state, all adding to 1
  # for (i in 1:n.states) {
  #  set.transition.params.qhmm(hmm, i, diff(c(0, sort(runif((rowSums(vtbl>0)[i])-1)), 1)))
  #}
  
  # Set emissions
  
  # Alternatively, set emission params randomly
  # (see notes above for what diff(c(0, sort(runif(n - 1)), 1)) is doing)
  nEmissionVals = 2  # number of possible emission values (0 and 1 in this case)
  for (i in 1:n.states) {
    for (slotNum in 1:nSlots) {
      set.emission.params.qhmm(hmm, i, diff(c(0, sort(runif(nEmissionVals - 1)), 1)), slot = slotNum)
      # set.emission.option.qhmm(hmm, i, "offset", 0, slot = slotNum)
    }
  }
  
  return(hmm)
}

# 4 states HMM, 1 slot, 3 covars
#
# state paths:
#    PNEG -> PG.B -> PG.D ->PNEG
#
#  state    #pr
#   PNEG     [1]
#   PG.B     [2]
#   PG.D     [3]
#   US       [4]
#
# #pro: [1] background level (gamma)
#       [2] gene body level (gamma)
#       [3] after gene level (gamma)
#       [4] unstable TSS (gamma)
#
hmm4states1slot3covarsContinuous <- function() {
  
  poisson.decay = FALSE
  
  # PNEG, PG.B, PG.D. US
  n.states = 4
  PNEG = 1
  PG.B = 2
  PG.D = 3
  US = 4
  
  # valid transitions
  nCovars = 3
  vtbl = NULL
  # tnames = c(rep("autocorr", n.states))
  tnames = c("autocorr_covar", "autocorr_covar", "autocorr", "autocorr")
  
  # valid transitions
  vtbl = matrix(data=0, nrow=n.states, ncol=n.states)
  
  vtbl[PNEG, PNEG] = 1
  vtbl[PNEG, PG.B] = 2
  vtbl[PNEG, US] = 3
  
  vtbl[PG.B, PG.B] = 1
  vtbl[PG.B, PG.D] = 2
  
  vtbl[PG.D, PG.D] = 1
  vtbl[PG.D, PNEG] = 2
  
  vtbl[US, US] = 1
  vtbl[US, PNEG] = 2
  
  # transition groups
  tgrps = NULL
  # tgrps = list(2:n.states)  # share distribution; it's a poor fit for all, but want to avoid getting a tiny 'T' state
  
  # emission groups
  nSlots = 1
  egrps = NULL
  if (!poisson.decay) {
    egrps = new.emission.groups(n.states, nSlots)
    egrps = add.emission.groups(egrps, states = c(2, 3), slots = c(1, 1)) # share GROseq over B and D [scaled]
  }
  
  tdist = "gamma"
  ddist = tdist
  if (poisson.decay) {
    ddist = "poisson"
    
    hmm = new.qhmm(list(rep(1, nSlots), rep(1, nCovars)),
                   vtbl,
                   tnames,
                   list(tdist,      # PNEG
                        tdist,      # B
                        ddist),     # D
                   transition.groups = tgrps,
                   emission.groups = egrps)
  } else {
    hmm <- new.qhmm(list(rep(1, nSlots), rep(1, nCovars)),
                    vtbl,
                    tnames,
                    rep(list(rep(tdist, nSlots)), n.states),
                    transition.groups = tgrps,
                    emission.groups = egrps)
  }
  
  # set initial parameters
  set.initial.probs.qhmm(hmm, c(1, rep(0, n.states - 1))) # start with background
  
  # set covariate slots
  set.transition.covars.qhmm(hmm, 1, 1)
  set.transition.covars.qhmm(hmm, 1, 2)
  set.transition.covars.qhmm(hmm, 2, 3)
  
  # set transitions
  set.transition.params.qhmm(hmm, 3:n.states, 0.99)
  # set.transition.params.qhmm(hmm, 2:n.states, 0.99)
  # set.transition.params.qhmm(hmm, 3, 0.5, fixed = TRUE)
  
  # Randomly generating n numbers that sum to 1 
  # (where n is the number of possible state transitions)
  # https://stackoverflow.com/questions/11003967/generate-3-random-number-that-sum-to-1-in-r
  # rowSums(tm>0)[i] gives the number of transitions for state i
  # diff(c(0, sort(runif(n-1)), 1)) gives n random numbers that sum to 1
  # (this gives an error if n < 1, so makes sure eash state always has at least one transition)
  # diff(c(0, sort(runif((rowSums(tm>0)[i])-1)), 1)) gives us a random number for each state, all adding to 1
  # for (i in 1:n.states) {
  #  set.transition.params.qhmm(hmm, i, diff(c(0, sort(runif((rowSums(vtbl>0)[i])-1)), 1)))
  #}
  
  # Set emissions
  
  # Alternatively, set emission params randomly
  # (see notes above for what diff(c(0, sort(runif(n - 1)), 1)) is doing)
  # //--- add , fixed = c(T, T) - see hmm.parse.R in grocap-tss
  nEmissionVals = 2  # number of possible emission values (0 and 1 in this case)
  if (poisson.decay) {
    set.emission.params.qhmm(hmm, 1, c(1, 1), slot = 1)
    set.emission.params.qhmm(hmm, 2, c(1, 1), slot = 1)
    set.emission.params.qhmm(hmm, 3, 0.1, slot = 1) # poisson params
  } else {
    for (i in 1:n.states) {
      for (slotNum in 1:nSlots) {
        # set.emission.params.qhmm(hmm, i, c(0.1, 0.1), slot = slotNum) # dgamma params
        set.emission.params.qhmm(hmm, i, diff(c(0, sort(runif(nEmissionVals - 1)), 1)), slot = slotNum)
        # set.emission.option.qhmm(hmm, i, "offset", 0, slot = slotNum)
      }
    }
  }
  
  return(hmm)
}

# 4 states HMM, 1 slot, 2 covars
#
# state paths:
#    PNEG -> PG.B -> PG.D ->PNEG
#
#  state    #pr
#   PNEG     [1]
#   PG.B     [2]
#   PG.D     [3]
#   US       [4]
#
# #pro: [1] background level (gamma)
#       [2] gene body level (gamma)
#       [3] after gene level (gamma)
#       [4] unstable TSS (gamma)
#
hmm4states1slot2covarsContinuous <- function() {
  
  poisson.decay = FALSE
  
  # PNEG, PG.B, PG.D. US
  n.states = 4
  PNEG = 1
  PG.B = 2
  PG.D = 3
  US = 4
  
  # valid transitions
  nCovars = 2
  vtbl = NULL
  # tnames = c(rep("autocorr", n.states))
  tnames = c("autocorr_covar", "autocorr_covar", "autocorr", "autocorr")
  
  # valid transitions
  vtbl = matrix(data=0, nrow=n.states, ncol=n.states)
  
  vtbl[PNEG, PNEG] = 1
  vtbl[PNEG, PG.B] = 2
  vtbl[PNEG, US] = 3
  
  vtbl[PG.B, PG.B] = 1
  vtbl[PG.B, PG.D] = 2
  
  vtbl[PG.D, PG.D] = 1
  vtbl[PG.D, PNEG] = 2
  
  vtbl[US, US] = 1
  vtbl[US, PNEG] = 2
  
  # transition groups
  tgrps = NULL
  # tgrps = list(2:n.states)  # share distribution; it's a poor fit for all, but want to avoid getting a tiny 'T' state
  
  # emission groups
  nSlots = 1
  egrps = NULL
  if (!poisson.decay) {
    egrps = new.emission.groups(n.states, nSlots)
    egrps = add.emission.groups(egrps, states = c(2, 3), slots = c(1, 1)) # share GROseq over B and D [scaled]
  }
  
  tdist = "gamma"
  ddist = tdist
  if (poisson.decay) {
    ddist = "poisson"
    
    hmm = new.qhmm(list(rep(1, nSlots), rep(1, nCovars)),
                   vtbl,
                   tnames,
                   list(tdist,      # PNEG
                        tdist,      # B
                        ddist),     # D
                   transition.groups = tgrps,
                   emission.groups = egrps)
  } else {
    hmm <- new.qhmm(list(rep(1, nSlots), rep(1, nCovars)),
                    vtbl,
                    tnames,
                    rep(list(rep(tdist, nSlots)), n.states),
                    transition.groups = tgrps,
                    emission.groups = egrps)
  }
  
  # set initial parameters
  set.initial.probs.qhmm(hmm, c(1, rep(0, n.states - 1))) # start with background
  
  # set covariate slots
  set.transition.covars.qhmm(hmm, 1, 1)
  set.transition.covars.qhmm(hmm, 2, 2)
  
  # set transitions
  set.transition.params.qhmm(hmm, 3:n.states, 0.99)
  # set.transition.params.qhmm(hmm, 2:n.states, 0.99)
  # set.transition.params.qhmm(hmm, 3, 0.5, fixed = TRUE)
  
  # Randomly generating n numbers that sum to 1 
  # (where n is the number of possible state transitions)
  # https://stackoverflow.com/questions/11003967/generate-3-random-number-that-sum-to-1-in-r
  # rowSums(tm>0)[i] gives the number of transitions for state i
  # diff(c(0, sort(runif(n-1)), 1)) gives n random numbers that sum to 1
  # (this gives an error if n < 1, so makes sure eash state always has at least one transition)
  # diff(c(0, sort(runif((rowSums(tm>0)[i])-1)), 1)) gives us a random number for each state, all adding to 1
  # for (i in 1:n.states) {
  #  set.transition.params.qhmm(hmm, i, diff(c(0, sort(runif((rowSums(vtbl>0)[i])-1)), 1)))
  #}
  
  # Set emissions
  
  # Alternatively, set emission params randomly
  # (see notes above for what diff(c(0, sort(runif(n - 1)), 1)) is doing)
  # //--- add , fixed = c(T, T) - see hmm.parse.R in grocap-tss
  nEmissionVals = 2  # number of possible emission values (0 and 1 in this case)
  if (poisson.decay) {
    set.emission.params.qhmm(hmm, 1, c(1, 1), slot = 1)
    set.emission.params.qhmm(hmm, 2, c(1, 1), slot = 1)
    set.emission.params.qhmm(hmm, 3, 0.1, slot = 1) # poisson params
  } else {
    for (i in 1:n.states) {
      for (slotNum in 1:nSlots) {
        # set.emission.params.qhmm(hmm, i, c(0.1, 0.1), slot = slotNum) # dgamma params
        set.emission.params.qhmm(hmm, i, diff(c(0, sort(runif(nEmissionVals - 1)), 1)), slot = slotNum)
        # set.emission.option.qhmm(hmm, i, "offset", 0, slot = slotNum)
      }
    }
  }
  
  return(hmm)
}

# 4 states HMM, 1 slot, 1 covar
#
# state paths:
#    PNEG -> PG.B -> PG.D ->PNEG
#
#  state    #pr
#   PNEG     [1]
#   PG.B     [2]
#   PG.D     [3]
#   US       [4]
#
# #pro: [1] background level (gamma)
#       [2] gene body level (gamma)
#       [3] after gene level (gamma)
#       [4] unstable TSS (gamma)
#
hmm4states1slot1covarContinuous <- function() {
  
  poisson.decay = FALSE
  
  # PNEG, PG.B, PG.D. US
  n.states = 4
  PNEG = 1
  PG.B = 2
  PG.D = 3
  US = 4
  
  # valid transitions
  nCovars = 1
  vtbl = NULL
  # tnames = c(rep("autocorr", n.states))
  tnames = c("autocorr_covar", "autocorr", "autocorr", "autocorr")
  
  # valid transitions
  vtbl = matrix(data=0, nrow=n.states, ncol=n.states)
  
  vtbl[PNEG, PNEG] = 1
  vtbl[PNEG, PG.B] = 2
  vtbl[PNEG, US] = 3
  
  vtbl[PG.B, PG.B] = 1
  vtbl[PG.B, PG.D] = 2
  
  vtbl[PG.D, PG.D] = 1
  vtbl[PG.D, PNEG] = 2
  
  vtbl[US, US] = 1
  vtbl[US, PNEG] = 2
  
  # transition groups
  tgrps = NULL
  # tgrps = list(2:n.states)  # share distribution; it's a poor fit for all, but want to avoid getting a tiny 'T' state
  
  # emission groups
  nSlots = 1
  egrps = NULL
  if (!poisson.decay) {
    egrps = new.emission.groups(n.states, nSlots)
    egrps = add.emission.groups(egrps, states = c(2, 3), slots = c(1, 1)) # share GROseq over B and D [scaled]
  }
  
  tdist = "gamma"
  ddist = tdist
  if (poisson.decay) {
    ddist = "poisson"
    
    hmm = new.qhmm(list(rep(1, nSlots), rep(1, nCovars)),
                   vtbl,
                   tnames,
                   list(tdist,      # PNEG
                        tdist,      # B
                        ddist),     # D
                   transition.groups = tgrps,
                   emission.groups = egrps)
  } else {
    hmm <- new.qhmm(list(rep(1, nSlots), rep(1, nCovars)),
                    vtbl,
                    tnames,
                    rep(list(rep(tdist, nSlots)), n.states),
                    transition.groups = tgrps,
                    emission.groups = egrps)
  }
  
  # set initial parameters
  set.initial.probs.qhmm(hmm, c(1, rep(0, n.states - 1))) # start with background
  
  # set covariate slots
  set.transition.covars.qhmm(hmm, 1, 1)

  # set transitions
  # set.transition.params.qhmm(hmm, 3:n.states, 0.99)
  set.transition.params.qhmm(hmm, 2:n.states, 0.99)
  # set.transition.params.qhmm(hmm, 3, 0.5, fixed = TRUE)
  
  # Randomly generating n numbers that sum to 1 
  # (where n is the number of possible state transitions)
  # https://stackoverflow.com/questions/11003967/generate-3-random-number-that-sum-to-1-in-r
  # rowSums(tm>0)[i] gives the number of transitions for state i
  # diff(c(0, sort(runif(n-1)), 1)) gives n random numbers that sum to 1
  # (this gives an error if n < 1, so makes sure eash state always has at least one transition)
  # diff(c(0, sort(runif((rowSums(tm>0)[i])-1)), 1)) gives us a random number for each state, all adding to 1
  # for (i in 1:n.states) {
  #  set.transition.params.qhmm(hmm, i, diff(c(0, sort(runif((rowSums(vtbl>0)[i])-1)), 1)))
  #}
  
  # Set emissions
  
  # Alternatively, set emission params randomly
  # (see notes above for what diff(c(0, sort(runif(n - 1)), 1)) is doing)
  # //--- add , fixed = c(T, T) - see hmm.parse.R in grocap-tss
  nEmissionVals = 2  # number of possible emission values (0 and 1 in this case)
  if (poisson.decay) {
    set.emission.params.qhmm(hmm, 1, c(1, 1), slot = 1)
    set.emission.params.qhmm(hmm, 2, c(1, 1), slot = 1)
    set.emission.params.qhmm(hmm, 3, 0.1, slot = 1) # poisson params
  } else {
    for (i in 1:n.states) {
      for (slotNum in 1:nSlots) {
        # set.emission.params.qhmm(hmm, i, c(0.1, 0.1), slot = slotNum) # dgamma params
        set.emission.params.qhmm(hmm, i, diff(c(0, sort(runif(nEmissionVals - 1)), 1)), slot = slotNum)
        # set.emission.option.qhmm(hmm, i, "offset", 0, slot = slotNum)
      }
    }
  }
  
  return(hmm)
}

# 3 state HMM, 1 slot, 3 covars
#
# state paths:
#    PNEG -> PG.B -> PG.D ->PNEG
#
#  state    #pro
#   PNEG     [1]
#   PG.B     [2]
#   PG.D     [3]
#
# #pro: [1] background level (gamma)
#       [2] gene body level (gamma)
#       [3] after gene level (gamma / poisson)
#
hmm3states1slot3covarsContinuous <- function() {
  
  poisson.decay = FALSE
  
  # PNEG, PG.B, PG.D
  n.states = 3
  PNEG = 1
  PG.B = 2
  PG.D = 3
  
  # valid transitions
  nCovars = 3
  vtbl = NULL
  # tnames = c(rep("autocorr", n.states))
  tnames = c("autocorr_covar", "autocorr_covar", "autocorr_covar")
  
  # valid transitions
  vtbl = matrix(data=0, nrow=n.states, ncol=n.states)
  
  vtbl[PNEG, PNEG] = 1
  vtbl[PNEG, PG.B] = 2
  
  vtbl[PG.B, PG.B] = 1
  vtbl[PG.B, PG.D] = 2
  
  vtbl[PG.D, PG.D] = 1
  vtbl[PG.D, PNEG] = 2
  
  # transition groups
  tgrps = NULL
  # tgrps = list(2:n.states)  # share distribution; it's a poor fit for all, but want to avoid getting a tiny 'T' state
  
  # emission groups
  nSlots = 1
  egrps = NULL
  if (!poisson.decay) {
    egrps = new.emission.groups(n.states, nSlots)
    egrps = add.emission.groups(egrps, states = c(2, 3), slots = c(1, 1)) # share GROseq over B and D [scaled]
  }
  
  tdist = "gamma"
  ddist = tdist
  if (poisson.decay) {
    ddist = "poisson"
    
    hmm = new.qhmm(list(rep(1, nSlots), rep(1, nCovars)),
                   vtbl,
                   tnames,
                   list(tdist,      # PNEG
                        tdist,      # B
                        ddist),     # D
                   transition.groups = tgrps,
                   emission.groups = egrps)
  } else {
    hmm <- new.qhmm(list(rep(1, nSlots), rep(1, nCovars)),
                    vtbl,
                    tnames,
                    rep(list(rep(tdist, nSlots)), n.states),
                    transition.groups = tgrps,
                    emission.groups = egrps)
  }
  
  # set initial parameters
  set.initial.probs.qhmm(hmm, c(1, rep(0, n.states - 1))) # start with background
  
  # set covariate slots
  set.transition.covars.qhmm(hmm, 1, 1)
  set.transition.covars.qhmm(hmm, 2, 2)
  set.transition.covars.qhmm(hmm, 3, 3)
  
  # set transitions
  # set.transition.params.qhmm(hmm, 3:n.states, 0.99)
  # set.transition.params.qhmm(hmm, 2:n.states, 0.99)
  # set.transition.params.qhmm(hmm, 3, 0.5, fixed = TRUE)
  
  # Randomly generating n numbers that sum to 1 
  # (where n is the number of possible state transitions)
  # https://stackoverflow.com/questions/11003967/generate-3-random-number-that-sum-to-1-in-r
  # rowSums(tm>0)[i] gives the number of transitions for state i
  # diff(c(0, sort(runif(n-1)), 1)) gives n random numbers that sum to 1
  # (this gives an error if n < 1, so makes sure eash state always has at least one transition)
  # diff(c(0, sort(runif((rowSums(tm>0)[i])-1)), 1)) gives us a random number for each state, all adding to 1
  # for (i in 1:n.states) {
  #  set.transition.params.qhmm(hmm, i, diff(c(0, sort(runif((rowSums(vtbl>0)[i])-1)), 1)))
  #}
  
  # Set emissions
  
  # Alternatively, set emission params randomly
  # (see notes above for what diff(c(0, sort(runif(n - 1)), 1)) is doing)
  # //--- add , fixed = c(T, T) - see hmm.parse.R in grocap-tss
  nEmissionVals = 2  # number of possible emission values (0 and 1 in this case)
  if (poisson.decay) {
    set.emission.params.qhmm(hmm, 1, c(1, 1), slot = 1)
    set.emission.params.qhmm(hmm, 2, c(1, 1), slot = 1)
    set.emission.params.qhmm(hmm, 3, 0.1, slot = 1) # poisson params
  } else {
    for (i in 1:n.states) {
      for (slotNum in 1:nSlots) {
        # set.emission.params.qhmm(hmm, i, c(0.1, 0.1), slot = slotNum) # dgamma params
        set.emission.params.qhmm(hmm, i, diff(c(0, sort(runif(nEmissionVals - 1)), 1)), slot = slotNum)
        # set.emission.option.qhmm(hmm, i, "offset", 0, slot = slotNum)
      }
    }
  }
  
  return(hmm)
}

# 3 state HMM, 1 slot, 2 covars
#
# state paths:
#    PNEG -> PG.B -> PG.D ->PNEG
#
#  state    #pro
#   PNEG     [1]
#   PG.B     [2]
#   PG.D     [3]
#
# #pro: [1] background level (gamma)
#       [2] gene body level (gamma)
#       [3] after gene level (gamma / poisson)
#
hmm3states1slot2covarsContinuous <- function() {
  
  poisson.decay = FALSE
  
  # PNEG, PG.B, PG.D
  n.states = 3
  PNEG = 1
  PG.B = 2
  PG.D = 3
  
  # valid transitions
  nCovars = 2
  vtbl = NULL
  # tnames = c(rep("autocorr", n.states))
  tnames = c("autocorr_covar", "autocorr_covar", "autocorr")
  
  # valid transitions
  vtbl = matrix(data=0, nrow=n.states, ncol=n.states)
  
  vtbl[PNEG, PNEG] = 1
  vtbl[PNEG, PG.B] = 2
  
  vtbl[PG.B, PG.B] = 1
  vtbl[PG.B, PG.D] = 2
  vtbl[PG.B, PNEG] = 3
  
  vtbl[PG.D, PG.D] = 1
  vtbl[PG.D, PNEG] = 2
  
  # transition groups
  tgrps = NULL
  # tgrps = list(2:n.states)  # share distribution; it's a poor fit for all, but want to avoid getting a tiny 'T' state
  
  # emission groups
  nSlots = 1
  egrps = NULL
  if (!poisson.decay) {
    egrps = new.emission.groups(n.states, nSlots)
    egrps = add.emission.groups(egrps, states = c(2, 3), slots = c(1, 1)) # share GROseq over B and D [scaled]
  }
  
  tdist = "gamma"
  ddist = tdist
  if (poisson.decay) {
    ddist = "poisson"
    
    hmm = new.qhmm(list(rep(1, nSlots), rep(1, nCovars)),
                   vtbl,
                   tnames,
                   list(tdist,      # PNEG
                        tdist,      # B
                        ddist),     # D
                   transition.groups = tgrps,
                   emission.groups = egrps)
  } else {
    hmm <- new.qhmm(list(rep(1, nSlots), rep(1, nCovars)),
                    vtbl,
                    tnames,
                    rep(list(rep(tdist, nSlots)), n.states),
                    transition.groups = tgrps,
                    emission.groups = egrps)
  }
  
  # set initial parameters
  set.initial.probs.qhmm(hmm, c(1, rep(0, n.states - 1))) # start with background
  
  # set covariate slots
  # set.transition.covars(hmm, state, covarIdx)
  set.transition.covars.qhmm(hmm, 1, 1)
  set.transition.covars.qhmm(hmm, 2, 2)
  
  # set transitions
  set.transition.params.qhmm(hmm, 3:n.states, 0.99)
  # set.transition.params.qhmm(hmm, 2:n.states, 0.99)
  # set.transition.params.qhmm(hmm, 3, 0.5, fixed = TRUE)
  
  # Randomly generating n numbers that sum to 1 
  # (where n is the number of possible state transitions)
  # https://stackoverflow.com/questions/11003967/generate-3-random-number-that-sum-to-1-in-r
  # rowSums(tm>0)[i] gives the number of transitions for state i
  # diff(c(0, sort(runif(n-1)), 1)) gives n random numbers that sum to 1
  # (this gives an error if n < 1, so makes sure eash state always has at least one transition)
  # diff(c(0, sort(runif((rowSums(tm>0)[i])-1)), 1)) gives us a random number for each state, all adding to 1
  # for (i in 1:n.states) {
  #  set.transition.params.qhmm(hmm, i, diff(c(0, sort(runif((rowSums(vtbl>0)[i])-1)), 1)))
  #}
  
  # Set emissions
  
  # Alternatively, set emission params randomly
  # (see notes above for what diff(c(0, sort(runif(n - 1)), 1)) is doing)
  # //--- add , fixed = c(T, T) - see hmm.parse.R in grocap-tss
  nEmissionVals = 2  # number of possible emission values (0 and 1 in this case)
  if (poisson.decay) {
    set.emission.params.qhmm(hmm, 1, c(1, 1), slot = 1)
    set.emission.params.qhmm(hmm, 2, c(1, 1), slot = 1)
    set.emission.params.qhmm(hmm, 3, 0.1, slot = 1) # poisson params
  } else {
    for (i in 1:n.states) {
      for (slotNum in 1:nSlots) {
        # set.emission.params.qhmm(hmm, i, c(0.1, 0.1), slot = slotNum) # dgamma params
        set.emission.params.qhmm(hmm, i, diff(c(0, sort(runif(nEmissionVals - 1)), 1)), slot = slotNum)
        # set.emission.option.qhmm(hmm, i, "offset", 0, slot = slotNum)
      }
    }
  }
  
  return(hmm)
}

# 11. 3 state HMM, 1 slot, 1 covar
#
# state paths:
#    PNEG -> PG.B -> PG.D ->PNEG
#
#  state    #pro
#   PNEG     [1]
#   PG.B     [2]
#   PG.D     [3]
#
# #pro: [1] background level (gamma)
#       [2] gene body level (gamma)
#       [3] after gene level (gamma / poisson)
#
hmm3states1slot1covarContinuous <- function() {
  
  poisson.decay = FALSE
  
  # PNEG, PG.B, PG.D
  n.states = 3
  PNEG = 1
  PG.B = 2
  PG.D = 3
  
  # valid transitions
  nCovars = 1
  vtbl = NULL
  # tnames = c(rep("autocorr", n.states))
  tnames = c("autocorr_covar", rep("autocorr", n.states - 1))
  
  # valid transitions
  vtbl = matrix(data=0, nrow=n.states, ncol=n.states)
  
  vtbl[PNEG, PNEG] = 1
  vtbl[PNEG, PG.B] = 2
  
  vtbl[PG.B, PG.B] = 1
  vtbl[PG.B, PG.D] = 2
  # vtbl[PG.B, PNEG] = 3
  
  vtbl[PG.D, PG.D] = 1
  vtbl[PG.D, PNEG] = 2
  
  # transition groups
  tgrps = NULL
  tgrps = list(2:n.states)  # share distribution; it's a poor fit for all, but want to avoid getting a tiny 'T' state
  
  # emission groups
  nSlots = 1
  egrps = NULL
  if (!poisson.decay) {
    egrps = new.emission.groups(n.states, nSlots)
    egrps = add.emission.groups(egrps, states = c(2, 3), slots = c(1, 1)) # share GROseq over B and D [scaled]
  }
  
  tdist = "gamma"
  ddist = tdist
  if (poisson.decay) {
    ddist = "poisson"
    
    hmm = new.qhmm(list(rep(1, nSlots), rep(1, nCovars)),
                   vtbl,
                   tnames,
                   list(tdist,      # PNEG
                        tdist,      # B
                        ddist),     # D
                   transition.groups = tgrps,
                   emission.groups = egrps)
  } else {
    hmm <- new.qhmm(list(rep(1, nSlots), rep(1, nCovars)),
                    vtbl,
                    tnames,
                    rep(list(rep(tdist, nSlots)), n.states),
                    transition.groups = tgrps,
                    emission.groups = egrps)
  }
  
  # set initial parameters
  set.initial.probs.qhmm(hmm, c(1, rep(0, n.states - 1))) # start with background
  
  # set transitions
  set.transition.params.qhmm(hmm, 2:n.states, 0.998)
  # set.transition.params.qhmm(hmm, 3, 0.5, fixed = TRUE)
  
  # Randomly generating n numbers that sum to 1 
  # (where n is the number of possible state transitions)
  # https://stackoverflow.com/questions/11003967/generate-3-random-number-that-sum-to-1-in-r
  # rowSums(tm>0)[i] gives the number of transitions for state i
  # diff(c(0, sort(runif(n-1)), 1)) gives n random numbers that sum to 1
  # (this gives an error if n < 1, so makes sure eash state always has at least one transition)
  # diff(c(0, sort(runif((rowSums(tm>0)[i])-1)), 1)) gives us a random number for each state, all adding to 1
  # for (i in 1:n.states) {
  #  set.transition.params.qhmm(hmm, i, diff(c(0, sort(runif((rowSums(vtbl>0)[i])-1)), 1)))
  #}
  
  # Set emissions
  
  # Alternatively, set emission params randomly
  # (see notes above for what diff(c(0, sort(runif(n - 1)), 1)) is doing)
  # //--- add , fixed = c(T, T) - see hmm.parse.R in grocap-tss
  nEmissionVals = 2  # number of possible emission values (0 and 1 in this case)
  if (poisson.decay) {
    set.emission.params.qhmm(hmm, 1, c(1, 1), slot = 1)
    set.emission.params.qhmm(hmm, 2, c(1, 1), slot = 1)
    set.emission.params.qhmm(hmm, 3, 0.1, slot = 1) # poisson params
  } else {
    for (i in 1:n.states) {
      for (slotNum in 1:nSlots) {
        # set.emission.params.qhmm(hmm, i, c(0.1, 0.1), slot = slotNum) # dgamma params
        set.emission.params.qhmm(hmm, i, diff(c(0, sort(runif(nEmissionVals - 1)), 1)), slot = slotNum)
        # set.emission.option.qhmm(hmm, i, "offset", 0, slot = slotNum)
      }
    }
  }
  
  return(hmm)
}

