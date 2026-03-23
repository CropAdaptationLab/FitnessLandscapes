# Title: MAPPING POPULATIONS
# Author: Ted Monyak
# Description: Create Mapping populations
# Assumes that CreateFounderPop and CreateIndepdendentPops have been called

# Creates a biparental recombinant inbred line (RIL) population with n.RILFams families
# pop1: the first population to sample an individual from
# pop2: the second population to sample an individual from.
# if FALSE, will cross two individuals from pop1
# Returns: a population where the first two individuals are parent1 and parent2,
# and the rest is the RIL
createRIL <- function(pop1, pop2) {
  # Select 1 random individual from pop1 and pop2
  parent1 <- pop1[sample.int(nInd(pop1),1)]
  parent2 <- pop2[sample.int(nInd(pop2),1)]

  # Cross parent A with parent B, and create n.RILs progeny
  F1 <- randCross2(parent1,
                   parent2,
                   nCrosses=n.RILs,
                   ignoreSexes = TRUE)

  # Create F10s of each RIL with SSD
  F2 <- self(F1, nProgeny = n.indPerRIL)
  F3 <- self(F2)
  F4 <- self(F3)
  F5 <- self(F4)
  F6 <- self(F5)
  F7 <- self(F6)
  F8 <- self(F7)
  F9 <- self(F8)
  RIL <- self(F9)
  return (c(parent1, parent2, RIL))
}

# Simulate the development of purelines
# landrace: One of the landrace populations
# Returns: a list of two pureline populations
makePurelinesBulk <- function(landrace) {
  # One round of crossing
  F1 <- randCross(landrace, nCrosses=n.F1)
  F2 <- self(F1, nProgeny = 4)
  F2 <- selectInd(F2, trait=selectionPheno, nInd=n.F2)
  F3 <- self(F2, nProgeny = 2)
  F3 <- selectInd(F3, trait=selectionPheno, nInd=n.F3)
  F4 <- self(F3, nProgeny = 2)
  F4 <- selectInd(F4, trait=selectionPheno, nInd=n.F4)
  F5 <- self(F4, nProgeny = 4)
  # Divide the F5 into 2, to create 2 purelines
  F5_1 <- F5[1:(nInd(F5)/2)]
  F5_2 <- F5[(nInd(F5)/2+1):nInd(F5)]
  
  # Create 1 pureline from each F5 population
  purifyF5 <- function(F5) {
    F5 <- selectInd(F5, trait=selectionPheno, nInd=n.F5)
    F6 <- self(F5, nProgeny = 2)
    F6 <- selectInd(F6, trait=selectionPheno, nInd=n.F6)
    F7 <- self(F6, nProgeny = 2)
    F7 <- selectInd(F7, trait=selectionPheno, nInd=n.F7)
    F8 <- self(F7, nProgeny = 4)
    F8 <- selectInd(F8, trait=selectionPheno, nInd=n.F8)
    F9 <- self(F8, nProgeny = 4)
    F9 <- selectInd(F9, trait=selectionPheno, nInd=n.F9)
    F10 <- self(F9, nProgeny = 4)
    F10 <- selectInd(F10, trait=selectionPheno, nInd=1)
    # Bulk the pureline
    pureline <- self(F10, nProgeny=50)
    pureline <- self(pureline, nProgeny=10)
    # Ensure the pureline is immortalized
    for (f in 1:8) {
      pureline <- self(pureline)
    }
    return (pureline)
  }
  
  # Create two purelines
  pureline_1 <- purifyF5(F5_1)
  pureline_2 <- purifyF5(F5_2)
  return (list(pureline_1, pureline_2))
}

# NOT BEING USED RIGHT NOW
# Creates a nested association mapping (NAM) population, to be used for GWAS
# Assumes that pops - a list of populations, exists
createNAM <- function() {
  # The "reference" population is the first population
  refPop <- pops[[1]]
  # Randomly select an individual from the reference population
  refLine <- refPop[sample.int(nInd(refPop),1)]
  # Make the reference line inbred by selfing it for 10 generations
  for (f in 1:10) {
    refLine <- self(refLine)
  }
  # The rest of the founder lines will come from the other populations in pops
  founderLines <- c()
  # Iterate through the rest of the populations and select a random individual
  # to become a founder line
  for (p in 2:n.nPops) {
    subPop <- pops[[p]]
    ind <- subPop[sample.int(nInd(subPop),1)]
    for (f in 1:10) {
      ind <- self(ind)
    }
    founderLines <- append(founderLines, ind)
  }
  
  # Cross each of the founder lines by the reference line
  crossPlan = matrix(c(rep(1:(n.nPops-1)), rep(1,(n.nPops-1))), nrow=n.nPops-1, ncol=2)
  
  F1 <- makeCross2(founderLines, refLine, crossPlan, simParam=SP)
  # Create F10s with SSD
  F2 <- self(F1, nProgeny=n.RILs)
  F3 <- self(F2)
  F4 <- self(F3)
  F5 <- self(F4)
  F6 <- self(F5)
  F7 <- self(F6)
  F8 <- self(F7)
  F9 <- self(F8)
  NAM <- self(F9)
  return (NAM)
}

