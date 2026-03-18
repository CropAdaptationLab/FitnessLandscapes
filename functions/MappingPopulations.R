# Title: MAPPING POPULATIONS
# Author: Ted Monyak
# Description: Create Mapping populations
# Assumes that CreateFounderPop and CreateIndepdendentPops have been called

# Creates a biparent1l recombinant inbred line (RIL) population with n.RILFams families
# pop1: the first population to sample an individual from
# pop2: the second population to sample an individual from.
# save_dir: directory to write plots t
# if FALSE, will cross two individuals from pop1
# Returns: a population where the first two individuals are parent1 and parent2,
# and the rest is the RIL
createRIL <- function(pop1, pop2, save_dir) {
  # Select 1 random individual from pop1 and pop2
  parent1 <- pop1[sample.int(nInd(pop1),1)]
  parent2 <- pop2[sample.int(nInd(pop2),1)]

  # Cross parent A with parent B, and create n.RILFams progeny
  F1 <- randCross2(parent1,
                   parent2,
                   nCrosses=n.RILFams,
                   ignoreSexes = TRUE)

  # Create F10s of each RIL Family with SSD
  F2 <- self(F1, nProgeny = n.indPerRILFam)
  F3 <- self(F2)
  F4 <- self(F3)
  F5 <- self(F4)
  F6 <- self(F5)
  F7 <- self(F6)
  F8 <- self(F7)
  F9 <- self(F8)
  
  # Each RIL Family has n.indPerRILFam replicates
  RIL <- self(F9)
  return (c(parent1, parent2, RIL))
}

# Runs population improvement for n.breedingGens generations
# pop: the initial population
# save_dir: the directory to save the genetic gain chart
# pop_id: the name of the population, for filename purposes
# returns: the improved population
populationImprovement <- function(pop, save_dir, pop_id) {
  # Store the initial mean yield
  yieldVals <- c(meanP(pop)[3])
  # Iterate through n.breedingGens
  for (gen in 1:n.breedingGens) {
    # Change this to select based on genomic selection instead to compare
    # Select the top n.purelines (using a 10% selection criteria)
    pop <- selectCross(pop, trait=fitFunc, nInd=n.purelines, nCrosses=n.purelines*10)
    yieldVals <- c(yieldVals, meanP(pop)[3])
  }

  # Create a line graph of the genetic gain
  if (saveTraitPlots) {
    yield.df <- data.frame(gen=c(1:(n.breedingGens+1)),
                           yield=yieldVals)
    g <- ggplot(yield.df, aes(x=gen, y=yield)) +
      geom_line() +
      labs(title=paste0("Genetic Gain ", pop_id), x="Generation", y="Mean Yield")
    ggplot2::ggsave(filename = paste0("GeneticGain", pop_id, ".pdf"),
                    path=save_dir,
                    device = "pdf",
                    width=10,
                    height=7)
  }

  return (pop)
}

# Simulate the development of purelines
# n.purelines are selected
makePurelinesBulk <- function(landrace) {
  # One round of crossing
  F1 <- randCross(landrace, nCrosses=n.F1)
  F2 <- self(F1, nProgeny = 4)
  F2 <- selectInd(F2, trait=fitFunc, nInd=n.F2)
  F3 <- self(F2, nProgeny = 2)
  F3 <- selectInd(F3, trait=fitFunc, nInd=n.F3)
  F4 <- self(F3, nProgeny = 2)
  F4 <- selectInd(F4, trait=fitFunc, nInd=n.F4)
  F5 <- self(F4, nProgeny = 4)
  F5_a <- F5[1:(nInd(F5)/2)]
  F5_b <- F5[(nInd(F5)/2+1):nInd(F5)]
  
  purifyF5 <- function(F5) {
    F5 <- selectInd(F5, trait=fitFunc, nInd=n.F5)
    F6 <- self(F5, nProgeny = 2)
    F6 <- selectInd(F6, trait=fitFunc, nInd=n.F6)
    F7 <- self(F6, nProgeny = 2)
    F7 <- selectInd(F7, trait=fitFunc, nInd=n.F7)
    F8 <- self(F7, nProgeny = 4)
    F8 <- selectInd(F8, trait=fitFunc, nInd=n.F8)
    F9 <- self(F8, nProgeny = 4)
    F9 <- selectInd(F9, trait=fitFunc, nInd=n.F9)
    F10 <- self(F9, nProgeny = 4)
    F10 <- selectInd(F10, trait=fitFunc, nInd=1)
    pureline <- self(F10, nProgeny=50)
    pureline <- self(pureline, nProgeny=10)
    for (f in 1:8) {
      pureline <- self(pureline)
    }
    return (pureline)
  }
  
  pureline_a <- purifyF5(F5_a)
  pureline_b <- purifyF5(F5_b)
  return (list(pureline_a, pureline_b))
}

# Simulate the development of purelines
# n.purelines are selected
makePurelinesPed <- function(landrace) {
  # One round of crossing
  F1 <- randCross(landrace, nCrosses=n.F1)
  F2 <- self(F1, nProgeny = 4)
  
  F3_seed <- selectFam(F2, trait=fitFunc, nFam=500)
  F3 <- self(F3_seed, nProgeny=2)
  
  # Select top F3 families, self to create F4 seed
  F4_seed <- selectFam(F3, trait=fitFunc, nFam=250)
  F4 <- self(F4_seed)
  
  # Select the top F4 families, self to create F5 seed
  F5_seed <- selectFam(F4, trait=fitFunc, nFam=125)
  F5 <- self(F5_seed)
  
  F6_seed <- selectFam(F5, trait=fitFunc, nFam=60)
  F6 <- self(F6_seed)
  
  F7_seed <- selectFam(F6, trait=fitFunc, nFam=30)
  F7 <- self(F7_seed)
  
  F8_seed <- selectFam(F7, trait=fitFunc, nFam=15)
  F8 <- self(F8_seed, nProgeny=2)
  
  F9_seed <- selectFam(F8, trait=fitFunc, nFam=8)
  F9 <- self(F9_seed, nProgeny=2)
  
  F10_seed <- selectFam(F9, trait=fitFunc, nFam=4)
  F10 <- self(F10_seed, nProgeny=2)
  
  purelines <- selectFam(F10, trait=fitFunc, nFam=2)
  purelines <- self(purelines, nProgeny=8)
  for (f in 1:8) {
    purelines <- self(purelines)
  }
  
  pureline_a <- purelines[1:(nInd(purelines)/2)]
  pureline_b <- purelines[(nInd(purelines)/2+1):nInd(purelines)]
  return (list(pureline_a, pureline_b))
}

# Simulate the development of purelines
# n.purelines are selected
makePurelinesHybrid <- function(landrace) {
  # One round of crossing
  F1 <- randCross(landrace, nCrosses=n.F1)
  F2 <- self(F1, nProgeny = 4)
  F2 <- selectInd(F2, trait=fitFunc, nInd=n.F2)
  F3 <- self(F2, nProgeny = 2)
  F3 <- selectInd(F3, trait=fitFunc, nInd=n.F3)
  F4 <- self(F3, nProgeny = 2)
  
  # Select the top F4 families, self to create F5 seed
  F5_seed <- selectFam(F4, trait=fitFunc, nFam=250)
  F5 <- self(F5_seed)
  
  F6_seed <- selectFam(F5, trait=fitFunc, nFam=125)
  F6 <- self(F6_seed)
  
  F7_seed <- selectFam(F6, trait=fitFunc, nFam=60)
  F7 <- self(F7_seed)
  
  F8_seed <- selectFam(F7, trait=fitFunc, nFam=30)
  F8 <- self(F8_seed, nProgeny=2)
  
  F9_seed <- selectFam(F8, trait=fitFunc, nFam=16)
  F9 <- self(F9_seed, nProgeny=2)
  
  F10_seed <- selectFam(F9, trait=fitFunc, nFam=4)
  F10 <- self(F10_seed, nProgeny=2)
  
  purelines <- selectFam(F10, trait=fitFunc, nFam=2)
  purelines <- self(purelines, nProgeny=16)
  for (f in 1:8) {
    purelines <- self(purelines)
  }
  
  pureline_a <- purelines[1:(nInd(purelines)/2)]
  pureline_b <- purelines[(nInd(purelines)/2+1):nInd(purelines)]
  return (list(pureline_a, pureline_b))
}

# Simulate a population going through a breeding program, under purifying selection
# First, the top landrace individuals are purified
# pop: the landrace
# Returns: an F8 population that has undergone selection and inbreeding
# NOT BEING USED RIGHT NOW
makeElite <- function(pop) {
  purifiedLandraces <- selectInd(pop, trait=fitFunc, nInd=n.purelines)
  # Purify each landrace by selfing it repeatedly
  for (f in 1:10) {
    purifiedLandraces <- self(purifiedLandraces)
  }
  # Create n.landrace x n.landrace F1s
  F1 <- randCross(purifiedLandraces, nCrosses=n.purelines)
  # Create enough F2s to select n.F2 to advance
  F2 <- self(F1, nProgeny=ceiling(n.F2/nInd(F1)))
  F2 <- selectInd(F2, trait=fitFunc, nInd=n.F2)
  F3 <- self(F2)
  F3 <- selectInd(F3, trait=fitFunc, nInd=n.F3)
  F4 <- self(F3)
  F4 <- selectInd(F4, trait=fitFunc, nInd=n.F4)
  F5 <- self(F4)
  F5 <- selectInd(F5, trait=fitFunc, nInd=n.F5)
  F6 <- self(F5)
  F6 <- selectInd(F6, trait=fitFunc, nInd=n.F6)
  F7 <- self(F6)
  F7 <- selectInd(F7, trait=fitFunc, nInd=n.F7)
  F8 <- self(F7)
  F8 <- selectInd(F8, trait=fitFunc, nInd=n.F8)
  
  pureline <- self(F8, nProgeny=20)
  return (pureline)
}

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
  F2 <- self(F1, nProgeny=250)
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

# Creates a diversity panel population by merging all of the subpopulations
# Returns: a population
# Assumes that pops - a list of populations, exists
createDP <- function() {
  DP <- c()
  # Iterate through all of the populations and merge them together
  for (p in 1:n.nPops) {
    DP <- append(DP, pops[[p]])
  }
  return (DP)
}

