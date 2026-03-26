# Title: GENOMIC PREDICTION
# Author: Ted Monyak
# Description: Contains functions for genomic prediction


# Train an RRBLUP model to predict one of the traits
# trainPop: the training population
# testPop: the test population
# trait: the AlphaSimR phenotype index
# Return correlation (r) between the EBVs and the actual genetic values in the test pop
evaluateGWP <- function(trainPop, testPop, trait) {
  # Update phenotype to have heritability associated with breeding programs
  trainPop <- setPheno(trainPop, h2=c(n.h2Breeding, n.h2Breeding, n.yieldH2Breeding))
  # Train the model
  # snpChip 2 is for genomic prediction
  model <- fastRRBLUP(trainPop, traits=trait, use="pheno", snpChip=2)
  # Set the estimated breeding values
  testPop <- setEBV(testPop, model)
  # Determine the correlation between genetic values and estimated breeding values
  # in the test population
  r <- cor(gv(testPop), ebv(testPop))[trait]
  return (r)
}

# Train an RRBLUP model to predict breeding fitness
# trainPop: the training population
# testPop: the test population
# trait: the AlphaSimR phenotype index
# Return correlation (r) between the EBVs and the actual genetic values in the test pop
evaluateGWP_W <- function(trainPop, testPop) {
  # Update phenotype to have heritabilities associated with breeding programs
  trainPop <- setPheno(trainPop, h2=c(n.h2Breeding, n.h2Breeding, n.yieldH2Breeding))
  # snpChip 2 is for genomic prediction
  model <- fastRRBLUP(trainPop, traits=calculateW_GWP, use="gv", snpChip=2)
  testPop <- setEBV(testPop, model)
  # Determine the correlation between genetic values and estimated breeding values
  # in the test population
  r <- cor(calculateW_GWP(gv(testPop)), ebv(testPop))[1]
  return (r)
}



# Run recurrent population improvement on the basePop
# The population will be improved according to phenotypic recurrent selection (PRS)
# and genomics-assisted recurrent selection (GARS)
# Store breeding fitness at each generation
# The GARS population has an extra selection each cycle, representative of a winter nursery
# Initial training population for RRBLUP is basePop. After that, the training population
# for cycle C is the top 20% of the population from the previous cycle based on EBV,
# plus the top 20% of the population from the previous cycle, based on phenotype
# Return a dataframe with n.C*2 rows (for both PRS and GARS)
recurrentSelection <- function(basePop) {
  # snpChip 2 is for genomic prediction
  model <- fastRRBLUP(basePop, traits=calculateW_GWP, use="gv", snpChip=2)
  
  # For storing results
  result <- data.frame(c=rep(seq(1:n.C), each=2),
                       sel=rep(c("GARS", "PRS"), times=n.C),
                       w=numeric(n.C*2))
  
  # Select the first cycle based on phenotype
  garsPop <- selectCross(basePop,
                         trait=breedingFitness,
                         nInd=nInd(basePop)*n.selInt,
                         nCrosses=nInd(basePop))
  # Replicate the population to conduct PRS
  prsPop <- garsPop
  
  # Iterate through each cycle
  for (c in 1:n.C) {
    # Calculate mean breeding fitness of GARS population
    result$w[(c-1)*2 + 1] <- as.data.frame(pheno(garsPop)) %>%
      dplyr::mutate(w=calculateBreedingFitness(Trait1, Trait2, Trait3)) %>%
      dplyr::summarize(meanW=mean(w)) %>%
      pull(meanW)
    
    # Mean breeding fitness of PRS population
    result$w[c*2] <- as.data.frame(pheno(prsPop)) %>%
      dplyr::mutate(w=calculateBreedingFitness(Trait1, Trait2, Trait3)) %>%
      dplyr::summarize(meanW=mean(w)) %>%
      pull(meanW)
  
    # Estimate BVs
    garsPop <- setEBV(garsPop, model)
    
    # Training population:
    # Top 20% of wGARS (based on EBV)
    # Top 20% of wGARS (based on pheno)
    topEBV <- selectInd(garsPop, nInd=0.2*nInd(garsPop), use="ebv")
    topPheno <- selectInd(garsPop, nInd=0.2*nInd(garsPop), trait=breedingFitness)
    trainPop <- c(topEBV, topPheno)
    
    # Winter nursery selection
    garsPop <- selectCross(garsPop, use="ebv", nInd=nInd(pop)*n.selInt, nCrosses=nInd(garsPop))
    
    # Update the model
    model <- fastRRBLUP(trainPop, traits=calculateW_GWP, use="gv", snpChip=2)
    
    # Estimate BVs
    garsPop <- setEBV(garsPop, model)
    
    # GS
    garsPop <- selectCross(garsPop, use="ebv", nInd=nInd(garsPop)*n.selInt, nCrosses=nInd(garsPop))
    
    # PS
    prsPop <- selectCross(prsPop, trait=breedingFitness, nInd=nInd(prsPop)*n.selInt, nCrosses=nInd(prsPop))
  }
  
  return (result)
}
