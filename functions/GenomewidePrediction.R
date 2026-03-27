# Title: GENOMEWIDE PREDICTION
# Author: Ted Monyak
# Description: Contains functions for genomewide prediction


# Calculates breeding fitness and removes a dimension to allow for easy 
# computation of correlation with EBVs
calculateW_GWP <- function(x, suitFunc=suitabilityGaussian) {
  realizedYield <- breedingFitness(x, suitFunc)
  # Remove the first dimension to enable correlation of values
  dim(realizedYield) <- c(length(realizedYield), 1)
  return (realizedYield)
}

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
  model <- fastRRBLUP(trainPop, traits=trait, use="gv", snpChip=2)
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
  # Set phenotypes for base population
  basePop <- setPheno(basePop, h2=c(n.h2Breeding, n.h2Breeding, n.yieldH2Breeding))
  
  # snpChip 2 is for genomic prediction
  model <- fastRRBLUP(basePop, traits=calculateW_GWP, use="gv", snpChip=2)
  
  # For storing results
  result <- data.frame(c=c(),
                       sel=c(),
                       w=c())
  
  # If the model does not fit any values, there is no genetic variance
  # in the population
  if (any(is.na(model@gv[[1]]@addEff))) {
    return(result)
  }
  
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
    meanWGARS <- as.data.frame(pheno(garsPop)) %>%
      dplyr::mutate(w=calculateBreedingFitness(Trait1, Trait2, Trait3)) %>%
      dplyr::summarize(meanW=mean(w)) %>%
      pull(meanW)
    result <- rbind(result,
                    data.frame(
                      c=c,
                      sel="GARS",
                      w=meanWGARS))
                    
    
    # Mean breeding fitness of PRS population
    meanWPRS <- as.data.frame(pheno(prsPop)) %>%
      dplyr::mutate(w=calculateBreedingFitness(Trait1, Trait2, Trait3)) %>%
      dplyr::summarize(meanW=mean(w)) %>%
      pull(meanW)
    
    result <- rbind(result,
                    data.frame(
                      c=c,
                      sel="PRS",
                      w=meanWPRS))

    # Estimate BVs
    garsPop <- setEBV(garsPop, model)
    
    # Update the model in even cycles
    if (c %% 2 == 0) {
      # Training population:
      # Top 20% of wGARS (based on EBV)
      topEBV <- selectInd(garsPop, nInd=0.2*nInd(garsPop), use="ebv")
      # Top 20% of wGARS (based on pheno)
      topPheno <- selectInd(garsPop, nInd=0.2*nInd(garsPop), trait=breedingFitness)
      trainPop <- c(topEBV, topPheno)
    
      # Retrain model
      model <- fastRRBLUP(trainPop, traits=calculateW_GWP, use="gv", snpChip=2)

      # If the model does not fit any values, there is no genetic variance
      # in the population
      if (any(is.na(model@gv[[1]]@addEff))) {
        return(result)
      }
    }

    garsPop <- selectCross(garsPop, use="ebv", nInd=nInd(garsPop)*n.selInt, nCrosses=nInd(garsPop))
    prsPop <- selectCross(prsPop, trait=breedingFitness, nInd=nInd(prsPop)*n.selInt, nCrosses=nInd(prsPop))
  }
  return (result)
}
