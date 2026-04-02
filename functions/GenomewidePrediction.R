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
# The population will be improved according to phenotypic recurrent selection (PS)
# and genomics-assisted recurrent selection (GS)
# Store breeding fitness at each generation
# The GS population has an extra selection each cycle, representative of a winter nursery
# Initial training population for RRBLUP is basePop. After that, the training population
# for cycle C is the top 20% of the population from the previous cycle based on EBV,
# plus the top 20% of the population from the previous cycle, based on phenotype
# Return a dataframe with n.C*2 rows (for both PS and GS)
recurrentSelection <- function(basePop) {
  # Set phenotypes for base population
  basePop <- setPheno(basePop, h2=c(n.h2Breeding, n.h2Breeding, n.yieldH2Breeding))
  
  # snpChip 2 is for genomic prediction
  model <- fastRRBLUP(basePop, traits=calculateW_GWP, use="gv", snpChip=2)
  
  # For storing results
  result <- data.frame(c=c(),
                       sel=c(),
                       w=c(),
                       r=c(),
                       genome_het=c(),
                       attained_het=c(),
                       desired_het=c(),
                       pop_ie=c())
  
  # If the model does not fit any values, there is no genetic variance
  # in the population
  if (any(is.na(model@gv[[1]]@addEff))) {
    return(result)
  }
  
  # Select the first cycle based on phenotype
  gsPop <- selectCross(basePop,
                         trait=breedingFitness,
                         nInd=nInd(basePop)*n.selInt,
                         nCrosses=nInd(basePop))
  # Replicate the population to conduct PS
  psPop <- gsPop
  
  # Iterate through each cycle
  for (c in 1:n.C) {
    
    # Estimate BVs
    gsPop <- setEBV(gsPop, model)
    
    # Calulate accuracy of predictions
    r <- cor(calculateW_GWP(gv(gsPop)), ebv(gsPop))[1]
    
    # Genome-wide heterozygosity
    genome_het <- meanHetLocus(pullSnpGeno(gsPop, snpChip=2))
    # Attained trait heterozygosity
    attained_het <- meanHetLocus(getUniqueQtl(gsPop))
    # Desired trait heterozygosity
    desired_het <- meanHetLocus(pullQtlGeno(gsPop, trait=3))
    # Get the population-level isoeliteness
    pop_ie <- popIsoeliteness(gsPop)
    
    # Calculate mean breeding fitness of GS population
    wGS <- as.data.frame(pheno(gsPop)) %>%
      dplyr::mutate(w=calculateBreedingFitness(Trait1, Trait2, Trait3)) %>%
      dplyr::summarize(meanW=mean(w)) %>%
      pull(meanW)

    result <- rbind(result,
                    data.frame(
                      c=c,
                      sel="GS",
                      w=wGS,
                      r=r,
                      genome_het=genome_het,
                      attained_het=attained_het,
                      desired_het=desired_het,
                      pop_ie=pop_ie,
                      gvar=varG(gsPop)[3,3]))
                    
    
    # Mean breeding fitness of PS population
    wPS <- as.data.frame(pheno(psPop)) %>%
      dplyr::mutate(w=calculateBreedingFitness(Trait1, Trait2, Trait3)) %>%
      dplyr::summarize(meanW=mean(w)) %>%
      pull(meanW)
    
    # Genome-wide heterozygosity
    genome_het <- meanHetLocus(pullSnpGeno(psPop, snpChip=2))
    # Attained trait heterozygosity
    attained_het <- meanHetLocus(getUniqueQtl(psPop))
    # Desired trait heterozygosity
    desired_het <- meanHetLocus(pullQtlGeno(psPop, trait=3))
    # Get the population-level isoeliteness
    pop_ie <- popIsoeliteness(psPop)
    
    result <- rbind(result,
                    data.frame(
                      c=c,
                      sel="PS",
                      w=wPS,
                      r=NA,
                      genome_het=genome_het,
                      attained_het=attained_het,
                      desired_het=desired_het,
                      pop_ie=pop_ie,
                      gvar=varG(psPop)[3,3]))
    
    # Update the model in even cycles
    if (c %% 2 == 0) {
      # Training population:
      # Top 20% of wGS (based on EBV)
      topEBV <- selectInd(gsPop, nInd=0.2*nInd(gsPop), use="ebv")
      # Top 20% of wGS (based on pheno)
      topPheno <- selectInd(gsPop, nInd=0.2*nInd(gsPop), trait=breedingFitness)
      trainPop <- c(topEBV, topPheno)
    
      # Retrain model
      model <- fastRRBLUP(trainPop, traits=calculateW_GWP, use="gv", snpChip=2)

      # If the model does not fit any values, there is no genetic variance
      # in the population
      if (any(is.na(model@gv[[1]]@addEff))) {
        return(result)
      }
    }

    gsPop <- selectCross(gsPop, use="ebv", nInd=nInd(gsPop)*n.selInt, nCrosses=nInd(gsPop))
    psPop <- selectCross(psPop, trait=breedingFitness, nInd=nInd(psPop)*n.selInt, nCrosses=nInd(psPop))
  }
  return (result)
}
