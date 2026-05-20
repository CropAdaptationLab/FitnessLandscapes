# Title: GENOMEWIDE PREDICTION
# Author: Ted Monyak
# Description: Contains functions for genomewide prediction


crossValAccuracy <- function(pop) {
  
  nSplits <- 5
  popSize <- pop@nInd
  
  results <- numeric(length=nSplits)
  
  for (c in seq(1:5)) {
    start <- (c-1) * popSize/nSplits + 1
    end <- c*popSize %/% nSplits
    testPop <- pop[start:end]
    trainPop <- pop[-(start:end)]
    results[c] <- evaluateGSModel(trainPop, testPop)
  }
  return(mean(results))
}

# Calculates breeding fitness and removes a dimension to allow for easy 
# computation of correlation with EBVs
calculateW_GWP <- function(x, suitFunc=suitabilityGaussian) {
  realizedYield <- breedingFitness(x, suitFunc)
  # Remove the first dimension to enable correlation of values
  dim(realizedYield) <- c(length(realizedYield), 1)
  return (realizedYield)
}

# Train an RRBLUP model to predict breeding fitness
# trainPop: the training population
# testPop: the test population
# trait: the AlphaSimR phenotype index
# Return correlation (r) between the EBVs and the actual genetic values in the test pop
evaluateGSModel <- function(trainPop, testPop) {
  # Update phenotype to have heritabilities associated with breeding programs
  trainPop <- setPheno(trainPop, h2=c(n.h2Breeding, n.h2Breeding, n.yieldH2Breeding))
  # snpChip 2 is for genomic prediction
  if (GS_MODEL == "RRBLUP") {
    model <- RRBLUP(trainPop, traits=calculateW_GWP, use=GS_PHENO, snpChip=2)
    testPop <- setEBV(testPop, model)
  } else if (GS_MODEL == "fastRRBLUP") {
    model <- fastRRBLUP(trainPop, traits=calculateW_GWP, use=GS_PHENO, snpChip=2)
    testPop <- setEBV(testPop, model)
  } else if (GS_MODEL == "GBLUP") {
    testPop <- estimateEBV_GBLUP(trainPop, testPop, subset=FALSE)
  }
  # Determine the correlation between genetic values and estimated breeding values
  # in the test population
  r <- cor(calculateW_GWP(gv(testPop)), ebv(testPop))[1]
  return (r)
}

estimateEBV_GBLUP <- function(trainPop, testPop, subset=FALSE) {
  # Update phenotypes
  trainPop <- setPheno(trainPop, h2=c(n.h2Breeding, n.h2Breeding, n.yieldH2Breeding))

  # Pull genotype matrices for both populations
  geno <- rbind(pullSnpGeno(trainPop, snpChip = 2),
                pullSnpGeno(testPop, snpChip = 2))
  
  # Get breeding fitness values for train pop
  w <- data.frame(pheno=trainPop@pheno) %>%
    dplyr::mutate(w=calculateBreedingFitness(pheno.Trait1, pheno.Trait2, pheno.Trait3)) %>%
    dplyr::pull(w)
  
  # Create a phenotype dataframe
  pheno <- data.frame(id=trainPop@id,
                      W=w)
  
  # Add the test pop ids with NAs as the phenotypes
  pheno <- dplyr::bind_rows(pheno,
                            data.frame(id=testPop@id))
  
  # Estimate GEBVS with GBLUP
  GEBV_W <- kin.blup(data=pheno,geno="id",pheno="W", K=A.mat(geno))$pred
  testPop@ebv <- matrix(tail(GEBV_W, testPop@nInd), nrow=testPop@nInd, ncol=1)
  return (testPop)
}

# Train an RR-BLUP model with the specified training population
# Returns: a model produced by fastRRBLUP
trainRRBLUPModel <- function(trainPop) {
  # Set phenotypes
  trainPop <- setPheno(trainPop, h2=c(n.h2Breeding, n.h2Breeding, n.yieldH2Breeding))
  
  # Retrain model
  if (GS_MODEL == "RRBLUP") {
    model <- RRBLUP(trainPop, traits=calculateW_GWP, use=GS_PHENO, snpChip=2)
  } else if (GS_MODEL == "fastRRBLUP") {
    model <- fastRRBLUP(trainPop, traits=calculateW_GWP, use=GS_PHENO, snpChip=2)
  }
  return (model)
}

# Create a training population with the top 20% of the population from the previous cycle
# based on EBV, plus the top 20% of the previous training population, based on phenotype
# Follows Muleta et al. 2019
# Returns: a training population
createTrainPop <- function(curPop, prevTrainPop, n=1000) {
  # Set phenotypes
  curPop <- selectInd(curPop, use="ebv", nInd=n/2)
  prevTrainPop <- setPheno(prevTrainPop, h2=c(n.h2Breeding, n.h2Breeding, n.yieldH2Breeding))
  prevTrainPop <- selectInd(prevTrainPop, trait=breedingFitness, nInd=n/2)
  return(c(curPop, prevTrainPop))
}

# Calculates metrics for a population at a particular cycle
# pop: An AlphaSimR population
# sel: The selection type (e.g. "GS", "PS", etc.)
# Returns: a dataframe tabulatireng all metrics
cycleMetrics <- function(pop, cycle, sel) {
  # Calulate accuracy of predictions
  r <- NA
  if (sel != "PS") {
    r <- cor(calculateW_GWP(gv(pop)), ebv(pop))[1]
  }
  pop <- setPheno(pop, h2=c(n.h2Breeding, n.h2Breeding, n.yieldH2Breeding))
  
  # Genome-wide heterozygosity
  genome_het <- meanHetLocus(pullSnpGeno(pop, snpChip=2))
  # Attained trait heterozygosity
  attained_het <- meanHetLocus(getUniqueQtl(pop))
  # Desired trait heterozygosity
  desired_het <- meanHetLocus(pullQtlGeno(pop, trait=3))
  # Get the population-level isoeliteness
  pop_ie <- popIsoeliteness(pop)
  
  # Calculate mean breeding fitness of GS population
  w <- as.data.frame(pheno(pop)) %>%
    dplyr::mutate(w=calculateBreedingFitness(Trait1, Trait2, Trait3)) %>%
    dplyr::summarize(meanW=mean(w)) %>%
    pull(meanW)
  
  # Calculate w based on genetic values
  wGV <- as.data.frame(gv(pop)) %>%
    dplyr::mutate(w=calculateBreedingFitness(Trait1, Trait2, Trait3)) %>%
    dplyr::summarize(meanW=mean(w)) %>%
    pull(meanW)
  
  return (data.frame(
    c=cycle,
    sel=sel,
    w=w,
    wGV=wGV,
    r=r,
    genome_het=genome_het,
    attained_het=attained_het,
    desired_het=desired_het,
    pop_ie=pop_ie,
    gvar=varG(pop)[3,3]))
}

# Run recurrent population improvement on the basePop
# Run epistatic linkage mapping to determine significant interactions,
# and run marker-assisted selection to recover parental haplotypes
# The population will be improved according to phenotypic recurrent selection (PS)
# and genomics-assisted recurrent selection (GS)
# Store breeding fitness at each generation
# The GS population has an extra selection each cycle, representative of a winter nursery
# Initial training population for RRBLUP is basePop. 
# Return a dataframe with n.C*2 rows (for both PS and GS)
recurrentSelection <- function(basePop, parent1, parent2) {
  # For storing results
  result <- data.frame(c=c(),
                       sel=c(),
                       w=c(),
                       wGV=c(),
                       r=c(),
                       genome_het=c(),
                       attained_het=c(),
                       desired_het=c(),
                       pop_ie=c(),
                       gvar=c())

  # Set phenotypes for base population
  basePop <- setPheno(basePop, h2=c(n.h2Breeding, n.h2Breeding, n.yieldH2Breeding))
  
  # Do marker-assisted selection

  # QTL-based haplotypes
  masHaplo_PERFECT <- data.frame(id=basePop@id)
  # Low resolution-based haplotypes
  masHaplo_LOW <- data.frame(id=basePop@id)
  # High resolution-based haplotypes
  masHaplo_HIGH <- data.frame(id=basePop@id)
  
  # 2D linkage mapping
  # Find all pairs of epistatic loci
  peaks <- epistaticLodPeaks(basePop, parent1, parent2, snpChip=2, trait=5)
  if (nrow(peaks) > 0) {
    masGeno_PERFECT <- getUniqueQtl(basePop)
    masGeno_HIGH <- as.data.frame(pullSnpGeno(basePop, snpChip=2))
    
    # Number of actual interactions
    nInt_PERFECT <- 0
    nInt_HIGH <- 0
    for (r in 1:nrow(peaks)) {
      # QTL interactions
      qtl1 <- peaks$qtl1[r]
      qtl2 <- peaks$qtl2[r]
      if (!is.na(qtl1) & !is.na(qtl2)) {
        if (segLocus(masGeno_PERFECT[,qtl1]) & segLocus(masGeno_PERFECT[,qtl2])) {
          haplos <- getHaplos(masGeno_PERFECT, qtl1, qtl2, parent1, parent2, snpChip=2, useQtls=TRUE)
          nInt_PERFECT <- nInt_PERFECT + 1
          masHaplo_PERFECT <- cbind(masHaplo_PERFECT, haplos)
        }
      }
      
      # High resolution
      m1 <- peaks$m1[r]
      m2 <- peaks$m2[r]
      if (!is.na(m1) & !is.na(m2)) {
        if (segLocus(masGeno_HIGH[,m1]) & segLocus(masGeno_HIGH[,m2])) {
          haplos <- getHaplos(masGeno_HIGH, m1, m2, parent1, parent2, snpChip=2, useQtls=FALSE)
          nInt_HIGH <- nInt_HIGH + 1
          masHaplo_HIGH <- cbind(masHaplo_HIGH, haplos)
        }
      }
    }
    # Number of QTL interactions
    if (nInt_PERFECT > 0) {
      colnames(masHaplo_PERFECT) <- c("id", paste0("INT_", 1:nInt_PERFECT))
    }
    
    # Number of high res marker interactions
    if (nInt_HIGH > 0) {
      colnames(masHaplo_HIGH) <- c("id", paste0("INT_", 1:nInt_HIGH))
    }
  }

  # Low resolution mapping
  low_res_peaks <- epistaticLodPeaks(basePop, parent1, parent2, snpChip=3, trait=5)
  if (nrow(low_res_peaks) > 0) {
    masGeno_LOW <- as.data.frame(pullSnpGeno(basePop, snpChip=3))
    nInt_LOW <- 0

    for (r in 1:nrow(low_res_peaks)) {
      # Low resolution
      m1 <- low_res_peaks$m1[r]
      m2 <- low_res_peaks$m2[r]
      if (!is.na(m1) & !is.na(m2)) {
        if (segLocus(masGeno_LOW[,m1]) & segLocus(masGeno_LOW[,m2])) {
          haplos <- getHaplos(masGeno_LOW, m1, m2, parent1, parent2, snpChip=3, useQtls=FALSE)
          nInt_LOW <- nInt_LOW + 1
          masHaplo_LOW <- cbind(masHaplo_LOW, haplos)
        }
      }
    }
    # Number of low res marker interactions
    if (nInt_LOW > 0) {
      colnames(masHaplo_LOW) <- c("id", paste0("INT_", 1:nInt_LOW))
    }
  }
  
  # Calculate sum of the favorable haplotypes and select individuals with
  # the maxmimum number
  masInds_PERFECT <- masHaplo_PERFECT %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      favHaplos = sum(c_across(starts_with("INT_")) == "P1")
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(favHaplos==max(favHaplos)) %>%
    dplyr::pull(id)
  
  masInds_HIGH <- masHaplo_HIGH %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      favHaplos = sum(c_across(starts_with("INT_")) == "P1")
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(favHaplos==max(favHaplos)) %>%
    dplyr::pull(id)

  masInds_LOW <- masHaplo_LOW %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      favHaplos = sum(c_across(starts_with("INT_")) == "P1")
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(favHaplos==max(favHaplos)) %>%
    dplyr::pull(id)
  
  # Initial training population
  gsTrainPop <- basePop
  # This will be the same training population throughout all cycles
  gsTrainPop_noUpdate <- gsTrainPop

  # Fit an RR-BLUP model to use initially for the 3 GS populations
  if (GS_MODEL %in% c("RRBLUP", "fastRRBLUP")) {
    # Train a model on the initial base population
    gsModel <- trainRRBLUPModel(gsTrainPop)
    # This model will not be updated throughout the recurrent selection
    gsModel_noUpdate <- gsModel
  }

  # Select the first cycle based on phenotype
  gsPop <- selectCross(basePop,
                       trait=breedingFitness,
                       nInd=nInd(basePop)*n.selInt,
                       nCrosses=nInd(basePop))
  # This population will use the initial GS model without updates
  gsPop_noUpdate <- gsPop
  
  # Use MAS to select individuals with the most favorable epistatic haplotypes
  # Train a GS model on this population
  masTrainPop_PERFECT <- basePop[masInds_PERFECT]
  if (GS_MODEL %in% c("RRBLUP", "fastRRBLUP")) {
    masModel_PERFECT <- trainRRBLUPModel(masTrainPop_PERFECT)
  }
  masPop_PERFECT <- selectCross(masTrainPop_PERFECT,
                                trait=breedingFitness,
                                nInd=min(length(masInds_PERFECT), nInd(basePop)*n.selInt),
                                nCrosses=nInd(basePop))

  # Store the initial training population, filtered by haplotypes
  masTrainPop_HIGH <- basePop[masInds_HIGH]
  # This will be the same training population for use throughout all cycles
  masTrainPop_HIGH_noUpdate <- masTrainPop_HIGH
  if (GS_MODEL %in% c("RRBLUP", "fastRRBLUP")) {
    masModel_HIGH <- trainRRBLUPModel(masTrainPop_HIGH)
    masModel_HIGH_noUpdate <- masModel_HIGH
  }
  
  masPop_HIGH <- selectCross(masTrainPop_HIGH,
                             trait=breedingFitness,
                             nInd=min(length(masInds_HIGH), nInd(basePop)*n.selInt),
                             nCrosses=nInd(basePop))
  # This population will use the same GS model without updates for all cycles
  masPop_HIGH_noUpdate <- masPop_HIGH

  masTrainPop_LOW <- basePop[masInds_LOW]
  if (GS_MODEL %in% c("RRBLUP", "fastRRBLUP")) {
    masModel_LOW <- trainRRBLUPModel(masTrainPop_LOW)
  }
  masPop_LOW <- selectCross(masTrainPop_LOW,
                            trait=breedingFitness,
                            nInd=min(length(masInds_LOW), nInd(basePop)*n.selInt),
                            nCrosses=nInd(basePop))
  
  # Phenotypic selection
  psPop <- gsPop
  # Conduct phenotypic selection on the population that has been filtered for MAS
  psIeMasPop <- masPop_HIGH
  
  # Iterate through each cycle
  for (cycle in 1:n.C) {
    
    if (GS_MODEL %in% c("RRBLUP", "fastRRBLUP")) {
      # If the model does not fit any values, there is no genetic variance
      # in the population
      if (any(is.na(gsModel@gv[[1]]@addEff)) |
          any(is.na(gsModel_noUpdate@gv[[1]]@addEff)) |
          any(is.na(masModel_PERFECT@gv[[1]]@addEff)) |
          any(is.na(masModel_LOW@gv[[1]]@addEff)) |
          any(is.na(masModel_HIGH@gv[[1]]@addEff)) |
          any(is.na(masModel_HIGH_noUpdate@gv[[1]]@addEff))) {
        return(result)
      }

      # NORMAL RRBLUP
      gsPop <- setEBV(gsPop, gsModel)
      gsPop_noUpdate <- setEBV(gsPop_noUpdate, gsModel_noUpdate)
      
      # PERFECT MAS
      masPop_PERFECT <- setEBV(masPop_PERFECT, masModel_PERFECT)
      
      # LOW RESOLUTION MAS
      masPop_LOW <- setEBV(masPop_LOW, masModel_LOW)
      
      # HIGH RESOLUTION MAS
      masPop_HIGH <- setEBV(masPop_HIGH, masModel_HIGH)
      masPop_HIGH_noUpdate <- setEBV(masPop_HIGH_noUpdate, masModel_HIGH_noUpdate)
    } else if (GS_MODEL == "GBLUP") {
      gsPop <- estimateEBV_GBLUP(gsTrainPop, gsPop)
      gsPop_noUpdate <- estimateEBV_GBLUP(gsTrainPop_noUpdate, gsPop_noUpdate)
      masPop_PERFECT <- estimateEBV_GBLUP(masTrainPop_PERFECT, masPop_PERFECT)
      masPop_LOW <- estimateEBV_GBLUP(masTrainPop_LOW, masPop_LOW)
      masPop_HIGH <- estimateEBV_GBLUP(masTrainPop_HIGH, masPop_HIGH)
      masPop_HIGH_noUpdate <- estimateEBV_GBLUP(masTrainPop_HIGH_noUpdate, masPop_HIGH_noUpdate)
    }
    result <- rbind(result, cycleMetrics(gsPop, cycle, "GS"))
    result <- rbind(result, cycleMetrics(gsPop_noUpdate, cycle, "GSnoUpdate"))
    result <- rbind(result, cycleMetrics(masPop_PERFECT, cycle, "ieMAS"))
    result <- rbind(result, cycleMetrics(masPop_LOW, cycle, "lowResMAS"))
    result <- rbind(result, cycleMetrics(masPop_HIGH, cycle, "highResMAS"))
    result <- rbind(result, cycleMetrics(masPop_HIGH_noUpdate, cycle, "highResMASnoUpdate"))

    # PHENOTYPIC SELECTION
    
    result <- rbind(result, cycleMetrics(psPop, cycle, "PS"))
    result <- rbind(result, cycleMetrics(psIeMasPop, cycle, "PSieMAS"))
    
    # Update the models in even cycles
    if (cycle %% 1 == 0) {
      # Retrain models
      gsTrainPop <- createTrainPop(gsPop, gsTrainPop)
      masTrainPop_PERFECT <- createTrainPop(masPop_PERFECT, masTrainPop_PERFECT)
      masTrainPop_LOW <- createTrainPop(masPop_LOW, masTrainPop_LOW)
      masTrainPop_HIGH <- createTrainPop(masPop_HIGH, masTrainPop_HIGH)
      if (GS_MODEL %in% c("RRBLUP", "fastRRBLUP")) {
        # TrainPop is top 200 w from prev train pop + top gebv from this pop 
        gsModel <- trainRRBLUPModel(gsTrainPop)
        masModel_PERFECT <- trainRRBLUPModel(masTrainPop_PERFECT)
        masModel_LOW <- trainRRBLUPModel(masTrainPop_LOW)
        masModel_HIGH <- trainRRBLUPModel(masTrainPop_HIGH)
      }
    }
    gsPop <- selectCross(gsPop, use="ebv", nInd=nInd(gsPop)*n.selInt, nCrosses=nInd(gsPop))
    gsPop_noUpdate <- selectCross(gsPop_noUpdate, use="ebv", nInd=nInd(gsPop_noUpdate)*n.selInt, nCrosses=nInd(gsPop_noUpdate))
    masPop_PERFECT <- selectCross(masPop_PERFECT, use="ebv", nInd=nInd(masPop_PERFECT)*n.selInt, nCrosses=nInd(masPop_PERFECT))
    masPop_LOW <- selectCross(masPop_LOW, use="ebv", nInd=nInd(masPop_LOW)*n.selInt, nCrosses=nInd(masPop_LOW))
    masPop_HIGH <- selectCross(masPop_HIGH, use="ebv", nInd=nInd(masPop_HIGH)*n.selInt, nCrosses=nInd(masPop_HIGH))
    masPop_HIGH_noUpdate <- selectCross(masPop_HIGH_noUpdate, use="ebv", nInd=nInd(masPop_HIGH_noUpdate)*n.selInt, nCrosses=nInd(masPop_HIGH_noUpdate))
    
    psPop <- setPheno(psPop, h2=c(n.h2Breeding, n.h2Breeding, n.yieldH2Breeding))
    psPop <- selectCross(psPop, trait=breedingFitness, nInd=nInd(psPop)*n.selInt, nCrosses=nInd(psPop))

    psIeMasPop <- setPheno(psIeMasPop, h2=c(n.h2Breeding, n.h2Breeding, n.yieldH2Breeding))
    psIeMasPop <- selectCross(psIeMasPop, trait=breedingFitness, nInd=nInd(psIeMasPop)*n.selInt, nCrosses=nInd(psIeMasPop))
  }
  return (result)
}
