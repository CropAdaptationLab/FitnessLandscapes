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
createTrainPop <- function(curPop, prevTrainPop, n=400) {
  # Set phenotypes
  curPop <- selectInd(curPop, use="ebv", nInd=n/2)
  prevTrainPop <- setPheno(prevTrainPop, h2=c(n.h2Breeding, n.h2Breeding, n.yieldH2Breeding))
  prevTrainPop <- selectInd(prevTrainPop, trait=breedingFitness, nInd=n/2)
  return(c(curPop, prevTrainPop))
}

# Calculates metrics for a population at a particular cycle
# pop: An AlphaSimR population
# sel: The selection type (e.g. "GS", "PS", etc.)
# Returns: a dataframe tabulating all metrics
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
  # Use the genetic values as phenotype, because these are supposed to be
  # replicated inbred lines
  basePop@pheno <- basePop@gv
  C0 <- selectInd(basePop, trait=breedingFitness, nInd=nInd(basePop)*n.selInt)
  
  result <- rbind(result, cycleMetrics(C0, 0, "GS"))
  result <- rbind(result, cycleMetrics(C0, 0, "GSnoUpdate"))
  result <- rbind(result, cycleMetrics(C0, 0, "PS"))
  
  gsTrainPop <- C0
  # This will be the same training population throughout all cycles
  gsTrainPop_noUpdate <- C0

  # Fit an RR-BLUP model to use initially for the 3 GS populations
  if (GS_MODEL %in% c("RRBLUP", "fastRRBLUP")) {
    # Train a model on the initial base population
    gsModel <- trainRRBLUPModel(gsTrainPop)
    # This model will not be updated throughout the recurrent selection
    gsModel_noUpdate <- gsModel
  }

  # Select the first cycle based on phenotype
  gs_S0 <- randCross(C0,
                     nCrosses=n.recurPopSize)
  # This population will use the initial GS model without updates
  gs_S0_noUpdate <- gs_S0
  # Phenotypic selection
  ps_S0 <- gs_S0
  
  # Use MAS to select individuals with the most favorable epistatic haplotypes
  # Train a GS model on this population
  gsTrainPop_PERFECT <- basePop[masInds_PERFECT]
  result <- rbind(result, cycleMetrics(gsTrainPop_PERFECT, 0, "ieMAS"))

  if (GS_MODEL %in% c("RRBLUP", "fastRRBLUP")) {
    gsModel_MAS_PERFECT <- trainRRBLUPModel(gsTrainPop_PERFECT)
  }
  mas_S0_PERFECT <- selectCross(gsTrainPop_PERFECT,
                                trait=breedingFitness,
                                nInd=min(length(masInds_PERFECT),
                                         nInd(basePop)*n.selInt),
                                nCrosses=n.recurPopSize)

  # Store the initial training population, filtered by haplotypes
  gsTrainPop_HIGH <- basePop[masInds_HIGH]
  # This will be the same training population for use throughout all cycles
  gsTrainPop_HIGH_noUpdate <- gsTrainPop_HIGH
  
  result <- rbind(result, cycleMetrics(gsTrainPop_HIGH, 0, "highResMAS"))
  result <- rbind(result, cycleMetrics(gsTrainPop_HIGH_noUpdate, 0, "highResMASnoUpdate"))
  result <- rbind(result, cycleMetrics(gsTrainPop_HIGH, 0, "PSieMAS"))
  
  if (GS_MODEL %in% c("RRBLUP", "fastRRBLUP")) {
    gsModel_MAS_HIGH <- trainRRBLUPModel(gsTrainPop_HIGH)
    gsModel_MAS_HIGH_noUpdate <- gsModel_MAS_HIGH
  }
  mas_S0_HIGH <- selectCross(gsTrainPop_HIGH,
                             trait=breedingFitness,
                             nInd=min(length(masInds_HIGH),
                                      nInd(basePop)*n.selInt),
                             nCrosses=n.recurPopSize)
  # This population will use the same GS model without updates for all cycles
  mas_S0_HIGH_noUpdate <- mas_S0_HIGH
  # Conduct phenotypic selection on the population that has been filtered for MAS
  ps_mas_S0 <- mas_S0_HIGH

  gsTrainPop_LOW <- basePop[masInds_LOW]
  
  result <- rbind(result, cycleMetrics(gsTrainPop_LOW, 0, "lowResMAS"))
  if (GS_MODEL %in% c("RRBLUP", "fastRRBLUP")) {
    gsModel_MAS_LOW <- trainRRBLUPModel(gsTrainPop_LOW)
  }
  mas_S0_LOW <- selectCross(gsTrainPop_LOW,
                            trait=breedingFitness,
                            nInd=min(length(masInds_LOW),
                                     nInd(basePop)*n.selInt),
                            nCrosses=n.recurPopSize)
  
  # Iterate through each cycle
  cycle <- 0
  for (year in 1:n.Y) {
    
    #-------SUMMER NURSERY-------#
    
    if (GS_MODEL %in% c("RRBLUP", "fastRRBLUP")) {
      # If the new models do not fit any values, there is no genetic variance
      # in the population
      if (any(is.na(gsModel@gv[[1]]@addEff)) |
          any(is.na(gsModel_noUpdate@gv[[1]]@addEff)) |
          any(is.na(gsModel_MAS_PERFECT@gv[[1]]@addEff)) |
          any(is.na(gsModel_MAS_LOW@gv[[1]]@addEff)) |
          any(is.na(gsModel_MAS_HIGH@gv[[1]]@addEff)) |
          any(is.na(gsModel_MAS_HIGH_noUpdate@gv[[1]]@addEff))) {
        return(result)
      }
      gs_S0 <- setEBV(gs_S0, gsModel)
      gs_S0_noUpdate <- setEBV(gs_S0_noUpdate, gsModel_noUpdate)
      mas_S0_PERFECT <- setEBV(mas_S0_PERFECT, gsModel_MAS_PERFECT)
      mas_S0_LOW <- setEBV(mas_S0_LOW, gsModel_MAS_LOW)
      mas_S0_HIGH <- setEBV(mas_S0_HIGH, gsModel_MAS_HIGH)
      mas_S0_HIGH_noUpdate <- setEBV(mas_S0_HIGH_noUpdate, gsModel_MAS_HIGH_noUpdate)
    } else if (GS_MODEL == "GBLUP") {
      gs_S0 <- estimateEBV_GBLUP(gsTrainPop, gs_S0)
      gs_S0_noUpdate <- estimateEBV_GBLUP(gsTrainPop_noUpdate, gs_S0_noUpdate)
      mas_S0_PERFECT <- estimateEBV_GBLUP(gsTrainPop_PERFECT, mas_S0_PERFECT)
      mas_S0_LOW <- estimateEBV_GBLUP(gsTrainPop_LOW, mas_S0_LOW)
      mas_S0_HIGH <- estimateEBV_GBLUP(gsTrainPop_HIGH, mas_S0_HIGH)
      mas_S0_HIGH_noUpdate <- estimateEBV_GBLUP(gsTrainPop_HIGH_noUpdate, mas_S0_HIGH_noUpdate)
    }

    cycle <- cycle + 0.5
    result <- rbind(result, cycleMetrics(gs_S0, cycle, "GS"))
    result <- rbind(result, cycleMetrics(gs_S0_noUpdate, cycle, "GSnoUpdate"))
    result <- rbind(result, cycleMetrics(mas_S0_PERFECT, cycle, "ieMAS"))
    result <- rbind(result, cycleMetrics(mas_S0_LOW, cycle, "lowResMAS"))
    result <- rbind(result, cycleMetrics(mas_S0_HIGH, cycle, "highResMAS"))
    result <- rbind(result, cycleMetrics(mas_S0_HIGH_noUpdate, cycle, "highResMASnoUpdate"))
    result <- rbind(result, cycleMetrics(ps_S0, cycle, "PS"))
    result <- rbind(result, cycleMetrics(ps_mas_S0, cycle, "PSieMAS"))

    # Establish new training population by selecting the top lines from the populations by gebv
    # Do this before selections and crossing because training happens after harvest
    gsTrainPop <- createTrainPop(gs_S0, gsTrainPop)
    gsTrainPop_PERFECT <- createTrainPop(mas_S0_PERFECT, gsTrainPop_PERFECT)
    gsTrainPop_LOW <- createTrainPop(mas_S0_LOW, gsTrainPop_LOW)
    gsTrainPop_HIGH <- createTrainPop(mas_S0_HIGH, gsTrainPop_HIGH)
    
    gs_S0 <- selectCross(gs_S0, use="ebv", nInd=nInd(gs_S0)*n.selInt, nCrosses=n.recurPopSize)
    gs_S0_noUpdate <- selectCross(gs_S0_noUpdate, use="ebv", nInd=nInd(gs_S0_noUpdate)*n.selInt, nCrosses=n.recurPopSize)
    mas_S0_PERFECT <- selectCross(mas_S0_PERFECT, use="ebv", nInd=nInd(mas_S0_PERFECT)*n.selInt, nCrosses=n.recurPopSize)
    mas_S0_LOW <- selectCross(mas_S0_LOW, use="ebv", nInd=nInd(mas_S0_LOW)*n.selInt, nCrosses=n.recurPopSize)
    mas_S0_HIGH <- selectCross(mas_S0_HIGH, use="ebv", nInd=nInd(mas_S0_HIGH)*n.selInt, nCrosses=n.recurPopSize)
    mas_S0_HIGH_noUpdate <- selectCross(mas_S0_HIGH_noUpdate, use="ebv", nInd=nInd(mas_S0_HIGH_noUpdate)*n.selInt, nCrosses=n.recurPopSize)

    # Phenotype the S0s
    ps_S0 <- setPheno(ps_S0, h2=c(n.h2Breeding, n.h2Breeding, n.yieldH2Breeding))
    # Select the top individuals by phenotype
    top_ps_S0 <- selectInd(ps_S0, trait=breedingFitness, nInd=nInd(ps_S0)*n.selInt)
    # In the simulation, only self the top lines, so they can be easily selected in
    # the next generation. In reality, all S1s would be planted
    ps_S1 <- self(top_ps_S0)
    
    ps_mas_S0 <- setPheno(ps_mas_S0, h2=c(n.h2Breeding, n.h2Breeding, n.yieldH2Breeding))
    top_ps_mas_S0 <- selectInd(ps_mas_S0, trait=breedingFitness, nInd=nInd(ps_mas_S0)*n.selInt)
    ps_mas_S1 <- self(top_ps_mas_S0)
    
    # RETRAIN MODELS
    if (GS_MODEL %in% c("RRBLUP", "fastRRBLUP")) {
      gsModel <- trainRRBLUPModel(gsTrainPop)
      gsModel_MAS_PERFECT <- trainRRBLUPModel(gsTrainPop_PERFECT)
      gsModel_MAS_LOW <- trainRRBLUPModel(gsTrainPop_LOW)
      gsModel_MAS_HIGH <- trainRRBLUPModel(gsTrainPop_HIGH)
    }
 
    # SET GEBVS
    if (GS_MODEL %in% c("RRBLUP", "fastRRBLUP")) {
      gs_S0 <- setEBV(gs_S0, gsModel)
      gs_S0_noUpdate <- setEBV(gs_S0_noUpdate, gsModel_noUpdate)
      mas_S0_PERFECT <- setEBV(mas_S0_PERFECT, gsModel_MAS_PERFECT)
      mas_S0_LOW <- setEBV(mas_S0_LOW, gsModel_MAS_LOW)
      mas_S0_HIGH <- setEBV(mas_S0_HIGH, gsModel_MAS_HIGH)
      mas_S0_HIGH_noUpdate <- setEBV(mas_S0_HIGH_noUpdate, gsModel_MAS_HIGH_noUpdate)
    } else if (GS_MODEL == "GBLUP") {
      gs_S0 <- estimateEBV_GBLUP(gsTrainPop, gs_S0)
      gs_S0_noUpdate <- estimateEBV_GBLUP(gsTrainPop_noUpdate, gs_S0_noUpdate)
      mas_S0_PERFECT <- estimateEBV_GBLUP(gsTrainPop_PERFECT, mas_S0_PERFECT)
      mas_S0_LOW <- estimateEBV_GBLUP(gsTrainPop_LOW, mas_S0_LOW)
      mas_S0_HIGH <- estimateEBV_GBLUP(gsTrainPop_HIGH, mas_S0_HIGH)
      mas_S0_HIGH_noUpdate <- estimateEBV_GBLUP(gsTrainPop_HIGH_noUpdate, mas_S0_HIGH_noUpdate)
    }
    
    #-------WINTER NURSERY-------#
    cycle <- cycle + 0.5
    result <- rbind(result, cycleMetrics(gs_S0, cycle, "GS"))
    result <- rbind(result, cycleMetrics(gs_S0_noUpdate, cycle, "GSnoUpdate"))
    result <- rbind(result, cycleMetrics(mas_S0_PERFECT, cycle, "ieMAS"))
    result <- rbind(result, cycleMetrics(mas_S0_LOW, cycle, "lowResMAS"))
    result <- rbind(result, cycleMetrics(mas_S0_HIGH, cycle, "highResMAS"))
    result <- rbind(result, cycleMetrics(mas_S0_HIGH_noUpdate, cycle, "highResMASnoUpdate"))
    result <- rbind(result, cycleMetrics(ps_S1, cycle, "PS"))
    result <- rbind(result, cycleMetrics(ps_mas_S1, cycle, "PSieMAS"))
    
    # ADVANCE LINES
    gs_S0 <- selectCross(gs_S0, use="ebv", nInd=nInd(gs_S0)*n.selInt, nCrosses=n.recurPopSize)
    gs_S0_noUpdate <- selectCross(gs_S0_noUpdate, use="ebv", nInd=nInd(gs_S0_noUpdate)*n.selInt, nCrosses=n.recurPopSize)
    mas_S0_PERFECT <- selectCross(mas_S0_PERFECT, use="ebv", nInd=nInd(mas_S0_PERFECT)*n.selInt, nCrosses=n.recurPopSize)
    mas_S0_LOW <- selectCross(mas_S0_LOW, use="ebv", nInd=nInd(mas_S0_LOW)*n.selInt, nCrosses=n.recurPopSize)
    mas_S0_HIGH <- selectCross(mas_S0_HIGH, use="ebv", nInd=nInd(mas_S0_HIGH)*n.selInt, nCrosses=n.recurPopSize)
    mas_S0_HIGH_noUpdate <- selectCross(mas_S0_HIGH_noUpdate, use="ebv", nInd=nInd(mas_S0_HIGH_noUpdate)*n.selInt, nCrosses=n.recurPopSize)
    ps_S0 <- randCross(ps_S1, nCrosses=n.recurPopSize)
    ps_mas_S0 <- randCross(ps_mas_S1, nCrosses=n.recurPopSize)
  }
  return (result)
}
