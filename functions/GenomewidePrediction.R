# Title: GENOMEWIDE PREDICTION
# Author: Ted Monyak
# Description: Contains functions for genomewide prediction

library(rrBLUP)

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
evaluateRRBLUP <- function(trainPop, testPop, trait) {
  # Update phenotype to have heritability associated with breeding programs
  trainPop <- setPheno(trainPop, h2=c(n.h2Breeding, n.h2Breeding, n.yieldH2Breeding))
  # Train the model
  # snpChip 2 is for genomic prediction
  model <- fastRRBLUP(trainPop, traits=trait, use=gsPheno, snpChip=2)
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
evaluateRRBLUP_W <- function(trainPop, testPop) {
  # Update phenotype to have heritabilities associated with breeding programs
  trainPop <- setPheno(trainPop, h2=c(n.h2Breeding, n.h2Breeding, n.yieldH2Breeding))
  # snpChip 2 is for genomic prediction
  model <- fastRRBLUP(trainPop, traits=calculateW_GWP, use=gsPheno, snpChip=2)
  testPop <- setEBV(testPop, model)
  # Determine the correlation between genetic values and estimated breeding values
  # in the test population
  r <- cor(calculateW_GWP(gv(testPop)), ebv(testPop))[1]
  return (r)
}

# Train a G-BLUP model from the rrBLUP package
# trainPop: the training population
# testPop: the test population
# Return correlation (r) between the EBVs and the actual genetic values in the test pop
evaluateGBLUP_W <- function(trainPop, testPop) {
  # Pull genotype matrices for both populations
  geno <- rbind(pullSnpGeno(trainPop, snpChip = 2),
                pullSnpGeno(testPop, snpChip = 2))
  
  # Update phenotypes
  trainPop <- setPheno(trainPop, h2=c(n.h2Breeding, n.h2Breeding, n.yieldH2Breeding))
  # Set W pheno
  
  # Get breeding fitness values for train pop
  w <- data.frame(gv=trainPop@gv) %>%
    dplyr::mutate(w=calculateBreedingFitness(gv.Trait1, gv.Trait2, gv.Trait3)) %>%
    dplyr::pull(w)
  
  # Create a phenotype dataframe
  pheno <- data.frame(id=trainPop@id,
                      W=w)
  
  # Add the test pop ids with NAs as the phenotypes
  pheno <- dplyr::bind_rows(pheno,
                            data.frame(id=testPop@id))
  
  # Estimate GEBVS with GBLUP
  GEBV_W <- kin.blup(data=pheno,geno="id",pheno="W", K=A.mat(geno))$pred

  r <- cor(calculateW_GWP(gv(testPop)), tail(GEBV_W, nInd(testPop)))[1]
  return (r)
}

# Train an RR-BLUP model with the top 20% of the population from the previous cycle
# based on EBV, plus the top 20% of the population from the previous cycle, based on phenotype
# Follows Muleta et al. 2019
# Returns: a model produced by fastRRBLUP
retrainModel <- function(pop) {
  # Top 20% based on EBV
  topEBV <- selectInd(pop, nInd=0.2*nInd(pop), use="ebv")
  # Top 20% based on pheno
  topPheno <- selectInd(pop, nInd=0.2*nInd(pop), trait=breedingFitness)
  trainPop <- c(topEBV, topPheno)
  
  # Retrain model
  newModel <- fastRRBLUP(trainPop, traits=calculateW_GWP, use=gsPheno, snpChip=2)
  return (newModel)
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
  
  return (data.frame(
    c=cycle,
    sel=sel,
    w=w,
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
  
  
  # Fit an RR-BLUP model to use initially for the 3 GS populations
  gsModel <- fastRRBLUP(basePop, traits=calculateW_GWP, use=gsPheno, snpChip=2)
  masModel_PERFECT <- gsModel
  masModel_HIGH <- gsModel
  masModel_LOW <- gsModel
  
  # For storing results
  result <- data.frame(c=c(),
                       sel=c(),
                       w=c(),
                       r=c(),
                       genome_het=c(),
                       attained_het=c(),
                       desired_het=c(),
                       pop_ie=c(),
                       gvar=c())
  
  # If the model does not fit any values, there is no genetic variance
  # in the population
  if (any(is.na(gsModel@gv[[1]]@addEff))) {
    return(result)
  }
  
  # Select the first cycle based on phenotype
  gsPop <- selectCross(basePop,
                         trait=breedingFitness,
                         nInd=nInd(basePop)*n.selInt,
                         nCrosses=nInd(basePop))
  # Replicate the population to conduct PS
  psPop <- gsPop
  
  # Use MAS to select individuals with the most favorable epistatic haplotypes
  masPop_PERFECT <- selectInd(basePop,
                      trait=breedingFitness,
                      nInd=min(nInd(basePop)*n.selInt, length(masInds_PERFECT)),
                      candidates=masInds_PERFECT)
  masPop_PERFECT <- randCross(masPop_PERFECT, nCrosses=nInd(basePop))
  
  masPop_HIGH <- selectInd(basePop,
                        trait=breedingFitness,
                        nInd=min(nInd(basePop)*n.selInt, length(masInds_HIGH)),
                        candidates=masInds_HIGH)
  masPop_HIGH <- randCross(masPop_HIGH, nCrosses=nInd(basePop))
  
  masPop_LOW <- selectInd(basePop,
                            trait=breedingFitness,
                            nInd=min(nInd(basePop)*n.selInt, length(masInds_LOW)),
                            candidates=masInds_LOW)
  masPop_LOW <- randCross(masPop_LOW, nCrosses=nInd(basePop))
  
  # Iterate through each cycle
  for (cycle in 1:n.C) {
    
    # GENOMIC SELECTION MAS
    gsPop <- setEBV(gsPop, gsModel)
    result <- rbind(result, cycleMetrics(gsPop, cycle, "GS"))
    
    # PHENOTYPIC SELECTION
    result <- rbind(result, cycleMetrics(psPop, cycle, "PS"))
    
    # PERFECT MAS
    masPop_PERFECT <- setEBV(masPop_PERFECT, masModel_PERFECT)
    result <- rbind(result, cycleMetrics(masPop_PERFECT, cycle, "ieMAS"))
    
    # LOW RESOLUTION MAS
    masPop_LOW <- setEBV(masPop_LOW, masModel_LOW)
    result <- rbind(result, cycleMetrics(masPop_LOW, cycle, "lowResMAS"))
    
    # HIGH RESOLUTION MAS
    masPop_HIGH <- setEBV(masPop_HIGH, masModel_HIGH)
    result <- rbind(result, cycleMetrics(masPop_HIGH, cycle, "highResMAS"))
    
    # Update the models in even cycles
    if (cycle %% 2 == 0) {
      # Retrain models
      gsModel <- retrainModel(gsPop)
      masModel_PERFECT <- retrainModel(masPop_PERFECT)
      masModel_LOW <- retrainModel(masPop_LOW)
      masModel_HIGH <- retrainModel(masPop_HIGH)
      
      # If the model does not fit any values, there is no genetic variance
      # in the population
      if (any(is.na(gsModel@gv[[1]]@addEff)) |
          any(is.na(masModel_PERFECT@gv[[1]]@addEff)) |
          any(is.na(masModel_LOW@gv[[1]]@addEff)) |
          any(is.na(masModel_HIGH@gv[[1]]@addEff))) {
        return(result)
      }
    }

    gsPop <- selectCross(gsPop, use="ebv", nInd=nInd(gsPop)*n.selInt, nCrosses=nInd(gsPop))
    psPop <- selectCross(psPop, trait=breedingFitness, nInd=nInd(psPop)*n.selInt, nCrosses=nInd(psPop))
    masPop_PERFECT <- selectCross(masPop_PERFECT, use="ebv", nInd=nInd(masPop_PERFECT)*n.selInt, nCrosses=nInd(masPop_PERFECT))
    masPop_LOW <- selectCross(masPop_LOW, use="ebv", nInd=nInd(masPop_LOW)*n.selInt, nCrosses=nInd(masPop_LOW))
    masPop_HIGH <- selectCross(masPop_HIGH, use="ebv", nInd=nInd(masPop_HIGH)*n.selInt, nCrosses=nInd(masPop_HIGH))
  }
  return (result)
}
