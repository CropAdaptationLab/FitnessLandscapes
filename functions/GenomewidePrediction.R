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

# Train an RRBLUP model to predict breeding fitness
# trainPop: the training population
# testPop: the test population
# trait: the AlphaSimR phenotype index
# Return correlation (r) between the EBVs and the actual genetic values in the test pop
evaluateGwpModel <- function(trainPop, testPop) {
  # Update phenotype to have heritability associated with breeding programs
  trainPop <- setPheno(trainPop,
                       h2=c(n.h2Breeding, n.h2Breeding, n.yieldH2Breeding),
                       reps=n.rilReps)
  # snpChip 2 is for genomic prediction
  if (GS_MODEL == "RRBLUP") {
    model <- RRBLUP(trainPop, traits=calculateW_GWP, use=GS_PHENO, snpChip=2)
    testPop <- setEBV(testPop, model)
  } else if (GS_MODEL == "fastRRBLUP") {
    model <- fastRRBLUP(trainPop, traits=calculateW_GWP, use=GS_PHENO, snpChip=2)
    testPop <- setEBV(testPop, model)
  } else if (GS_MODEL == "GBLUP") {
    testPop <- calculateEBV_GBLUP(trainPop=trainPop, testPop=testPop)
  }
  # Determine the correlation between genetic values and estimated breeding values
  # in the test population
  r <- cor(calculateW_GWP(gv(testPop)), ebv(testPop))[1]
  return (r)
}

# Calculates the 5-fold cross validation accuracy for a GWP model trained on 'pop'
crossValAccuracy <- function(pop) {
  
  # Number of folks
  nSplits <- 5
  # Number of individuals
  popSize <- pop@nInd
  
  # Stores the results
  results <- numeric(length=nSplits)
  
  # For a population of size 1000, the test pops are:
  # [1,200], [201,400], etc.
  for (c in 1:nSplits) {
    start <- (c-1) * (popSize %/% nSplits) + 1
    end <- c* (popSize %/% nSplits)
    # Determine the sequence of individuals in the test pop
    testPop <- pop[start:end]
    # Train on the remaining individuals
    trainPop <- pop[-(start:end)]
    # Evaluate the GWP model and save the result
    results[c] <- evaluateGwpModel(trainPop, testPop)
  }
  
  # Cross-validation accuracy is the mean
  return(mean(results))
}

# Calculates estimated breeding values (EBV) with GBLUP
calculateEBV_GBLUP <- function(trainPop, testPop) {
  # Update phenotypes
  trainPop <- setPheno(trainPop, h2=c(n.h2Breeding, n.h2Breeding, n.yieldH2Breeding))

  # Pull genotype matrices for both populations
  geno <- rbind(pullSnpGeno(trainPop, snpChip = 2),
                pullSnpGeno(testPop, snpChip = 2))
  
  # Get breeding fitness values for train pop
  w <- data.frame(pheno(trainPop)) %>%
    dplyr::mutate(w=calculateBreedingFitness(Trait1, Trait2, Trait3)) %>%
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
createTrainPop <- function(curPop,
                           prevTrainPop,
                           n=n.trainPopSize) {
  topFams <- selectFam(curPop, nFam=n.topFamilies, use="ebv")
  topInds <- selectWithinFam(topFams, use="ebv", nInd=1)
  prevTrainPop <- setPheno(prevTrainPop, h2=c(n.h2Breeding, n.h2Breeding, n.yieldH2Breeding))
  prevTrainPop <- selectInd(prevTrainPop, trait=breedingFitness, nInd=n-nInd(topInds))
  return(c(topInds, prevTrainPop))
}

# Calculates metrics for a population at a particular cycle
# pop: An AlphaSimR population
# sel: The selection type (e.g. "GS", "PS", etc.)
# Returns: a dataframe tabulating all metrics
cycleMetrics <- function(pop, season, cycle, sel) {
  # Calulate accuracy of predictions
  r <- NA
  if (!(sel %in% c("PS", "PSieMAS"))) {
    r <- cor(calculateW_GWP(gv(pop)), ebv(pop))[1]
  }
  # Update the phenotypes
  #pop <- setPheno(pop, h2=c(n.h2Breeding, n.h2Breeding, n.yieldH2Breeding))
  
  # Genome-wide heterozygosity
  genome_het <- meanHetLocus(pullSnpGeno(pop, snpChip=2))
  # Attained trait heterozygosity
  attained_het <- meanHetLocus(getUniqueQtl(pop))
  # Desired trait heterozygosity
  desired_het <- meanHetLocus(pullQtlGeno(pop, trait=3))
  
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
    season=season,
    c=cycle,
    sel=sel,
    w=w,
    wGV=wGV,
    r=r,
    genome_het=genome_het,
    attained_het=attained_het,
    desired_het=desired_het,
    gvar=varG(pop)[3,3]))
}

# Carry out marker-assisted selection
getMasInds <- function(pop, peaks, masSelInt, snpChip, useQtls) {
  masHaplo <- data.frame(id=pop@id)
  if (nrow(peaks) > 0) {
    if (useQtls) {
      masGeno <- getUniqueQtl(pop)
    } else {
      masGeno <- as.data.frame(pullSnpGeno(pop, snpChip))
    }
    
    # Interaction size
    LOD <- c()
    for (r in 1:nrow(peaks)) {
      # Interactions
      if (useQtls) {
        m1 <- peaks$qtl1[r]
        m2 <- peaks$qtl2[r]
      } else {
        m1 <- peaks$m1[r]
        m2 <- peaks$m2[r]
      }
      if (!is.na(m1) & !is.na(m2)) {
        if (segLocus(masGeno[,m1]) & segLocus(masGeno[,m2])) {
          haplos <- getHaplos(masGeno, m1, m2, parent1, parent2, snpChip, useQtls)
          LOD <- c(LOD, peaks$lod.int[r])
          masHaplo <- cbind(masHaplo, haplos)
        }
      }
    }
    # Number of QTL interactions
    if (length(LOD) > 0) {
      colnames(masHaplo) <- c("id", paste0("INT_", LOD))
      # Calculate the majority haplotype between P1 and P2 in each column
      majority_haplo <- masHaplo %>%
        dplyr::select(starts_with("INT_")) %>%
        dplyr::summarise(across(everything(), ~ {
          p1_count <- sum(. == "P1", na.rm = TRUE)
          p2_count <- sum(. == "P2", na.rm = TRUE)
          if (p1_count >= p2_count) "P1" else "P2"
        })) %>%
        as.character()
      
      # Pull out the names of the columsn
      int_cols <- names(masHaplo)[startsWith(names(masHaplo), "INT_")]
      
      # If the haplotype matches the majority haplotype, set it as a '1', otherwise as a '0'
      masHaplo[int_cols] <- (masHaplo[int_cols] == rep(majority_haplo, each = nrow(masHaplo))) * 1
      
      # Ensure the columns are ordered by LOD score (descending)
      sorted_int_cols <- int_cols[order(LOD, decreasing = TRUE)]
      masHaplo <- masHaplo[, c("id", sorted_int_cols)]
    }
  }
  
  # Order by whether individuals have each favorable haplotype, then by 'w'
  orderedInds <- as.data.frame(pheno(pop)) %>%
    dplyr::mutate(id=pop@id,
                  w=calculateBreedingFitness(Trait1, Trait2, Trait3)) %>%
    dplyr::select(-c("Trait1", "Trait2", "Trait3")) %>%
    dplyr::left_join(masHaplo, by="id") %>%
    dplyr::arrange(
      across(starts_with("INT_"), desc),
      desc(w)
    ) %>%
    dplyr::slice_head(n=nInd(pop)*masSelInt) %>%
    dplyr::pull(id)
  return(pop[orderedInds])
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
  results_list <- list(cycleMetrics(basePop, 0, 0, "PS"),
                       cycleMetrics(basePop, 0, 0, "PS_ieMAS"),
                       cycleMetrics(basePop, 0, 0, "GS"),
                       cycleMetrics(basePop, 0, 0, "ieMAS"))
                       #cycleMetrics(basePop, 0, 0, "ieMAS_perfect"),
                       #cycleMetrics(basePop, 0, 0, "ieMAS_low"),
                       #cycleMetrics(basePop, 0, 0, "GS_noUpdate"),
                       #cycleMetrics(basePop, 0, 0, "ieMAS_noUpdate"))

  # Set phenotypes for base population
  basePop <- setPheno(basePop, h2=c(n.h2Breeding, n.h2Breeding, n.yieldH2Breeding), reps=n.rilReps)

  # 2D linkage mapping
  # Find all pairs of epistatic loci
  peaks <- epistaticLodPeaks(basePop, parent1, parent2, snpChip=2, trait=5)
  # Low resolution mapping
  #low_res_peaks <- epistaticLodPeaks(basePop, parent1, parent2, snpChip=3, trait=5)
  
  # Select the first cycle based on phenotype
  topRILs <- selectInd(basePop,
                       nInd=n.topFamilies,
                       trait=breedingFitness)
  ps_S0 <- randCross(topRILs,
                     nCrosses=n.families,
                     nProgeny=20)
  ps_S0 <- setPheno(ps_S0, h2=c(n.h2Breeding, n.h2Breeding, n.yieldH2Breeding))
  psMAS_S0 <- ps_S0
  gs_S0 <- randCross(topRILs,
                     nCrosses=n.families,
                     nProgeny=10)
  gs_S0 <- setPheno(gs_S0, h2=c(n.h2Breeding, n.h2Breeding, n.yieldH2Breeding))
  gsMAS_high_S0 <- gs_S0
  gsMAS_perfect_S0 <- gs_S0
  gsMAS_low_S0 <- gs_S0
  gs_noUpdate_S0 <- gs_S0
  gsMAS_high_noUpdate_S0 <- gs_S0
  
  # Initial training population is the top 40 RILs, selected by phenotype, plus
  # 360 randomly chosen RILs from the remaining 960 lines
  randLines <- selectInd(basePop[-match(topRILs@id, basePop@id)],
                         use="rand",
                         nInd=n.trainPopSize-n.topFamilies)
  gsTrainPop <- c(topRILs, randLines)
  # This will be the same training population throughout all cycles
  gsTrainPop_HIGH <- gsTrainPop
  #gsTrainPop_PERFECT <- gsTrainPop
  #gsTrainPop_LOW <- gsTrainPop
  #gsTrainPop_noUpdate <- gsTrainPop
  #gsTrainPop_HIGH_noUpdate <- gsTrainPop
  
  # Fit an RR-BLUP model to use initially for the GS populations
  if (GS_MODEL %in% c("RRBLUP", "fastRRBLUP")) {
    # Train a model on the initial base population
    gsModel <- trainRRBLUPModel(gsTrainPop)
    gsModel_MAS_HIGH <- gsModel
    #gsModel_MAS_PERFECT <- gsModel
    #gsModel_MAS_LOW <- gsModel
    #gsModel_MAS_HIGH_noUpdate <- gsModel
    # This model will not be updated throughout the recurrent selection
    #gsModel_noUpdate <- gsModel
  }

  # There are 3 years per phenotypic selection cycle
  PHENO_CYCLES <- n.Y / 3
  
  # PHENOTYPIC SELECTION
  season <- 1
  for (cycle in 1:PHENO_CYCLES) {

    # --------- SEASON 1 (WINTER) ---------
    
    results_list[[length(results_list) + 1]] <- cycleMetrics(ps_S0, season, cycle, "PS")
    results_list[[length(results_list) + 1]] <- cycleMetrics(psMAS_S0, season, cycle, "PS_ieMAS")
    ps_S1 <- self(ps_S0)
    psMAS_S1 <- self(psMAS_S0)
    season <- season + 1
    
    # --------- SEASON 2 (SUMMER) ---------
    
    ps_S1 <- setPheno(ps_S1, h2=c(n.h2Breeding, n.h2Breeding, n.yieldH2Breeding))
    psMAS_S1 <- setPheno(psMAS_S1, h2=c(n.h2Breeding, n.h2Breeding, n.yieldH2Breeding))
    results_list[[length(results_list) + 1]] <- cycleMetrics(ps_S1, season, cycle, "PS")
    results_list[[length(results_list) + 1]] <- cycleMetrics(psMAS_S1, season, cycle, "PS_ieMAS")
    ps_S1 <- selectWithinFam(ps_S1, nInd=10, trait=suitability)
    ps_S2 <- self(ps_S1)
    # MAS
    if (cycle == 1) {
      psMAS_S1 <- getMasInds(pop=psMAS_S1, peaks=peaks, masSelInt=n.masSelInt, snpChip=2, useQtls=FALSE)
    }
    psMAS_S1 <- selectWithinFam(psMAS_S1, nInd=10, trait=suitability)
    psMAS_S2 <- self(psMAS_S1)
    
    season <- season + 1

    # --------- SEASON 3 (WINTER) ---------
    
    ps_S2 <- setPheno(ps_S2, h2=c(n.h2Breeding, n.h2Breeding, n.yieldH2Breeding))
    psMAS_S2 <- setPheno(psMAS_S2, h2=c(n.h2Breeding, n.h2Breeding, n.yieldH2Breeding))
    results_list[[length(results_list) + 1]] <- cycleMetrics(ps_S2, season, cycle, "PS")
    results_list[[length(results_list) + 1]] <- cycleMetrics(psMAS_S2, season, cycle, "PS_ieMAS")
    ps_S3 <- self(ps_S2, nProgeny=2)
    psMAS_S3 <- self(psMAS_S2, nProgeny=2)
    season <- season + 1
    
    # --------- SEASON 4 (SUMMER) ---------
    
    ps_S3 <- setPheno(ps_S3, h2=c(n.h2Breeding, n.h2Breeding, n.yieldH2Breeding))
    psMAS_S3 <- setPheno(psMAS_S3, h2=c(n.h2Breeding, n.h2Breeding, n.yieldH2Breeding))
    results_list[[length(results_list) + 1]] <- cycleMetrics(ps_S3, season, cycle, "PS")
    results_list[[length(results_list) + 1]] <- cycleMetrics(psMAS_S3, season, cycle, "PS_ieMAS")
    # Make selections based on realized yield
    ps_S3 <- selectFam(ps_S3, trait=breedingFitness, nFam=n.topFamilies)
    ps_S4 <- self(ps_S3)
    psMAS_S3 <- selectFam(psMAS_S3, trait=breedingFitness, nFam=n.topFamilies)
    psMAS_S4 <- self(psMAS_S3)
    season <- season + 1
    
    # --------- SEASON 5 (WINTER) ---------
    
    ps_S4 <- setPheno(ps_S4, h2=c(n.h2Breeding, n.h2Breeding, n.yieldH2Breeding))
    psMAS_S4 <- setPheno(psMAS_S4, h2=c(n.h2Breeding, n.h2Breeding, n.yieldH2Breeding))
    results_list[[length(results_list) + 1]] <- cycleMetrics(ps_S4, season, cycle, "PS")
    results_list[[length(results_list) + 1]] <- cycleMetrics(psMAS_S4, season, cycle, "PS_ieMAS")
    ps_S5 <- self(ps_S4)
    psMAS_S5 <- self(psMAS_S4)
    season <- season + 1
    
    # --------- SEASON 6 (SUMMER) ---------
    ps_S5 <- setPheno(ps_S5, h2=c(n.h2Breeding, n.h2Breeding, n.yieldH2Breeding))
    psMAS_S5 <- setPheno(psMAS_S5, h2=c(n.h2Breeding, n.h2Breeding, n.yieldH2Breeding))
    results_list[[length(results_list) + 1]] <- cycleMetrics(ps_S5, season, cycle, "PS")
    results_list[[length(results_list) + 1]] <- cycleMetrics(psMAS_S5, season, cycle, "PS_ieMAS")
    ps_topLines <- selectWithinFam(ps_S5, nInd=1, trait=breedingFitness)
    psMAS_topLines <- selectWithinFam(psMAS_S5, nInd=1, trait=breedingFitness)
    season <- season + 1
    
    ps_S0 <- randCross(ps_topLines,
                       nCrosses=n.families,
                       nProgeny=20)
    psMAS_S0 <- randCross(psMAS_topLines,
                           nCrosses=n.families,
                           nProgeny=20)
  }
  
  # GENOMIC SELECTION
  season <- 1
  for (cycle in 1:n.Y) {
    
    # --------- SEASON 1 (WINTER) ---------
    results_list[[length(results_list) + 1]] <- cycleMetrics(gs_S0, season, cycle, "GS")
    results_list[[length(results_list) + 1]] <- cycleMetrics(gsMAS_high_S0, season, cycle, "ieMAS")
    #results_list[[length(results_list) + 1]] <- cycleMetrics(gsMAS_perfect_S0, season, cycle, "ieMAS_perfect")
    #results_list[[length(results_list) + 1]] <- cycleMetrics(gsMAS_low_S0, season, cycle, "ieMAS_low")
    #results_list[[length(results_list) + 1]] <- cycleMetrics(gs_noUpdate_S0, season, cycle, "GS_noUpdate")
    #results_list[[length(results_list) + 1]] <- cycleMetrics(gsMAS_high_noUpdate_S0, season, cycle, "ieMAS_noUpdate")
    
    gs_S1 <- self(gs_S0)
    gsMAS_high_S1 <- self(gsMAS_high_S0)
    #gsMAS_perfect_S1 <- self(gsMAS_perfect_S0)
    #gsMAS_low_S1 <- self(gsMAS_low_S0)
    #gs_noUpdate_S1 <- self(gs_noUpdate_S0)
    #gsMAS_high_noUpdate_S1 <- self(gsMAS_high_noUpdate_S0)
    
    season <- season + 1
    
    # --------- SEASON 2 (SUMMER) ---------

    if (GS_MODEL %in% c("RRBLUP", "fastRRBLUP")) {
      # If the new models do not fit any values, there is no genetic variance
      # in the population
      if (any(is.na(gsModel@gv[[1]]@addEff)) |
          any(is.na(gsModel_MAS_HIGH@gv[[1]]@addEff))) {
          #any(is.na(gsModel_MAS_PERFECT@gv[[1]]@addEff)) |
          #any(is.na(gsModel_MAS_LOW@gv[[1]]@addEff)) |
          #any(is.na(gsModel_noUpdate@gv[[1]]@addEff)) |
          #any(is.na(gsModel_MAS_HIGH_noUpdate@gv[[1]]@addEff))) {
        return(do.call(rbind, results_list))
      }
      gs_S1 <- setEBV(gs_S1, gsModel)
      gsMAS_high_S1 <- setEBV(gsMAS_high_S1, gsModel_MAS_HIGH)
      #gsMAS_perfect_S1 <- setEBV(gsMAS_perfect_S1, gsModel_MAS_PERFECT)
      #gsMAS_low_S1 <- setEBV(gsMAS_low_S1, gsModel_MAS_LOW)
      #gs_noUpdate_S1 <- setEBV(gs_noUpdate_S1, gsModel_noUpdate)
      #gsMAS_high_noUpdate_S1 <- setEBV(gsMAS_high_noUpdate_S1, gsModel_noUpdate)
    } else if (GS_MODEL == "GBLUP") {
      gs_S1 <- calculateEBV_GBLUP(trainPop=gsTrainPop, testPop=gs_S1)
      gsMAS_high_S1 <- calculateEBV_GBLUP(trainPop=gsTrainPop_HIGH, testPop=gsMAS_high_S1)
      #gsMAS_perfect_S1 <- calculateEBV_GBLUP(trainPop=gsTrainPop_PERFECT, testPop=gsMAS_perfect_S1)
      #gsMAS_low_S1 <- calculateEBV_GBLUP(trainPop=gsTrainPop_LOW, testPop=gsMAS_low_S1)
      #gs_noUpdate_S1 <- calculateEBV_GBLUP(trainPop=gsTrainPop_noUpdate, testPop=gs_noUpdate_S1)
      #gsMAS_high_noUpdate_S1 <- calculateEBV_GBLUP(trainPop=gsTrainPop_HIGH_noUpdate, testPop=gsMAS_high_noUpdate_S1)
    }
    
    results_list[[length(results_list) + 1]] <- cycleMetrics(gs_S1, season, cycle, "GS")
    results_list[[length(results_list) + 1]] <- cycleMetrics(gsMAS_high_S1, season, cycle, "ieMAS")
    #results_list[[length(results_list) + 1]] <- cycleMetrics(gsMAS_perfect_S1, season, cycle, "ieMAS_perfect")
    #results_list[[length(results_list) + 1]] <- cycleMetrics(gsMAS_low_S1, season, cycle, "ieMAS_low")
    #results_list[[length(results_list) + 1]] <- cycleMetrics(gs_noUpdate_S1, season, cycle, "GS_noUpdate")
    #results_list[[length(results_list) + 1]] <- cycleMetrics(gsMAS_high_noUpdate_S1, season, cycle, "ieMAS_noUpdate")

    # Establish new training population by selecting the top lines from the populations by gebv
    gsTrainPop <- createTrainPop(curPop=gs_S1, prevTrainPop=gsTrainPop)
    gsTrainPop_HIGH <- createTrainPop(curPop=gsMAS_high_S1, prevTrainPop=gsTrainPop_HIGH)
    #gsTrainPop_PERFECT <- createTrainPop(curPop=gsMAS_perfect_S1, prevTrainPop=gsTrainPop_PERFECT)
    #gsTrainPop_LOW <- createTrainPop(curPop=gsMAS_low_S1, prevTrainPop=gsTrainPop_LOW)
    # RETRAIN MODELS
    if (GS_MODEL %in% c("RRBLUP", "fastRRBLUP")) {
      gsModel <- trainRRBLUPModel(gsTrainPop)
      gsModel_MAS_HIGH <- trainRRBLUPModel(gsTrainPop_HIGH)
      #gsModel_MAS_PERFECT <- trainRRBLUPModel(gsTrainPop_PERFECT)
      #gsModel_MAS_LOW <- trainRRBLUPModel(gsTrainPop_LOW)
    }
    
    # MARKER-ASSISTED SELECTION
    if (cycle == 1) {
      gsMAS_high_S1 <- getMasInds(pop=gsMAS_high_S1, peaks=peaks, masSelInt=n.masSelInt, snpChip=2, useQtls=FALSE)
      #gsMAS_perfect_S1 <- getMasInds(pop=gsMAS_perfect_S1, peaks=peaks, masSelInt=n.masSelInt, snpChip=2, useQtls=TRUE)
      #gsMAS_low_S1 <- getMasInds(pop=gsMAS_low_S1, peaks=low_res_peaks, masSelInt=n.masSelInt, snpChip=3, useQtls=FALSE)
      #gsMAS_high_noUpdate_S1 <- getMasInds(pop=gsMAS_high_noUpdate_S1, peaks=peaks, masSelInt=n.masSelInt, snpChip=2, useQtls=FALSE)
    }
    
    # Select the top families, then the top individual out of those selected families,
    # and randomly intermate those lines
    advanceGsPop <- function(pop) {
      topFams <- selectFam(pop, nFam=n.topFamilies, use="ebv")
      topLines <- selectWithinFam(topFams, nInd=1, use="ebv")
      S0 <- randCross(topLines,
                      nCrosses=n.families,
                      nProgeny=10)
      return(S0)
    }
    gs_S0 <- advanceGsPop(gs_S1)
    gsMAS_high_S0 <- advanceGsPop(gsMAS_high_S1)
    #gsMAS_perfect_S0 <- advanceGsPop(gsMAS_perfect_S1)
    #gsMAS_low_S0 <- advanceGsPop(gsMAS_low_S1)
    #gs_noUpdate_S0 <- advanceGsPop(gs_noUpdate_S1)
    #gsMAS_high_noUpdate_S0 <- advanceGsPop(gsMAS_high_noUpdate_S1)
  
    season <- season + 1
  }
  return (do.call(rbind, results_list))
}
