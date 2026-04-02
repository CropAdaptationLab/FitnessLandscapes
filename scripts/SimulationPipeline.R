# Title: SIMULATION PIPELINE
# Author: Ted Monyak
# Description:
# This will run n.popResets * n.reps simulations, and for each replication,
# create two subpopulations, then two biparental RIL populations:
# 1 "admixed" (parents come from different subpopulations),
# 1 "unadmixed" (parents come from the same subpopulation)

library(AlphaSimR)
library(devtools)
library(dplyr)
library(ggpubr)
library(ggplot2)
library(ggpmisc)
library(patchwork)
library(plotly)
library(purrr)
library(qtl)
library(qtl2)
library(qtl2convert)
library(RColorBrewer)
library(reshape2)
library(tibble)
library(tidyr)
library(viridis)

rm(list = ls())

set.seed(123)

# Set to false if running on Alpine
runLocal = TRUE

if (runLocal) {
  setwd("~/Documents/CSU/FitnessLandscapes")
  n.cores <- 8
} else {
  setwd("/pl/active/Morris_CSU/Ted_Monyak/FitnessLandscapes")
  n.cores <- 4
}
base_dir <- file.path(getwd(), "output/AggregatedResults")
if (!dir.exists(base_dir)) dir.create(base_dir)

source("functions/Fitness.R")
source("functions/GenoConversions.R")
source("functions/MappingPopulations.R")
source("functions/PopGen.R")
source("functions/QtlMapping.R")
source("functions/QuantGen.R")
source("functions/TraitArchitecture.R")
source("figures/G2FLandscapes.R")
source("figures/LinkageMaps.R")
source("figures/ReactionNorm.R")
source("figures/T2FLandscapes.R")
source("scripts/GlobalParameters.R")

# Number of founder populations to simulate
n.popResets <- 1
# Number of adaptive walk replications per pair of subpopulations
n.reps <- 10

SAMPLING <- 'geometric'
n.minMAF <- 0.3

# FUNCTIONS
qtlMapping <- TRUE
twoQtlMapping <- TRUE
compareUnadmixed <- TRUE

# SAVE ALLELE FREQUENCIES
saveFixationOrder <- TRUE
saveEffectSizes <- TRUE

# PLOTTING
saveQtlPlots <- FALSE
saveTraitPlots <- FALSE
saveAllelePlots <- FALSE
saveFitnessPlots <- FALSE

# Set to true for doing G > F landscape creation (not currently working)
sampleInds <- FALSE

# All the parameter combinations to iterate through
qtl_vec <- c(10,20,50)

output_dir <- file.path(base_dir, paste0("Sim_", format(Sys.time(), "%F_%H_%M")))
if (!dir.exists(output_dir)) dir.create(output_dir)

# Result dataframe
res.df <- data.frame(var=c(), a=c(), popSize=c(), qtl=c(), pop=c(), rep=c(), type=c(), fst=c(),
                     nLod_T1=c(), nLod_T2=c(), nLod_T3=c(), nLod_Suit=c(), nLod_W=c(),
                     nLod_Int=c(),
                     ev_T1=c(), ev_T2=c(), ev_T3=c(), ev_Suit=c(), ev_W=c(),
                     isoElite_T1=c(), isoElite_T2=c(), isoElite_T3=c(),
                     hamm_T1=c(), hamm_T2=c(),
                     initA_T1=c(), initA_T2=c(), initVar_T1=c(), initVar_T2=c(),
                     rilA_T1=c(), rilA_T2=c(), rilVar_T1=c(), rilVar_T2=c())

# Effect size dataframe
effectSizes.df <- data.frame(
  rank=c(),
  id=c(),
  eff_size=c(),
  qtl=c())

# Fixation order dataframe
fixedAlleles.df <- data.frame(
  id=c(),
  eff_size=c(),
  rank=c(),
  subpop=c(),
  qtl=c())

for (lx in 1:length(qtl_vec)) {
  n.L <- qtl_vec[lx]
  print(paste0("QTL: ", n.L))
  
  # The directory for this number of QTL
  sim_dir <- file.path(output_dir, paste0("QTL_", n.L))
  if (!dir.exists(sim_dir)) dir.create(sim_dir)
  
  # Reset the founder population n.popResets times
  for (f in 1:n.popResets) {
    pop_dir <- file.path(sim_dir, paste0("FounderPopulation", f))
    if (!dir.exists(pop_dir)) dir.create(pop_dir)
    print(paste0("Founder Reset ", f))
    source("scripts/CreateFounderPop.R")
    # Get the largest effect size and variance of the two attained traits
    t1_arch <- singleTraitArchitecture(founderPop,1)$eff_size
    t2_arch <- singleTraitArchitecture(founderPop,2)$eff_size
    initA_T1 <- max(t1_arch)
    initA_T2 <- max(t2_arch)
    initVar_T1 <- var(t1_arch)
    initVar_T2 <- var(t2_arch)
    
    # For each founder population, create independent populations n.sims times
    for (r in 1:n.reps) {
      print(paste0("Rep ", r))
      rep_dir <- file.path(pop_dir, paste0("Rep", r))
      if (!dir.exists(rep_dir)) dir.create(rep_dir)
      
      suitFunc <- suitabilityGaussian
      selectionPheno <- breedingFitness
      
      # Landrace heritability
      SP$setVarE(h2=c(n.h2, n.h2, n.yieldH2))
      source("Scripts/CreateIndependentPops.R")
      # Calculate FST
      fst <- FST(pops)
      
      # Update to use the suitability selection function
      selectionPheno <- suitability
      SP$setVarE(h2=c(n.h2Breeding, n.h2Breeding, n.yieldH2Breeding))
      
      ril_dir <- file.path(rep_dir, "Admixed")
      if (!dir.exists(ril_dir)) dir.create(ril_dir)

      # Develop purelines from each population
      purelines1 <- makePurelinesBulk(pops[[1]])
      pop1 <- purelines1[[1]]
      purelines2 <- makePurelinesBulk(pops[[2]])
      pop2 <- purelines2[[1]]
      
      # Create a biparental RIL population by sampling one individual from each subpopulation
      res <- createRIL(pop1, pop2)
      # Select the parents from the result, and the RIL population
      RIL <- res[-(1:2)]
      parent1 <- res[1]
      parent2 <- res[2]
      
      # Get effect sizes of each attained trait in the RIL
      effSizes_t1 <- singleTraitArchitecture(RIL, 1)
      if (nrow(effSizes_t1) > 0) {
        effSizes_t1$qtl <- n.L
        effectSizes.df <- rbind(effectSizes.df, effSizes_t1)
        rilA_T1 <- max(effSizes_t1$eff_size)
        rilVar_T1 <- var(effSizes_t1$eff_size)
      } else {
        rilA_T1 <- 0
        rilVar_T1 <- 0
      }

      effSizes_t2 <- singleTraitArchitecture(RIL, 2)
      if (nrow(effSizes_t2) > 0) {
        effSizes_t2$qtl <- n.L
        effectSizes.df <- rbind(effectSizes.df, effSizes_t2)
        rilA_T2 <- max(effSizes_t2$eff_size)
        rilVar_T2 <- var(effSizes_t2$eff_size)
      } else {
        rilA_T2 <- 0
        rilVar_T2 <- 0
      }

      # Calculate isoeliteness for each trait
      isoElite_T1 <- isoEliteness(parent1, parent2, founderPop, 1)
      isoElite_T2 <- isoEliteness(parent1, parent2, founderPop, 2)
      isoElite_T3 <- isoEliteness(parent1, parent2, founderPop, 3)
      
      # Calculate Hamming distance for each attained trait
      hamm_T1 <- hammingDistance(parent1, parent2, 1)
      hamm_T2 <- hammingDistance(parent1, parent2, 2)
      source("figures/GenerateSubpopPlots.R")
      
      # Get the phenotypes of the RIL family and the purelines
      RIL_pheno <- as.data.frame(pheno(RIL)) %>%
        dplyr::mutate(Suitability=suitFunc(Trait1, Trait2),
                      W=calculateBreedingFitness(Trait1, Trait2, Trait3),
                      Population="RIL Family")
      pureline1_pheno <- as.data.frame(pheno(pop1)) %>%
        dplyr::mutate(Suitability=suitFunc(Trait1, Trait2),
                      W=calculateBreedingFitness(Trait1, Trait2, Trait3),
                      Population="Pureline 1")
      pureline2_pheno <- as.data.frame(pheno(pop2)) %>%
        dplyr::mutate(Suitability=suitFunc(Trait1, Trait2),
                      W=calculateBreedingFitness(Trait1, Trait2, Trait3),
                      Population="Pureline 2")

      pureline_pheno <- rbind(pureline1_pheno, pureline2_pheno)
      
      # Calculate the increase in RIL variance compared to the pureline variance
      ev_T1 <- excessVariance(RIL_pheno[,1], pureline_pheno[,1])
      ev_T2 <- excessVariance(RIL_pheno[,2], pureline_pheno[,2])
      ev_T3 <- excessVariance(RIL_pheno[,3], pureline_pheno[,3])
      ev_Suit <- excessVariance(RIL_pheno[,4], pureline_pheno[,4])
      ev_W <- excessVariance(RIL_pheno[,5], pureline_pheno[,5])
      source("figures/TransgressiveSegregation.R")

      # Determine the number of significant LOD peaks from single linkage mapping
      if (qtlMapping) {
        nLOD <- getLodPeaks(RIL, parent1, parent2, ril_dir)
      } else {
        nLOD <- c(0,0,0,0,0)
      }
      # Get number of significant interaction effects
      if (twoQtlMapping) {
        # Get the peak locations
        peaks <- epistaticLodPeaks(RIL, parent1, parent2, trait=5, ril_dir)
        if (nrow(peaks) > 0) {
          maxLodIdx <- which.max(peaks$lod.int)
          qtl1 <- peaks$qtl1[maxLodIdx]
          qtl2 <- peaks$qtl2[maxLodIdx]
          if ((qtl1 != qtl2) & (saveQtlPlots)) {
            plotReactionNorm(RIL, qtl1, qtl2, parent1, parent2, suitFunc, ril_dir)
            plot1DLandscape(RIL, pops[[1]], pops[[2]], parent1, parent2, qtl1, qtl2, ril_dir)
          }
        }
        nIntPeaks <- nrow(peaks)
      } else {
        nIntPeaks <- 0
      }

      # Update the results table
      res.df <- rbind(res.df, data.frame(var=n.var, a=n.a, popSize=n.subPopSize, qtl=n.L, pop=f, rep=r, type="Admixed", fst=fst,
                                         nLod_T1=nLOD[1], nLod_T2=nLOD[2], nLod_T3=nLOD[3], nLod_Suit=nLOD[4],
                                         nLod_W=nLOD[5], nLod_Int=nIntPeaks,
                                         ev_T1=ev_T1, ev_T2=ev_T2, ev_T3=ev_T3, ev_Suit=ev_Suit, ev_W=ev_W,
                                         isoElite_T1=isoElite_T1, isoElite_T2=isoElite_T2, isoElite_T3=isoElite_T3,
                                         hamm_T1=hamm_T1, hamm_T2=hamm_T2,
                                         initA_T1=initA_T1,
                                         initA_T2=initA_T2,
                                         initVar_T1=initVar_T1,
                                         initVar_T2=initVar_T2,
                                         rilA_T1=rilA_T1,
                                         rilA_T2=rilA_T2,
                                         rilVar_T1=rilVar_T1,
                                         rilVar_T2=rilVar_T2))
      
      # Create a biparental RIL population with purelines derived from the same subpopulation
      if (compareUnadmixed) {
        ril_dir <- file.path(rep_dir, "Unadmixed")
        if (!dir.exists(ril_dir)) dir.create(ril_dir)
        
        pop1 <- purelines1[[1]]
        pop2 <- purelines1[[2]]
        res <- createRIL(pop1, pop2)
        RIL <- res[-(1:2)]
        parent1 <- res[1]
        parent2 <- res[2]
        isoElite_T1 <- isoEliteness(parent1, parent2, founderPop, 1)
        isoElite_T2 <- isoEliteness(parent1, parent2, founderPop, 2)
        isoElite_T3 <- isoEliteness(parent1, parent2, founderPop, 3)
        hamm_T1 <- hammingDistance(parent1, parent2, 1)
        hamm_T2 <- hammingDistance(parent1, parent2, 2)
        if (qtlMapping) {
          nLOD <- getLodPeaks(RIL, parent1, parent2, ril_dir)
        } else {
          nLOD <- c(0,0,0,0,0)
        }
        
        if (twoQtlMapping) {
          peaks <- epistaticLodPeaks(RIL, parent1, parent2, trait=4, ril_dir)
          nIntPeaks <- nrow(peaks)
        } else {
          nIntPeaks <- 0
        }
        RIL_pheno <- as.data.frame(pheno(RIL)) %>%
          dplyr::mutate(Suitability=suitFunc(Trait1, Trait2),
                        W=calculateBreedingFitness(Trait1, Trait2, Trait3))
        pureline_pheno <- as.data.frame(pheno(c(pop1, pop2))) %>%
          dplyr::mutate(Suitability=suitFunc(Trait1, Trait2),
                        W=calculateBreedingFitness(Trait1, Trait2, Trait3))
        
        ev_T1 <- excessVariance(RIL_pheno[,1], pureline_pheno[,1])
        ev_T2 <- excessVariance(RIL_pheno[,2], pureline_pheno[,2])
        ev_T3 <- excessVariance(RIL_pheno[,3], pureline_pheno[,3])
        ev_Suit <- excessVariance(RIL_pheno[,4], pureline_pheno[,4])
        ev_W <- excessVariance(RIL_pheno[,5], pureline_pheno[,5])

        # Update the result dataframe
        res.df <- rbind(res.df, data.frame(var=n.var, a=n.a, popSize=n.subPopSize, qtl=n.L, pop=f, rep=r, type="Unadmixed", fst=NA,
                                           nLod_T1=nLOD[1], nLod_T2=nLOD[2], nLod_T3=nLOD[3], nLod_Suit=nLOD[4],
                                           nLod_W=nLOD[5], nLod_Int=nIntPeaks,
                                           ev_T1=ev_T1, ev_T2=ev_T2, ev_T3=ev_T3, ev_Suit=ev_Suit, ev_W=ev_W,
                                           isoElite_T1=isoElite_T1, isoElite_T2=isoElite_T2, isoElite_T3=isoElite_T3,
                                           hamm_T1=hamm_T1, hamm_T2=hamm_T2,
                                           initA_T1=NA, initA_T2=NA, initVar_T1=NA, initVar_T2=NA,
                                           rilA_T1=NA, rilA_T2=NA, rilVar_T1=NA, rilVar_T2=NA))
      } # end compareIntra
    } # end n.reps
  } # end n.popResets
  
  if (saveEffectSizes) {
    write.table((effectSizes.df %>% dplyr::filter(qtl==n.L)),
                file.path(sim_dir, "effect_size.csv"), col.names=TRUE, quote=FALSE, sep=",")
    
  }
  if (saveFixationOrder) {
    write.table((fixedAlleles.df %>% dplyr::filter(qtl==n.L)),
                file.path(sim_dir, "fixed_alleles.csv"), col.names=TRUE, quote=FALSE, sep=",")
    
  }
  write.table((res.df %>% dplyr::filter(qtl==n.L)),
              file.path(sim_dir, "sim_results.csv"), col.names=TRUE, quote=FALSE, sep=",")
  write.table(getParams(), file.path(sim_dir, "params.txt"), col.names=FALSE, quote=FALSE, sep=":\t")
} # end qtl_vec
source("figures/aggregation/AggregateFigures.R")
