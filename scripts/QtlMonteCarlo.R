# Title: QTL Monte Carlo
# Author Ted Monyak
# Description:
# This will run n.popResets * n.sims simulations, and for each simulation,
# create two subpopulations, then two biparental RIL populations:
# 1 "admixed" (parents come from different subpopulations),
# 1 "unadmixed" (parents come from the same subpopulation)
# For each RIL, calculate and store the number of significant LOD peaks in the linkage map
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
library(rrBLUP)
library(tibble)
library(tidyr)
library(viridis)

rm(list = ls())

set.seed(123)

# Set to false if running on Alpine
runLocal = TRUE

if (runLocal) {
  setwd("~/Documents/CSU/Lab-Notebooks/members/Ted/breeding_sims")
  output_dir <- file.path(getwd(), "output")
  n.cores <- 8
} else {
  setwd("/pl/active/Morris_CSU/Ted_Monyak/BreedingSims")
  output_dir <- ("/scratch/alpine/c837220672@colostate.edu/output")
  n.cores <- 4
}

if (!dir.exists(output_dir)) dir.create(output_dir)

source("Functions/Epistasis.R")
source("Functions/Fitness.R")
source("Functions/GenoConversions.R")
source("Functions/MappingPopulations.R")
source("Functions/PopGen.R")
source("Functions/QtlMapping.R")
source("Functions/QuantGen.R")
source("Functions/TraitArchitecture.R")
source("Plotting/LinkageMaps.R")
source("Plotting/PCA.R")
source("Plotting/ReactionNorm.R")
source("Scripts/GlobalParameters.R")

# Number of founder populations to simulate
n.popResets <- 1
# Number of adaptive walk simulations per pair of subpopulations
n.sims <- 1

SAMPLING <- 'geometric'
n.minMAF <- 0.3

# FUNCTIONS
qtlMapping <- TRUE
twoQtlMapping <- TRUE
genomicPrediction <- FALSE
compareUnadmixed <- FALSE

# KEEP THESE TRUE
saveFixationOrder <- FALSE
saveEffectSizes <- FALSE

# PLOTTING
saveQtlPlots <- TRUE
saveTraitPlots <- TRUE
saveAllelePlots <- TRUE
saveFitnessPlots <- TRUE

# Set to true for doing fitness landscape creation
sampleInds <- FALSE

# All the parameter combinations to iterate through
qtl_vec <- c(1)
subpop_vec = c(1000)
var_vec <- c(0.4)

for (qx in 1:length(qtl_vec)) {
  n.qtlPerChr <- qtl_vec[qx]
  for (px in 1:length(subpop_vec)){
    n.subPopSize <- subpop_vec[px]
    for (vx in 1:length(var_vec)) {
      n.var <- var_vec[vx]
      print(paste0("QTL: ", n.qtlPerChr, " | POP: ", n.subPopSize, " | v: ", n.var))
      
      # Result dataframe
      res.df <- data.frame(var=c(), a=c(), popSize=c(), qtl=c(), pop=c(), sim=c(), type=c(), fst=c(),
                           nLod_T1=c(), nLod_T2=c(), nLod_T3=c(), nLod_Suit=c(), nLod_W=c(),
                           nLod_Int=c(),
                           ev_T1=c(), ev_T2=c(), ev_T3=c(), ev_Suit=c(), ev_W=c(),
                           isoElite_T1=c(), isoElite_T2=c(), isoElite_T3=c(),
                           hamm_T1=c(), hamm_T2=c(),
                           gwpR_T1=c(), gwpR_T2=c(), gwpR_W=c(),
                           gwpR_T1_Pop2=c(), gwpR_T2_Pop2=c(), gwpR_W_Pop2=c(),
                           initA_T1=c(), initA_T2=c(), initVar_T1=c(), initVar_T2=c(),
                           rilA_T1=c(), rilA_T2=c(), rilVar_T1=c(), rilVar_T2=c())

      # Effect size dataframe
      effectSizes.df <- data.frame(
        id=c(),
        eff_size=c(),
        rank=c(),
        var=c(),
        qtl=c())
      
      # Tidy dataframe to store the order in which each allele is fixed, and the effect size
      fixedAlleles.df <- data.frame(
        id=c(),
        eff_size=c(),
        order_fixed=c(),
        subpop=c(),
        var=c(),
        qtl=c())
      
      base_dir <- file.path(output_dir, "QtlMonteCarlo")
      if (!dir.exists(base_dir)) dir.create(base_dir)
      base_fname <- paste0(paste0("qtl_", n.qtlPerChr, "_Ne_", n.subPopSize, "_var_", n.var, "_"))
      base_dir <- file.path(base_dir, paste0(base_fname, format(Sys.time(), "%F_%H_%M")))
      if (!dir.exists(base_dir)) dir.create(base_dir)
      
      # Reset the founder population n.popResets times
      for (r in 1:n.popResets) {
        # Update selProp to be a random value if using random parameters
        pop_dir <- file.path(base_dir, paste0("FounderPopulation", r))
        if (!dir.exists(pop_dir)) dir.create(pop_dir)
        print(paste0("Pop Reset ", r))
        source("Scripts/CreateFounderPop.R")
        # Get the largest effect size of the  traits
        t1_arch <- singleTraitArchitecture(founderPop,1)$eff_size
        t2_arch <- singleTraitArchitecture(founderPop,2)$eff_size
        initA_T1 <- max(t1_arch)
        initA_T2 <- max(t2_arch)
        initVar_T1 <- var(t1_arch)
        initVar_T2 <- var(t2_arch)
        
        # For each founder population, create independent populations n.sims times
        for (s in 1:n.sims) {
          print(paste0("Sim ", s))
          sim_dir <- file.path(pop_dir, paste0("Sim", s))
          if (!dir.exists(sim_dir)) dir.create(sim_dir)
          
          fitFunc <- gaussianLandraceFitFunc
          fitCalc <- calculateFitnessGaussian
          
          # Landrace heritability
          SP$setVarE(h2=c(n.h2, n.h2, n.yieldH2))
          source("Scripts/CreateIndependentPops.R")
          # Update to use the breeding fitness function
          fitFunc <- gaussianFitFunc
          SP$setVarE(h2=c(n.h2Breeding, n.h2Breeding, n.yieldH2Breeding))
          
          save_dir <- file.path(sim_dir, "Admixed")
          if (!dir.exists(save_dir)) dir.create(save_dir)

          # Develop purelines
          fst <- FST(pops)
          purelines1 <- makePurelinesBulk(pops[[1]])
          pop1 <- purelines1[[1]]
          purelines2 <- makePurelinesBulk(pops[[2]])
          pop2 <- purelines2[[1]]
          
          # Create a biparental RIL population by sampling one individual from each subpopulation
          n.RILFams <- 250
          res <- createRIL(pop1, pop2, save_dir)
          # Select the parents from the result, and the RIL population
          RIL <- res[-(1:2)]
          parent1 <- res[1]
          parent2 <- res[2]
          # Get effect sizes of RIL
          
          effSizes_t1 <- sortedEffectSizes(RIL, c(1))
          if (nrow(effSizes_t1) > 0) {
            effSizes_t1$var <- n.var
            effSizes_t1$qtl <- n.L
            effectSizes.df <- rbind(effectSizes.df, effSizes_t1)
            rilA_T1 <- max(effSizes_t1$eff_size)
            rilVar_T1 <- var(effSizes_t1$eff_size)
          } else {
            rilA_T1 <- 0
            rilVar_T1 <- 0
          }

          effSizes_t2 <- sortedEffectSizes(RIL, c(2))
          if (nrow(effSizes_t2) > 0) {
            effSizes_t2$var <- n.var
            effSizes_t2$qtl <- n.L
            effectSizes.df <- rbind(effectSizes.df, effSizes_t2)
            rilA_T2 <- max(effSizes_t2$eff_size)
            rilVar_T2 <- var(effSizes_t2$eff_size)
          } else {
            rilA_T2 <- 0
            rilVar_T2 <- 0
          }

          isoElite_T1 <- isoEliteness(parent1, parent2, founderPop, 1)
          isoElite_T2 <- isoEliteness(parent1, parent2, founderPop, 2)
          isoElite_T3 <- isoEliteness(parent1, parent2, founderPop, 3)
          
          hamm_T1 <- hammingDistance(parent1, parent2, 1)
          hamm_T2 <- hammingDistance(parent1, parent2, 2)
          source("Plotting/GenerateSubpopPlots.R")
          
          # Calculate the ratio of RIL variance to pureline variance for each of the traits
          RIL_pheno <- as.data.frame(pheno(RIL)) %>%
            dplyr::mutate(Suitability=fitCalc(Trait1, Trait2),
                          W=calculateBreedingFitness(Trait1, Trait2, Trait3),
                          Population="RIL Family")
          pureline1_pheno <- as.data.frame(pheno(pop1)) %>%
            dplyr::mutate(Suitability=fitCalc(Trait1, Trait2),
                          W=calculateBreedingFitness(Trait1, Trait2, Trait3),
                          Population="Pureline 1")
          pureline2_pheno <- as.data.frame(pheno(pop2)) %>%
            dplyr::mutate(Suitability=fitCalc(Trait1, Trait2),
                          W=calculateBreedingFitness(Trait1, Trait2, Trait3),
                          Population="Pureline 2")
  
          pureline_pheno <- rbind(pureline1_pheno, pureline2_pheno)

          ev_T1 <- excessVariance(RIL_pheno[,1], pureline_pheno[,1])
          ev_T2 <- excessVariance(RIL_pheno[,2], pureline_pheno[,2])
          ev_T3 <- excessVariance(RIL_pheno[,3], pureline_pheno[,3])
          ev_Suit <- excessVariance(RIL_pheno[,4], pureline_pheno[,4])
          ev_W <- excessVariance(RIL_pheno[,5], pureline_pheno[,5])
          source("Plotting/TransgressiveSegregation.R")

          # Determine the number of significant QTL
          if (qtlMapping) {
            nLOD <- getLodPeaks(RIL, parent1, parent2, save_dir)
          } else {
            nLOD <- c(0,0,0,0,0)
          }
          if (twoQtlMapping) {
            # CHANGE THIS FUNCTION SO THAT IT GETS MARKERS NOT QTLS
            peaks <- epistaticLodPeaks(RIL, parent1, parent2, trait=5, save_dir)
            if (nrow(peaks) > 0) {
              maxLodIdx <- which.max(peaks$lod.int)
              qtl1 <- peaks$qtl1[maxLodIdx]
              qtl2 <- peaks$qtl2[maxLodIdx]
              if ((qtl1 != qtl2) & (saveQtlPlots)) {
                plot_reaction_norm(RIL, qtl1, qtl2, parent1, parent2, calculateBreedingFitness, save_dir)
                plot_1D_PCA(RIL, pops[[1]], pops[[2]], parent1, parent2, qtl1, qtl2, save_dir)
              }
            }
            nIntPeaks <- nrow(peaks)
          } else {
            nIntPeaks <- 0
          }
          
          if (genomicPrediction) {
            gwpR_T1 <- evaluateGWP(trainPop=c(pops[[1]],pops[[2]]), testPop=RIL, trait=1)
            gwpR_T2 <- evaluateGWP(trainPop=c(pops[[1]],pops[[2]]), testPop=RIL, trait=2)
            gwpR_W <- evaluateGWP_W(trainPop=c(pops[[1]],pops[[2]]), testPop=RIL)
          } else {
            gwpR_T1 <- 0
            gwpR_T2 <- 0
            gwpR_W <- 0
          }

          res.df <- rbind(res.df, data.frame(var=n.var, a=n.a, popSize=n.subPopSize, qtl=n.L, pop=r, sim=s, type="Admixed", fst=fst,
                                             nLod_T1=nLOD[1], nLod_T2=nLOD[2], nLod_T3=nLOD[3], nLod_Suit=nLOD[4],
                                             nLod_W=nLOD[5], nLod_Int=nIntPeaks,
                                             ev_T1=ev_T1, ev_T2=ev_T2, ev_T3=ev_T3, ev_Suit=ev_Suit, ev_W=ev_W,
                                             isoElite_T1=isoElite_T1, isoElite_T2=isoElite_T2, isoElite_T3=isoElite_T3,
                                             hamm_T1=hamm_T1, hamm_T2=hamm_T2,
                                             gwpR_T1=gwpR_T1, gwpR_T2=gwpR_T2, gwpR_W=gwpR_W,
                                             gwpR_T1_Pop2=NA, gwpR_T2_Pop2=NA, gwpR_W_Pop2=NA,
                                             initA_T1=initA_T1,
                                             initA_T2=initA_T2,
                                             initVar_T1=initVar_T1,
                                             initVar_T2=initVar_T2,
                                             rilA_T1=rilA_T1,
                                             rilA_T2=rilA_T2,
                                             rilVar_T1=rilVar_T1,
                                             rilVar_T2=rilVar_T2))
          
          if (compareUnadmixed) {
            # Create a biparental RIL population by sampling two individuals from the same subpopulation
            pop1 <- purelines1[[1]]
            pop2 <- purelines1[[2]]
            res <- createRIL(pop1, pop2, save_dir)
            RIL <- res[-(1:2)]
            parent1 <- res[1]
            parent2 <- res[2]
            isoElite_T1 <- isoEliteness(parent1, parent2, founderPop, 1)
            isoElite_T2 <- isoEliteness(parent1, parent2, founderPop, 2)
            isoElite_T3 <- isoEliteness(parent1, parent2, founderPop, 3)
            hamm_T1 <- hammingDistance(parent1, parent2, 1)
            hamm_T2 <- hammingDistance(parent1, parent2, 2)
            if (qtlMapping) {
              nLOD <- getLodPeaks(RIL, parent1, parent2, save_dir)
            } else {
              nLOD <- c(0,0,0,0,0)
            }
            
            if (twoQtlMapping) {
              peaks <- epistaticLodPeaks(RIL, parent1, parent2, trait=4, save_dir)
              nIntPeaks <- nrow(peaks)
            } else {
              nIntPeaks <- 0
            }
            RIL_pheno <- as.data.frame(pheno(RIL)) %>%
              dplyr::mutate(Fitness=fitCalc(Trait1, Trait2),
                            W=calculateBreedingFitness(Trait1, Trait2, Trait3))
            pureline_pheno <- as.data.frame(pheno(c(pop1, pop2))) %>%
              dplyr::mutate(Fitness=fitCalc(Trait1, Trait2),
                            W=calculateBreedingFitness(Trait1, Trait2, Trait3))
            
            ev_T1 <- excessVariance(RIL_pheno[,1], pureline_pheno[,1])
            ev_T2 <- excessVariance(RIL_pheno[,2], pureline_pheno[,2])
            ev_T3 <- excessVariance(RIL_pheno[,3], pureline_pheno[,3])
            ev_Suit <- excessVariance(RIL_pheno[,4], pureline_pheno[,4])
            ev_W <- excessVariance(RIL_pheno[,5], pureline_pheno[,5])
            if (genomicPrediction) {
              gwpR_T1 <- evaluateGWP(trainPop=pops[[1]], testPop=RIL, trait=1)
              gwpR_T2 <- evaluateGWP(trainPop=pops[[1]], testPop=RIL, trait=2)
              gwpR_W <- evaluateGWP_W(trainPop=pops[[1]], testPop=RIL)
              
              # Do an intra RIL for population 2
              pop1 <- purelines2[[1]]
              pop2 <- purelines2[[2]]
              res <- createRIL(pop1, pop2, save_dir)
              RIL <- res[-(1:2)]
              gwpR_T1_Pop2 <- evaluateGWP(trainPop=pops[[2]], testPop=RIL, trait=1)
              gwpR_T2_Pop2 <- evaluateGWP(trainPop=pops[[2]], testPop=RIL, trait=2)
              gwpR_W_Pop2 <- evaluateGWP_W(trainPop=pops[[2]], testPop=RIL)
            } else {
              gwpR_T1 <- 0
              gwpR_T2 <- 0
              gwpR_W <- 0
              gwpR_T1_Pop2 <- 0
              gwpR_T2_Pop2 <- 0
              gwpR_W_Pop2 <- 0
            }
            # Update the result dataframe
            res.df <- rbind(res.df, data.frame(var=n.var, a=n.a, popSize=n.subPopSize, qtl=n.L, pop=r, sim=s, type="Unadmixed", fst=NA,
                                               nLod_T1=nLOD[1], nLod_T2=nLOD[2], nLod_T3=nLOD[3], nLod_Suit=nLOD[4],
                                               nLod_W=nLOD[5], nLod_Int=nIntPeaks,
                                               ev_T1=ev_T1, ev_T2=ev_T2, ev_T3=ev_T3, ev_Suit=ev_Suit, ev_W=ev_W,
                                               isoElite_T1=isoElite_T1, isoElite_T2=isoElite_T2, isoElite_T3=isoElite_T3,
                                               hamm_T1=hamm_T1, hamm_T2=hamm_T2,
                                               gwpR_T1=gwpR_T1, gwpR_T2=gwpR_T2, gwpR_W=gwpR_W,
                                               gwpR_T1_Pop2=gwpR_T1_Pop2, gwpR_T2_Pop2=gwpR_T2_Pop2, gwpR_W_Pop2=gwpR_W_Pop2,
                                               initA_T1=NA, initA_T2=NA, initVar_T1=NA, initVar_T2=NA,
                                               rilA_T1=NA, rilA_T2=NA, rilVar_T1=NA, rilVar_T2=NA))
          } # end compareIntra
        } # end n.sims
      } # end n.popResets
      #source("Plotting/Aggregation/Isoeliteness.R")
      #source("Plotting/Aggregation/ExcessVariance.R")
      #source("Plotting/Aggregation/AllogenicRank.R")
      #source("Plotting/Aggregation/GWP.R")
      #source("Plotting/Aggregation/LODPeaks.R")
      source("Plotting/GenerateSeriesPlots.R")
      write.table(res.df, file.path(base_dir, "sim_results.csv"), col.names=TRUE, quote=FALSE, sep=",")
      write.table(getParams(), file.path(base_dir, "params.txt"), col.names=FALSE, quote=FALSE, sep=":\t")
    } 
  } # end pop_vec
} # end qtl_vec
