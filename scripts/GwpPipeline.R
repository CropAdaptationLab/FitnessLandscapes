# Title: GWP PIPELINE
# Author: Ted Monyak
# Description:
# This will run n.popResets * n.reps simulations, and for each replication,
# create two subpopulations, then two biparental RIL populations:
# 1 "admixed" (parents come from different subpopulations),
# 1 "unadmixed" (parents come from the same subpopulation)
# This script is similar to SimulationPipeline.R, but is focused on
# genomewide prediction

library(AlphaSimR)
library(devtools)
library(dplyr)
library(ggpubr)
library(ggplot2)
library(ggpmisc)
library(patchwork)
library(qtl2convert)
library(RColorBrewer)
library(reshape2)
library(tibble)
library(tidyr)

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
base_dir <- file.path(getwd(), "output/GWP")
if (!dir.exists(base_dir)) dir.create(base_dir)

source("functions/Fitness.R")
source("functions/GenomewidePrediction.R")
source("functions/MappingPopulations.R")
source("functions/QuantGen.R")
source("functions/TraitArchitecture.R")
source("scripts/GlobalParameters.R")

# Number of founder populations to simulate
n.popResets <- 100
# Number of adaptive walk replications per pair of subpopulations
n.reps <- 2

n.C <- 20

compareUnadmixed = TRUE

# Store the results of the GARS
res.df <- data.frame(
  qtl=c(), # L
  founder=c(), # f
  rep=c(), # r
  type=c(), # Admixed / Unadmixed
  isoElite=c(), # mean isoeliteness of parents
  gwpR=c(), # accuracy of genomic prediction
  gwpR_P2=c(), # (only for unadmixed population) GWP accuracy on second subpop
  c=c(), # cycle
  sel=c(), # GARS or PRS
  w=c() # breeding fitness
)

# All the parameter combinations to iterate through
qtl_vec <- c(1,2,5)

output_dir <- file.path(base_dir, paste0("Sim_", format(Sys.time(), "%F_%H_%M")))
if (!dir.exists(output_dir)) dir.create(output_dir)

for (qx in 1:length(qtl_vec)) {
  n.qtlPerChr <- qtl_vec[qx]
  print(paste0("QTL: ", n.qtlPerChr))
  
  # The directory for this number of QTL
  sim_dir <- file.path(output_dir, paste0("QTL_", n.qtlPerChr))
  if (!dir.exists(sim_dir)) dir.create(sim_dir)
  
  # Reset the founder population n.popResets times
  for (f in 1:n.popResets) {
    pop_dir <- file.path(sim_dir, paste0("FounderPopulation", f))
    if (!dir.exists(pop_dir)) dir.create(pop_dir)
    print(paste0("Founder Reset ", f))
    source("scripts/CreateFounderPop.R")

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
      
      # Update to use the suitability selection function
      selectionPheno <- suitability
      SP$setVarE(h2=c(n.h2Breeding, n.h2Breeding, n.yieldH2Breeding))

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
  
      # Calculate isoeliteness for each trait
      isoElite_T1 <- isoEliteness(parent1, parent2, founderPop, 1)
      isoElite_T2 <- isoEliteness(parent1, parent2, founderPop, 2)
      isoElite <- mean(c(isoElite_T1, isoElite_T2))
      
      # RRBLUP
      gwpR <- evaluateGWP_W(trainPop=c(pops[[1]],pops[[2]]), testPop=RIL)
      
      # Run recurrent selection to improve the RIL
      # Run recurrent selection to improve the RIL
      newRow <- recurrentSelection(RIL) %>%
        dplyr::mutate(qtl=n.L,
                      founder=f,
                      rep=r,
                      type="Admixed",
                      isoElite=isoElite,
                      gwpR=gwpR,
                      gwpR_P2=NA,
                      .before=1)
      res.df <- rbind(res.df, newRow)
      
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
        isoElite <- mean(c(isoElite_T1, isoElite_T2))
        
        gwpR <- evaluateGWP_W(trainPop=pops[[1]], testPop=RIL)
        
        # Do an unadmixed RIL for population 2
        pop1 <- purelines2[[1]]
        pop2 <- purelines2[[2]]
        res <- createRIL(pop1, pop2)
        RIL <- res[-(1:2)]
        gwpR_P2 <- evaluateGWP_W(trainPop=pops[[2]], testPop=RIL)
        
        # Run recurrent selection to improve the RIL
        newRow <- recurrentSelection(RIL) %>%
          dplyr::mutate(qtl=n.L,
                        founder=f,
                        rep=r,
                        type="Unadmixed",
                        isoElite=isoElite,
                        gwpR=gwpR,
                        gwpR_P2=gwpR_P2,
                        .before=1)
        res.df <- rbind(res.df, newRow)
      } # end compareIntra
    } # end n.reps
  } # end n.popResets
  write.table((res.df %>% dplyr::filter(qtl==n.L)),
              file.path(sim_dir, "sim_results.csv"), col.names=TRUE, quote=FALSE, sep=",")
  write.table(getParams(), file.path(sim_dir, "params.txt"), col.names=FALSE, quote=FALSE, sep=":\t")
} # end qtl_vec

source("figures/aggregation/GwpFigures.R")
