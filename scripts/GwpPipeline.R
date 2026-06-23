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
library(purrr)
library(qtl)
library(RColorBrewer)
library(rrBLUP)
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
source("functions/GenoConversions.R")
source("functions/GenomewidePrediction.R")
source("functions/MappingPopulations.R")
source("functions/QtlMapping.R")
source("functions/QuantGen.R")
source("functions/TraitArchitecture.R")
source("scripts/GlobalParameters.R")

# Number of founder populations to simulate
n.popResets <- 2000
# Number of adaptive walk replications per pair of subpopulations
n.reps <- 1

# Recurrent selection years
n.Y <- 18

# Phenotype to use for genomic selection
GS_PHENO <- "pheno" # gv

n.trainPopSize <- 400

# Store the results of GWP from landrace into the RIL family
RIL.list <- list()

# Store the results of the recurrent selection 
RS.list <- list()

# All the parameter combinations to iterate through
model_vec <- c("RRBLUP")

output_dir <- file.path(base_dir, paste0("Sim_", format(Sys.time(), "%F_%H_%M")))
if (!dir.exists(output_dir)) dir.create(output_dir)

for (GS_MODEL in model_vec) {
  print(paste0("Model: ", GS_MODEL))
  
  # The directory for this number of QTL
  sim_dir <- file.path(output_dir, GS_MODEL)
  if (!dir.exists(sim_dir)) dir.create(sim_dir)
  
  # Reset the founder population n.popResets times
  for (f in 1238:n.popResets) {
    pop_dir <- file.path(sim_dir, paste0("FounderPopulation", f))
    if (!dir.exists(pop_dir)) dir.create(pop_dir)
    print(paste0("Founder Reset ", f))
    source("scripts/CreateFounderPop.R")

    # For each founder population, create independent populations n.sims times
    for (rep in 1:n.reps) {
      print(paste0("Rep ", rep))
      rep_dir <- file.path(pop_dir, paste0("Rep", rep))
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
      purelines2 <- makePurelinesBulk(pops[[2]])
      
      # ADMIXED
      # Create a biparental RIL population by sampling one individual from each subpopulation
      res_admixed <- createRIL(purelines1[[1]], purelines2[[1]])
      # Select the parents from the result, and the RIL population
      RIL_admixed <- res_admixed[-(1:2)]
      parent1 <- res_admixed[1]
      parent2 <- res_admixed[2]
      
      # POP1
      res_pop1 <- createRIL(purelines1[[1]], purelines1[[2]])
      RIL_pop1 <- res_pop1[-(1:2)]
      
      # POP2
      res_pop2 <- createRIL(purelines2[[1]], purelines2[[2]])
      RIL_pop2 <- res_pop2[-(1:2)]
      
      rAdmixed=evaluateGwpModel(trainPop=c(pops[[1]],pops[[2]]), testPop=RIL_admixed)
      rP1=evaluateGwpModel(trainPop=c(pops[[1]]), testPop=RIL_pop1)
      rP2=evaluateGwpModel(trainPop=c(pops[[2]]), testPop=RIL_pop2)
      rP1P2=evaluateGwpModel(trainPop=c(pops[[1]]), testPop=RIL_pop2)
      rP2P1=evaluateGwpModel(trainPop=c(pops[[2]]), testPop=RIL_pop1)
      rCvP1P2=crossValAccuracy(pop=c(pops[[1]],pops[[2]]))
      rCvP1=crossValAccuracy(pop=pops[[1]])
      rCvP2=crossValAccuracy(pop=pops[[2]])

      # Calculate isoeliteness for each trait
      isoElite_T1 <- isoEliteness(parent1, parent2, founderPop, 1)
      isoElite_T2 <- isoEliteness(parent1, parent2, founderPop, 2)
      isoEliteAtt <- mean(c(isoElite_T1, isoElite_T2))
      isoElite_T3 <- isoEliteness(parent1, parent2, founderPop, 3)
      
      

      # Store admixed, within, and cross-population prediction accuracies
      RIL.list[[length(RIL.list) + 1]] <- data.frame(
                        qtl=n.L,
                        founder=f,
                        rep=rep,
                        model=GS_MODEL,
                        isoElite=isoEliteAtt,
                        isoEliteDes=isoElite_T3,
                        rAdmixed=rAdmixed,
                        rP1=rP1,
                        rP2=rP2,
                        rP1P2=rP1P2,
                        rP2P1=rP2P1,
                        rCvP1P2=rCvP1P2,
                        rCvP1=rCvP1,
                        rCvP2=rCvP2
                      )
      
      # Run recurrent selection to improve the admixed RIL
      rs_result <- recurrentSelection(RIL_admixed, parent1, parent2)
      if (length(rs_result) > 0) {
        RS.list[[length(RS.list) + 1]] <- rs_result %>%
          dplyr::mutate(qtl=n.L,
                        founder=f,
                        rep=rep,
                        model=GS_MODEL,
                        type="Admixed",
                        isoElite=isoEliteAtt,
                        isoEliteDes=isoElite_T3,
                        .before=1)
      }

      # Just calculate isoeliteness for the P1 unadmixed family
      isoElite_T1 <- isoEliteness(res_pop1[1], res_pop1[2], founderPop, 1)
      isoElite_T2 <- isoEliteness(res_pop1[1], res_pop1[2], founderPop, 2)
      isoEliteAtt <- mean(c(isoElite_T1, isoElite_T2))
      isoElite_T3 <- isoEliteness(res_pop1[1], res_pop1[2], founderPop, 3)

      # Run recurrent selection to improve the unadmixed RIL
      # Only do this in even pop resets (because it's only being used as a negative control)
      # Only use P1
      #if (f %% 2 == 0) {
      rs_result <- recurrentSelection(RIL_pop1, res_pop1[1], res_pop1[2])
      if (length(rs_result) > 0) {
        RS.list[[length(RS.list) + 1]] <- rs_result %>%
          dplyr::mutate(qtl=n.L,
                        founder=f,
                        rep=rep,
                        model=GS_MODEL,
                        type="Unadmixed",
                        isoElite=isoEliteAtt,
                        isoEliteDes=isoElite_T3,
                        .before=1)
      }
      
    } # end n.reps
  } # end n.popResets
  RIL.df <- do.call(rbind, RIL.list)
  write.table((RIL.df %>% dplyr::filter(model == GS_MODEL)),
              file.path(sim_dir, "ril_results.csv"), col.names=TRUE, quote=FALSE, sep=",")
  RS.df <- do.call(rbind, RS.list)
  write.table((RS.df %>% dplyr::filter(model == GS_MODEL)),
              file.path(sim_dir, "rs_results.csv"), col.names=TRUE, quote=FALSE, sep=",")
  write.table(getParams(), file.path(sim_dir, "params.txt"), col.names=FALSE, quote=FALSE, sep=":\t")
} # end model_vec

source("figures/aggregation/GwpFigures.R")
