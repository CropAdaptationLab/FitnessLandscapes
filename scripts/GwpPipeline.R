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
n.popResets <- 150
# Number of adaptive walk replications per pair of subpopulations
n.reps <- 2

# Recurrent selection cycles
n.C <- 20

# Store the results of GWP from landrace into the RIL family
RIL.df <- data.frame(
  qtl=c(), # L
  founder=c(), # f
  rep=c(), # r
  isoElite=c(), # mean isoeliteness of parents
  rAdmixed=c(), # GWP accuracy in admixed family
  rP1=c(), # GWP accuracy within pop 1
  rP2=c(), # GWP accuracy within pop 2
  rP1P2=c(), # GWP cross-pop accuracy (pop1 to pop1)
  rP2P1=c() # GWP cross-pop accuracy (pop2 to pop1)
)

# Store the results of the recurrent selection 
RS.df <- data.frame(
  qtl=c(), # L
  founder=c(), # f
  rep=c(), # r
  type=c(), # Admixed / Unadmixed
  isoElite=c(), # mean isoeliteness of parents
  c=c(), # cycle
  sel=c(), # GARS or PRS
  w=c(), # breeding fitness
  r=c(), # GWP accuracy per cycle
  genome_het=c(), # genome-wide heterozygosity
  attained_het=c(), # attained trait heterozygosity
  desired_het=c(), # desired trait heterozygosity
  pop_ie=c() # population-level isoeliteness
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
      # Admixed GWP
      rAdmixed <- evaluateGWP_W(trainPop=c(pops[[1]],pops[[2]]), testPop=RIL_admixed)
      
      # POP1
      res_pop1 <- createRIL(purelines1[[1]], purelines1[[2]])
      RIL_pop1 <- res_pop1[-(1:2)]
      
      # POP2
      res_pop2 <- createRIL(purelines2[[1]], purelines2[[2]])
      RIL_pop2 <- res_pop2[-(1:2)]
      
      # Calculate isoeliteness for each trait
      isoElite_T1 <- isoEliteness(res_admixed[1], res_admixed[2], founderPop, 1)
      isoElite_T2 <- isoEliteness(res_admixed[1], res_admixed[2], founderPop, 2)
      isoEliteAdmixed <- mean(c(isoElite_T1, isoElite_T2))

      # Store admixed, within, and cross-population prediction accuracies
      RIL.df <- rbind(RIL.df,
                      data.frame(
                        qtl=n.L,
                        founder=f,
                        rep=rep,
                        isoElite=isoEliteAdmixed,
                        rAdmixed=evaluateGWP_W(trainPop=c(pops[[1]],pops[[2]]), testPop=RIL_admixed),
                        rP1=evaluateGWP_W(trainPop=c(pops[[1]]), testPop=RIL_pop1),
                        rP2=evaluateGWP_W(trainPop=c(pops[[2]]), testPop=RIL_pop2),
                        rP1P2=evaluateGWP_W(trainPop=c(pops[[1]]), testPop=RIL_pop2),
                        rP2P1=evaluateGWP_W(trainPop=c(pops[[2]]), testPop=RIL_pop1)
                      ))
      
      # Run recurrent selection to improve the admixed RIL
      newRow <- recurrentSelection(RIL_admixed) %>%
        dplyr::mutate(qtl=n.L,
                      founder=f,
                      rep=rep,
                      type="Admixed",
                      isoElite=isoEliteAdmixed,
                      .before=1)
      RS.df <- rbind(RS.df, newRow)
      
      
      isoElite_T1 <- isoEliteness(res_pop1[1], res_pop1[2], founderPop, 1)
      isoElite_T2 <- isoEliteness(res_pop1[1], res_pop1[2], founderPop, 2)
      isoEliteUnadmixed <- mean(c(isoElite_T1, isoElite_T2))

      # Run recurrent selection to improve the unadmixed RIL
      newRow <- recurrentSelection(RIL_pop1) %>%
        dplyr::mutate(qtl=n.L,
                      founder=f,
                      rep=rep,
                      type="Unadmixed",
                      isoElite=isoEliteUnadmixed,
                      .before=1)
      RS.df <- rbind(RS.df, newRow)
    } # end n.reps
  } # end n.popResets
  write.table((RIL.df %>% dplyr::filter(qtl==n.L)),
              file.path(sim_dir, "ril_results.csv"), col.names=TRUE, quote=FALSE, sep=",")
  write.table((RS.df %>% dplyr::filter(qtl==n.L)),
              file.path(sim_dir, "rs_results.csv"), col.names=TRUE, quote=FALSE, sep=",")
  write.table(getParams(), file.path(sim_dir, "params.txt"), col.names=FALSE, quote=FALSE, sep=":\t")
} # end qtl_vec

source("figures/aggregation/GwpFigures.R")
