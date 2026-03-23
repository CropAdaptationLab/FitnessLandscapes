# Title: TRAIT ARCHITECTURE FOUNDER
# Author: Ted Monyak
# Description: This script will simulate several founder populations and
# determine the mean series of allele substituion effect sizes

library(AlphaSimR)
library(dplyr)
library(ggpubr)
library(ggplot2)
library(ggpmisc)
library(grid)
library(patchwork)
library(tibble)

rm(list = ls())

set.seed(123)

setwd("~/Documents/CSU/FitnessLandscapes")
output_dir <- file.path(getwd(), "output/Geometric/3_15")
if (!dir.exists(output_dir)) dir.create(output_dir)

source("functions/TraitArchitecture.R")
source("scripts/GlobalParameters.R")

# Number of simulations per combination of parameters
n.sims <- 1

# QTL per chromosome
q_vec <- c(1,2,5)

# To store the output data
init_df <- data.frame(rank=numeric(),
                   id=character(),
                   eff_size=numeric(),
                   qtl=numeric(),
                   scaled=numeric())

# Iterate through each QTL per chromosome
for (qx in 1:length(q_vec)) {
  n.qtlPerChr <- q_vec[qx]
  print(paste0("QTL: ", n.qtlPerChr))
  for (s in 1:n.sims) {
    print(paste0("Sim: ", s))
    source("Scripts/CreateFounderPop.R")
    # Order by eff_size, and add 'rank' column
    eff_sizes_1 <- singleTraitArchitecture(founderPop, 1)
    alpha_1 <- eff_sizes_1[1,3]
    # Approximation of the percent variance explained by each QTL
    # Based on Lande and Thompson 1990
    eff_sizes_1 <- eff_sizes_1 %>%
      dplyr::mutate(scaled=((1-n.a)/alpha_1)*eff_size)

    eff_sizes_2 <- singleTraitArchitecture(founderPop, 2)
    alpha_2 <- eff_sizes_2[1,3]
    
    # Approximation of the percent variance explained by each QTL
    eff_sizes_2 <- eff_sizes_2 %>%
      dplyr::mutate(scaled=((1-n.a)/alpha_2)*eff_size)

    eff_sizes <- rbind(eff_sizes_1,
                       eff_sizes_2)
    eff_sizes$qtl <- n.qtlPerChr
    
    init_df <- rbind(init_df,
                    eff_sizes)
  } # end n.sims
} # end q_vec

write.table(init_df, file.path(output_dir, "initial_architecture.csv"), col.names=TRUE, quote=FALSE, sep=",")
