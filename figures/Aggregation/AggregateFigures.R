# Title: AGGREGATE FIGURES
# Author: Ted Monyak
# Description: This is called at the end of SimulationPipeline.R to aggregate
# the reults into a multi-panel, but can also be overridden by manually specifying paths
# to datasets

library(dplyr)
library(ggplot2)
library(ggpubr)
library(grid)
library(patchwork)

# Manually override this to load old data
if (FALSE) {
  setwd("~/Documents/CSU/FitnessLandscapes/output/AggregatedResults")
  output_dir <- file.path(getwd(), "Sim_manuscript2")
  if (!dir.exists(output_dir)) dir.create(output_dir)
  
  res.df <- readRDS("Sim_manuscript/figure_data/sim_results.rds")
  effectSizes.df <- readRDS("Sim_manuscript/figure_data/effect_sizes.rds")
  fixedAlleles.df <- readRDS("Sim_manuscript/figure_data/fixed_alleles.rds")
}

res.df$qtl <- as.factor(res.df$qtl)
res.df$var <- as.factor(res.df$var)
res.df$type <- factor(res.df$type)

# Post-processing of data for plotting and tables
res.df <- res.df %>%
  dplyr::mutate(
    relRank_T1 = as.numeric(log(rilA_T1 / initA_T1) / log(a)), # Solve for k from the founder architecture
    relRank_T2 = as.numeric(log(rilA_T2 / initA_T2) / log(a)),
    emergentLod_W = nLod_Int-nLod_W, # The difference between interaction LOD peaks and single LOD peaks
    hamm_Att=hamm_T1+hamm_T2,
  ) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    relRank_T1 = replace(relRank_T1, is.infinite(relRank_T1), NA),
    relRank_T2 = replace(relRank_T1, is.infinite(relRank_T2), NA),
  ) %>%
  dplyr::mutate(
    isoElite_Att = mean(c(isoElite_T1, isoElite_T2), na.rm=TRUE),
    attR = mean(c(gwpR_T1,gwpR_T2), na.rm=TRUE),
    attR_Pop2 = mean(c(gwpR_T1_Pop2, gwpR_T2_Pop2), na.rm=TRUE),
    ev_Att = mean(c(ev_T1,ev_T2), na.rm=TRUE),
    nLod_Att = mean(c(nLod_T1,nLod_T2), na.rm=TRUE),
    relRankMean = mean(c(relRank_T1, relRank_T2), na.rm=TRUE),
    rilA = mean(c(rilA_T1,rilA_T2), na.rm=TRUE),
    rilVar = mean(c(rilVar_T1,rilVar_T2), na.rm=TRUE),
    initA = mean(c(initA_T1,initA_T2), na.rm=TRUE),
    initVar = mean(c(initVar_T1,initVar_T2), na.rm=TRUE),
  )

# Common theme elements for each plot
scale_fill <- scale_fill_manual(name = "RIL Family",
                  values = c("Admixed" = "gold2",
                             "Unadmixed" = "#FFFDD9"),
                  labels = c("Admixed", "Unadmixed"))

scale_color <- scale_color_manual(name = "QTL per\nAttained Trait",
                                  values = c("10" = "#4A1A6B",
                                             "20" = "#9B59B6",
                                             "50" = "#D7B8F3"))

scale_fill_qtl <- scale_fill_manual(name = "QTL per\nAttained Trait",
                                  values = c("10" = "#4A1A6B",
                                             "20" = "#9B59B6",
                                             "50" = "#D7B8F3"),
                                  guide="none")

# The maximum allelic effect size, to have a consistent y-axis
a <- 0.6

# Generate all plots
setwd("~/Documents/CSU/FitnessLandscapes/figures/aggregation")
source("Isoeliteness.R")
source("ExcessVariance.R")
source("LODPeaks.R")
source("AllogenicRank.R")
source("GWP.R")
source("DataTables.R")
source("AlleleFixationOrder.R")
source("TraitArchitectureRIL.R")

