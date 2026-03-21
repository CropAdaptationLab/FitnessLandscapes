library(dplyr)
library(ggplot2)
library(ggpubr)
library(grid)
library(patchwork)

setwd("~/Documents/CSU/Lab-Notebooks/members/Ted/breeding_sims/output/QtlMonteCarlo/figure_data")
output_dir <- file.path(getwd(), "../../AggregateQtl/3_15")
if (!dir.exists(output_dir)) dir.create(output_dir)

# 400 sims
qtl1 <- read.csv("qtl_1_Ne_1000_var_0.4_2026-03-14_00_04/sim_results.csv")
qtl2 <- read.csv("qtl_2_Ne_1000_var_0.4_2026-03-14_19_10/sim_results.csv")
qtl5 <- read.csv("qtl_5_Ne_1000_var_0.4_2026-03-15_05_51/sim_results.csv")

res.df <- rbind(qtl1,
                qtl2,
                qtl5)

res.df$qtl <- as.factor(res.df$qtl)
res.df$var <- as.factor(res.df$var)
res.df$type <- factor(res.df$type)

## DATA ANALYSIS FOR TABLES
res.df <- res.df %>%
  dplyr::mutate(
    relRank_T1 = as.numeric(log(rilA_T1 / initA_T1) / log(a)),
    relRank_T2 = as.numeric(log(rilA_T2 / initA_T2) / log(a)),
    emergentLod_Suit = nLod_Int-nLod_Suit,
    emergentLod_W = nLod_Int-nLod_W,
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

a <- 0.6

setwd("~/Documents/CSU/Lab-Notebooks/members/Ted/breeding_sims/Plotting/Aggregation")

source("Isoeliteness.R")
source("ExcessVariance.R")
source("LODPeaks.R")
#source("AllogenicRank.R")
#source("GWP.R")
source("DataTables.R")

