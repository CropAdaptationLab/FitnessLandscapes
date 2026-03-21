library(AlphaSimR)
library(dplyr)
library(ggpubr)
library(ggplot2)
library(ggpmisc)
library(grid)
library(patchwork)
library(tibble)

#rm(list = ls())

set.seed(123)

setwd("~/Documents/CSU/Lab-Notebooks/members/Ted/breeding_sims")
output_dir <- file.path(getwd(), "output/Geometric/2_27")
if (!dir.exists(output_dir)) dir.create(output_dir)

source("Functions/TraitArchitecture.R")
source("Scripts/GlobalParameters.R")

generateNewData <- TRUE

if (generateNewData) {
  nSims <- 10
  v_vec <- c(0.2, 0.4)
  q_vec <- c(1,2,5)
  
  init_df <- data.frame(rank=numeric(),
                       id=character(),
                       eff_size=numeric(),
                       qtl=numeric(),
                       var=numeric(),
                       scaled=numeric())
  
  for (qx in 1:length(q_vec)) {
    n.qtlPerChr <- q_vec[qx]
    print(paste0("QTL: ", n.qtlPerChr))
    for (vx in 1:length(v_vec)) {
      n.var <- v_vec[vx]
      print(paste0("Var: ", n.var))
      for (s in 1:nSims) {
        print(paste0("Sim: ", s))
        source("Scripts/CreateFounderPop.R")
        # Order by eff_size, and add 'rank' column
        eff_sizes_1 <- singleTraitArchitecture(founderPop, 1)
        alpha_1 <- eff_sizes_1[1,3]
        eff_sizes_1 <- eff_sizes_1 %>%
          dplyr::mutate(scaled=((1-n.R)/alpha_1)*eff_size)

        eff_sizes_2 <- singleTraitArchitecture(founderPop, 2)
        alpha_2 <- eff_sizes_2[1,3]
        eff_sizes_2 <- eff_sizes_2 %>%
          dplyr::mutate(scaled=((1-n.R)/alpha_2)*eff_size)

        eff_sizes <- rbind(eff_sizes_1,
                           eff_sizes_2)
        eff_sizes$qtl <- n.qtlPerChr
        eff_sizes$var <- n.var
        
        init_df <- rbind(init_df,
                        eff_sizes)
      }
    }
  }
} else {
  init_df <- read.csv(file.path(output_dir, "initial_architecture.csv"))
}

write.table(init_df, file.path(output_dir, "initial_architecture.csv"), col.names=TRUE, quote=FALSE, sep=",")
