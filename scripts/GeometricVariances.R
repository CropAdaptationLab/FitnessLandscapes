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

setwd("~/Documents/CSU/Lab-Notebooks/members/Ted/breeding_sims")
output_dir <- file.path(getwd(), "output/TraitArchitecture/geometric")
if (!dir.exists(output_dir)) dir.create(output_dir)

n.cores <- 8
source("Functions/TraitArchitecture.R")
source("Scripts/GlobalParameters.R")

SAMPLING <- 'geometric'
n.popSize <- 1000

nSims <- 10
r_vec <- c(0.7,0.9)
a_vec <- c(0.4,0.6,0.8)

all_pop_vals <- data.frame(A=numeric(),
                           R=numeric(),
                           var=numeric(),
                           mean=numeric())

for (ax in 1:length(a_vec)) {
  n.A <- a_vec[ax]
  print(paste0("A: ", n.A))
  for (rx in 1:length(r_vec)) {
    n.R <- r_vec[rx]
    print(paste0("R: ", n.R))
    for (s in 1:nSims) {
      print(paste0("Sim: ", s))
      source("Scripts/CreateFounderPop.R")
      all_pop_vals <- rbind(all_pop_vals, data.frame(A=c(n.A, n.A),
                                                     R=c(n.R, n.R),
                                                     var=c(varG(founderPop)[1,1], varG(founderPop)[2,2]),
                                                     mean=c(meanG(founderPop)[1], meanG(founderPop)[2])))
    }
  }
}

all_pop_vals$R <- as.factor(all_pop_vals$R)
all_pop_vals$A <- as.factor(all_pop_vals$A)


p1 <- ggplot(all_pop_vals, aes(x = factor(A), y = var, fill = factor(A))) +
  geom_boxplot() +
  labs(title = "Var by A", x = "Category A", y = "Value") +
  theme_minimal() +
  theme(legend.position = "none")

# 2. mean across A
p2 <- ggplot(all_pop_vals, aes(x = factor(A), y = mean, fill = factor(A))) +
  geom_boxplot() +
  labs(title = "Mean by A", x = "Category A", y = "Value") +
  theme_minimal() +
  theme(legend.position = "none")

# 3. var across R
p3 <- ggplot(all_pop_vals, aes(x = factor(R), y = var, fill = factor(R))) +
  geom_boxplot() +
  labs(title = "Var by R", x = "Category R", y = "Value") +
  theme_minimal() +
  theme(legend.position = "none")

# 4. mean across R
p4 <- ggplot(all_pop_vals, aes(x = factor(R), y = mean, fill = factor(R))) +
  geom_boxplot() +
  labs(title = "Mean by R", x = "Category R", y = "Value") +
  theme_minimal() +
  theme(legend.position = "none")

# Display them all together
(p1 | p2) / (p3 | p4)
fname <- file.path(output_dir, "vars_and_means.jpg")
ggplot2::ggsave(filename = fname,
                device = "jpg",
                height=6,
                width=6,
                units="in")

write.table(all_pop_vals, file.path(output_dir, "geometric_data.csv"), col.names=TRUE, quote=FALSE, sep=",")
