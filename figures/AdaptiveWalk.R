# Title: ADAPTIVE WALK
# Author: Ted Monyak
# Description: Creates a plot of the adaptive walk of a single population to a fitness peak

library(AlphaSimR)
library(dplyr)
library(ggplot2)
library(grDevices)
library(plotly)
library(RColorBrewer)
library(tibble)
library(tidyr)
library(viridis)

rm(list = ls())

setwd("~/Documents/CSU/FitnessLandscapes")
output_dir <- file.path(getwd(), "output/AdaptiveWalks")
n.cores <- 8

source("figures/T2FLandscapes.R")
source("functions/Fitness.R")
source("scripts/GlobalParameters.R")

saveFitnessPlots <- FALSE
suitFunc <- suitabilityGaussian

source("scripts/CreateFounderPop.R")

# Run for fewer generations to create a cleaner plot
n.gens <- 50

# For storing the results at each generation
res.df <- data.frame(Trait1=c(),
                     Trait2=c(),
                     Suit=c(),
                     Gen=c())


# Select 1000 random individuals
pop <- selectInd(founderPop, nInd=1000, use="rand")
large_res.df <- data.frame(Trait1=c(),
                         Trait2=c(),
                         Suit=c(),
                         Gen=c())

# Specify the "middle" window of the large population from which to
# sample individuals, based on the two attained traits
# This is to make the adaptive walk look less jittery, by selecting
# the individuals most representative of the population
window <- 100

# Select a "medium" sized population of 200 individuals
medPop <- selectInd(founderPop, nInd=500, use="rand")
med_res.df <- data.frame(Trait1=c(),
                         Trait2=c(),
                         Suit=c(),
                         Gen=c())

# Select a "small" sized population of 200 individuals
smallPop <- selectInd(founderPop, nInd=200, use="rand")
small_res.df <- data.frame(Trait1=c(),
                           Trait2=c(),
                           Suit=c(),
                           Gen=c())

for (gen in 1:n.gens) {
  pheno <- as.data.frame(pheno(pop)) %>%
    rownames_to_column("idx")
  
  # The midpoint of the number of individuals in the population
  mid <- nInd(pop)/2
  
  # Sort phenotypes by attained trait 1
  trait1 <- pheno %>%
    arrange(Trait1)
  # Get the middle individuals w.r.t. trait 1
  middle_t1 <- trait1[(mid-(window/2)):(mid+(window/2)),]  
  trait2 <- pheno %>%
    arrange(Trait2)
  # Get the middle individuals w.r.t. trait 2
  middle_t2 <- trait2[(mid-(window/2)):(mid+(window/2)),]
  
  # Find the individuals with middle phenotypes for both of the traits
  common_ids <- intersect(middle_t1$idx,
                          middle_t2$idx)
  
  # Get the phenotypes of the selected individuals
  pheno_filt <- pheno[pheno$idx %in% common_ids[1:5],]
  
  data_pts <- pheno_filt %>%
    dplyr::mutate(Suit=suitFunc(Trait1, Trait2)) %>%
    dplyr::mutate(Gen=gen) %>%
    dplyr::select(-c(Trait3, idx))
  
  res.df <- rbind(res.df,
                  data_pts)
  
  # Store mean phenotypes of the large population
  large_res.df <- rbind(large_res.df,
                      data.frame(Trait1=mean(pheno$Trait1),
                                 Trait2=mean(pheno$Trait2),
                                 Suit=suitFunc(mean(pheno$Trait1), mean(pheno$Trait2)),
                                 Gen=gen)
                      )
  
  # Store mean phenotypes of the medium population
  med_pheno <- as.data.frame(pheno(medPop))
  med_res.df <- rbind(med_res.df,
                      data.frame(Trait1=mean(med_pheno$Trait1),
                                 Trait2=mean(med_pheno$Trait2),
                                 Suit=suitFunc(mean(med_pheno$Trait1), mean(med_pheno$Trait2)),
                                 Gen=gen))
  
  
  # Store mean phenotypes of the small population
  small_pheno <- as.data.frame(pheno(smallPop))
  small_res.df <- rbind(small_res.df,
                        data.frame(Trait1=mean(small_pheno$Trait1),
                                   Trait2=mean(small_pheno$Trait2),
                                   Suit=suitFunc(mean(small_pheno$Trait1), mean(small_pheno$Trait2)),
                                   Gen=gen)
  )
  
  # Select and advance each population
  pop <- selectCross(pop, trait=breedingFitness, nInd=nInd(pop)*n.selProp, nCrosses=nInd(pop))
  medPop <- selectCross(medPop, trait=breedingFitness, nInd=nInd(medPop)*n.selProp, nCrosses=nInd(medPop))
  smallPop <- selectCross(smallPop, trait=breedingFitness, nInd=nInd(smallPop)*n.selProp, nCrosses=nInd(smallPop))
}

surface_fig <- plotAdaptiveWalkWithIndividuals(mean_res.df=large_res.df,
                                               sampled_inds.df=res.df)
surface_fig
htmlwidgets::saveWidget(as_widget(surface_fig), file.path(output_dir, "figure1.html"))

pts_fig <- plotSampledInds(res.df)
pts_fig
htmlwidgets::saveWidget(as_widget(pts_fig), file.path(output_dir, "sampled_individuals.html"))

three_pop_fig <- plotAdaptiveWalk3Populations(large_res.df,
                                              med_res.df,
                                              small_res.df)
three_pop_fig
htmlwidgets::saveWidget(as_widget(three_pop_fig), file.path(output_dir, "three_populations.html"))
