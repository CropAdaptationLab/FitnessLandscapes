# Title: FITNESS
# Author: Ted Monyak
# Description: This file contains functions for calculating fitness and
# plotting fitness landscapes

library(ggplot2)
library(plotly)

# Calculates the fitness with a constant slope
# Renders in 3d space as a cone
# w = -sqrt(t1^2 + t2^2) + 1
suitabilityConstantSlope <- function(t1,t2) {
  res <- -sqrt(t1^2 + t2^2) + 1
  return (res)
}

# Calculates fitness based on an optimum value of zero for each trait
# w = -(t1^2) + -(t2^2) + 1
suitabilityQuadratic <- function(t1,t2) {
  # Add 1 so the range is [0,1]
  res <- -(t1^2) - (t2^2) + 1
  return (res)
}

# Calculate fitness with a bivariate normal peak
suitabilityBivariateNormal <- function(t1, t2) {
  sigma <- 2
  res <- (1/(2*pi*sqrt(1-sigma^2))) * exp(-(1/(2*(1-sigma^2)))*(t1^2 - 2*sigma*t1*t2 + t2^2))
  return (res)
}

# Calculate fitness with a gaussian peak
suitabilityGaussian <- function(t1, t2) {
  res <- exp(-(t1^2 + t2^2) / 2)
  return (res)
}

# Calculate breeding fitness based on the suitability of the two attained traits (t1 and t2)
# And the yield potential (t3)
calculateBreedingFitness <- function(t1, t2, t3, suitFunc=suitabilityGaussian) {
  suit <- suitFunc(t1, t2)
  realizedYield <- suit * t3
  return (realizedYield)
}

# Breeding fitness phenotype
breedingFitness <- function(x, suitFunc=suitabilityGaussian) {
  return (calculateBreedingFitness(x[,1], x[,2], x[,3], suitFunc))
}

# Calculates suitability alone
suitability <- function(x, suitFunc=suitabilityGaussian) {
  return (suitFunc(x[,1], x[,2]))
}

# CURRENTLY NOT BEING USED
# Calculate a decaying selection ratio based on the distance from the fitness optimum
# Uses a geometric series to determine the result, where a=(1-n.selProp),
# r is set as an initial parameter (n.r), and n is a function of the distance from the initial fitness
# w: the current fitness
# Returns: a ratio between 0 and 1 which determines what percentage of individuals to advance
selectionRatio <- function(w, suitFunc) {
  # If at fitness optimum, return n.selProp
  if (w == 0) {
    return (n.selProp)
  }
  # If not using a geometric decay
  if (n.selR == 1) {
    return (n.selProp)
  }
  # Based on the simulation parameters, this is the starting fitness value
  initSuit <- suitFunc(n.initTraitVal,n.initTraitVal)
  # The initial selection ratio
  a <- 1-n.selProp
  # The "n" term in the geometric series increases as the the distance from the initial fitness increases
  n <- abs(initFit/w)
  # Return a geometrically increasing value (which increases with n)
  return (1-(a * n.selR^(n-1)))
}

# This section will return 1s for all of the indices in randVec that match 'idx',
# specifying which individuals to select. This ensures there are no overlaps.
# idx: The index of the subpopulation being selected
# randVec: should be created previously with this line: sample(rep(c(1:n.nPops), times=n.popSize/n.nPops))
# Returns: a vector of size randVec where all the values matching idx are 1, and all others are 0
selectSubPop <- function(x, idx, randVec) {
  as.numeric(lapply(randVec, idx, FUN=setequal))
}