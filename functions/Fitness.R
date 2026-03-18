# Title: FITNESS
# Author: Ted Monyak
# Description: This file contains functions for calculating fitness and
# plotting fitness landscapes

library(ggplot2)
library(plotly)

# ONE TRAIT FITNESS FUNCTION

# Calculates fitness based on an optimum value of zero for one trait
# w = -(x^2)
oneTraitFitFunc <- function(x) {
  res <- -(x)^2
  return (res)
}

# SIMPLE QUADRATIC FITNESS FUNCTIONS


# Calculates fitness based on an optimum value of zero for each trait
# w = -(x^2) + -(y^2)
calculateFitnessTwoTrait <- function(x,y) {
  res <- -(x^2) - (y^2)
  return (res)
}

# Calculates fitness based on an optimum value of zero for each trait
# w = -(x^2) + -(y^2)
twoTraitFitFunc <- function(x) {
  return (calculateFitnessTwoTrait((x[,1]), x[,2]))
}

# For plotting purposes only - adds a small amount to each calculated value
# to be able to overlay it onto the surface
calculateFitnessTwoTraitModified <- function(x,y) {
  res <- calculateFitnessTwoTrait(x,y)
  res <- res + 0.01
  return (res)
}

# CONSTANT SLOPE FITNESS FUNCTION

# Calculates the fitness with a constant slope
# Renders in 3d space as a cone
# w = -sqrt(x^2 + y^2)
calculateFitnessConstantSlope <- function(x,y) {
  res <- -sqrt(x^2 + y^2)
  return (res)
}

# Calculates the fitness with a constant slope
# Renders in 3d space as a cone
# w = -sqrt(x^2 + y^2)
constantSlopeFitFunc <- function(x) {
  return (calculateFitnessConstantSlope(x[,1], x[,2]))
}

# SELECTION INDEX FITNESS FUNCTIONS

# Calculates fitness as a realized yield, as measured by the yield potential
# adjusted for the distance from optimal values for the 2 acquired traits
# w = -(x^2) + -(y^2)
# Return: (1-abs(w) * z)
landraceFitFunc <- function(x) {
  w <- twoTraitFitFunc(x)

  # Determine realized yield so that it is penalized by a lower fitness (as determined by traits A and B)
  realizedYield <- (1-abs(w)) * x[,3]
  #print(paste0("W: ", w, " P: ", x[,3], " Y: ", realizedYield))
  return (realizedYield)
}

# Calculates fitness as a weighted average between fitness of the two acquired traits
# and yield
# w = -(x^2) + -(y^2)
# yieldPotential = (1-abs(w) * z)
# Return: w*(1-yieldProp) + yield*yieldProp
breedingFitFunc <- function(x, yieldProp=n.yieldProp) {
  w <- twoTraitFitFunc(x)
  
  # Determine yieldPotential so that it is negatively penalized by a lower fitness (as determined by traits A and B)
  realizedYield <- (1-abs(w)) * x[,3]
  
  # Return a weighted average of fitness and yield
  return (w*(1-yieldProp) + realizedYield*yieldProp)
}

# GAUSSIAN FITNESS FUNCTIONS


calculateFitnessGaussian <- function(x,y) {
  sigma <- 1
  res <- exp(-(x^2 + y^2) / (2*sigma^2))
  return (res)
}

calculateBreedingFitness <- function(x,y,z) {
  w <- calculateFitnessGaussian(x,y)
  realizedYield <- w * z
  return (realizedYield)
}

gaussianFitFunc <- function(x) {
  return (calculateFitnessGaussian(x[,1], x[,2]))
}

gaussianLandraceFitFunc <- function(x) {
  w <- gaussianFitFunc(x)
  realizedYield <- w * x[,3]
  return (realizedYield)
}


calculateW_GWP <- function(x) {
  w <- gaussianFitFunc(x)
  realizedYield <- w * x[,3]
  dim(realizedYield) <- c(length(realizedYield), 1)
  return (realizedYield)
}

calculateFitnessGaussianModified <- function(x,y) {
  res <- calculateFitnessGaussian(x,y)
  res <- res + 0.01
  return (res)
}

bivariateNormal <- function(x,y) {
  sigma <- 2
  res <- (1/(2*pi*sqrt(1-sigma^2))) * exp(-(1/(2*(1-sigma^2)))*(x^2 - 2*sigma*x*y + y^2))
  return (res)
}

# Calculate a decaying selection ratio based on the distance from the fitness optimum
# Uses a geometric series to determine the result, where a=(1-n.selProp),
# r is set as an initial parameter (n.r), and n is a function of the distance from the initial fitness
# w: the current fitness
# Returns: a ratio between 0 and 1 which determines what percentage of individuals to advance
selectionRatio <- function(w) {
  # If at fitness optimum, return n.selProp
  if (w == 0) {
    return (n.selProp)
  }
  # If not using a geometric decay
  if (n.selR == 1) {
    return (n.selProp)
  }
  # Based on the simulation parameters, this is the starting fitness value
  initFit <- calculateFitnessTwoTrait(n.initTraitVal,n.initTraitVal)
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


# Set theme elements
theme <- theme(plot.background = ggplot2::element_blank(),
               panel.background = ggplot2::element_blank(),
               axis.line = ggplot2::element_line(linewidth = 0.2),
               plot.title = ggplot2::element_text(hjust = 0.5,
                                                  face = 'bold',
                                                  size = 12),
               axis.text = ggplot2::element_text(size  = 12,
                                                 color = 'black'),
               axis.title = ggplot2::element_text(size  = 12,
                                                  face = 'bold',
                                                  color = 'black'))



# Create a plot_ly figure of a fitness landscape, with an adaptive walk
# overlaid on it.
# pop1_df: a dataframe with 2 columns: traitValA, traitValB (for type "CONTOUR"),
# with a 3rd column (fitness) (for type "SURFACE")
# pop2_df: another dataframe
# type: one of 'CONTOUR' (for a 2D landscape) or 'SURFACE' (for a 3D landscape)
# fit calc: the function for determining fitness based on two trait values
# Also, supply the min and max trait values for the fitness landscape
# TODO: add a fixed value to df, because phenotypic data is 'under' the fitness curve
overlayWalkOnLandscape <- function(pop1_df,
                                   pop2_df,
                                   type="CONTOUR",
                                   fitCalc=calculateFitnessTwoTrait,
                                   traitMin=-1,
                                   traitMax=1,
                                   popId_1="Subpopulation 1",
                                   popId_2="Subpopulation 2") {
  
  # Create a matrix of fitness values, with a small increment along the x and y axes.
  fitness_x = seq(traitMin,traitMax, by=(traitMax-traitMin)/40)
  fitness_y = seq(traitMin,traitMax, by=(traitMax-traitMin)/40)
  fitness_z = outer(fitness_x,fitness_y,fitCalc)
  
  # Create the contours
  contour_list <- contourLines(x=fitness_x, y=fitness_y, z=fitness_z, nlevels=10)
  cx <- unlist(lapply(contour_list, function(cl) c(cl$x, NA)))
  cy <- unlist(lapply(contour_list, function(cl) c(cl$y, NA)))
  cz <- unlist(lapply(contour_list, function(cl) c(rep(cl$level + 0.005, length(cl$x)), NA)))
  
  if (type == "CONTOUR") {
    f <- list(family="Helvetica", size=20)
    fig <- plot_ly() %>%
      layout(font=list(
              family="Helvetica",
              size=24,
              color="black"),
            legend = list(title=list(text="Subpopulation")),
             xaxis = list(title = "Attained Trait 1", constrain = "domain"),
             yaxis = list(title = "Attained Trait 2", scaleanchor="x")) %>%
      add_trace(
        fig,
        x=fitness_x,
        y=fitness_y,
        z=fitness_z,
        type='contour',
        colors = 'Greys', # viridis(n=10)
        reversescale=TRUE,
        colorbar=list(title = "Suitability"),
        line = list(color = 'black', width = 1),
        opacity=1) %>%
      add_trace(
        fig,
        df,
        name = popId_1,
        x = pop1_df$traitValA,
        y = pop1_df$traitValB,
        type='scatter',
        mode = 'lines',
        line = list(color = '#CC0000', width = 4, dash = 'solid'),
        opacity = 1) %>%
      add_trace(
        fig,
        df,
        name = popId_2,
        x = pop2_df$traitValA,
        y = pop2_df$traitValB,
        type='scatter',
        mode = 'lines',
        line = list(color = '#3C78D8', width = 4, dash = 'solid'),
        opacity = 1)
    return (fig)
  } else if (type == "SURFACE"){
    f <- list(family="Helvetica", size=20)
    fig <- plot_ly() %>%
      layout(
        font=list(
          family="Helvetica",
          size=14,
          color="black"
        ),
        legend = list(title=list(text="Subpopulation")),
        scene = list(xaxis = list(title = "Attained Trait 1"),
                     yaxis = list(title = "Attained Trait 2"),
                     zaxis = list(title = "Suitability"),
                     aspectmode='cube'))
    fig <- fig %>%
      add_trace(
        fig,
        df,
        name = popId_1,
        x = pop1_df$traitValA,
        y = pop1_df$traitValB,
        z = pop1_df$fitness+0.01,
        type = 'scatter3d',
        mode = 'lines',
        opacity = 1,
        line = list(autocolorscale=FALSE, color="#CC0000", which=2, width = 10)
      ) %>%
      add_trace(
        fig,
        df,
        name = popId_2,
        x = pop2_df$traitValA,
        y = pop2_df$traitValB,
        z = pop2_df$fitness+0.01,
        type = 'scatter3d',
        mode = 'lines',
        opacity = 1,
        line = list(autocolorscale=FALSE, color="#3C78D8", which=2, width = 10)
      ) %>%
      add_trace(
        fig,
        x=fitness_x,
        y=fitness_y,
        z=fitness_z,
        type='surface',
        colorbar=list(title = "Suitability"),
        #colors = viridis(n=10),
        colorscale = 'Greys',
        opacity=1.0,
        contours = list(
          x = list(show = TRUE, start = min(fitness_x), 
                   end = max(fitness_x), size = 0.2, usecolormap = FALSE,
                   color = "black", width = 1),
          y = list(show = TRUE, start = min(fitness_y), 
                   end = max(fitness_y), size = 0.2, usecolormap = FALSE,
                   color = "black", width = 1)
        ))
    return (fig)
  } else {
    print("Type Not Supported")
  }
}

# Plot a population on a a 3D fitness surface (only works for populations with 2 traits)
# pop: The population to plot
# fitCalc: The function for calculating fitness based on two traits
plot3dPopulationFitness <- function(pop, fitCalc=calculateFitnessTwoTrait) {
  popSize <- nInd(pop)
  df <- data.frame(trait1=pheno(pop)[,1],
                   trait2=pheno(pop)[,2],
                   fitness=numeric(popSize))
  # Calculate the fitness for each individual in the population
  df <- df %>% mutate(fitness=fitCalc(traitA, traitB))
  fig <- plot_ly()
  fig <- plot_ly() %>%
    layout(scene = list(xaxis = list(title = "Attained Trait 1"),
                        yaxis = list(title = "Attained Trait 2"),
                        zaxis = list(title = "Suitability"),
                        aspectmode='cube'),
           annotations=list(text="Suitability", showarrow=FALSE, x=1.18, y=1.03)) %>%
    add_trace(
      fig,
      df,
      name = popSize,
      x = df$trait1,
      y = df$trait2,
      z = df$fitness,
      type = 'scatter3d',
      mode = 'markers',
      color=df$fitness)
  return (fig)
}

# Plot 2 populations on a 3D fitness surface (only works for populations with 2 traits)
# pop: The population to plot
# fitCalc: The function for calculating fitness based on two traits
plot3dPopulationFitnessTwoPops <- function(pop1, pop2, fitCalc) {
  df.a <- data.frame(trait1=pheno(pop1)[,1],
                     trait2=pheno(pop1)[,2],
                     fitness=numeric(nInd(pop1)),
                     pop=rep("pop1", times=nInd(pop1)))
  df.b <- data.frame(trait1=pheno(pop2)[,1],
                     trait2=pheno(pop2)[,2],
                     fitness=numeric(nInd(pop2)),
                     pop=rep("pop2", times=nInd(pop2)))
  df <- rbind(df.a, df.b)
  # Calculate the fitness for each individual in the population
  df <- df %>% mutate(fitness=fitCalc(trait1, trait2))
  fig <- plot_ly()
  fig <- plot_ly() %>%
    layout(scene = list(xaxis = list(title = "Attained Trait 1"),
                        yaxis = list(title = "Attained Trait 2"),
                        zaxis = list(title = "Suitability"),
                        aspectmode='cube')) %>%
    add_trace(
      df,
      name = df$pop,
      x = df$trait1,
      y = df$trait2,
      z = df$fitness,
      type = 'scatter3d',
      mode = 'markers',
      opacity=0.9,
      color=df$pop,
      colors=c("#CC0000", "#3C78D8"))
  return (fig)
}

# Plots a 2D density plot of the fitness values of a population
plot2DFitness <- function(pop, fitFunc=calculateFitnessTwoTrait) {
  df <- data.frame(trait1=pheno(pop)[,1],
                   trait2=pheno(pop)[,2],
                   fitness=numeric(nInd(pop)))
  df <- df %>% mutate(fitness=fitCalc(traitA, traitB))
  g <- ggplot(df, aes(x=fitness)) +
    geom_density()
  return(g)
}

# Creates a 3D fitness landscape and returns a plot_ly rendering of it.
# Valid colors arguments:
# colors = "PuBuGn"
# colors = colorRampPalette(c("blue", "orange"))(15)
# colors = magma(50, alpha = 1, begin = 0, end = 1, direction = 1) (viridis, plasma, magma, inferno)
plotFitnessLandscape <- function(fitFunc=calculateFitnessTwoTrait) {
  fitness_x = seq(-3,3, by=0.05)
  fitness_y = seq(-3,3, by=0.05)
  fitness_z = outer(fitness_x,fitness_y,fitFunc)
  
  p <- plot_ly(x=fitness_x,
               y=fitness_y,
               z=fitness_z,
               type='surface',
               colors = viridis(50, alpha = 1, begin = 0, end = 1, direction = 1),
               opacity=1) %>%
    layout(scene = list(xaxis = list(title = "Attained Trait 1"),
                        yaxis = list(title = "Attained Trait 2"),
                        zaxis = list(title = "Suitability"),
                        aspectmode='cube'))
  return (p) 
}



# Show a 2D curve of a population's fitness over generations
plotFitness <- function(df) {
  g <- ggplot(data=df, aes(x=gen, y=fitness)) +
    geom_line() +
    geom_point() +
    labs(x = "Generation", y = "Fitness")
  return (g)
}

# Plot Genetic Values for Two Traits
plotHist <- function(pop) {
  df <- as.data.frame(cbind(gv(pop)[,1], gv(pop)[,2]))
  colnames(df) <- c("trait1", "trait2")
  df <- df %>%
    pivot_longer(c("trait1", "trait2"), names_to="trait", values_to="gv")
  g <- ggplot(df, aes(gv, fill=trait, color=trait)) +
    geom_density(alpha=0.1) +
    labs(title="Genetic Values")
  return (g)
}

