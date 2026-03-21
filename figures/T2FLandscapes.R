
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
                                   suitFunc=suitabilityGaussian,
                                   traitMin=-1,
                                   traitMax=1,
                                   popId_1="Subpopulation 1",
                                   popId_2="Subpopulation 2") {
  
  # Create a matrix of fitness values, with a small increment along the x and y axes.
  fitness_x = seq(traitMin,traitMax, by=(traitMax-traitMin)/40)
  fitness_y = seq(traitMin,traitMax, by=(traitMax-traitMin)/40)
  fitness_z = outer(fitness_x,fitness_y,suitFunc)
  
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
        x = pop1_df$traitVal1,
        y = pop1_df$traitVal2,
        type='scatter',
        mode = 'lines',
        line = list(color = '#CC0000', width = 4, dash = 'solid'),
        opacity = 1) %>%
      add_trace(
        fig,
        df,
        name = popId_2,
        x = pop2_df$traitVal1,
        y = pop2_df$traitVal2,
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
        x = pop1_df$traitVal1,
        y = pop1_df$traitVal2,
        z = pop1_df$suit+0.01,
        type = 'scatter3d',
        mode = 'lines',
        opacity = 1,
        line = list(autocolorscale=FALSE, color="#CC0000", which=2, width = 10)
      ) %>%
      add_trace(
        fig,
        df,
        name = popId_2,
        x = pop2_df$traitVal1,
        y = pop2_df$traitVal2,
        z = pop2_df$suit+0.01,
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
# suitFunc: The function for calculating fitness based on two traits
plot3dPopulationFitness <- function(pop, suitFunc=suitabilityGaussian) {
  popSize <- nInd(pop)
  df <- data.frame(trait1=pheno(pop)[,1],
                   trait2=pheno(pop)[,2],
                   suit=numeric(popSize))
  # Calculate the suitability for each individual in the population
  df <- df %>% mutate(suit=suitFunc(trait1, trait2))
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
      z = df$suit,
      type = 'scatter3d',
      mode = 'markers',
      color=df$suit)
  return (fig)
}

# Plot 2 populations on a 3D fitness surface (only works for populations with 2 traits)
# pop: The population to plot
# suitFunc: The function for calculating suitability based on two traits
plot3dPopulationFitnessTwoPops <- function(pop1, pop2, suitFunc) {
  df.a <- data.frame(trait1=pheno(pop1)[,1],
                     trait2=pheno(pop1)[,2],
                     suit=numeric(nInd(pop1)),
                     pop=rep("pop1", times=nInd(pop1)))
  df.b <- data.frame(trait1=pheno(pop2)[,1],
                     trait2=pheno(pop2)[,2],
                     suit=numeric(nInd(pop2)),
                     pop=rep("pop2", times=nInd(pop2)))
  df <- rbind(df.a, df.b)
  # Calculate the suitability for each individual in the population
  df <- df %>% mutate(suit=suitFunc(trait1, trait2))
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
      z = df$suit,
      type = 'scatter3d',
      mode = 'markers',
      opacity=0.9,
      color=df$pop,
      colors=c("#CC0000", "#3C78D8"))
  return (fig)
}

# Plots a 2D density plot of the suitability values of a population
plot2DFitness <- function(pop, suitFunc=suitabilityGaussian) {
  df <- data.frame(trait1=pheno(pop)[,1],
                   trait2=pheno(pop)[,2],
                   suit=numeric(nInd(pop)))
  df <- df %>% mutate(suit=suitFunc(trait1, trait2))
  g <- ggplot(df, aes(x=suit)) +
    geom_density()
  return(g)
}

# Creates a 3D fitness landscape and returns a plot_ly rendering of it.
# Valid colors arguments:
# colors = "PuBuGn"
# colors = colorRampPalette(c("blue", "orange"))(15)
# colors = magma(50, alpha = 1, begin = 0, end = 1, direction = 1) (viridis, plasma, magma, inferno)
plotFitnessLandscape <- function(suitFunc=suitabilityGaussian) {
  fitness_x = seq(-3,3, by=0.05)
  fitness_y = seq(-3,3, by=0.05)
  fitness_z = outer(fitness_x,fitness_y,suitFunc)
  
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