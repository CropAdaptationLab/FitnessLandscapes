# Title: Trait to Fitness Landscapes
# Author: Ted Monyak
# Description: Contains functions for plotting trait-to-fitness surfaces and overlaying adaptive walks

# Create a plot_ly figure of a fitness landscape, with two adaptive walks
# overlaid on it, one for each subpopulation.
# Suitability values are calculated based on the mean values of each of the attained traits
# pop1_df: a dataframe with 2 columns: traitVal1, traitVal2 (for type "CONTOUR"),
# with a 3rd column (fitness) (for type "SURFACE")
# pop2_df: another dataframe
# type: one of 'CONTOUR' (for a 2D landscape) or 'SURFACE' (for a 3D landscape)
# suitFunc: the function for determining suitability based on two trait values
# Also, supply the min and max trait values for the fitness landscape
overlayWalkOnLandscape <- function(pop1_df,
                                   pop2_df,
                                   type="CONTOUR",
                                   suitFunc=suitabilityGaussian,
                                   traitMin=-n.initTraitVal-1,
                                   traitMax=n.initTraitVal+1,
                                   popId_1="Subpopulation 1",
                                   popId_2="Subpopulation 2") {
  
  # Create a matrix of fitness values, with a small increment along the x and y axes.
  fitness_x = seq(traitMin,traitMax, by=(traitMax-traitMin)/40)
  fitness_y = seq(traitMin,traitMax, by=(traitMax-traitMin)/40)
  fitness_z = outer(fitness_x,fitness_y,suitFunc)
  
  f <- list(family="Helvetica", size=20)
  fig <- plot_ly()
  if (type == "CONTOUR") {
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

# Plots the adaptive walk of a population, based on the mean values of the attained
# traits, not the mean suitability value at each generation
plotAdaptiveWalk <- function(df) {
  fig <- plot_ly()
  fig <- add_trace(
    fig,
    df,
    name = p,
    x = df$traitVal1,
    y = df$traitVal2,
    z = df$suit,
    type = 'scatter3d',
    mode = 'lines',
    opacity = 1,
    color = p,
    line = list(width = 5))
  fig <- fig %>% layout(legend=list(title=list(text='Population Size')),
                        scene = list(xaxis = list(title = "Trait 1"),
                                     yaxis = list(title = "Trait 2"),
                                     zaxis = list(title = "Suitability"),
                                     aspectmode='cube')) %>% hide_colorbar()
}


# Creates a 3D fitness landscape and returns a plot_ly rendering of it.
# Valid colors arguments:
# colors = "PuBuGn"
# colors = colorRampPalette(c("blue", "orange"))(15)
# colors = magma(50, alpha = 1, begin = 0, end = 1, direction = 1) (viridis, plasma, magma, inferno)
plotFitnessLandscape <- function(traitMin=-n.initTraitVal-1,
                                 traitMax=n.initTraitVal+1,
                                 suitFunc=suitabilityGaussian) {
  fitness_x = seq(traitMin,traitMax, by=(traitMax-traitMin)/40)
  fitness_y = seq(traitMin,traitMax, by=(traitMax-traitMin)/40)
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

# Plots an adaptive walk of a single population, along with a sampling of the
# most representative individuals from that population
# Run this on data generated from AdaptiveWalk.R
# mean_res.df: Contains the mean trait values from the adaptive walk
# sampled_inds.df: Contains the sampled individuals from each generation of
# the adaptive walk (based on vicinity to the mean attained trait values)
plotAdaptiveWalkWithIndividuals <- function(mean_res.df,
                                            sampled_inds.df,
                                            suitFunc=suitabilityGaussian,
                                            traitMin=-n.initTraitVal-1,
                                            traitMax=n.initTraitVal+1) {
  fitness_x = seq(traitMin,traitMax, by=(traitMax-traitMin)/40)
  fitness_y = seq(traitMin,traitMax, by=(traitMax-traitMin)/40)
  fitness_z = outer(fitness_x,fitness_y,suitFunc)
  
  surface_fig <- plot_ly() %>%
    add_trace(
      data=sampled_inds.df,
      name = ~Gen,
      x = ~Trait1,
      y = ~Trait2,
      z = ~Suit+0.01,
      type = 'scatter3d',
      mode = 'markers',
      marker=list(
        autocolorscale=FALSE,
        size=5,
        colorscale="Electric", # Other colors: Blackbody, Hot, Portland, Greens, Reds, Blues, Inferno, Cividis, Jet, Earth, Picnic
        reversescale=TRUE,
        color=sampled_inds.df$Gen,
        cmin=min(sampled_inds.df$Gen),
        cmax=max(sampled_inds.df$Gen),
        colorbar = list(
          title = "Generation",
          len=0.5
        )
      ),
      showlegend=FALSE) %>%
    add_trace(
      data=mean_res.df,
      x = ~Trait1,
      y = ~Trait2,
      z = ~Suit+0.01,
      type = 'scatter3d',
      mode = 'lines',
      opacity = 1,
      line = list(autocolorscale=FALSE, color="white", which=2, width = 10, dash="solid")
    ) %>%
    add_surface(
      x=fitness_x,
      y=fitness_y,
      z=fitness_z,
      showscale=FALSE,
      colorscale='Greys',
      opacity=1.0,
      contours = list(
        x = list(show = TRUE, start = min(fitness_x), 
                 end = max(fitness_x), size = 0.2, usecolormap = FALSE,
                 color = "black", width = 1),
        y = list(show = TRUE, start = min(fitness_y), 
                 end = max(fitness_y), size = 0.2, usecolormap = FALSE,
                 color = "black", width = 1)
      )) %>%
    layout(
      font=list(
        family="Helvetica",
        size=16,
        color="black"
      ),
      scene = list(xaxis = list(title = "Attained Trait 1"),
                   yaxis = list(title = "Attained Trait 2"),
                   zaxis = list(title = "Suitability"),
                   aspectmode='cube'))
  return (surface_fig)
}

# Plots the mean adaptive walks of 3 populations of different sizes
# To visualize the effect of genetic drift
# Should be run after AdaptiveWalk.R
# large_res.df: Stores the adaptive walk of the large population (N=1000)
# med_res.df: (N=500)
# small_res.df (N=200)
plotAdaptiveWalk3Populations <- function(large_res.df,
                                         med_res.df,
                                         small_res.df,
                                         suitFunc=suitabilityGaussian,
                                         traitMin=-n.initTraitVal-1,
                                         traitMax=n.initTraitVal+1) {
  fitness_x = seq(traitMin,traitMax, by=(traitMax-traitMin)/40)
  fitness_y = seq(traitMin,traitMax, by=(traitMax-traitMin)/40)
  fitness_z = outer(fitness_x, fitness_y, suitFunc)

  three_pop_fig <- plot_ly(
    x=fitness_x,
    y=fitness_y,
    z=fitness_z,
    type="surface",
    colorscale='Greys',
    opacity=1.0,
    hidesurface = FALSE,
    opacity=1,
    colorbar=list(title="Suitability"),
    contours = list(
      x = list(show = TRUE, start = min(fitness_x), 
               end = max(fitness_x), size = 0.2, usecolormap = FALSE,
               color = "black", width = 1),
      y = list(show = TRUE, start = min(fitness_y), 
               end = max(fitness_y), size = 0.2, usecolormap = FALSE,
               color = "black", width = 1)
    ))
  three_pop_fig <- add_trace(p=three_pop_fig,
                           data=large_res.df,
                           name="1000",
                           x = ~Trait1,
                           y = ~Trait2,
                           z = ~Suit+0.03,
                           type = 'scatter3d',
                           mode = 'lines',
                           opacity = 1,
                           line = list(autocolorscale=FALSE, color="#228B22", which=2, width = 10)
  )
  three_pop_fig <- add_trace(
    p=three_pop_fig,
    data=med_res.df,
    name="500",
    x = ~Trait1,
    y = ~Trait2,
    z = ~Suit+0.03,
    type = 'scatter3d',
    mode = 'lines',
    opacity = 1,
    line = list(autocolorscale=FALSE, color="#F4A7B9", which=2, width = 10)
  )
  three_pop_fig <- add_trace(
    p=three_pop_fig,
    data=small_res.df,
    name="200",
    x = ~Trait1,
    y = ~Trait2,
    z = ~Suit+0.03,
    type = 'scatter3d',
    mode = 'lines',
    opacity = 1,
    line = list(autocolorscale=FALSE, color="#F4A460", which=2, width = 10)
  )
  three_pop_fig <- three_pop_fig %>%
    layout(
      font=list(
        family="Helvetica",
        size=16,
        color="black"
      ),
      legend = list(title=list(text="Population Size")),
      scene = list(xaxis = list(title = "Attained Trait 1"),
                   yaxis = list(title = "Attained Trait 2"),
                   zaxis = list(title = "Suitability"),
                   aspectmode='cube'))
  return (three_pop_fig)
}

# Plot the sampled individuals only (with no surface)
# res.df has the attained trait and suitability values of 5 sampled individuals
# per generation
# Run this on data generated from AdaptiveWalk.R
plotSampledInds <- function(res.df) {
  pts_fig <- plot_ly() %>%
    add_trace(
      data=res.df,
      name = ~Gen,
      x = ~Trait1,
      y = ~Trait2,
      z = ~Suit,
      type = 'scatter3d',
      mode = 'markers',
      color=res.df$Gen,
      colors=colorRampPalette(brewer.pal(9, "Greens"))(n.gens),
      marker=list(
        size=6
      ),
      showlegend=FALSE) %>%
    hide_colorbar() %>%
    layout(
      font=list(
        family="Helvetica",
        size=16,
        color="black"
      ),
      scene = list(xaxis = list(title = "Attained Trait 1"),
                   yaxis = list(title = "Attained Trait 2"),
                   zaxis = list(title = "Suitability"),
                   aspectmode='cube'))
  return (pts_fig)
}