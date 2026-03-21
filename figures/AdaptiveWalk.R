# 

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

setwd("~/Documents/CSU/Lab-Notebooks/members/Ted/breeding_sims")
output_dir <- file.path(getwd(), "Output")
save_dir <- file.path(output_dir, "Wright")
n.cores <- 8

source("Functions/Fitness.R")
source("Scripts/GlobalParameters.R")
source("Scripts/CreateFounderPop.R")

suitFunc <- suitabilityGaussian
n.gens <- 50

# For storing the results at each generation
res.df <- data.frame(Trait1=c(),
                     Trait2=c(),
                     Suit=c(),
                     Gen=c())


# Select 1000 random individuals
pop <- selectInd(founderPop, nInd=1000, use="rand")
res_agg.df <- data.frame(Trait1=c(),
                         Trait2=c(),
                         Suit=c(),
                         Gen=c())

# Specify the "middle" range of the large population from which to
# sample individuals
mid <- 500
int <- 50

smallPop <- selectInd(founderPop, nInd=500, use="rand")
small_res.df <- data.frame(Trait1=c(),
                         Trait2=c(),
                         Fitness=c(),
                         Gen=c())

smallerPop <- selectInd(founderPop, nInd=200, use="rand")
smaller_res.df <- data.frame(Trait1=c(),
                           Trait2=c(),
                           Fitness=c(),
                           Gen=c())

for (gen in 1:n.gens) {
  pheno <- as.data.frame(pheno(pop)) %>%
    rownames_to_column("idx")
  
  trait1 <- pheno %>%
    arrange(Trait1)
  middle_t1 <- trait1[(mid-int):(mid+int),]  
  trait2 <- pheno %>%
    arrange(Trait2)
  middle_t2 <- trait2[(mid-int):(mid+int),]  
  
  common_ids <- intersect(middle_t1$idx,
                          middle_t2$idx)
  
  pheno_filt <- pheno[pheno$idx %in% common_ids[1:5],]
  
  data_pts <- pheno_filt %>%
    dplyr::mutate(Fitness=fitCalc(Trait1, Trait2)) %>%
    dplyr::mutate(Gen=gen) %>%
    dplyr::select(-c(Trait3, idx))
  
  res.df <- rbind(res.df,
                  data_pts)
  
  res_agg.df <- rbind(res_agg.df,
                      data.frame(Trait1=mean(pheno$Trait1),
                                 Trait2=mean(pheno$Trait2),
                                 Fitness=fitCalc(mean(pheno$Trait1), mean(pheno$Trait2)),
                                 Gen=gen)
                      )
  
  small_pheno <- as.data.frame(pheno(smallPop))

  small_res.df <- rbind(small_res.df,
                      data.frame(Trait1=mean(small_pheno$Trait1),
                                 Trait2=mean(small_pheno$Trait2),
                                 Fitness=fitCalc(mean(small_pheno$Trait1), mean(small_pheno$Trait2)),
                                 Gen=gen)
  )
  
  smaller_pheno <- as.data.frame(pheno(smallerPop))
  
  smaller_res.df <- rbind(smaller_res.df,
                        data.frame(Trait1=mean(smaller_pheno$Trait1),
                                   Trait2=mean(smaller_pheno$Trait2),
                                   Fitness=fitCalc(mean(smaller_pheno$Trait1), mean(smaller_pheno$Trait2)),
                                   Gen=gen)
  )
  
  pop <- selectCross(pop, trait=fitFunc, nInd=nInd(pop)*0.1, nCrosses=nInd(pop))
  smallPop <- selectCross(smallPop, trait=fitFunc, nInd=nInd(smallPop)*0.1, nCrosses=nInd(smallPop))
  smallerPop <- selectCross(smallerPop, trait=fitFunc, nInd=nInd(smallerPop)*0.1, nCrosses=nInd(smallerPop))
}


traitMin <- -3
traitMax <- 3
fitness_x = seq(traitMin,traitMax, by=(traitMax-traitMin)/40)
fitness_y = seq(traitMin,traitMax, by=(traitMax-traitMin)/40)
fitness_z = outer(fitness_x,fitness_y,fitCalc)

# Create the contours
contour_list <- contourLines(x=fitness_x, y=fitness_y, z=fitness_z, nlevels=10)
cx <- unlist(lapply(contour_list, function(cl) c(cl$x, NA)))
cy <- unlist(lapply(contour_list, function(cl) c(cl$y, NA)))
cz <- unlist(lapply(contour_list, function(cl) c(rep(cl$level + 0.005, length(cl$x)), NA)))

surface_fig <- plot_ly() %>%
  add_trace(
    data=res.df,
    name = ~Gen,
    x = ~Trait1,
    y = ~Trait2,
    z = ~Fitness+0.01,
    type = 'scatter3d',
    mode = 'markers',
    marker=list(
      autocolorscale=FALSE,
      size=5,
      colorscale="Electric", # Other colors: Blackbody, Hot, Portland, Greens, Reds, Blues, Inferno, Cividis, Jet, Earth, Picnic
      reversescale=TRUE,
      color=res.df$Gen,
      cmin=min(res.df$Gen),
      cmax=max(res.df$Gen),
      colorbar = list(
        title = "Generation",
        len=0.5
      )
    ),
    showlegend=FALSE) %>%
  add_trace(
    data=res_agg.df,
    x = ~Trait1,
    y = ~Trait2,
    z = ~Fitness+0.01,
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
surface_fig
htmlwidgets::saveWidget(as_widget(surface_fig), file.path(save_dir, "wright_surface_trajectory_landrace.html"))

two_pop_fig <- plot_ly(
  x=fitness_x,
  y=fitness_y,
  z=fitness_z,
  type="surface",
  #showscale=FALSE,
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
two_pop_fig <- add_trace(p=two_pop_fig,
    data=res_agg.df,
    name="1000",
    x = ~Trait1,
    y = ~Trait2,
    z = ~Fitness+0.03,
    type = 'scatter3d',
    mode = 'lines',
    opacity = 1,
    line = list(autocolorscale=FALSE, color="#228B22", which=2, width = 10)
  )
two_pop_fig <- add_trace(
  p=two_pop_fig,
    data=small_res.df,
    name="500",
    x = ~Trait1,
    y = ~Trait2,
    z = ~Fitness+0.03,
    type = 'scatter3d',
    mode = 'lines',
    opacity = 1,
    line = list(autocolorscale=FALSE, color="#F4A7B9", which=2, width = 10)
  )
two_pop_fig <- add_trace(
  p=two_pop_fig,
    data=smaller_res.df,
    name="200",
    x = ~Trait1,
    y = ~Trait2,
    z = ~Fitness+0.03,
    type = 'scatter3d',
    mode = 'lines',
    opacity = 1,
    line = list(autocolorscale=FALSE, color="#F4A460", which=2, width = 10)
  )
two_pop_fig <- two_pop_fig %>%
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
two_pop_fig
htmlwidgets::saveWidget(as_widget(two_pop_fig), file.path(save_dir, "two_populations.html"))

pts_fig <- plot_ly() %>%
  add_trace(
    data=res.df,
    name = ~Gen,
    x = ~Trait1,
    y = ~Trait2,
    z = ~Fitness+0.01,
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
pts_fig
htmlwidgets::saveWidget(as_widget(pts_fig), file.path(save_dir, "wright_trajectory_landrace.html"))
