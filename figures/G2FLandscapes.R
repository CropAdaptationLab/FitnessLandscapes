# Title: Genotype-to-Fitness Landscapes
# Author: Ted Monyak
# Description: This file contains plotting functions for creating genotype to
# fitness landscapes

library(fields)
library(ggdensity)
library(reshape2)
library(sp)
library(stats)
library(terra)

# Plots a one-dimensional fitness landscape where PC1 is plotted on the x-axis
# and breeding fitness is on the y-axis
# PCA is run on the oligogenic component of the landrace populations and the
# admixed RIL family is projected onto PC1
# Samples are color-coded by their haplotype w.r.t. the most significant
# epistatic effect for breeding fitness in the RIL family
# A trendline is drawn as the average fitness across PC1
# RIL: an admixed RIL family developed from a cross between parent1 and parent2
# pop1: The landrace subpopulation from which parent1 was derived
# pop2: The landrace subpopulation from which parent2 was derived
# qtl1: One of the QTL from a significant pairwise interaction
# qtl2: The other QTL from a signficant pairwise interaction
# ril_dir: The directory in which to save the resulting plot
plot1DLandscape <- function(RIL, pop1, pop2, parent1, parent2, qtl1, qtl2, snpChip, ril_dir) {
  # Make a larger RIL to reduce noise
  RIL <- self(RIL, nProgeny=5)
  
  # Get the genotypes for the attained traits only
  RIL_geno <- getUniqueQtl(RIL, traits=c(1,2))
  subpops_geno <- rbind(getUniqueQtl(pop1, traits=c(1,2)),
                        getUniqueQtl(pop2, traits=c(1,2)))
  
  haplos <- getHaplos(RIL_geno, qtl1, qtl2, parent1, parent2, snpChip)
  
  # Calculate the breeding fitness of each RIL sample
  RIL_pheno.df <- data.frame(fitness=breedingFitness(pheno(RIL)),
                             haplo=haplos)
  
  # Run PCA on the landrace subpopulations
  pca <- prcomp(subpops_geno)
  VAF <- summary(pca)$importance[2,1] * 100
  
  # Project the admixed RIL family onto the principal components
  pc_RIL <- predict(pca, RIL_geno)
  
  # Create a dataframe for plotting
  df.PCA = data.frame(
    "Fitness" = RIL_pheno.df$fitness,
    "PC1" = pc_RIL[, 1],
    "PC2" = pc_RIL[, 2],
    "Haplotype" = RIL_pheno.df$haplo)
  
  # Plot the [haplotype] color-coded samples on PC1
  ggplot(df.PCA, aes(x = PC1, y = Fitness)) +
    geom_hdr(aes(fill = Haplotype), alpha = 0.3, color = NA, probs = seq(0.1, 0.9, by = 0.1)) +
    geom_point(aes(color = Haplotype), alpha = 0.1, size = 0.8) +
    geom_smooth(
      #aes(group = Haplotype, color = Haplotype),
      color="black", # comment this and uncomment above to get haplo-specific lines
      method = "loess",
      se = FALSE,
      linewidth = 0.8
    ) +
    scale_color_manual(values = c("#CC0000", "#3C78D8", "gold"),,
                       labels = c("Parental Type 1", "Parental Type 2", "Recombinant"),
                       name="Haplotype at\nEpistatic Loci") +
    scale_fill_manual(values = c("#CC0000", "#3C78D8", "gold"),
                      labels = c("Parental Type 1", "Parental Type 2", "Recombinant"),
                      name="Haplotype at\nEpistatic Loci") +
    guides(alpha="none") +
    xlab("PC1 of Oligogenic Component") +
    ylab("Breeding Fitness") +
    theme_minimal(base_size = 14,
                  base_family="Helvetica") +
    theme(
      plot.title = element_blank(),
      axis.title.y = element_text(margin=margin(t=0, r=10, b=0, l=10, unit="pt")),
      axis.text.x = element_blank(),
      legend.position = "right",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white")
    )
  ggplot2::ggsave(filename = file.path(ril_dir, "PCA_Bridge.jpg"),
                  device = "jpg",
                  width=5.5,
                  height=3.5,
                  dpi=600)
  
  ggplot2::ggsave(filename = file.path(ril_dir, "PCA_Bridge.pdf"),
                  device = "pdf",
                  width=5.5,
                  height=3.5)
}

# Plots individuals along an adaptive walk on a genotype-to-fitness landscape
# Assumes the following objects exist: sampled_inds, RIL, pops[[1]], pops[[2]]
# sampled_inds: an AlphaSimR population made up of one or more subpopulations
# traversing to a fitness peak
# Returns: a plot_ly object with the fitness surface and the individuals projected onto it
plotAdaptiveWalkLandscape <- function() {
  # Run PCA on the RIL and both subpopulations
  pca_geno <- pullSnpGeno(c(RIL,
                            pops[[1]],
                            pops[[2]]), snpChip=2)
  pca <- prcomp(pca_geno)
  # Variance explained
  VAF <- summary(pca)$importance[2,] * 100
  
  # Project the admixed RIL family onto the principal components
  sampled_inds_geno <- pullSnpGeno(sampled_inds, snpChip=2)
  pc_pred <- predict(pca, sampled_inds_geno)
  
  # Determine suitability based on genotypic values
  # Heritability is low, and this is a G > F landscape, so genotypic values
  # as phenotypes are appropriate
  pheno <- as.data.frame(gv(sampled_inds)) %>%
    dplyr::mutate(Suit=suitFunc(Trait1, Trait2)) %>%
    dplyr::pull(Suit)
  
  # Create a dataframe for plotting
  suit_df = data.frame(
    "PC1" = pc_pred[, 1],
    "PC2" = pc_pred[, 2],
    "Suit" = pheno)
  
  # Plot color-coded PCA
  #ggplot(suit_df, aes(x=PC1, y=PC2, color=Suit)) +
  #  geom_point()
  
  # Add in the metadata (subpop and generation)
  suit_df <- cbind(suit_df, sampled_inds_metadata)

  # Add color-coding
  suit_df <- suit_df %>%
    dplyr::mutate(
      family_color= case_when(
        subpop == "Founder" ~ "#B04EDE",
        subpop == "Subpop. 1" ~ "#CC0000",
        subpop == "Subpop. 2" ~ "#3C78D8"
      )) %>%
    dplyr::rename(
      "Family" = "subpop"
    )

  # Downsample so that there is just one individual per subpop per generation
  # Take the "mean" indvidual based on euclidean distance from the mean PC1 and PC2
  # from each generation
  adaptive_walk_df <- suit_df %>%
    dplyr::group_by(gen, Family) %>%
    dplyr::mutate(
      mean_PC1 = mean(PC1, na.rm = TRUE),
      mean_PC2 = mean(PC2, na.rm = TRUE),
      # Euclidean distance from the centroid
      dist_to_mean = sqrt((PC1 - mean_PC1)^2 + (PC2 - mean_PC2)^2)
    ) %>%
    # Get closest individual
    dplyr::slice_min(dist_to_mean, n = 1, with_ties = FALSE) %>%
    dplyr::select(-mean_PC1, -mean_PC2, -dist_to_mean) %>%
    ungroup() %>%
    dplyr::filter(Family != "Founder")

  # Create a smoothed landscape
  landscape <- generate_landscape(df=suit_df,
                                  pred_df=adaptive_walk_df,
                                  pcx=1,
                                  pcy=2)

  # Render the new smoothed landscape
  # render_2d_landscape(landscape[[2]], VAF, 1, 2)
  
  p <- render_3d_landscape(smoothed=landscape[[3]],
                      df=landscape[[1]],
                      plot_inds=FALSE,
                      families=c("Founder", "Subpop. 1", "Subpop. 2"),
                      pcx=1,
                      pcy=2)
  return (p)
}

# Created a smoothed landscape based on a color-coded PCA result
# Uses a Thin Plate Spline from the fields() package to estimate a smooth surface
# df should have as its first columns the principal components, in order
# Suit should be a column appearing after the PCs
# pcx and pcy are the principal coordinates to use as the axes
# pred_df contains a set of points to project onto the landscape
# Returns a list:
#   - pred_df: The original pred_df, with the points projected onto the estimated landscape
#   - grid: A dataframe with a 'Suitability' value for each coordinate
#   - smoothed: A list with x (pcx), y (pcy), and z (a wide matrix format of the suitability values)
generate_landscape <- function(df, pred_df, pcx, pcy) {
  # Fit a thin plate spline to the data
  tps_fit <- fields::Tps(x = cbind(df[[pcx]], df[[pcy]]), Y = df$Suit)
  # Predict values on a full 100x100 square grid
  n_grid <- 100
  grid <- expand.grid(
    PCX = seq(min(df[[pcx]]), max(df[[pcx]]), length = n_grid),
    PCY = seq(min(df[[pcy]]), max(df[[pcy]]), length = n_grid)
  )
  grid$Suitability <- predict(tps_fit, x = as.matrix(grid[, c("PCX", "PCY")]))
  
  # Transform the long-format grid into an array
  suit_matrix <- reshape2::acast(grid,  PCX ~ PCY, value.var = "Suitability")
  smoothed <- list(
    x = as.numeric(colnames(suit_matrix)),  # PCX
    y = as.numeric(rownames(suit_matrix)),  # PCY
    z = suit_matrix
  )
  
  # Project the PCA coordinates onto the plane estimated by thin plate spline
  # Add a slight offset to get the dots to render above the surface
  pred_df$Suitability <- predict(tps_fit, 
                                 x = cbind(pred_df[[pcx]],
                                           pred_df[[pcy]])) + 0.005
  
  return (list(pred_df, grid, smoothed))
}

# Render a 2d image of the smoothed G > F landscape
# grid: a dataframe with PCX, PCY, and Suitability values
# PCX and PCY should span an entire square grid coordinate system
# VAF: the variance explained by each PC
# pcx: The principal component for the x-axis
# pcy: The principal component for the y-axis
render_2d_landscape <- function(grid, VAF, pcx=1, pcy=2) {
  ggplot(grid, aes(x = PCX, y = PCY, fill = Suitability)) +
    geom_tile() +
    scale_fill_viridis_c(name = "Suitability") +
    labs(x = paste0("PC", pcx),
         y = paste0("PC", pcy)) +
    #labs(x = paste0("PC", pcx, ": ", round(VAF[pcx]), "%"),
    #     y = paste0("PC", pcy, ": ", round(VAF[pcy]), "%")) +
    theme_minimal(base_family="Helvetica", base_size=8) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 1),
          axis.title.x = element_text(hjust = 0.5),
          axis.title.y = element_text(hjust = 0.5),
          axis.text = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "white", color = "black"),
          aspect.ratio=1)
}

# Recreate Wright's canonical fitness landscape (Wright 1932)
# smoothed: The result of a smooth.2d() call
# pcx: The principal component for the x-axis
# pcy: The principal component for the y-axis
# window: the size of the focal window to use to locate peaks and valleys
wright_landscape <- function(smoothed, pcx, pcy, window=5) {
  # Transpose the z matrix
  z_mat <- t(smoothed[["z"]])
  
  # Find peaks and valleys
  find_local_extrema <- function(z_mat, x_vec, y_vec, type = "peak", window = 11) {
    # Convert to a raster
    r <- rast(z_mat)
    # Create a square window
    w <- matrix(1, nrow = window, ncol = window)
    
    if (type == "peak") {
      # Slide a window over every cell
      ref <- focal(r, w = w, fun = max, na.policy = "omit")
      # Use a boolean mask to determine if the original raster value is
      # also the neighborhood max
      is_extr <- (r == ref)
    } else {
      ref <- focal(r, w = w, fun = min, na.policy = "omit")
      is_extr  <- (r == ref)
    }
    
    # Get coordinates of extreme values
    idx <- which(as.matrix(is_extr, wide = TRUE), arr.ind = TRUE)
    
    # Clamp to valid range
    row_idx <- pmax(1, pmin(length(y_vec), idx[, 1]))
    col_idx <- pmax(1, pmin(length(x_vec), idx[, 2]))
    
    # Return the coordinates
    data.frame(
      x = x_vec[col_idx],
      y = y_vec[row_idx]
    )
  }
  
  # Find local peaks and valleys
  peaks   <- find_local_extrema(z_mat, smoothed[["x"]], smoothed[["y"]], type = "peak", window)
  valleys <- find_local_extrema(z_mat, smoothed[["x"]], smoothed[["y"]], type = "valley", window)
  
  p <- plot_ly(
    x = smoothed[["x"]],
    y = smoothed[["y"]],
    z = z_mat,
    type = 'contour',
    contours = list(
      coloring = 'none',
      showlabels = FALSE,
      start = min(z_mat, na.rm = TRUE),
      end = max(z_mat, na.rm = TRUE),
      size = (max(z_mat, na.rm = TRUE) - min(z_mat, na.rm = TRUE)) / 10
    ),
    line = list(color = 'black', width = 1, dash='dot'),
    showscale = FALSE 
  ) 
  # Add the peaks as 'plusses'
  if (nrow(peaks) > 0) {
    for (k in 1:nrow(peaks)) {
      p <- add_trace(p,
                     inherit=FALSE,
                     x = peaks$x[k], y = peaks$y[k],
                     type = 'scatter', mode = 'text',
                     text = '+',
                     textfont = list(size = 24, color = 'black', family = 'Helvetica'),
                     line = list(width = 0),
                     showlegend = FALSE,
                     hoverinfo = 'none',
                     z = NULL,
                     contours = NULL,
                     showscale = NULL
      )
    }
  }

  # Add the valleys as 'minuses'
  if (nrow(valleys) > 0) {
    for (k in 1:nrow(valleys)) {
      p <- add_trace(p,
                     inherit=FALSE,
                     x = valleys$x[k], y = valleys$y[k],
                     type = 'scatter', mode = 'text',
                     text = '\u2212',
                     textfont = list(size = 24, color = 'black', family = 'Helvetica'),
                     line = list(width = 0),
                     showlegend = FALSE,
                     hoverinfo = 'none',
                     z          = NULL,
                     contours   = NULL,
                     showscale  = NULL
      )
    }
  }

  p <- p %>% layout(
    font=list(
      family="Helvetica",
      size=20,
      color="black"),
    xaxis = list(
      title = paste0("PC", pcx),
      showticklabels = FALSE,
      ticks = '',
      showline = TRUE,
      mirror = TRUE,
      zeroline = FALSE,
      showgrid = FALSE
    ),
    yaxis = list(
      title = paste0("PC", pcy),
      showticklabels = FALSE,
      ticks = '',
      showline = TRUE,
      mirror = TRUE,
      zeroline = FALSE,
      showgrid = FALSE
    ),
    plot_bgcolor  = 'white',
    paper_bgcolor = 'white'
  )
  return (p)
}

# Create a 3D rendering of the fitness landscape,
# with samples projected onto the surface
# smoothed: A list with x, y and z values
# pca_plot_df: Contains the PC coordinates of each sample
# families: a list of family names in the plot
# plot_inds: if TRUE, plots samples. If FALSE, plots adaptive walks
# founders: a list of each of the founder line names (only relevant for NAM)
# pcx: The principal component for the x-axis
# pcy: The principal component for the y-axis
render_3d_landscape <- function(smoothed,
                                df,
                                plot_inds=TRUE,
                                families,
                                pcx=1,
                                pcy=2,
                                founders=c()) {
  p <- plot_ly(
    x=smoothed[["x"]],
    y=smoothed[["y"]],
    z=smoothed[["z"]], # Transpose because smooth.2d() produces a transposed matrix
    type='surface',
    colorscale = 'Greys',
    #colors = viridis(50, alpha = 1, begin = 0, end = 1, direction = 1),
    opacity=1,
    colorbar=list(title="Suitability"),
    contours = list(
      x = list(show = TRUE, start = min(smoothed[["x"]]), 
               end = max(smoothed[["x"]]), size = 8,
               color = "black", width = 1, highlightcolor = "black"),
      y = list(show = TRUE, start = min(smoothed[["y"]]), 
               end = max(smoothed[["y"]]), size = 8,
               color = "black", width = 1, highlightcolor = "black")
    ))
  # Project each of the samples onto the surface (hovering slightly above)
  for (fam in families) {
    df_sub <- df[df$Family == fam, ]
    
    # Founders should be marked with a large diamond
    is_founder <- FALSE
    if (length(founders) > 0) {
      is_founder <- !grepl("_RIL$", fam)
    }
    
    if (plot_inds) {
      p <- add_trace(p,
                     data = df_sub,
                     x = df_sub[[pcy]], y = df_sub[[pcx]], z = df_sub$Suitability,
                     type = 'scatter3d', mode = 'markers',
                     name = fam,
                     marker = list(size = if (is_founder) 12 else 5,
                                   symbol = if (is_founder) "diamond" else "circle",
                                   color = df_sub$family_color
                                   #color = ~Suit,
                                   #colorscale = "Viridis",
                                   #showscale = TRUE,
                                   #colorbar = list(title = "Suitability", len = 0.5)
                     ),
                     showlegend = TRUE) # set to false if doing suitability colorscale
    } else { # Plot adaptive walk
      p <- add_trace(p,
                     data = df_sub,
                     x = df_sub[[pcy]], y = df_sub[[pcx]], z = df_sub$Suitability[,1],
                     type = 'scatter3d', mode = 'lines',
                     name = fam,
                     line = list(autocolorscale=FALSE, color=df_sub$family_color, which=2, width = 10),
                     showlegend = TRUE) 
    }
    p
  }
  return(p)
}

# Create a 3D rendering of the each sample, with its height determined by its suitability
# This is only tested with the sorghum NAM
# smoothed: The output of smooth.2d()
# df: Contains the PC coordinates of each sample
# families: a list of family names in the plot
# founders: a list of each of the founder line names
# pcx: The principal component for the x-axis
# pcy: The principal component for the y-axis
render_individuals <- function(df, families, founders, pcx=1, pcy=2) {
  p <- plot_ly()
  
  for (fam in families) {
    df_sub <- df[df$Family == fam, ]
    
    # Founders should be marked with a large diamond
    is_founder <- !grepl("_RIL$", fam)
    marker_symbol <- if (is_founder) "diamond" else "circle"
    marker_size   <- if (is_founder) 12 else 5
    
    p <- add_trace(p,
                   data = df_sub,
                   x = df_sub[[pcx]], y = df_sub[[pcy]], z = df_sub[["Suit"]],
                   type = 'scatter3d', mode = 'markers', # set to 'lines' for adaptive walk
                   name = fam,
                   marker = list(size = marker_size,
                                 symbol = marker_symbol,
                                 color = df_sub[["family_color"]]
                                 #color = ~Suit,
                                 #colorscale = "Viridis",
                                 #showscale = TRUE,
                                 #colorbar = list(title = "Suitability", len = 0.5)
                   ),
                   # UNCOMMENT THIS FOR DOING ADAPTIVE WALK
                   # line = list(width = 4, dash = 'solid'),
                   showlegend = TRUE) # set to false if doing suitabillity colorscale
  }
  
  p <- p %>% layout(legend = 
                      list(
                        font=list(
                          family="Helvetica",
                          size=20,
                          color="black"),
                        title=list(text="Family")),
                    scene = list(xaxis = list(title = paste0("PC", pcx)),
                                 yaxis = list(title = paste0("PC", pcy)),
                                 zaxis = list(title = "Suitability", showgrid=FALSE, zeroline=FALSE),
                                 aspectmode = 'cube'))
  return (p)
}


