# Title: Genotype-to-Fitness Landscapes
# Author: Ted Monyak
# Description: This file contains plotting functions for creating genotype to
# fitness landscapes

library(fields)
library(ggdensity)
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
plot1DLandscape <- function(RIL, pop1, pop2, parent1, parent2, qtl1, qtl2, ril_dir) {
  # Make a larger RIL to reduce noise
  RIL <- self(RIL, nProgeny=5)
  
  # Get the genotypes for the attained traits only
  RIL_geno <- getUniqueQtl(RIL, traits=c(1,2))
  subpops_geno <- rbind(getUniqueQtl(pop1, traits=c(1,2)),
                        getUniqueQtl(pop2, traits=c(1,2)))
  
  
  # Determine the haplotypes of parent1 and parent 2 w.r.t. the two QTL markers
  haplo1 <- getUniqueQtl(parent1) %>%
    dplyr::select(all_of(c(qtl1, qtl2)))
  haplo2 <- getUniqueQtl(parent2) %>%
    dplyr::select(all_of(c(qtl1, qtl2)))
  
  # Determine whether a haplotype corresponds to parent 1, parent 2, or neither (recombinant)
  # allele1: The allele at the first locus
  # allele2: The allele at the second locus
  categorize <- function(allele1, allele2) {
    if (all(c(allele1, allele2) == haplo1)) {
      return("P1")
    } else if (all(c(allele1, allele2) == haplo2)) {
      return("P2")
    }
    return("R")
  }

  # Determine the parental haplotypes of each individual in the RIL family  
  haplo <- RIL_geno %>%
    dplyr::select(all_of(c(qtl1,qtl2))) %>%
    setNames(c("qtl1", "qtl2")) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(cat=categorize(qtl1, qtl2)) %>%
    dplyr::ungroup()
  
  # Calculate the breeding fitness of each RIL sample
  RIL_pheno.df <- data.frame(fitness=breedingFitness(pheno(RIL)),
                             haplo=haplo$cat)
  
  # Run PCA on the landrace subpopulations
  pca <- prcomp(subpops_geno)
  VAF <- summary(pca)$importance[2,1] * 100
  
  # Project the admixed RIL family onto the principal components
  pc_RIL <- predict(pca, RIL_geno)
  
  # Creat a dataframe for plotting
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

# Create a smoothed landscape by iteratively running a gaussian filter
# suit_df: A datafame containing the principal component coordinates
# and 'Suit', suitability values for each sample
# pcx: The principal component for the x-axis
# pcy: The principal component for the y-axis
# thetas: A list of kernel sizes for each iteration of gaussian smoothing
generate_landscape <- function(suit_df, pcx=1, pcy=2, thetas=c(1:14)) {
  # Obtain the initial smoothed image
  z <- suit_df$Suit
  xy <- cbind(suit_df[,pcx], suit_df[,pcy])
  smoothed <- smooth.2d(Y = z, x = xy, theta = thetas[1])
  
  # Iterate through each iteration of the gaussian smoothing
  for (i in 2:(length(thetas))) {
    
    # Convert the smooth.2d result into a datafame to filter out Suit values
    # that are infinite, or fall outside the 0 to 1 range
    smoothed_df <- as.data.frame(smoothed$z) %>%
      setNames(smoothed$y) %>%
      mutate(x = smoothed$x) %>%
      pivot_longer(cols = -x, 
                   names_to = "y", 
                   values_to = "Suit") %>%
      mutate(y = as.numeric(y)) %>%
      rename(PCX = x, PCY = y) %>%
      drop_na(Suit) %>%
      filter(is.finite(Suit), Suit <= 1, Suit >= 0)
    
    smoothed <- smooth.2d(Y = smoothed_df$Suit, 
                          x = cbind(smoothed_df$PCX, smoothed_df$PCY),
                          theta = thetas[i])
  }
  return(smoothed)
}

# Render a 2d image of the smoothed G > F landscape
# smoothed: The result of a smooth.2d() call
# VAF: the variance explained by each PC
# pcx: The principal component for the x-axis
# pcy: The principal component for the y-axis
render_2d_landscape <- function(smoothed, VAF, pcx=1, pcy=1) {
  # Clip the smoothed result between 0 and 1
  smoothed_df <- as.data.frame(smoothed$z) %>%
    setNames(smoothed$y) %>%
    mutate(x = smoothed$x) %>%
    pivot_longer(cols = -x, names_to = "y", values_to = "Suit") %>%
    mutate(y = as.numeric(y)) %>%
    rename(PCX = x, PCY = y) %>%
    filter(Suit >= 0, Suit <= 1)
  
  # Plot a 2-dimensional version of the genotype-to-fitness landscape
  ggplot(smoothed_df, aes(x = PCX, y = PCY)) +
    geom_point(aes(colour = Suit)) +
    scale_color_viridis(alpha=0.9, name="Suitability") +
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
      size = (max(z_mat, na.rm = TRUE) - min(z_mat, na.rm = TRUE)) / 8
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
                     textfont = list(size = 16, color = 'black', family = 'Arial'),
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
                     textfont = list(size = 16, color = 'black', family = 'Arial'),
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
    xaxis = list(
      title = '',
      showticklabels = FALSE,
      ticks = '',
      showline = TRUE,
      mirror = TRUE,
      zeroline = FALSE,
      showgrid = FALSE
    ),
    yaxis = list(
      title = '',
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

# Create a 3D rendering of the fitness landscape, with samples projected onto the surface
# This is only tested with the sorghum NAM
# smoothed: The output of smooth.2d()
# pca_plot_df: Contains the PC coordinates of each sample
# families: a list of family names in the plot
# founders: a list of each of the founder line names
# pcx: The principal component for the x-axis
# pcy: The principal component for the y-axis
render_3d_landscape <- function(smoothed,
                                pca_plot_df,
                                families,
                                founders,
                                pcx=1,
                                pcy=2) {
  p <- plot_ly(
    x=smoothed[["x"]],
    y=smoothed[["y"]],
    z=t(smoothed[["z"]]), # Transpose because smooth.2d() produces a transposed matrix
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
    df_sub <- pca_plot_df[pca_plot_df$Family == fam, ]
    
    # Founders should be marked with a large diamond
    is_founder <- !grepl("_RIL$", fam)
    marker_symbol <- if (is_founder) "diamond" else "circle"
    marker_size   <- if (is_founder) 12 else 5

    p <- add_trace(p,
                   data = df_sub,
                   x = df_sub[[pcx]], y = df_sub[[pcy]], z = ~Suitability,
                   type = 'scatter3d', mode = 'markers', # set to 'lines' for adaptive walk
                   name = fam,
                   marker = list(size = marker_size,
                                 symbol = marker_symbol,
                                 color = ~family_color
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
                scene = list(xaxis = list(title = paste0("PC", pcx), 
                                          showgrid = FALSE,
                                          zeroline = FALSE,
                                          showticklabels = FALSE,
                                          showline = FALSE),
                             yaxis = list(title = paste0("PC", pcy), 
                                          showgrid = FALSE,
                                          zeroline = FALSE,
                                          showticklabels = FALSE,
                                          showline = FALSE),
                             zaxis = list(title = "Suitability", showgrid=FALSE, zeroline=FALSE),
                             aspectmode = 'cube'))
  return (p)
}

# Create a 3D rendering of the each sample, with its height determined by its suitability
# This is only tested with the sorghum NAM
# smoothed: The output of smooth.2d()
# pca_plot_df: Contains the PC coordinates of each sample
# families: a list of family names in the plot
# founders: a list of each of the founder line names
# pcx: The principal component for the x-axis
# pcy: The principal component for the y-axis
render_individuals <- function(pca_plot_df, families, founders, pcx=1, pcy=2) {
  p <- plot_ly()
  
  for (fam in families) {
    df_sub <- pca_plot_df[pca_plot_df$Family == fam, ]
    
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
