library(fields)
library(terra)

plot_NAM_suit_distribution <- function(pca_plot_df) {
  theme <- theme_minimal(base_family = "Helvetica", base_size = 8) +
    theme(
      axis.ticks.x     = element_blank(),
      axis.text.y = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", color = "black"),
    )
  
  family_colors <- setNames(
    pca_plot_df$family_color,
    pca_plot_df$Family
  )
  pca_plot_df %>%
    dplyr::filter(grepl("_RIL", Family)) %>%
    ggplot(aes(x=Suit, color=Family)) +
    geom_density(size=0.4) +
    scale_color_manual(values=c(family_colors),
                       aesthetics=c("color")) +
    labs(x = "Suitability",
         y = "Density") +
    theme
}

plot_NAM_trait_distribution <- function(pca_plot_df) {
  df_long <- pca_plot_df %>%
    dplyr::filter(grepl("_RIL", Family)) %>%
    dplyr::select(Family, family_color, FT, PH) %>%
    dplyr::mutate(Family = sub("_.*", "", Family))
  
  family_colors <- setNames(
    df_long$family_color,
    df_long$Family
  )
  
  theme <- theme_minimal(base_family = "Helvetica", base_size = 8) +
    theme(
      axis.text.x      = element_blank(),
      axis.ticks.x     = element_blank(),
      axis.title.x     = element_blank(),panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", color = "black"),
      legend.position  = "none"
    )

  p_ft <- df_long %>%
    ggplot(aes(x = Family, y = FT, fill = Family)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_hline(aes(yintercept = 72.5, linetype = 'RTx430'),
               color = 'black', linewidth = 0.4) +
    scale_fill_manual(values = family_colors) +
    scale_linetype_manual(name = NULL, values = c("RTx430" = "dashed")) +
    labs(y = "Flowering Time (days)", fill = "Family") +
    theme

  p_ph <- df_long %>%
    ggplot(aes(x = Family, y = PH, fill = Family)) +
    geom_violin(trim = FALSE, alpha = 0.7) +
    geom_hline(aes(yintercept = 106.5, linetype = 'RTx430'),
               color = 'black', linewidth = 0.4) +
    scale_fill_manual(values = family_colors) +
    scale_linetype_manual(name = NULL, values = c("RTx430" = "dashed")) +
    labs(y = "Plant Height (cm)", fill = "Family") +
    theme +
    theme(legend.position = "right") +
    guides(
      fill     = guide_legend(order = 1),   # Family on top
      linetype = guide_legend(order = 2)    # Optimal below
    )
  
  p_ft + p_ph
}

generate_landscape <- function(suit_df, pcx=1, pcy=2, thetas=c(1:10)) {
  z <- suit_df$Suit
  xy <- cbind(suit_df[,pcx], suit_df[,pcy])
  
  smoothed <- smooth.2d(Y = z, x = xy, theta = thetas[1])
  
  for (i in 2:(length(thetas)-1)) {
    
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

render_2d_landscape <- function(smoothed, VAF, pcx=1, pcy=1) {
  smoothed_df <- as.data.frame(smoothed$z) %>%
    setNames(smoothed$y) %>%
    mutate(x = smoothed$x) %>%
    pivot_longer(cols = -x, names_to = "y", values_to = "Suit") %>%
    mutate(y = as.numeric(y)) %>%
    rename(PCX = x, PCY = y) %>%
    filter(Suit >= 0, Suit <= 1)
  
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

wright_landscape <- function(smoothed, pcx, pcy, window=5) {
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

render_3d_landscape <- function(smoothed, pca_plot_df, families, founders, pcx=1, pcy=2) {
  p <- plot_ly(
    x=smoothed[["x"]],
    y=smoothed[["y"]],
    z=t(smoothed[["z"]]),
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
  
  for (fam in families) {
    df_sub <- pca_plot_df[pca_plot_df$Family == fam, ]
    
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

render_individuals <- function(smoothed, pca_plot_df, families, founders, pcx=1, pcy=2) {
  p <- plot_ly()
  
  for (fam in families) {
    df_sub <- pca_plot_df[pca_plot_df$Family == fam, ]
    
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
