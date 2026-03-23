# Title: NAM PCA
# Author: Ted Monyak
# Description: Contains plotting functions for the NAM phenotypes and PCA projections

# Plots the distribution of suitability values for the NAM
plot_NAM_suit_distribution <- function(pca_plot_df) {
  theme <- theme_minimal(base_family = "Helvetica", base_size = 8) +
    theme(
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", color = "black"),
    )
  
  # Get the family-specific colors
  family_colors <- setNames(
    pca_plot_df$family_color,
    pca_plot_df$Family
  )
  # Create a density plot
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

# Plot the distribution of trait values for the NAM
# Creates a separate plot for each phenotype (flowering time and plant height)
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
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", color = "black"),
      legend.position = "none"
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


# Plot the NAM samples, color-coded by family
# pca_plot_df: contains the principal component for each sample, family_color, and Family
# VAF: the variance explained by each PC
# pcx: The principal component for the x-axis
# pcy: The principal component for the y-axis
plot_NAM_family <- function(pca_plot_df, VAF, pcx=1, pcy=2) {
  # Get the column names of each PC
  xcol <- colnames(pca_plot_df)[pcx]
  ycol <- colnames(pca_plot_df)[pcy]
  
  ggplot(pca_plot_df, aes(x = .data[[xcol]], y = .data[[ycol]], color = family_color)) +
    geom_point(data = subset(pca_plot_df, grepl("_RIL$", Family)),
               size = 0.5, shape = 16) +
    geom_point(data = subset(pca_plot_df, !grepl("_RIL$", Family)),
               size = 3, shape = 18) +
    scale_color_identity(guide="legend",
                         breaks = pca_plot_df$family_color,
                         labels = pca_plot_df$Family,
                         name="Family") +
    guides(color = guide_legend(override.aes = list(size = 2))) +
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

# Plot the NAM samples, color-coded by suitability
# pca_plot_df: contains the principal component for each sample, Family, and 'Suit'
# VAF: the variance explained by each PC
# pcx: The principal component for the x-axis
# pcy: The principal component for the y-axis
plot_NAM_suit <- function(pca_plot_df, VAF, pcx=1, pcy=2) {
  # Plot PCA with colors by group
  xcol <- colnames(pca_plot_df)[pcx]
  ycol <- colnames(pca_plot_df)[pcy]
  
  ggplot(pca_plot_df, aes(x = .data[[xcol]], y = .data[[ycol]])) +
    geom_point(data = subset(pca_plot_df, grepl("_RIL$", Family)),
               size = 0.5, shape = 16,
               aes(colour = Suit, alpha=0.9), show.legend = c(color = TRUE, size = FALSE, alpha = FALSE)) +
    geom_point(data = subset(pca_plot_df, !grepl("_RIL$", Family)),
               size = 3, shape = 18,
               aes(colour = Suit, alpha=0.9), show.legend = c(color = TRUE, size = FALSE, alpha = FALSE)) +
    scale_color_viridis(name="Suitability") +
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