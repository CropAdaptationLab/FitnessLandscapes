# Title: GENERATE SUBPOP PLOTS
# Author: Ted Monyak
# Description: After CreateIndependentPops.R, this script creates figures 
# reflecting the independent adaptive walks and allele frequencies

# Create allele frequenciy plots
if (saveAllelePlots) {
  for (p in 1:n.nPops) {
    pop <- pops[[p]]
    fit.df <- fit_dfs[[p]]
    subpop_dir <- subpop_dirs[[p]]
    
    # Make the dataframe tidy
    freq.df <- fit.df[-c(2:6)]
    freq.df<- melt(freq.df, id="gen", variable.name="id", value.name="freq")
    
    # Add the qtl effect size data to the dataframe
    freq.df <- merge(freq.df, qtlEff.df, by="id", all.x=TRUE)
    
    freq.df <- freq.df %>%
      separate(id, into = c("chr", "loc"), sep = "_", remove = FALSE) %>%
      mutate(
        chr = as.numeric(chr),
        loc = as.numeric(loc)
      ) %>%
      arrange(chr, loc)
    freq.df$id <- factor(freq.df$id, levels = unique(freq.df$id))
    
    theme <- theme_minimal(base_size = 10,
                           base_family="Helvetica") +
      theme(
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0, hjust=1),
        legend.text = element_text(),
        legend.title = element_text(),
        legend.key = element_rect(linewidth=0.05),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black"),
        panel.background = element_rect(fill = "white", color = "black"))
    
    
    # Create a line plot for the change in frequency of alleles over time
    # Each line's color is a function of its effect size
    if (p == 1) {
      colors_palette <- "Reds"
    } else {
      colors_palette <- "Blues"
    }
    ggplot(freq.df, aes(x=gen, y=freq, group=id)) +
      geom_line(aes(color=eff_size), linewidth=0.4, show.legend=TRUE) +
      scale_color_distiller(palette=colors_palette, direction=1, "Effect Size") +
      labs(x="Generation", y="Allele Frequency", title=paste0("Subpopulation ", p)) +
      theme
    ggplot2::ggsave(filename = "allelefrequencies.jpg",
                    path=subpop_dir,
                    device = "jpg",
                    height=2,
                    width=4,
                    dpi=600)
    ggplot2::ggsave(filename = "allelefrequencies.pdf",
                    path=subpop_dir,
                    device = "pdf",
                    height=2,
                    width=4)

    # Show the allele frequencies in heatmap form, in a constant order
    ggplot(freq.df, aes(x=gen, y=id, fill=freq)) +
      geom_tile() +
      scale_fill_distiller(palette=colors_palette, direction=1) +
      scale_y_discrete(limits=rev) +
      labs(x = "Generation", y = "QTL", fill="Frequency of\n'1' Allele") +
      theme_minimal()
    
    ggplot2::ggsave(filename = "allelefrequencies_heatmap.jpg",
                    path=subpop_dir,
                    device = "jpg",
                    width=10,
                    height=7,
                    dpi=600)

    # Join the allele frequencies of both subpopulations
    freq_1.df <- fit_dfs[[1]][-c(2:6)]
    freq_1.df <- freq_1.df[freq_1.df$gen==n.gens,-1]
    freq_1.df <- pivot_longer(freq_1.df, cols=everything(), names_to="id", values_to="freq")
    freq_1.df$pop <- 1
    freq_2.df <- fit_dfs[[2]][-c(2:6)]
    freq_2.df <- freq_2.df[freq_2.df$gen==n.gens,-1]
    freq_2.df <- pivot_longer(freq_2.df, cols=everything(), names_to="id", values_to="freq")
    freq_2.df$pop <- 2
    freq.df <- rbind(freq_1.df, freq_2.df)
    
    # Add a border color column
    freq.df$border_col <- "black"
    freq.df$border_col[freq.df$pop == 1]  <- "#CC000080"
    freq.df$border_col[freq.df$pop == 2] <- "#3C78D880"
    
    # Create a heatmap of the final allele frequencies of each subpopulation, comparing each QTL
    ggplot(freq.df, aes(x=id, y=pop, fill=freq)) +
      geom_tile(aes(color=border_col), linewidth=0.5) +
      scale_color_identity() +
      coord_fixed() +
      scale_fill_distiller(palette='Greys', direction=1) +
      theme +
      labs(x="QTLs", y="Subpopulation", fill="Frequency of\n'1' Allele",
           subtitle = paste0("Mean Isoeliteness = ", round(mean(isoElite_T1, isoElite_T2), 2))) +
      theme_minimal(base_size = 8,
                    base_family="Helvetica") +
      theme(axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            axis.title.x = element_text(),
            axis.title.y = element_text(angle=90, vjust=2),
            legend.position="bottom",
            legend.direction="horizontal",
            legend.text=element_text(),
            legend.title=element_text(),
            legend.key.height = unit(0.15, "cm"),
            legend.key.width = unit(0.4,  "cm"),
            legend.margin = margin(0, 0, 0, 0),
            plot.margin = margin(2, 2, 2, 2, "pt"),
            plot.subtitle=element_text(size=8, hjust=1))
            
    
    ggplot2::ggsave(filename = "allelefrequencies_contrast.jpg",
                    path=rep_dir,
                    device = "jpg",
                    height=1.3,
                    width=4,
                    dpi=600)
    ggplot2::ggsave(filename = "allelefrequencies_contrast.pdf",
                    path=rep_dir,
                    device = "pdf",
                    height=1.3,
                    width=4)
  }
}

# Create trait-to-fitness landscape plots
if (saveFitnessPlots) {
  for (p in 1:n.nPops) {
    pop <- pops[[p]]
    fit.df <- fit_dfs[[p]]
    subpop_dir <- subpop_dirs[[p]]
  
    adaptiveWalk <- plotAdaptiveWalk(fit.df)
    fname <- file.path(subpop_dir, "adaptivewalk.html")
    htmlwidgets::saveWidget(as_widget(adaptiveWalk), fname)
    
    # Plot the yield potential gain over generations
    g <- ggplot(fit.df, aes(x=gen, y=yieldPotential)) +
      geom_line()
    ggplot2::ggsave(filename = "yieldPotential.jpg",
                    path=subpop_dir,
                    device = "jpg",
                    width=10,
                    height=7)
    
    # Plot the suitability gain over generations
    g <- ggplot(fit.df, aes(x=gen, y=meanSuit)) +
      geom_line()
    ggplot2::ggsave(filename = "meanSuit.jpg",
                    path=subpop_dir,
                    device = "jpg",
                    width=10,
                    height=7)
    
    write.table(fit.df, file.path(subpop_dir, "fitness.csv"), col.names=NA, quote=FALSE, sep=",")
  }
  
  contour <- overlayWalkOnLandscape(fit_dfs[[1]],
                                    fit_dfs[[2]],
                                    type="CONTOUR",
                                    traitMin=-n.initTraitVal*1.1,
                                    traitMax=n.initTraitVal*1.1,
                                    popId_1="1",
                                    popId_2="2")
  fname <- file.path(rep_dir, "adaptivewalk_contour.html")
  htmlwidgets::saveWidget(as_widget(contour), fname)
  
  surface <- overlayWalkOnLandscape(fit_dfs[[1]],
                                    fit_dfs[[2]],
                                    type="SURFACE",
                                    suitFunc=suitFunc,
                                    traitMin=-n.initTraitVal-1,
                                    traitMax=n.initTraitVal+1,
                                    popId_1="1",
                                    popId_2="2")
  fname <- file.path(rep_dir, "adaptivewalk_surface.html")
  htmlwidgets::saveWidget(as_widget(surface), fname)
 
  # Plot the suitability of both subpopulations
  fig <- plot3dPopulationFitnessTwoPops(pops[[1]], pops[[2]], suitFunc)
  fname <- file.path(rep_dir, "3DFitness.html")
  htmlwidgets::saveWidget(as_widget(fig), fname)
}

# Plot the adaptive walks on a G > F landscape
if (sampleInds) {
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
  red_pal  <- scales::col_numeric(palette = c("#ffc9c9","#cc0000", "#a30000", "#7a0000"), domain = range(suit_df$gen)) 
  blue_pal <- scales::col_numeric(palette = c("#accbfc", "#0957d6", "#0748b2", "#063a8f"), domain = range(suit_df$gen))
  
  # Define a custom palette
  # stops: named numeric vector, e.g. c(0, 0.1, 0.5, 1)
  # colors: hex colors of same length
  custom_pal <- function(stops, colors) {
    ramp <- colorRamp(colors)
    function(x) {
      # Remap x_norm through your custom stop positions
      x_remapped <- approx(x = stops, y = seq(0, 1, length.out = length(stops)), xout = x)$y
      rgb_vals <- ramp(x_remapped)
      rgb(rgb_vals[,1], rgb_vals[,2], rgb_vals[,3], maxColorValue = 255)
    }
  }
  
  red_pal <- custom_pal(
    stops  = c(0, 0.5, 1)*n.gens,   # custom positions
    colors = c("#ffc9c9", "#cc0000", "#7a0000")
  )
  
  blue_pal <- custom_pal(
    stops  = c(0, 0.5, 1)*n.gens,   # custom positions
    colors = c("#accbfc", "#0957d6", "#063a8f")
  )
  
  suit_df <- suit_df %>%
    dplyr::rename(
      "Family" = "subpop"
    ) %>%
    dplyr::mutate(
      family_color= case_when(
        Family == "Founder" ~ "#B04EDE",
        Family == "Subpop. 1" ~ red_pal(gen),
        Family == "Subpop. 2" ~ blue_pal(gen)
      ))
  
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
  p
  htmlwidgets::saveWidget(as_widget(p), file.path(rep_dir, "adaptive_walk_genotype.html"))
}