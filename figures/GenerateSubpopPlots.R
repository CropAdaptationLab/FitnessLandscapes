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