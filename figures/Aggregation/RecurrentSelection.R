# Title: RECURRENT SELECTION
# Author: Ted Monyak
# Description: Aggregate results from the recurrent selection simulations

library(ggplot2)

# ASTHETIC ELEMENTS
theme <- theme_minimal(base_size = 8,
                       base_family="Helvetica") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(),
    panel.background = element_rect(fill = "white", color = "black"),
    plot.margin= unit(c(10,10,10,10), unit="pt"))

colors <- c(
  "Admixed GARS" = "gold2",
  "Admixed PRS" = "gold2",
  "Unadmixed GARS" = "#fa8034",
  "Unadmixed PRS" = "#fa8034"
)

shapes <- c(
  "Admixed GARS" = 16,
  "Admixed PRS" = 17,
  "Unadmixed GARS" = 16,
  "Unadmixed PRS" = 17
)

scale_color <- scale_color_manual(name = "QTL per\nAttained Trait",
                                  values = c("10" = "#4A1A6B",
                                             "20" = "#9B59B6",
                                             "50" = "#D7B8F3"))

# Create a unique identifier for the RILTYPE_SELECTION TYPE
res.df <- res.df %>%
  dplyr::mutate(pop=paste0(type, " ", sel))
res.df$pop <- factor(res.df$pop,
                     levels=c("Admixed GARS",
                              "Admixed PRS",
                              "Unadmixed GARS",
                              "Unadmixed PRS"))
res.df$qtl <- as.factor(as.character(res.df$qtl))

# Plot the average breeding fitness per cycle as line plots
meanWPerCycle <- function(df, nQtl) {
  df %>%
    dplyr::filter(qtl==nQtl) %>%
    ggplot(aes(x = c, y = meanW, color = pop, shape = pop, group = pop)) +
      geom_line() +
      geom_point() +
      scale_color_manual(
        name = "Population",
        values = colors,
      ) +
      scale_shape_manual(
        name = "Population",
        values = shapes,
      ) +
      labs(
        title  = paste("QTL =", nQtl),
        x = "Cycle",
        y = "Breeding Fitness",
        color = "Family"
      ) + 
      theme
}

# Calculate average breeding fitness per pop per cycle
meanW.df <-  res.df %>%
  dplyr::group_by(qtl, pop, c) %>%
  dplyr::summarize(meanW = mean(w))

meanW10 <- meanWPerCycle(meanW.df, "10")
meanW20 <- meanWPerCycle(meanW.df, "20")
meanW50 <- meanWPerCycle(meanW.df, "50")

(meanW10 | meanW20 | meanW50) + plot_layout(guides = "collect", axes = "collect")
ggplot2::ggsave(filename = "breedingFitness.jpg",
                path=output_dir,
                device = "jpg",
                width=6.5,
                height=2.5,
                dpi=600)
ggplot2::ggsave(filename = "breedingFitness.pdf",
                path=output_dir,
                device = "pdf",
                width=6.5,
                height=2.5)


# Plot the genetic gain per replicate as a boxplot
plotGeneticGain <- function(df, nQtl) {
  df %>%
    dplyr::filter(qtl==nQtl) %>%
    ggplot(aes(x = pop, y = geneticGain, fill=pop)) +
    geom_boxplot() +
    labs(
      title  = paste("QTL =", nQtl),
      x = "Group",
      y = "Genetic Gain"
    ) + 
    scale_fill_manual(
      values = colors
    ) +
    theme +
    theme(
      axis.text.x = element_blank(),
      axis.title.y = element_blank()
    )
}

# Calculate genetic gain per rep
geneticGain.df <- res.df %>%
  dplyr::group_by(qtl, pop, founder, rep) %>%
  dplyr::summarize(
    geneticGain = w[c == max(c)] - w[c == 1],
    isoElite = mean(isoElite),
    .groups = "drop"
  )

gain10 <- plotGeneticGain(geneticGain.df, 10)
gain20 <- plotGeneticGain(geneticGain.df, 20)
gain50 <- plotGeneticGain(geneticGain.df, 50)

(gain10 | gain20 | gain50) + plot_layout(guides = "collect", axes = "collect")
ggplot2::ggsave(filename = "geneticGain.jpg",
                path=output_dir,
                device = "jpg",
                width=6.5,
                height=2.5,
                dpi=600)
ggplot2::ggsave(filename = "geneticGain.pdf",
                path=output_dir,
                device = "pdf",
                width=6.5,
                height=2.5)

# Plot isoeliteness against genetic gain, color-coded by number of QTL
plotIeGain <- function(df) {
  df %>%
    dplyr::filter(pop=="Admixed GARS") %>%
    ggplot(aes(x = isoElite, y = geneticGain)) +
      geom_point(aes(color=qtl), alpha=0.6) +
      geom_smooth(method="lm", se=FALSE, color="black", linewidth=0.4) +
      stat_cor(method="pearson", 
               aes(label=paste(
                 sub("R", "r", after_stat(r.label)), 
                 after_stat(p.label), 
                 sep="~`,`~"
               ))) +
      guides(color = guide_legend(override.aes = list(size = 2))) +
      labs(
        x = "Isoeliteness of Parents",
        y = "Genetic Gain"
      ) + 
      scale_color +
      theme
}

plotIeGain(geneticGain.df)
ggplot2::ggsave(filename = "ieGain.jpg",
                path=output_dir,
                device = "jpg",
                width=6,
                height=4,
                dpi=600)
ggplot2::ggsave(filename = "ieGain.pdf",
                path=output_dir,
                device = "pdf",
                width=6,
                height=4)
