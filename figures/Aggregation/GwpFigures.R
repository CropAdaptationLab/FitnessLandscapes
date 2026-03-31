# Title: GWP Figures
# Author: Ted Monyak
# Description: This script produces summary plots for cross population prediction
# and recurrent selection with data generated in GwpPipeline.R

library(ggpubr)
library(grid)
library(patchwork)

theme <- theme_minimal(base_size = 8,
                       base_family="Helvetica") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(),
    panel.background = element_rect(fill = "white", color = "black"),
    plot.margin= unit(c(5,5,5,5), unit="pt"),
    legend.position = "right",
    legend.title.position = "top",
    legend.direction="vertical",
    aspect.ratio=1)

scale_fill_gwp <- scale_fill_manual(name = "Test Family",
                                values = c("Admixed" = "gold2",
                                           "Within-population" = "#CC0000",
                                           "Cross-population" = "#9000de"),
                                breaks = c("Within-population", "Admixed", "Cross-population"))

scale_color <- scale_color_manual(name = "QTL per\nAttained Trait",
                                  values = c("10" = "#4A1A6B",
                                             "20" = "#9B59B6",
                                             "50" = "#D7B8F3"))

type_colors <- c(
  "Admixed" = "gold2",
  "Unadmixed" = "#CC0000"
)

colors <- c(
  "Admixed GARS" = unname(type_colors[1]),
  "Admixed PRS" = unname(type_colors[1]),
  "Unadmixed GARS" = unname(type_colors[2]),
  "Unadmixed PRS" = unname(type_colors[2])
)

shapes <- c(
  "Admixed GARS" = 16,
  "Admixed PRS" = 17,
  "Unadmixed GARS" = 16,
  "Unadmixed PRS" = 17
)


# Determine a string representation of the correlation and significance
sig_cor <- stat_cor(method="pearson",
                    aes(label=paste(
                      sub("R", "r", after_stat(r.label)),
                      ifelse(after_stat(p) < 0.001, '"***"',
                             ifelse(after_stat(p) < 0.01,  '"**"',
                                    ifelse(after_stat(p) < 0.05,  '"*"', '"ns"'))),
                      sep="~"
                    )))

# Merge the within and across cases
gwp.df <- RIL.df %>%
  pivot_longer(
    cols = c(rAdmixed, rP1, rP2, rP1P2, rP2P1),
    names_to = "type",
    values_to = "r"
  ) %>%
  mutate(type = case_when(
    type %in% c("rP1", "rP2") ~ "Within-population",
    type %in% c("rP1P2", "rP2P1") ~ "Cross-population",
    type == "rAdmixed" ~ "Admixed",
    TRUE ~ type
  ))


gwp.df$type <- factor(gwp.df$type,
                      levels=c("Within-population", "Cross-population", "Admixed"))
gwp.df$qtl <- as.factor(gwp.df$qtl)

# Plot the GWP results by RIL family type
gwp.df %>%
  ggplot(aes(x=qtl, y=r, fill=type)) +
  geom_boxplot(
    outlier.shape=NA,
    position = position_dodge2(width = 1, padding = 0, preserve = "single"),
    linewidth = 0.1,
  ) +
  stat_compare_means(
    aes(group = type),
    label = "p.signif",
    hide.ns = FALSE,
  ) +
  scale_fill_gwp +
  scale_y_continuous(expand=c(0,0.1)) +
  labs(x="QTL per Attained Trait", y="GWP accuracy for breeding fitness (r)") +
  theme

ggplot2::ggsave(filename = "cross_pop.jpg",
                path=output_dir,
                device = "jpg",
                width=3.5,
                height=3,
                dpi=600)
ggplot2::ggsave(filename = "cross_pop.pdf",
                path=output_dir,
                device = "pdf",
                width=3.5,
                height=3)

minIe <- min(gwp.df$isoElite)
minR <- min(gwp.df$r)
maxR <- max(gwp.df$r)

# Plot isoeliteness against GWP accuracy
plot_RIL_IE <- function(df, nQtl, ylabel=TRUE) {
  df %>%
    dplyr::filter(type=="Admixed",
                  qtl==nQtl) %>%
    ggplot(aes(x=isoElite, y=r)) +
    geom_point(size=0.5) +
    geom_smooth(method="lm", se=FALSE, color="black", linewidth=0.4) +
    sig_cor +
    labs(x="Mean Isoeliteness",
         y="GWP accuracy for breeding fitness (r)",
         title=paste0("QTL: ", nQtl)) +
    scale_x_continuous(limits=c(minIe, 1)) +
    scale_y_continuous(limits=c(minR, maxR)) +
    theme +
    theme(
      axis.title.y = if (!ylabel) element_blank() else element_text(),
      axis.text.y = if (!ylabel) element_blank() else element_text()
    )
}


rilIe10 <- plot_RIL_IE(gwp.df, 10)
rilIe20 <- plot_RIL_IE(gwp.df, 20, ylabel=FALSE)
rilIe50 <- plot_RIL_IE(gwp.df, 50, ylabel=FALSE)
(rilIe10 | rilIe20 | rilIe50) + plot_layout(guides = "collect", axes = "collect")

ggplot2::ggsave(filename = "gwp_ie.jpg",
                path=output_dir,
                device = "jpg",
                width=6.5,
                height=2.5,
                dpi=600)
ggplot2::ggsave(filename = "gwp_ie.pdf",
                path=output_dir,
                device = "pdf",
                width=6.5,
                height=2.5)


# Create a unique identifier for the RILTYPE_SELECTION TYPE
gs.df <- RS.df %>%
  dplyr::mutate(pop=paste0(type, " ", sel))
gs.df$pop <- factor(gs.df$pop,
                     levels=c("Admixed GARS",
                              "Admixed PRS",
                              "Unadmixed GARS",
                              "Unadmixed PRS"))
gs.df$qtl <- as.factor(as.character(gs.df$qtl))

# Calculate average breeding fitness per pop per cycle
cycleMean.df <-  gs.df %>%
  dplyr::group_by(qtl, pop, c) %>%
  dplyr::summarize(w = mean(w),
                   genHt = mean(genome_het),
                   attHt = mean(attained_het),
                   desHt = mean(desired_het),
                   r = mean(r),
                   ie = mean(pop_ie))

minW <- min(cycleMean.df$w)
maxW <- max(cycleMean.df$w)

# Plot the average breeding fitness per cycle as line plots
meanWPerCycle <- function(df, nQtl, ylabel=TRUE) {
  df %>%
    dplyr::filter(qtl==nQtl) %>%
    ggplot(aes(x = c, y = w, color = pop, shape = pop, group = pop)) +
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
      title  = paste("QTL: ", nQtl),
      x = "Cycle",
      y = if (ylabel) "Breeding Fitness" else NULL,
      color = "Family"
    ) + 
    scale_y_continuous(limits=c(minW, maxW)) +
    theme +
    theme(
      axis.text.y = if (!ylabel) element_blank() else element_text()
    )
}

meanW10 <- meanWPerCycle(cycleMean.df, 10)
meanW20 <- meanWPerCycle(cycleMean.df, 20, ylabel=FALSE)
meanW50 <- meanWPerCycle(cycleMean.df, 50, ylabel=FALSE)

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

minHt <- min(c(cycleMean.df$genHt, cycleMean.df$attHt, cycleMean.df$desHt))
maxHt <- max(c(cycleMean.df$genHt, cycleMean.df$attHt, cycleMean.df$desHt))

ht.df <- cycleMean.df %>%
  pivot_longer(cols=c(genHt, attHt, desHt),
               names_to="qtlType",
               values_to="het")
ht.df$qtlType <- factor(ht.df$qtlType,
                        levels=c("attHt", "desHt", "genHt"))

# Plot the average heterozygosity per cycle as line plots
meanHetPerCycle <- function(df, nQtl, ylabel=TRUE) {
  df %>%
    dplyr::filter(pop=="Admixed GARS") %>%
    dplyr::filter(qtl==nQtl) %>%
    ggplot(aes(x = c, y = het, color = qtlType)) +
    geom_line() +
    geom_point() +
    labs(
      title  = paste("QTL: ", nQtl),
      x = "Cycle",
      y = if (ylabel) "Heterozygosity" else NULL
    ) + 
    scale_color_manual(name = "QTL",
                       values = c("genHt" = '#264653',
                                  "attHt" = '#2a9d8f',
                                  "desHt" = '#e9c46a'),
                       labels = c("genHt" = "Genomewide",
                                "attHt" = "Attained Trait",
                                "desHt" = "Desired Trait")) +
    scale_y_continuous(limits=c(minHt, maxHt)) +
    theme +
    theme(
      axis.text.y = if (!ylabel) element_blank() else element_text()
    )
}

ht10 <- meanHetPerCycle(ht.df, 10)
ht20 <- meanHetPerCycle(ht.df, 20, ylabel=FALSE)
ht50 <- meanHetPerCycle(ht.df, 50, ylabel=FALSE)
(ht10 | ht20 | ht50) + plot_layout(guides = "collect", axes = "collect")
ggplot2::ggsave(filename = "heterozygosity.jpg",
                path=output_dir,
                device = "jpg",
                width=6.5,
                height=2.5,
                dpi=600)
ggplot2::ggsave(filename = "heterozygosity.pdf",
                path=output_dir,
                device = "pdf",
                width=6.5,
                height=2.5)


minPopIe <- min(cycleMean.df$ie)
maxPopIe <- max(cycleMean.df$ie)

# Plot population isoeliteness over cycles
meanIePerCycle <- function(df, nQtl, ylabel=TRUE) {
  df %>%
    dplyr::filter(qtl==nQtl) %>%
    ggplot(aes(x = c, y = ie, color = pop, shape = pop, group = pop)) +
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
      title  = paste("QTL: ", nQtl),
      x = "Cycle",
      y = if (ylabel) "Population Isoeliteness" else NULL,
      color = "Family"
    ) + 
    scale_y_continuous(limits=c(minPopIe, maxPopIe)) +
    theme +
    theme(
      axis.text.y = if (!ylabel) element_blank() else element_text()
    )
}

meanIe10 <- meanIePerCycle(cycleMean.df, 10)
meanIe20 <- meanIePerCycle(cycleMean.df, 20, ylabel=FALSE)
meanIe50 <- meanIePerCycle(cycleMean.df, 50, ylabel=FALSE)

(meanIe10 | meanIe20 | meanIe50) + plot_layout(guides = "collect", axes = "collect")
ggplot2::ggsave(filename = "isoeliteness.jpg",
                path=output_dir,
                device = "jpg",
                width=6.5,
                height=2.5,
                dpi=600)
ggplot2::ggsave(filename = "isoeliteness.pdf",
                path=output_dir,
                device = "pdf",
                width=6.5,
                height=2.5)

minRsR <- min(gs.df$r, na.rm=TRUE)
maxRsR <- max(gs.df$r, na.rm=TRUE)

# Plot R over cycles
cycleMean.df %>%
  dplyr::filter(pop=="Admixed GARS") %>%
  ggplot(aes(x = c, y = r, color = qtl)) +
  geom_line() +
  geom_point() +
  labs(
    title  = paste("QTL: ", nQtl),
    x = "Cycle",
    y = "GWP accuracy for breeding fitness (r)"
  ) +
  scale_color +
  theme
ggplot2::ggsave(filename = "gwpAccuracy.jpg",
                path=output_dir,
                device = "jpg",
                width=3.5,
                height=3,
                dpi=600)
ggplot2::ggsave(filename = "gwpAccuracy.pdf",
                path=output_dir,
                device = "pdf",
                width=3.5,
                height=3)

# Plot population isoeliteness against GWP R
popIeR <- function(df, nQtl, ylabel=TRUE) {
  df %>%
    dplyr::filter(pop=="Admixed GARS",
                  qtl==nQtl) %>%
    ggplot(aes(x=pop_ie, y=r)) +
    geom_point(size=0.5) +
    geom_smooth(method="lm", se=FALSE, color="black", linewidth=0.4) +
    sig_cor +
    labs(x="Population Isoeliteness",
         y="GWP accuracy for breeding fitness (r)",
         title=paste0("QTL: ", nQtl)) +
    scale_x_continuous(limits=c(minPopIe, 1)) +
    scale_y_continuous(limits=c(minRsR, maxRsR)) +
    theme +
    theme(
      axis.title.y = if (!ylabel) element_blank() else element_text(),
      axis.text.y = if (!ylabel) element_blank() else element_text()
    )
}

popIeR10 <- popIeR(gs.df, 10)
popIeR20 <- popIeR(gs.df, 20, ylabel=FALSE)
popIeR50 <- popIeR(gs.df, 50, ylabel=FALSE)

(popIeR10 | popIeR20 | popIeR50) + plot_layout(guides = "collect", axes = "collect")
ggplot2::ggsave(filename = "popIeR.jpg",
                path=output_dir,
                device = "jpg",
                width=6.5,
                height=2.5,
                dpi=600)
ggplot2::ggsave(filename = "popIeR.pdf",
                path=output_dir,
                device = "pdf",
                width=6.5,
                height=2.5)

minAttHet <- min(gs.df$attained_het, na.rm=TRUE)
maxAttHet <- max(gs.df$attained_het, na.rm=TRUE)

# Plot attained trait heterozygosity against GWP R
attHetR <- function(df, nQtl, ylabel=TRUE) {
  df %>%
    dplyr::filter(pop=="Admixed GARS",
                  qtl==nQtl) %>%
    ggplot(aes(x=attained_het, y=r)) +
    geom_point(size=0.5) +
    geom_smooth(method="lm", se=FALSE, color="black", linewidth=0.4) +
    sig_cor +
    labs(x="Attained Trait Heterozygosity",
         y="GWP accuracy for breeding fitness (r)",
         title=paste0("QTL: ", nQtl)) +
    scale_x_continuous(limits=c(minAttHet, maxAttHet)) +
    scale_y_continuous(limits=c(minRsR, maxRsR)) +
    theme +
    theme(
      axis.title.y = if (!ylabel) element_blank() else element_text(),
      axis.text.y = if (!ylabel) element_blank() else element_text()
    )
}

hetR10 <- attHetR(gs.df, 10)
hetR20 <- attHetR(gs.df, 20, ylabel=FALSE)
hetR50 <- attHetR(gs.df, 50, ylabel=FALSE)

(hetR10 | hetR20 | hetR50) + plot_layout(guides = "collect", axes = "collect")
ggplot2::ggsave(filename = "hetR.jpg",
                path=output_dir,
                device = "jpg",
                width=6.5,
                height=2.5,
                dpi=600)
ggplot2::ggsave(filename = "hetR.pdf",
                path=output_dir,
                device = "pdf",
                width=6.5,
                height=2.5)

# Calculate genetic gain per rep
geneticGain.df <- gs.df %>%
  dplyr::group_by(qtl, founder, rep, type, sel) %>%
  dplyr::summarize(
    gain = w[c == max(c)] - w[c == 1],
    isoElite = mean(isoElite),
    .groups = "drop"
  )

maxGain <- max(geneticGain.df$gain)

# Plot the genetic gain per replicate as a boxplot
plotGeneticGain <- function(df, nQtl, ylabel=TRUE) {
  df %>%
    dplyr::filter(qtl==nQtl) %>%
    ggplot(aes(x = sel, y = gain, fill=type)) +
    geom_boxplot(
      outlier.shape=NA,
      position = position_dodge2(width = 1, padding = 0, preserve = "single"),
      linewidth = 0.1,
    ) +
    labs(
      title  = paste("QTL: ", nQtl),
      x = "Group",
      y = "Genetic Gain"
    ) + 
    scale_fill_manual(
      name = "RIL Family",
      values = type_colors
    ) +
    scale_y_continuous(limits=c(0, maxGain)) +
    theme +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = if (!ylabel) element_blank() else element_text(),
      axis.text.y = if (!ylabel) element_blank() else element_text()
    )
}


gain10 <- plotGeneticGain(geneticGain.df, 10)
gain20 <- plotGeneticGain(geneticGain.df, 20, ylabel=FALSE)
gain50 <- plotGeneticGain(geneticGain.df, 50, ylabel=FALSE)

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

minGain <- min(geneticGain.df$gain)
maxGain <- max(geneticGain.df$gain)

# Plot isoeliteness against genetic gain
plotIeGain <- function(df, nQtl, ylabel=TRUE) {
  df %>%
    dplyr::filter(qtl==nQtl) %>%
    ggplot(aes(x = isoElite, y = gain)) +
    geom_point(size=0.5) +
    geom_smooth(method="lm", se=FALSE, color="black", linewidth=0.4) +
    sig_cor +
    labs(
      title  = paste("QTL: ", nQtl),
      x = "Biparental Cross Isoeliteness",
      y = "Genetic Gain"
    ) + 
    scale_x_continuous(limits=c(minIe, 1)) +
    scale_y_continuous(limits=c(minGain, maxGain)) +
    theme +
    theme(
      axis.title.y = if (!ylabel) element_blank() else element_text(),
      axis.text.y = if (!ylabel) element_blank() else element_text()
    )
}

garsIe10 <- plotIeGain(geneticGain.df, 10)
garsIe20 <- plotIeGain(geneticGain.df, 20, FALSE)
garsIe50 <- plotIeGain(geneticGain.df, 50, FALSE)

(garsIe10 | garsIe20 | garsIe50) + plot_layout(guides = "collect", axes = "collect")

ggplot2::ggsave(filename = "ieGain.jpg",
                path=output_dir,
                device = "jpg",
                width=6.5,
                height=2.5,
                dpi=600)
ggplot2::ggsave(filename = "ieGain.pdf",
                path=output_dir,
                device = "pdf",
                width=6.5,
                height=2.5)
