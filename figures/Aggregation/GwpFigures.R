# Title: GWP Figures
# Author: Ted Monyak
# Description: This script produces summary plots for cross population prediction
# and recurrent selection with data generated in GwpPipeline.R

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
                                           "Unadmixed 1" = "#CC0000",
                                           "Unadmixed 2" = "#3C78D8"),
                                breaks = c("Unadmixed 1", "Unadmixed 2", "Admixed"))

scale_color <- scale_color_manual(name = "QTL per\nAttained Trait",
                                  values = c("10" = "#4A1A6B",
                                             "20" = "#9B59B6",
                                             "50" = "#D7B8F3"))

type_colors <- c(
  "Admixed" = "gold2",
  "Unadmixed" = "#fa8034"
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


# Pull out the admixed RIL family results
gwp_admixed <- res.df %>%
  dplyr::filter(type=="Admixed") %>%
  dplyr::select(isoElite, gwpR, qtl, type)

# Ungroup the data from both unadmixed families and designate them by the
# subpopulation from which they were derived
gwp_p1 <- res.df %>%
  dplyr::filter(type=="Unadmixed") %>%
  dplyr::select(isoElite, gwpR, qtl) %>%
  dplyr::mutate(type="Unadmixed 1")

gwp_p2 <- res.df %>%
  dplyr::filter(type=="Unadmixed") %>%
  dplyr::select(isoElite, gwpR_P2, qtl) %>%
  dplyr::rename(
    "gwpR"=gwpR_P2
  ) %>%
  dplyr::mutate(type="Unadmixed 2")

# Merge data from all RIL families
gwp.df <- dplyr::bind_rows(gwp_admixed, gwp_p1, gwp_p2) %>%
  drop_na()

gwp.df$type <- factor(gwp.df$type,
                      levels=c("Unadmixed 1", "Unadmixed 2", "Admixed"))
gwp.df$qtl <- as.factor(gwp.df$qtl)

# Test population should only have 2 categories
gwp.df <- gwp.df %>%
  dplyr::mutate(test_pop = ifelse(type=="Admixed", "Admixed", "Unadmixed"))

# Plot the GWP results by RIL family type
gwp.df %>%
  ggplot(aes(x=qtl, y=gwpR, fill=type)) +
  geom_boxplot(
    outlier.shape=NA,
    position = position_dodge2(width = 1, padding = 0, preserve = "single"),
    linewidth = 0.1,
  ) +
  stat_compare_means(
    aes(group = test_pop),
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

# Plot isoeliteness against GWP accuracy
gIE <- gwp.df %>%
  dplyr::filter(type=="Admixed") %>%
  ggplot(aes(x=isoElite, y=gwpR)) +
  geom_point(aes(color=qtl), alpha=0.6, size=0.5) +
  geom_smooth(method="lm", se=FALSE, color="black", linewidth=0.4) +
  sig_cor +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  scale_color +
  labs(x="Mean Isoeliteness", y="GWP accuracy for breeding fitness (r)") +
  theme
gIE

ggplot2::ggsave(filename = "gwp_ie.jpg",
                path=output_dir,
                device = "jpg",
                width=3.5,
                height=3,
                dpi=600)
ggplot2::ggsave(filename = "gwp_ie.pdf",
                path=output_dir,
                device = "pdf",
                width=3.5,
                height=3)


# Create a unique identifier for the RILTYPE_SELECTION TYPE
gs.df <- res.df %>%
  dplyr::mutate(pop=paste0(type, " ", sel))
gs.df$pop <- factor(gs.df$pop,
                     levels=c("Admixed GARS",
                              "Admixed PRS",
                              "Unadmixed GARS",
                              "Unadmixed PRS"))
gs.df$qtl <- as.factor(as.character(gs.df$qtl))

# Calculate average breeding fitness per pop per cycle
meanW.df <-  gs.df %>%
  dplyr::group_by(qtl, pop, c) %>%
  dplyr::summarize(meanW = mean(w))

minW <- min(meanW.df$meanW)
maxW <- max(meanW.df$meanW)

# Plot the average breeding fitness per cycle as line plots
meanWPerCycle <- function(df, nQtl, ylabel=TRUE) {
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
      y = if (ylabel) "Breeding Fitness" else NULL,
      color = "Family"
    ) + 
    scale_y_continuous(limits=c(minW, maxW)) +
    theme +
    theme(
      axis.text.y = if (!ylabel) element_blank() else element_text()
    )
}

meanW10 <- meanWPerCycle(meanW.df, "10")
meanW20 <- meanWPerCycle(meanW.df, "20", ylabel=FALSE)
meanW50 <- meanWPerCycle(meanW.df, "50", ylabel=FALSE)

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
      title  = paste("QTL =", nQtl),
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

# Plot isoeliteness against genetic gain, color-coded by number of QTL
plotIeGain <- function(df, ril, selType) {
  df %>%
    dplyr::filter(type==ril,
                  sel==selType) %>%
    ggplot(aes(x = isoElite, y = gain)) +
    geom_point(aes(color=qtl), alpha=0.6, size=0.5) +
    guides(color = guide_legend(override.aes = list(size = 2))) +
    labs(
      x = "Biparental Cross Isoeliteness",
      y = "Genetic Gain"
    ) + 
    scale_color +
    theme
}

garsIe <- plotIeGain(geneticGain.df, "Admixed", "GARS")
garsIe
ggplot2::ggsave(filename = "ieGain.jpg",
                path=output_dir,
                device = "jpg",
                width=3.5,
                height=3,
                dpi=600)
ggplot2::ggsave(filename = "ieGain.pdf",
                path=output_dir,
                device = "pdf",
                width=3.5,
                height=3)
