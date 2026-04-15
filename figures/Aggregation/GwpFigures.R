# Title: GWP Figures
# Author: Ted Monyak
# Description: This script produces summary plots for cross population prediction
# and recurrent selection with data generated in GwpPipeline.R

library(dplyr)
library(ggplot2)
library(ggpubr)
library(grid)
library(patchwork)
library(tibble)
library(tidyr)

#setwd("~/Documents/CSU/FitnessLandscapes/output/GWP/current_best")
#RIL.df <- rbind(read.csv("QTL_10/ril_results.csv"),
#                read.csv("QTL_20/ril_results.csv"),
#                read.csv("QTL_50/ril_results.csv"))
#GS.df <- rbind(read.csv("QTL_10/rs_results.csv"),
#               read.csv("QTL_20/rs_results.csv"),
#               read.csv("QTL_50/rs_results.csv"))

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
                                values = c("Admixed" = "#E8B84B",
                                           "Within-population" = "#d5a6bd",
                                           "Cross-population" = "#8e7cc3"),
                                breaks = c("Within-population", "Admixed", "Cross-population"))

scale_color <- scale_color_manual(name = "QTL per\nAttained Trait",
                                  values = c("2" = "black",
                                             "10" = "#8F7511",
                                             "20" = "#E69900",
                                             "50" = "gold"))

# Determine a string representation of the correlation and significance
sig_cor <- stat_cor(method="pearson",
                    aes(label=paste(
                      sub("R", "r", after_stat(r.label)),
                      ifelse(after_stat(p) < 0.001, '"***"',
                             ifelse(after_stat(p) < 0.01,  '"**"',
                                    ifelse(after_stat(p) < 0.05,  '"*"', '"ns"'))),
                      sep="~"
                    )),
                    show.legend=FALSE)

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
                      levels=c("Within-population", "Admixed", "Cross-population"))
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

# Plot the GWP results by RIL family type for QTL == 10
gwp.df %>%
  dplyr::filter(qtl==10) %>%
  ggplot(aes(x=type, y=r, fill=type)) +
  geom_boxplot(
    outlier.shape=NA,
    linewidth = 0.1,
  ) +
  stat_compare_means(
    label = "p.signif",
    hide.ns = FALSE,
  ) +
  scale_fill_gwp +
  scale_y_continuous(expand=c(0,0.1)) +
  labs(x="QTL per Attained Trait", y="GWP accuracy for breeding fitness (r)") +
  theme +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())

ggplot2::ggsave(filename = "cross_pop_qtl10.jpg",
                path=output_dir,
                device = "jpg",
                width=3.5,
                height=3,
                dpi=600)
ggplot2::ggsave(filename = "cross_pop_qtl10.pdf",
                path=output_dir,
                device = "pdf",
                width=3.5,
                height=3)

minIe <- min(gwp.df$isoElite)
minMaxR <- gwp.df %>%
  dplyr::filter(type=="Admixed") %>%
  dplyr::summarize(min=min(r),
                   max=max(r))
minR <- minMaxR$min
maxR <- minMaxR$max

# Plot isoeliteness against GWP accuracy
gwp.df %>%
  dplyr::filter(type=="Admixed") %>%
  ggplot(aes(x=isoElite, y=r, color=qtl)) +
  geom_point(size=0.5) + #, aes(color=isoEliteDes)) +
  geom_smooth(method="lm", se=FALSE, aes(color=qtl), linewidth=0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey40") +
  guides(color = guide_legend(override.aes = list(shape = 16, linetype = 1, size=2))) +
  sig_cor +
  labs(x="Mean Isoeliteness",
       y="GWP accuracy for breeding fitness (r)") +
  scale_x_continuous(limits=c(minIe, 1)) +
  scale_y_continuous(limits=c(minR, maxR)) +
  scale_color +
  theme +
  theme(
    legend.position = "none",
  )

plot_RIL_IE(gwp.df)
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

# Plot isoeliteness against GWP accuracy for QTL = 10
gwp.df %>%
  dplyr::filter(type=="Admixed",
                qtl==10) %>%
  ggplot(aes(x=isoElite, y=r)) +
  geom_point(size=0.5) +  #, aes(color=isoEliteDes)) +
  geom_smooth(method="lm", se=FALSE, linewidth=0.4, color= "black") +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey40") +
  guides(color = guide_legend(override.aes = list(shape = 16, linetype = 1, size=2))) +
  sig_cor +
  labs(x="Mean Isoeliteness",
       y="GWP accuracy for breeding fitness (r)") +
  theme +
  theme(
    legend.position = "none",
  )

ggplot2::ggsave(filename = "gwp_ie_qtl10.jpg",
                path=output_dir,
                device = "jpg",
                width=3.5,
                height=3,
                dpi=600)
ggplot2::ggsave(filename = "gwp_ie_qtl10.pdf",
                path=output_dir,
                device = "pdf",
                width=3.5,
                height=3)


# Create a unique identifier for the RILTYPE_SELECTION TYPE
gs.df <- RS.df %>%
  dplyr::mutate(pop=paste0(type, " ", sel))

gs.df$pop <- factor(gs.df$pop,
                     levels=c("Admixed GS",
                              "Admixed PS",
                              "Unadmixed GS",
                              "Unadmixed PS"))
gs.df$qtl <- as.factor(as.character(gs.df$qtl))

# Calculate average breeding fitness per pop per cycle
cycleMean.df <-  gs.df %>%
  dplyr::group_by(qtl, pop, type, sel, c) %>%
  dplyr::summarize(w_se  = sd(w) / sqrt(n()),
                   w = mean(w),
                   genHt = mean(genome_het),
                   attHt = mean(attained_het),
                   desHt = mean(desired_het),
                   r = mean(r),
                   ie = mean(pop_ie)) %>%
                   #gvar = mean(gvar)) %>%
  dplyr::mutate(
    pt_fill = case_when(
      pop == "Admixed GS"   ~ "gold2",
      pop == "Admixed PS"    ~ "white",
      pop == "Unadmixed GS" ~ "#CC0000",
      pop == "Unadmixed PS"  ~ "white"
    )
  )

minW <- min(cycleMean.df$w)
maxW <- max(cycleMean.df$w)

# Plot the average breeding fitness per cycle as line plots
meanWPerCycle <- function(df, nQtl, ylabel=TRUE) {
  df %>%
    dplyr::filter(qtl==nQtl) %>%
    ggplot(aes(x = c, y = w, group = pop)) +
    #geom_errorbar(aes(ymin = w - w_se, ymax = w + w_se), width = 0.3) +
    geom_line(aes(color = pop)) +
    geom_point(aes(fill = pt_fill, color = pop),
               shape = 21, stroke = 0.5, size = 1) +
    scale_color_manual(
      name = "Population",
      values = c(
        "Admixed GS"   = "gold2",
        "Admixed PS"    = "gold2",
        "Unadmixed GS" = "#CC0000",
        "Unadmixed PS"  = "#CC0000"
      ),
      guide = guide_legend(override.aes = list(
        fill  = c("gold2", "white", "#CC0000", "white"),
        shape = 21,
        size  = 3
      ))
    ) +
    scale_fill_identity() +
    labs(
      title = paste("QTL: ", nQtl),
      x = "Cycle",
      y = if (ylabel) "Breeding Fitness" else NULL
    ) +
    scale_y_continuous(limits = c(minW, maxW)) +
    theme +
    theme(
      axis.text.y = if (!ylabel) element_blank() else element_text()
    )
}

#meanW2 <- meanWPerCycle(cycleMean.df, 2)
meanW10 <- meanWPerCycle(cycleMean.df, 10, ylabel=FALSE)
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

# Plot for QTL == 10 only
cycleMean.df %>%
  dplyr::filter(qtl==10) %>%
  ggplot(aes(x = c, y = w, group = pop)) +
  #geom_errorbar(aes(ymin = w - w_se, ymax = w + w_se), width = 0.3) +
  geom_line(aes(color = pop)) +
  geom_point(aes(fill = pt_fill, color = pop),
             shape = 21, stroke = 0.5, size = 1) +
  scale_color_manual(
    name = "Population",
    values = c(
      "Admixed GS"   = "gold2",
      "Admixed PS"    = "gold2",
      "Unadmixed GS" = "#CC0000",
      "Unadmixed PS"  = "#CC0000"
    ),
    guide = guide_legend(override.aes = list(
      fill  = c("gold2", "white", "#CC0000", "white"),
      shape = 21,
      size  = 3
    ))
  ) +
  scale_fill_identity() +
  labs(
    x = "Cycle",
    y = "Breeding Fitness"
  ) +
  theme
ggplot2::ggsave(filename = "breedingFitness_qtl10.jpg",
                path=output_dir,
                device = "jpg",
                width=3.5,
                height=3,
                dpi=600)
ggplot2::ggsave(filename = "breedingFitness_qtl10.pdf",
                path=output_dir,
                device = "pdf",
                width=3.5,
                height=3)

# PLOT EACH REPLICATE CURVE

# Group by replicate
reps.df <- gs.df %>%
  dplyr::filter(pop == "Admixed GS",
                qtl == 10) %>%
  dplyr::group_by(founder, rep, c) %>%
  dplyr::summarize(meanIe = mean(isoElite),
                   w = mean(w),
                   .groups = "drop")

# Create discrete categories of isoeliteness
ie_means.df <- plot_data %>%
  dplyr::mutate(ie_cat = round(meanIe, 1)) %>%
  dplyr::group_by(ie_cat, c) %>%
  dplyr::summarize(w = mean(w), .groups = "drop")

ggplot() +
  geom_line(data = reps.df,
            aes(x = c, y = w, group = interaction(founder, rep), color = meanIe),
            alpha = 0.2) +
  geom_line(data = ie_means.df,
            aes(x = c, y = w, group = ie_cat, color = ie_cat),
            linewidth = 1.2) +
  labs(
    x = "Cycle",
    y = "Breeding Fitness"
  ) +
  scale_color_viridis_c(name = "Isoeliteness") +
  #scale_y_continuous(limits = c(minW, maxW)) +
  theme

ggplot2::ggsave(filename = "allFitnessCurves_qtl10.jpg",
                path=output_dir,
                device = "jpg",
                width=3.5,
                height=3,
                dpi=600)
ggplot2::ggsave(filename = "allFitnessCurves_qtl10.pdf",
                path=output_dir,
                device = "pdf",
                width=3.5,
                height=3)

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
    dplyr::filter(pop=="Admixed GS") %>%
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

ht2 <- meanHetPerCycle(ht.df, 2)
ht10 <- meanHetPerCycle(ht.df, 10, ylabel=FALSE)
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

ht.df %>%
  dplyr::filter(pop=="Admixed GS") %>%
  dplyr::filter(qtl==10) %>%
  ggplot(aes(x = c, y = het, color = qtlType)) +
  geom_line() +
  geom_point() +
  labs(
    x = "Cycle",
    y = "Heterozygosity"
  ) + 
  scale_color_manual(name = "QTL",
                     values = c("genHt" = '#264653',
                                "attHt" = '#2a9d8f',
                                "desHt" = '#e9c46a'),
                     labels = c("genHt" = "Genomewide",
                                "attHt" = "Attained Trait",
                                "desHt" = "Desired Trait")) +
  theme

ggplot2::ggsave(filename = "heterozygosity_qtl10.jpg",
                path=output_dir,
                device = "jpg",
                width=3.5,
                height=3,
                dpi=600)
ggplot2::ggsave(filename = "heterozygosity_qtl10.pdf",
                path=output_dir,
                device = "pdf",
                width=3.5,
                height=3)

minPopIe <- min(cycleMean.df$ie)
maxPopIe <- max(cycleMean.df$ie)

# Plot population isoeliteness over cycles
meanIePerCycle <- function(df, nQtl, ylabel=TRUE) {
  df %>%
    dplyr::filter(qtl==nQtl) %>%
    ggplot(aes(x = c, y = ie, group = pop)) +
    geom_line(aes(color = pop)) +
    geom_point(aes(fill = pt_fill, color = pop), shape = 21, stroke = 0.5, size = 1) +
    scale_color_manual(
      name = "Population",
      values = c(
        "Admixed GS"   = "gold2",
        "Admixed PS"    = "gold2",
        "Unadmixed GS" = "#CC0000",
        "Unadmixed PS"  = "#CC0000"
      ),
      guide = guide_legend(override.aes = list(
        fill  = c("gold2", "white", "#CC0000", "white"),
        shape = 21,
        size  = 3
      ))
    ) +
    scale_fill_identity() +
    labs(
      title  = paste("QTL: ", nQtl),
      x = "Cycle",
      y = if (ylabel) "Population Isoeliteness" else NULL,
    ) + 
    scale_y_continuous(limits=c(minPopIe, maxPopIe)) +
    theme +
    theme(
      axis.text.y = if (!ylabel) element_blank() else element_text()
    )
}

meanIe2 <- meanIePerCycle(cycleMean.df, 2)
meanIe10 <- meanIePerCycle(cycleMean.df, 10, ylabel=FALSE)
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

# Plot for QTL=10 only
cycleMean.df %>%
  dplyr::filter(qtl==10) %>%
  ggplot(aes(x = c, y = ie, group = pop)) +
  geom_line(aes(color = pop)) +
  geom_point(aes(fill = pt_fill, color = pop), shape = 21, stroke = 0.5, size = 1) +
  scale_color_manual(
    name = "Population",
    values = c(
      "Admixed GS"   = "gold2",
      "Admixed PS"    = "gold2",
      "Unadmixed GS" = "#CC0000",
      "Unadmixed PS"  = "#CC0000"
    ),
    guide = guide_legend(override.aes = list(
      fill  = c("gold2", "white", "#CC0000", "white"),
      shape = 21,
      size  = 3
    ))
  ) +
  scale_fill_identity() +
  labs(
    x = "Cycle",
    y = "Population Isoeliteness"
  ) +
  theme

ggplot2::ggsave(filename = "isoeliteness_qtl10.jpg",
                path=output_dir,
                device = "jpg",
                width=3.5,
                height=3,
                dpi=600)
ggplot2::ggsave(filename = "isoeliteness_qtl10.pdf",
                path=output_dir,
                device = "pdf",
                width=3.5,
                height=3)

minRsR <- min(gs.df$r, na.rm=TRUE)
maxRsR <- max(gs.df$r, na.rm=TRUE)

# Plot R over cycles
cycleMean.df %>%
  dplyr::filter(pop=="Admixed GS") %>%
  dplyr::filter(c%%2 == 1) %>%
  ggplot(aes(x = c, y = r, color = qtl)) +
  geom_line() +
  geom_point() +
  labs(
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
    dplyr::filter(pop=="Admixed GS",
                  qtl==nQtl) %>%
    dplyr::filter(c%%2 == 1) %>%
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
    dplyr::filter(pop=="Admixed GS",
                  qtl==nQtl) %>%
    dplyr::filter(c%%2 == 1) %>%
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

hetR10
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
    "5" =  w[c == 5] - w[c == 1],
    "10" =  w[c == 10] - w[c == 1],
    "20" = w[c == max(c)] - w[c == 1],
    isoElite = mean(isoElite),
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(cols=c("5", "10", "20"),
                      names_to="cycles",
                      values_to = "gain") %>%
  dplyr::filter(qtl==10)

geneticGain.df$cycles <- factor(geneticGain.df$cycles, levels=c("5", "10", "20"))

maxGain <- max(geneticGain.df$gain)

# Plot the genetic gain per replicate as a boxplot
plotGeneticGain <- function(df, cycle, ylabel=TRUE) {
  df %>%
    dplyr::filter(cycles==cycle) %>%
    ggplot(aes(x = sel, y = gain, fill=type)) +
    geom_boxplot(
      outlier.shape=NA,
      position = position_dodge2(width = 1, padding = 0, preserve = "single"),
      linewidth = 0.1,
    ) +
    labs(
      title  = paste("Cycles: ", cycle),
      x = "Group",
      y = "Genetic Gain"
    ) + 
    scale_fill_manual(
      name = "RIL Family",
      values = c(
          "Admixed" = "gold2",
          "Unadmixed" = "#CC0000"
      )
    ) +
    scale_y_continuous(limits=c(0, maxGain)) +
    theme +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = if (!ylabel) element_blank() else element_text(),
      axis.text.y = if (!ylabel) element_blank() else element_text()
    )
}


#gain2 <- plotGeneticGain(geneticGain.df, 2)
gain5 <- plotGeneticGain(geneticGain.df, "5", ylabel=FALSE)
gain10 <- plotGeneticGain(geneticGain.df, "10", ylabel=FALSE)
gain20 <- plotGeneticGain(geneticGain.df, "20", ylabel=FALSE)

(gain5 |gain10 | gain20) + plot_layout(guides = "collect", axes = "collect")
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

minMaxGain <- geneticGain.df %>%
  dplyr::filter(type=="Admixed", sel=="GS") %>%
  summarize(min=min(gain),
            max=max(gain))
minGain <- minMaxGain$min
maxGain <- minMaxGain$max

# Plot isoeliteness against genetic gain after varying numbers of cycles
geneticGain.df %>%
  dplyr::filter(type=="Admixed", sel=="GS") %>%
  ggplot(aes(x = isoElite, y = gain, color = cycles)) +
  geom_point(size=0.5) +
  geom_smooth(method="lm", se=FALSE, aes(color=cycles), linewidth=0.4) +
  sig_cor +
  labs(
    #title  = paste("QTL: ", nQtl),
    x = "Mean Isoeliteness",
    y = "Genetic Gain"
  ) + 
  scale_color_manual(
    name = "Number of Cycles",
    values = c("#4A1A6B", "#9B59B6", "#D7B8F3")
  ) +
  #scale_x_continuous(limits=c(minIe, 1)) +
  scale_y_continuous(limits=c(minGain, maxGain)) +
  theme

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
