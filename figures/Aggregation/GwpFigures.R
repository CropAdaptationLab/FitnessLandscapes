# Title: GWP Figures
# Author: Ted Monyak
# Description: This script produces summary plots for cross population prediction
# and recurrent selection with data generated in GwpPipeline.R

library(dplyr)
library(ggnewscale)
library(ggplot2)
library(ggpubr)
library(grid)
library(LSD)
library(patchwork)
library(tibble)
library(tidyr)

#setwd("~/Documents/CSU/FitnessLandscapes/output/GWP/Sim_")
#output_dir <- getwd()
#RIL.df <- rbind(read.csv("QTL_10/ril_results.csv"))
#RS.df <- rbind(read.csv("QTL_10/rs_results.csv"))

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
  ),
  pt_fill = dplyr::case_when(
    r > 0 ~ "pos",
    r <= 0 ~ "neg"
  ))

gwp.df$type <- factor(gwp.df$type,
                      levels=c("Within-population", "Admixed", "Cross-population"))
gwp.df$qtl <- as.factor(gwp.df$qtl)

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
  labs(x="Test Family", y="GWP accuracy for breeding fitness (r)") +
  theme +
  theme(legend.position = "none")

ggplot2::ggsave(filename = "gwp_accuracy.jpg",
                path=output_dir,
                device = "jpg",
                width=3.5,
                height=3,
                dpi=600)
ggplot2::ggsave(filename = "gwp_accuracy.pdf",
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

# Plot isoeliteness against GWP accuracy for QTL = 10
gwp.df %>%
  dplyr::filter(type == "Admixed", qtl == 10) %>%
  ggplot(aes(x = isoElite, y = r)) +
  stat_density_2d(
    aes(fill = after_stat(level)),
    geom = "polygon",
    contour = TRUE,
    bins = 4,                          # number of contour levels
    alpha = 0.7
  ) +
  stat_density_2d(
    color = "grey30",                  # contour outlines
    linewidth = 0.2,
    contour = TRUE,
    bins = 4
  ) +
  scale_fill_gradientn(
    colors = LSD::colorpalette("heat"),
    name = "Density"
  ) +
  geom_point(size = 0.5, alpha = 0.5, color = "black") +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.4, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey40") +
  sig_cor +
  labs(
    x = "Attained Trait Isoeliteness",
    y = "GWP accuracy for breeding fitness (r)"
  ) +
  theme +
  theme(legend.position = "none")


ggplot2::ggsave(filename = "gwp_ie_scatterheat.jpg",
                path=output_dir,
                device = "jpg",
                width=3.5,
                height=3,
                dpi=600)
ggplot2::ggsave(filename = "gwp_ie_scatterheat.pdf",
                path=output_dir,
                device = "pdf",
                width=3.5,
                height=3)

gwp.df %>%
  dplyr::filter(type=="Admixed",
                qtl==10) %>%
  ggplot(aes(x=isoElite, y=r)) +
  geom_point(size=2, aes(color=isoEliteDes)) +
  geom_smooth(method="lm", se=FALSE, linewidth=0.4, color= "black") +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey40") +
  guides(color = guide_legend(override.aes = list(shape = 16, linetype = 1, size=2))) +
  sig_cor +
  labs(x="Attained Trait Isoeliteness",
       y="GWP accuracy for breeding fitness (r)") +
  theme

ggplot2::ggsave(filename = "gwp_ieAtt.jpg",
                path=output_dir,
                device = "jpg",
                width=3.5,
                height=3,
                dpi=600)
ggplot2::ggsave(filename = "gwp_ieAtt.pdf",
                path=output_dir,
                device = "pdf",
                width=3.5,
                height=3)


# Plot desired trait isoeliteness against GWP accuracy for QTL = 10
gwp.df %>%
  dplyr::filter(type=="Admixed",
                qtl==10) %>%
  ggplot(aes(x=isoEliteDes, y=r)) +
  geom_point(size=0.5) +
  geom_smooth(method="lm", se=FALSE, linewidth=0.4, color= "black") +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey40") +
  guides(color = guide_legend(override.aes = list(shape = 16, linetype = 1, size=2))) +
  sig_cor +
  labs(x="Desired Trait Isoeliteness",
       y="GWP accuracy for breeding fitness (r)") +
  theme +
  theme(
    legend.position = "none",
  )

ggplot2::ggsave(filename = "gwp_ieDes.jpg",
                path=output_dir,
                device = "jpg",
                width=3.5,
                height=3,
                dpi=600)
ggplot2::ggsave(filename = "gwp_ieDes.pdf",
                path=output_dir,
                device = "pdf",
                width=3.5,
                height=3)

gwp.df %>%
  dplyr::filter(type=="Admixed",
                qtl==10) %>%
  ggplot(aes(x=isoEliteDes, y=isoElite)) +
  geom_point(size=0.5) +
  geom_smooth(method="lm", se=FALSE, linewidth=0.4, color= "black") +
  guides(color = guide_legend(override.aes = list(shape = 16, linetype = 1, size=2))) +
  sig_cor +
  labs(x="Desired Trait Isoeliteness",
       y="Attained Trait Isoeliteness") +
  theme +
  theme(
    legend.position = "none",
  )

ggplot2::ggsave(filename = "des_att.jpg",
                path=output_dir,
                device = "jpg",
                width=3.5,
                height=3,
                dpi=600)
ggplot2::ggsave(filename = "des_att.pdf",
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
                              "Admixed MAS",
                              "Unadmixed GS",
                              "Unadmixed PS",
                              "Unadmixed MAS"))
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
                   ie = mean(pop_ie),
                   gvar = mean(gvar)) %>%
  dplyr::mutate(
    pt_fill = case_when(
      pop == "Admixed GS" ~ "black",
      pop == "Admixed PS" ~ "white",
      pop == "Unadmixed GS" ~ "gray",
      pop == "Unadmixed PS"  ~ "white"
    )
  )

minCycleW <- min(cycleMean.df$w)
maxCycleW <- max(cycleMean.df$w)

# Plot for QTL == 10 only
cycleMean.df %>%
  dplyr::filter(qtl==10) %>%
  dplyr::filter(sel != "MAS") %>%
  ggplot(aes(x = c, y = w, group = pop)) +
  #geom_errorbar(aes(ymin = w - w_se, ymax = w + w_se), width = 0.3) +
  geom_line(aes(color = pop)) +
  geom_point(aes(fill = pt_fill, color = pop),
             shape = 21, stroke = 0.5, size = 1) +
  scale_color_manual(
    name = "Population",
    values = c(
      "Admixed GS"   = "black",
      "Admixed PS"    = "black",
      "Unadmixed GS" = "gray",
      "Unadmixed PS"  = "gray"
    ),
    guide = guide_legend(override.aes = list(
      fill  = c("black", "white", "gray", "white"),
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
ggplot2::ggsave(filename = "fitness_per_cycle.jpg",
                path=output_dir,
                device = "jpg",
                width=3.5,
                height=3,
                dpi=600)
ggplot2::ggsave(filename = "fitness_per_cycle.pdf",
                path=output_dir,
                device = "pdf",
                width=3.5,
                height=3)

# PLOT EACH REPLICATE CURVE
# Group by replicate
reps.df <- gs.df %>%
  dplyr::filter(type == "Admixed",
                qtl == 10) %>%
  dplyr::group_by(founder, rep, sel, c) %>%
  dplyr::summarize(meanIe = mean(isoElite),
                   w = mean(w),
                   .groups = "drop")

# Create discrete categories of isoeliteness
ie_means.df <- reps.df %>%
  dplyr::mutate(ie_cat = dplyr::case_when(
    meanIe > 0.9 ~ "High",
    meanIe >= 0.8 & meanIe <= 0.9 ~ "Moderate",
    meanIe <  0.8 ~ "Low"
  ))

ie_cats.df <- ie_means.df %>%
  dplyr::group_by(ie_cat, sel, c) %>%
  dplyr::summarize(
    n  = n(),
    se = sd(w) / sqrt(n()),
    w  = mean(w),
    .groups = "drop"
  )
ie_cats.df$ie_cat <- factor(ie_cats.df$ie_cat,
                             levels=c("Low", "Moderate", "High"))

minW <- min(reps.df$w)
maxW <- max(reps.df$w)

plotAllCurves <- function(selType) {
  ggplot() +
    geom_line(data = reps.df %>% dplyr::filter(sel == selType),
              aes(x = c, y = w, group = interaction(founder, rep)),
              color = "grey",
              alpha = 0.1) +
    geom_errorbar(data = ie_cats.df %>% dplyr::filter(sel == selType),
                  aes(x = c, ymin = w - se, ymax = w + se, group = ie_cat, color = ie_cat),
                  width = 0.2,
                  linewidth = 0.3,
                  alpha = 0.6) +
    geom_line(data = ie_cats.df %>% dplyr::filter(sel == selType),
              aes(x = c, y = w, group = ie_cat, color = ie_cat),
              linewidth = 0.8) +
    scale_color_manual(
      name = "Attained Trait\nIsoeliteness",
      values = c("High"    = "steelblue",
                 "Moderate" = "grey40",
                 "Low"  = "firebrick")
    ) +
    labs(
      x = "Cycle",
      y = "Breeding Fitness",
      title = selType
    ) +
    #scale_y_continuous(limits = c(minW, maxW)) +
    scale_y_continuous(limits = c(120, 175)) +
    theme
}

(plotAllCurves("GS") | plotAllCurves("MAS")) + plot_layout(guides = "collect", axes = "collect")

ggplot2::ggsave(filename = "fitness_by_ieCat.jpg",
                path=output_dir,
                device = "jpg",
                width=6.5,
                height=3,
                dpi=600)
ggplot2::ggsave(filename = "fitness_by_ieCat.pdf",
                path=output_dir,
                device = "pdf",
                width=6.5,
                height=3)

# Calculate genetic gain per rep

# Safe lookup of w at c == 1
# Will return NA if there is no w value at that cycle
# Otherwise, returns the value of w at that cycle
safe_w <- function(w, c, target_c) {
  val <- w[c == target_c]
  if (length(val) == 0) NA else val[1]
}

geneticGain.df <- gs.df %>%
  dplyr::group_by(qtl, founder, rep, type, sel) %>%
  dplyr::summarize(
    "5"  = safe_w(w, c, 5) - safe_w(w, c, 1),
    "10" = safe_w(w, c, 10) - safe_w(w, c, 1),
    "20" = safe_w(w, c, max(c)) - safe_w(w, c, 1),
    meanIe = mean(isoElite),
    .groups = "drop"
  ) %>%
  drop_na() %>%
  dplyr::mutate(ie_cat = dplyr::case_when(
    meanIe > 0.9 ~ "High",
    meanIe >= 0.8 ~ "Moderate",
    meanIe <  0.8 ~ "Low"
  )) %>%
  tidyr::pivot_longer(cols=c("5", "10", "20"),
                      names_to="cycles",
                      values_to = "gain") %>%
  dplyr::filter(qtl==10)

geneticGain.df$cycles <- factor(geneticGain.df$cycles, levels=c("5", "10", "20"))
geneticGain.df$ie_cat <- factor(geneticGain.df$ie_cat,
                            levels=c("Low", "Moderate", "High"))

maxGain <- max(geneticGain.df$gain)

# Plot the genetic gain per replicate as a boxplot
plotGeneticGain <- function(df, cycle, ylabel=TRUE) {
  df %>%
    dplyr::filter(cycles==cycle, type=="Admixed") %>%
    ggplot(aes(x = sel, y = gain, fill=ie_cat)) +
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
      name = "Attained Trait\nIsoeliteness",
      values = c("High"    = "steelblue",
                 "Moderate" = "grey40",
                 "Low"  = "firebrick")
    ) +
    scale_y_continuous(limits=c(0,maxGain)) +
    stat_compare_means(
      label = "p.signif",
      hide.ns = TRUE,
      label.y = c(35,35),
      bracket.size = 0
    ) +
    theme +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = if (!ylabel) element_blank() else element_text(),
      axis.text.y = if (!ylabel) element_blank() else element_text()
    )
}


gain5 <- plotGeneticGain(geneticGain.df, "5")
gain10 <- plotGeneticGain(geneticGain.df, "10", ylabel=FALSE)
gain20 <- plotGeneticGain(geneticGain.df, "20", ylabel=FALSE)

(gain5 |gain10 | gain20) + plot_layout(guides = "collect", axes = "collect")
ggplot2::ggsave(filename = "genetic_gain_boxplot.jpg",
                path=output_dir,
                device = "jpg",
                width=6.5,
                height=2.5,
                dpi=600)
ggplot2::ggsave(filename = "genetic_gain_boxplot.pdf",
                path=output_dir,
                device = "pdf",
                width=6.5,
                height=2.5)

minMaxGain <- geneticGain.df %>%
  dplyr::filter(sel=="GS") %>%
  summarize(min=min(gain),
            max=max(gain))
minGain <- minMaxGain$min
maxGain <- minMaxGain$max

# Plot isoeliteness against genetic gain after varying numbers of cycles
geneticGain.df %>%
  dplyr::filter(type == "Admixed", sel=="GS") %>%
  ggplot(aes(x = meanIe, y = gain, color = cycles)) +
  geom_point(size=0.5) +
  geom_smooth(method="lm", se=FALSE, aes(color=cycles), linewidth=0.4) +
  sig_cor +
  labs(
    #title  = paste("QTL: ", nQtl),
    x = "Attained Trait Isoeliteness",
    y = "Genetic Gain"
  ) + 
  scale_color_manual(
    name = "Number of Cycles",
    values = c("#4A1A6B", "#9B59B6", "#D7B8F3")
  ) +
  #scale_x_continuous(limits=c(minIe, 1)) +
  scale_y_continuous(limits=c(minGain, maxGain)) +
  theme

ggplot2::ggsave(filename = "genetic_gain_ie.jpg",
                path=output_dir,
                device = "jpg",
                width=3.5,
                height=3,
                dpi=600)
ggplot2::ggsave(filename = "genetic_gain_ie.pdf",
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

ggplot2::ggsave(filename = "het.jpg",
                path=output_dir,
                device = "jpg",
                width=3.5,
                height=3,
                dpi=600)
ggplot2::ggsave(filename = "het.pdf",
                path=output_dir,
                device = "pdf",
                width=3.5,
                height=3)

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
    theme +
    theme(
      axis.title.y = if (!ylabel) element_blank() else element_text(),
      axis.text.y = if (!ylabel) element_blank() else element_text()
    )
}

attHetR(gs.df, 10)
ggplot2::ggsave(filename = "het_r.jpg",
                path=output_dir,
                device = "jpg",
                width=6.5,
                height=2.5,
                dpi=600)
ggplot2::ggsave(filename = "het_r.pdf",
                path=output_dir,
                device = "pdf",
                width=6.5,
                height=2.5)


# Plot R over cycles
cycleMean.df %>%
  dplyr::filter(pop=="Admixed GS") %>%
  dplyr::filter(c%%2 == 1) %>%
  ggplot(aes(x = c, y = r)) +
  geom_line() +
  geom_point() +
  labs(
    x = "Cycle",
    y = "GWP accuracy for breeding fitness (r)"
  ) +
  theme
ggplot2::ggsave(filename = "r_per_cycle.jpg",
                path=output_dir,
                device = "jpg",
                width=3.5,
                height=3,
                dpi=600)
ggplot2::ggsave(filename = "r_per_cycle.pdf",
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
    theme +
    theme(
      axis.title.y = if (!ylabel) element_blank() else element_text(),
      axis.text.y = if (!ylabel) element_blank() else element_text()
    )
}

popIeR(gs.df, 10)
ggplot2::ggsave(filename = "popIe_r.jpg",
                path=output_dir,
                device = "jpg",
                width=6.5,
                height=2.5,
                dpi=600)
ggplot2::ggsave(filename = "popIe_r.pdf",
                path=output_dir,
                device = "pdf",
                width=6.5,
                height=2.5)

minPopIe <- min(cycleMean.df$ie)
maxPopIe <- max(cycleMean.df$ie)

# Plot for QTL=10 only
cycleMean.df %>%
  dplyr::filter(qtl==10) %>%
  dplyr::filter(sel != "MAS") %>%
  ggplot(aes(x = c, y = ie, group = pop)) +
  geom_line(aes(color = pop)) +
  geom_point(aes(fill = pt_fill, color = pop), shape = 21, stroke = 0.5, size = 1) +
  scale_color_manual(
    name = "Population",
    values = c(
      "Admixed GS"   = "black",
      "Admixed PS"    = "black",
      "Unadmixed GS" = "gray",
      "Unadmixed PS"  = "gray"
    ),
    guide = guide_legend(override.aes = list(
      fill  = c("black", "white", "gray", "white"),
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

ggplot2::ggsave(filename = "popIe_by_cycle.jpg",
                path=output_dir,
                device = "jpg",
                width=3.5,
                height=3,
                dpi=600)
ggplot2::ggsave(filename = "popIe_by_cycle.pdf",
                path=output_dir,
                device = "pdf",
                width=3.5,
                height=3)

minRsR <- min(gs.df$r, na.rm=TRUE)
maxRsR <- max(gs.df$r, na.rm=TRUE)

# PLOT TABLES
gwp.df %>%
  dplyr::group_by(qtl, type) %>%
  dplyr::summarize(meanR=mean(r)) %>%
  write.csv(file.path(output_dir, "gwp_accuracy.csv"),quote=FALSE)

geneticGain.df %>%
  dplyr::filter(type=="Admixed") %>%
  dplyr::group_by(qtl, sel, cycles, ie_cat) %>%
  dplyr::summarize(meanGain=mean(gain),
                   varGain=var(gain))%>%
  write.csv(file.path(output_dir, "genetic_gain.csv"),quote=FALSE)


# ANOVA
sig.df <- geneticGain.df %>%
  dplyr::filter(type=="Admixed",
                qtl==10,
                sel %in% c("GS", "MAS"))


# Test whether there is more variance in GS than MAS
var.test(gain ~ sel, data=sig.df %>% filter(cycles==5,
                                                 ie_cat=="Low"))

# SIGNIFICANT
anova(lm(gain ~ sel, data=sig.df %>% filter(cycles==5,
                                                 ie_cat == "Low")))
anova(lm(gain ~ sel, data=sig.df %>% filter(cycles==10,
                                                 ie_cat == "Low")))
anova(lm(gain ~ sel, data=sig.df %>% filter(cycles==20,
                                                 ie_cat == "Low")))
anova(lm(gain ~ sel, data=sig.df %>% filter(cycles==5,
                                                 ie_cat == "Moderate")))
anova(lm(gain ~ sel, data=sig.df %>% filter(cycles==10,
                                                 ie_cat == "Moderate")))
anova(lm(gain ~ sel, data=sig.df %>% filter(cycles==20,
                                                 ie_cat == "Moderate")))
