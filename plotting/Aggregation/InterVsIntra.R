# Title: INTERVSINTRA
# Author: Ted Monyak
# Description: This script aggregates several sim_resultss from QTL simulation runs
# and compiles them into one plot.

library("ggplot2")
library("ggpubr")
library("grid")
library("patchwork")

setwd("~/Documents/CSU/Lab-Notebooks/members/Ted/breeding_sims/Output/QtlMonteCarlo")
output_dir <- file.path(getwd(), "../InterVsIntra")
if (!dir.exists(output_dir)) dir.create(output_dir)

ne50qtl2h10 <- read.csv("Ne_50_qtl_2_selProp_0.1_h2_0.1_gens_200_var_0.05_2025-02-06_00_11/sim_results.csv") %>%
  mutate(h2=0.1,qtl=2,ne=50)
ne50qtl2h20 <- read.csv("Ne_50_qtl_2_selProp_0.1_h2_0.2_gens_200_var_0.05_2025-02-06_04_00/sim_results.csv") %>%
  mutate(h2=0.2,qtl=2,ne=50)
ne50qtl2h30 <- read.csv("Ne_50_qtl_2_selProp_0.1_h2_0.3_gens_200_var_0.05_2025-02-06_08_13/sim_results.csv") %>%
  mutate(h2=0.3,qtl=2,ne=50)
olig_small.df <- rbind(ne50qtl2h10,ne50qtl2h20, ne50qtl2h30)

ne50qtl20h10 <- read.csv("Ne_50_qtl_20_selProp_0.1_h2_0.1_gens_200_var_0.05_2025-02-06_01_02/sim_results.csv") %>%
  mutate(h2=0.1,qtl=20,ne=50)
ne50qtl20h20 <- read.csv("Ne_50_qtl_20_selProp_0.1_h2_0.2_gens_200_var_0.05_2025-02-06_04_40/sim_results.csv") %>%
  mutate(h2=0.2,qtl=20,ne=50)
ne50qtl20h30 <- read.csv("Ne_50_qtl_20_selProp_0.1_h2_0.3_gens_200_var_0.05_2025-02-06_09_38/sim_results.csv") %>%
  mutate(h2=0.3,qtl=20,ne=50)
poly_small.df <- rbind(ne50qtl20h10,ne50qtl20h20, ne50qtl20h30)

ne500qtl2h10 <- read.csv("Ne_500_qtl_2_selProp_0.1_h2_0.1_gens_200_var_0.05_2025-02-06_01_43/sim_results.csv") %>%
  mutate(h2=0.1,qtl=2,ne=500)
ne500qtl2h20 <- read.csv("Ne_500_qtl_2_selProp_0.1_h2_0.2_gens_200_var_0.05_2025-02-06_05_19/sim_results.csv") %>%
  mutate(h2=0.2,qtl=2,ne=500)
ne500qtl2h30 <- read.csv("Ne_500_qtl_2_selProp_0.1_h2_0.3_gens_200_var_0.05_2025-02-06_10_13/sim_results.csv") %>%
  mutate(h2=0.3,qtl=2,ne=500)
olig_large.df <- rbind(ne500qtl2h10,ne500qtl2h20, ne500qtl2h30)

ne500qtl20h10 <- read.csv("Ne_500_qtl_20_selProp_0.1_h2_0.1_gens_200_var_0.05_2025-02-06_02_52/sim_results.csv") %>%
  mutate(h2=0.1,qtl=2,ne=500)
ne500qtl20h20 <- read.csv("Ne_500_qtl_20_selProp_0.1_h2_0.2_gens_200_var_0.05_2025-02-06_06_23/sim_results.csv") %>%
  mutate(h2=0.2,qtl=2,ne=500)
ne500qtl20h30 <- read.csv("Ne_500_qtl_20_selProp_0.1_h2_0.3_gens_200_var_0.05_2025-02-06_11_20/sim_results.csv") %>%
  mutate(h2=0.3,qtl=2,ne=500)
poly_large.df <- rbind(ne500qtl20h10,ne500qtl20h20, ne500qtl20h30)

ymax = max(c(olig_small.df$nSigQtl,
           olig_large.df$nSigQtl,
           poly_small.df$nSigQtl,
           poly_large.df$n.SigQtl))

theme <- theme(
  axis.title.x = element_text(family="Helvetica", size=10),
  axis.text.x = element_text(angle = 0, hjust=1, size=8),
  axis.title.y = element_text(family="Helvetica", size=10),
  plot.title = element_text(family="Helvetica", size=14, hjust = 0.5),
  legend.text = element_text(family="Helvetica", size=10),
  legend.title = element_text(family="Helvetica", size=12),
  legend.key = element_rect(linewidth=0.05),
  plot.caption = element_text(family="Helvetica", size=10, hjust = 0.5),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_rect(fill = "white", color = "black"),
  plot.margin= unit(c(0,0,0,0), unit="pt"),
  aspect.ratio = 1)

colors <- c("darkorange2", "darkorchid2")

olig_small <- ggplot(data=olig_small.df, aes(x=as.factor(h2), y=nSigQtl, fill=type)) +
  geom_boxplot(outlier.shape = NA) +
  ylim(0,ymax) +
  scale_fill_manual(name = "Parental Cross Type",
                    values = c("Inter" = colors[1],
                               "Intra" = colors[2]),
                    labels = c("Inter-subpopulation", "Intra-subpopulation")) +
  theme +
  labs(x="h2",
       y="Trait QTL Detected")

olig_large <- ggplot(data=olig_large.df, aes(x=as.factor(h2), y=nSigQtl, fill=type)) +
  geom_boxplot(outlier.shape = NA) +
  ylim(0,ymax) +
  scale_fill_manual(name = "Parental Cross Type",
                    values = c("Inter" = colors[1],
                               "Intra" = colors[2]),
                    labels = c("Inter-subpopulation", "Intra-subpopulation")) +
  theme +
  labs(x="h2",
       y="Trait QTL Detected")

poly_small <- ggplot(data=poly_small.df, aes(x=as.factor(h2), y=nSigQtl, fill=type)) +
  geom_boxplot(outlier.shape = NA) +
  ylim(0,ymax) +
  scale_fill_manual(name = "Parental Cross Type",
                    values = c("Inter" = colors[1],
                               "Intra" = colors[2]),
                    labels = c("Inter-subpopulation", "Intra-subpopulation")) +
  theme +
  labs(x="h2",
       y="Trait QTL Detected")

poly_large <- ggplot(data=poly_large.df, aes(x=as.factor(h2), y=nSigQtl, fill=type)) +
  geom_boxplot(outlier.shape = NA) +
  ylim(0,ymax) +
  scale_fill_manual(name = "Parental Cross Type",
                    values = c("Inter" = colors[1],
                               "Intra" = colors[2]),
                    labels = c("Inter-subpopulation", "Intra-subpopulation")) +
  theme +
  labs(x="h2",
       y="Trait QTL Detected")
olig_label <- wrap_elements(panel = textGrob('Oligogenic\nTrait Architecture\n(QTL per Trait = 20',
                                             rot=90,
                                             gp=gpar(fontsize=12)))
poly_label <- wrap_elements(panel = textGrob('Polygenic\nTrait Architecture\n(QTL per Trait = 200)',
                                             rot=90,
                                             gp=gpar(fontsize=12))) 
small_label <- wrap_elements(panel = textGrob('Small Subopulation\n(Ne = 50)',
                                              gp=gpar(fontsize=12),
                                              vjust=2))
large_label <- wrap_elements(panel = textGrob('Large Subpopulation\n(Ne = 500)',
                                              gp=gpar(fontsize=12),
                                              vjust=2))

p <- plot_spacer() + small_label + large_label +
  olig_label + olig_small + olig_large +
  poly_label + poly_small + poly_large +
  plot_layout(guides='collect', axes='collect', widths=c(1,2,2), heights=c(0.5,0.5,0.5))& theme(legend.position="right")

fname <- file.path(output_dir, "aggregate_qtl_02_06_25.pdf")
ggplot2::ggsave(filename = fname,
                device = "pdf",
                height=7,
                width=8)
ymin = min(c(olig_small.df$fst,
             olig_large.df$fst,
             poly_small.df$fst,
             poly_large.df$fst))

ymax = max(c(olig_small.df$fst,
             olig_large.df$fst,
             poly_small.df$fst,
             poly_large.df$fst))

olig_small_fst <- ggplot(data=olig_small.df, aes(x=as.factor(h2), y=fst, fill=type)) +
  geom_boxplot(outlier.shape = NA) +
  ylim(ymin,ymax) +
  scale_fill_manual(name = "Parental Cross Type",
                    values = c("Inter" = colors[1],
                               "Intra" = colors[2]),
                    labels = c("Inter-subpopulation", "Intra-subpopulation")) +
  theme +
  labs(x="h2",
       y="FST")

olig_large_fst <- ggplot(data=olig_large.df, aes(x=as.factor(h2), y=fst, fill=type)) +
  geom_boxplot(outlier.shape = NA) +
  ylim(ymin,ymax) +
  scale_fill_manual(name = "Parental Cross Type",
                    values = c("Inter" = colors[1],
                               "Intra" = colors[2]),
                    labels = c("Inter-subpopulation", "Intra-subpopulation")) +
  theme +
  labs(x="h2",
       y="FST")

poly_small_fst <- ggplot(data=poly_small.df, aes(x=as.factor(h2), y=fst, fill=type)) +
  geom_boxplot(outlier.shape = NA) +
  ylim(ymin,ymax) +
  scale_fill_manual(name = "Parental Cross Type",
                    values = c("Inter" = colors[1],
                               "Intra" = colors[2]),
                    labels = c("Inter-subpopulation", "Intra-subpopulation")) +
  theme +
  labs(x="h2",
       y="FST")

poly_large_fst <- ggplot(data=poly_large.df, aes(x=as.factor(h2), y=fst, fill=type)) +
  geom_boxplot(outlier.shape = NA) +
  ylim(ymin,ymax) +
  scale_fill_manual(name = "Parental Cross Type",
                    values = c("Inter" = colors[1],
                               "Intra" = colors[2]),
                    labels = c("Inter-subpopulation", "Intra-subpopulation")) +
  theme +
  labs(x="h2",
       y="FST")

p <- plot_spacer() + small_label + large_label +
  olig_label + olig_small_fst + olig_large_fst +
  poly_label + poly_small_fst + poly_large_fst +
  plot_layout(guides='collect', axes='collect', widths=c(1,2,2), heights=c(0.5,0.5,0.5))& theme(legend.position="right")

fname <- file.path(output_dir, "aggregate_fst_02_06_25.pdf")
ggplot2::ggsave(filename = fname,
                device = "pdf",
                height=7,
                width=8)
