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

#setwd("~/Documents/CSU/FitnessLandscapes/output/GWP/Sim_2026-05-14_22_53")
#output_dir <- getwd()
#RIL.df <- rbind(read.csv("QTL_10/ril_results.csv"))
#RS.df <- rbind(read.csv("QTL_10/rs_results.csv"))

theme <- theme_minimal(base_size = 10,
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
                                values = c("Within-population\n(Positive Control)" = "#7EBD4F",
                                           "Admixed" = "goldenrod1",
                                           "Cross-population\n(Negative Control)" = "#A66A28"),
                                breaks = c("Within-population\n(Positive Control)",
                                           "Admixed", 
                                           "Cross-population\n(Negative Control)"),
                                labels = c("Within-population" = "Within-population\n(Positive Control)",
                                           "Admixed" = "Admixed",
                                           "Cross-population" = "Cross-population\n(Negative Control)"))

scale_color_sel <- scale_fill_manual(name = "Selection",
                                      values = c("PS" = "#8E44C6",
                                                 "GS" = "#2EBD65",
                                                 "ieMAS (Low Density) + GS" = "#FFD4A8",
                                                 "ieMAS (High Density) + GS" = "#F06A00",
                                                 "ieMAS (Perfect Markers) + GS" = "#C14B00"))

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
    cols = c(rAdmixed, rP1, rP2, rP1P2, rP2P1,
             rCvP1P2, rCvP1, rCvP2),
    names_to = "type",
    values_to = "r"
  ) %>%
  dplyr::mutate(type = case_when(
    type %in% c("rP1", "rP2") ~ "Within-population\n(Positive Control)",
    type %in% c("rP1P2", "rP2P1") ~ "Cross-population\n(Negative Control)",
    type %in% c("rCvP1", "rCvP2") ~ "Cross-validation\n(Single Population)",
    type == "rCvP1P2" ~ "Cross-validation\n(Both Populations)",
    type == "rAdmixed" ~ "Admixed",
    TRUE ~ type
  ),
  pt_fill = dplyr::case_when(
    r > 0 ~ "pos",
    r <= 0 ~ "neg"
  ))

gwp.df$type <- factor(gwp.df$type,
                      levels=c(
                        "Within-population\n(Positive Control)",
                        "Cross-validation\n(Single Population)",
                        "Cross-validation\n(Both Populations)",
                        "Admixed",
                        "Cross-population\n(Negative Control)"))
gwp.df$qtl <- as.factor(gwp.df$qtl)

# Plot the GWP results by RIL family type for QTL == 10
gwp.df %>%
  dplyr::filter(qtl==10) %>%
  ggplot(aes(x=type, y=r, fill=model)) +
  geom_boxplot(
    outlier.shape=NA,
    linewidth = 0.1,
  ) +
  stat_compare_means(
    label = "p.signif",
    hide.ns = FALSE,
  ) +
  #scale_fill_gwp +
  scale_y_continuous(expand=c(0,0.1)) +
  labs(x="Test Family", y=expression("GWP accuracy for breeding fitness (r"[GWP]*")")) +
  theme

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
for (GS_MODEL in c("RRBLUP", "GBLUP")) {
  # Plot isoeliteness against GWP accuracy for QTL = 10
  gwp.df %>%
    dplyr::filter(type == "Admixed", model == GS_MODEL) %>%
    ggplot(aes(x = isoElite, y = r)) +
    stat_density_2d(
      aes(fill = after_stat(level)),
      geom = "polygon",
      contour = TRUE,
      bins = 5,                          # number of contour levels
      alpha = 0.7
    ) +
    stat_density_2d(
      color = "grey30",                  # contour outlines
      linewidth = 0.2,
      contour = TRUE,
      bins = 5
    ) +
    scale_fill_gradientn(
      colours = c("#FFF9C4", "#FFF176", "#FFEB3B", "#FDD835", "#FBC02D"),
      name = "Density"
    ) +
    geom_point(size = 0.5, alpha = 0.4, color = "black") +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.4, color = "black") +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey40") +
    sig_cor +
    labs(
      x = "Mean Isoeliteness",
      y = expression("GWP accuracy for breeding fitness (r"[GWP]*")"),
    ) +
    theme +
    theme(legend.position = "none")
  
  
  ggplot2::ggsave(filename = paste0("gwp_ie_scatterheat_", GS_MODEL, ".jpg"),
                  path=output_dir,
                  device = "jpg",
                  width=3.5,
                  height=3,
                  dpi=600)
  ggplot2::ggsave(filename = paste0("gwp_ie_scatterheat_", GS_MODEL, ".pdf"),
                  path=output_dir,
                  device = "pdf",
                  width=3.5,
                  height=3)
  
  # Plot desired trait heterogeneity against GWP accuracy for QTL = 10
  ie_des <- gwp.df %>%
    dplyr::filter(type=="Admixed",
                  qtl==10,
                  model==GS_MODEL) %>%
    ggplot(aes(x=isoEliteDes, y=r)) +
    geom_point(size=0.5) +
    geom_smooth(method="lm", se=FALSE, linewidth=0.4, color= "black") +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, color = "grey40") +
    guides(color = guide_legend(override.aes = list(shape = 16, linetype = 1, size=2))) +
    sig_cor +
    labs(x="Desired Trait Heterogeneity",
         y="GWP accuracy for breeding fitness (r)") +
    theme +
    theme(
      legend.position = "none",
    )
  
  att_des <- gwp.df %>%
    dplyr::filter(type=="Admixed",
                  qtl==10,
                  model==GS_MODEL) %>%
    ggplot(aes(x=isoEliteDes, y=isoElite)) +
    geom_point(size=0.5) +
    geom_smooth(method="lm", se=FALSE, linewidth=0.4, color= "black") +
    guides(color = guide_legend(override.aes = list(shape = 16, linetype = 1, size=2))) +
    sig_cor +
    labs(x="Desired Trait Heterogeneity",
         y="Mean Isoeliteness") +
    theme +
    theme(
      legend.position = "none",
    )
  
  (ie_des | att_des)
  
  ggplot2::ggsave(filename = paste0("gwp_ieDes_", GS_MODEL, ".jpg"),
                  path=output_dir,
                  device = "jpg",
                  width=6.5,
                  height=3,
                  dpi=600)
  ggplot2::ggsave(filename = paste0("gwp_ieDes_", GS_MODEL, ".pdf"),
                  path=output_dir,
                  device = "pdf",
                  width=6.5,
                  height=3)
  
  # Create a unique identifier for the RILTYPE_SELECTION TYPE
  gs.df <- RS.df %>%
    dplyr::filter(model == GS_MODEL) %>%
    dplyr::mutate(sel = case_when(
      sel == "highResMAS" ~ "ieMAS (High Density) + GS",
      sel == "lowResMAS" ~ "ieMAS (Low Density) + GS",
      sel == "ieMAS" ~ "ieMAS (Perfect Markers) + GS",
      TRUE ~ sel
    )) %>%
    dplyr::mutate(pop=paste0(type, " ", sel)) %>%
    dplyr::mutate(pop=factor(pop, levels=c("Admixed GS",
                                           "Admixed PS",
                                           "Admixed ieMAS",
                                           "Admixed ieMAS (Low Density) + GS",
                                           "Admixed ieMAS (High Density) + GS",
                                           "Admixed ieMAS (Perfect Markers) + GS",
                                           "Unadmixed GS",
                                           "Unadmixed PS",
                                           "Unadmixed ieMAS",
                                           "Unadmixed ieMAS (Low Density) + GS",
                                           "Unadmixed ieMAS (High Density) + GS",
                                           "Unadmixed ieMAS (Perfect Markers) + GS")))
  
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
        pop == "Admixed GS" ~ "#2EBD65",
        pop == "Admixed PS" ~ "#8E44C6",
        pop == "Admixed ieMAS (High Density) + GS" ~ "#F06A00",
        pop == "Unadmixed GS" ~ "#2EBD65",
        pop == "Unadmixed PS"  ~ "#8E44C6",
        pop == "Unadmixed ieMAS (High Density) + GS"  ~ "#F06A00"
      )
    )
  
  
  # PLOT EACH REPLICATE CURVE
  # Group by replicate
  reps.df <- gs.df %>%
    dplyr::filter(type == "Admixed",
                  qtl == 10,
                  sel %in% c("GS", "ieMAS (High Density) + GS")) %>%
    dplyr::group_by(founder, rep, sel, c) %>%
    dplyr::summarize(meanIe = mean(isoElite),
                     w = mean(w),
                     .groups = "drop") %>%
    dplyr::mutate(ie_cat = dplyr::case_when(
      meanIe > 0.9 ~ "High",
      meanIe >= 0.8 & meanIe <= 0.9 ~ "Moderate",
      meanIe <  0.8 ~ "Low"
    ))
  
  ie_cats.df <- reps.df %>%
    dplyr::group_by(ie_cat, sel, c) %>%
    dplyr::summarize(
      n  = n(),
      se = sd(w) / sqrt(n()),
      w = mean(w),
      .groups = "drop"
    ) %>%
    # Update the name here
    dplyr::mutate(sel = factor(sel, levels = c("GS", "ieMAS (High Density) + GS"))) %>%
    dplyr::mutate(ie_cat = factor(ie_cat, levels = c("Low", "Moderate", "High"))) %>%
    dplyr::mutate(pop = paste0(sel, " ", ie_cat)) %>%
    dplyr::mutate(pop = factor(pop, levels=c("GS Low",
                                             "GS Moderate",
                                             "GS High",
                                             "ieMAS (High Density) + GS Low",
                                             "ieMAS (High Density) + GS Moderate",
                                             "ieMAS (High Density) + GS High"))) %>%
    dplyr::mutate(
      pt_fill = case_when(
        ie_cat == "Low" ~ "#cc0000",
        ie_cat == "Moderate" ~ "grey40",
        ie_cat == "High" ~ "#0957d6"
      )
    )
  
  minCycleW <- min(ie_cats.df$w)
  maxCycleW <- max(ie_cats.df$w)
  
  # Plot for QTL == 10 only
  cycleMean.df %>%
    dplyr::filter(qtl==10) %>%
    dplyr::filter(sel %in% c("GS", "PS")) %>%
    ggplot(aes(x = c, y = w, group = pop, color=pt_fill)) +
    geom_line() +
    geom_point(aes(fill = ifelse(type == "Unadmixed", "white", pt_fill), shape = pop), stroke = 0.5, size = 1) +
    scale_shape_manual(
      name = "Population",
      values = c(
        "Admixed GS"   = 21,
        "Admixed PS"    = 21,
        "Unadmixed GS" = 21,
        "Unadmixed PS"  = 21
      ),
      guide = guide_legend(override.aes = list(
        fill  = c("#2EBD65", "#8E44C6", "white", "white"),
        color = c("#2EBD65", "#8E44C6", "#2EBD65", "#8E44C6"),
        size  = 3
      ))
    ) +
    scale_fill_identity() +
    scale_color_identity() +
    scale_y_continuous(limits = c(minCycleW, maxCycleW)) +
    labs(
      x = "Cycle",
      y = "Realized Yield\n(Breeding Fitness)"
    ) +
    theme
  
  ggplot2::ggsave(filename = paste0("GS_vs_PS_", GS_MODEL, ".jpg"),
                  path=output_dir,
                  device = "jpg",
                  width=3.5,
                  height=3,
                  dpi=600)
  ggplot2::ggsave(filename = paste0("GS_vs_PS_", GS_MODEL, ".pdf"),
                  path=output_dir,
                  device = "pdf",
                  width=3.5,
                  height=3)

  plotAllCurves <- function(selType, ylabel = TRUE) {
    ggplot() +
      geom_line(data = reps.df %>% dplyr::filter(sel == selType),
                aes(x = c, y = w, group = interaction(founder, rep)),
                color = "grey",
                alpha = 0.1) +
      geom_errorbar(data = ie_cats.df %>% dplyr::filter(sel == selType),
                    aes(x = c, ymin = w - se, ymax = w + se, group = ie_cat, color = pt_fill),
                    width = 0.2,
                    linewidth = 0.3,
                    alpha = 0.6) +
      geom_line(data = ie_cats.df %>% dplyr::filter(sel == selType),
                aes(x = c, y = w, group = ie_cat, color = pt_fill),
                linewidth = 0.8) +
      scale_color_identity(
        name = "Mean Isoeliteness",
        guide = "legend",
        breaks = c("#0957d6", "grey40", "#cc0000"),
        labels = c("High", "Moderate", "Low")
      ) +
      labs(
        x = "Cycle",
        y = "Realized Yield\n(Breeding Fitness)",
        title = selType
      ) +
      scale_y_continuous(limits = c(minCycleW, maxCycleW)) +
      theme +
      theme(
        axis.title.y = if (!ylabel) element_blank() else element_text(),
        axis.text.y = if (!ylabel) element_blank() else element_text()
      )
  }
  
  (plotAllCurves("GS") | plotAllCurves("ieMAS (High Density) + GS", FALSE)) + plot_layout(guides = "collect", axes = "collect")
  
  ggplot2::ggsave(filename = paste0("GS_vs_MAS_Curves_", GS_MODEL, ".jpg"),
                  path=output_dir,
                  device = "jpg",
                  width=6.5,
                  height=3,
                  dpi=600)
  ggplot2::ggsave(filename = paste0("GS_vs_MAS_Curves_", GS_MODEL, ".pdf"),
                  path=output_dir,
                  device = "pdf",
                  width=6.5,
                  height=3)
  
  ie_cats.df %>%
    dplyr::filter(sel %in% c("GS", "ieMAS (High Density) + GS")) %>%
    ggplot(aes(x = c, y = w, group = pop, color = sel)) +
    # Main trend lines
    geom_line(linewidth = 0.8) +
    # Points with shapes mapped to isoeliteness
    geom_point(aes(shape = ie_cat), size = 2.5) +
    # Selection colors (Green vs Orange)
    scale_color_manual(
      name = "Selection Scheme",
      values = c(
        "GS" = "#2EBD65", 
        "ieMAS (High Density) + GS" = "#F06A00"
      )
    ) +
    # Different shapes for isoeliteness levels
    scale_shape_manual(
      name = "Mean Isoeliteness",
      values = c(
        "High"     = 16, # Solid Circle
        "Moderate" = 17, # Solid Triangle
        "Low"      = 15  # Solid Square
      )
    ) +
    labs(
      x = "Cycle", 
      y = "Realized Yield\n(Breeding Fitness)"
    ) +
    scale_y_continuous(limits = c(minCycleW, maxCycleW)) +
    theme
  
  ggplot2::ggsave(filename = paste0("GS_vs_MAS_Overlaid_", GS_MODEL, ".jpg"),
                  path=output_dir,
                  device = "jpg",
                  width=3.5,
                  height=3,
                  dpi=600)
  ggplot2::ggsave(filename = paste0("GS_vs_MAS_Overlaid_", GS_MODEL, ".pdf"),
                  path=output_dir,
                  device = "pdf",
                  width=6.5,
                  height=3)
  
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
    dplyr::filter(qtl==10) %>%
    dplyr::mutate(cycles=factor(cycles, levels=c("5", "10", "20")),
                  ie_cat=factor(ie_cat, levels=c("Low", "Moderate", "High")),
                  sel=factor(sel, levels=c("PS", "GS", "ieMAS (Low Density) + GS", "ieMAS (High Density) + GS", "ieMAS (Perfect Markers) + GS")))
  
  # Plot the genetic gain per replicate as a boxplot
  plotGeneticGain <- function(df, cycle, selTypes, ylabel=TRUE) {
    df %>%
      dplyr::filter(cycles==cycle,
                    type=="Admixed",
                    sel %in% selTypes) %>%
      ggplot(aes(x = ie_cat, y = gain, fill=sel)) +
      geom_boxplot(
        outlier.shape=NA,
        position = position_dodge2(width = 1, padding = 0, preserve = "single"),
        linewidth = 0.1,
      ) +
      labs(
        title  = paste("Cycles: ", cycle),
        x = "Mean Isoeliteness",
        y = "Genetic Gain"
      ) + 
      scale_color_sel +
      #scale_y_continuous(limits=c(0,maxGain)) +
      stat_compare_means(
        label = "p.signif",
        hide.ns = TRUE,
        label.y = c(35,35),
        bracket.size = 0
      ) +
      theme +
      theme(
        axis.title.y = if (!ylabel) element_blank() else element_text(),
        axis.text.y = if (!ylabel) element_blank() else element_text()
      )
  }
  
  gain5 <- plotGeneticGain(geneticGain.df, "5", c("GS", "ieMAS (High Density) + GS"))
  gain10 <- plotGeneticGain(geneticGain.df, "10", c("GS", "ieMAS (High Density) + GS"), ylabel=FALSE)
  gain20 <- plotGeneticGain(geneticGain.df, "20", c("GS", "ieMAS (High Density) + GS"), ylabel=FALSE)
  
  (gain5 | gain10 | gain20) + plot_layout(guides = "collect", axes = "collect")
  ggplot2::ggsave(filename = paste0("GS_vs_MAS_Boxplot_", GS_MODEL, ".jpg"),
                  path=output_dir,
                  device = "jpg",
                  width=6.5,
                  height=2.5,
                  dpi=600)
  ggplot2::ggsave(filename = paste0("GS_vs_MAS_Boxplot_", GS_MODEL, ".pdf"),
                  path=output_dir,
                  device = "pdf",
                  width=6.5,
                  height=2.5)
  
  gain5_MAS <- plotGeneticGain(geneticGain.df, "5", c("GS", "PS", "ieMAS (Low Density) + GS", "ieMAS (High Density) + GS", "ieMAS (Perfect Markers) + GS"))
  gain10_MAS <- plotGeneticGain(geneticGain.df, "10", c("GS", "PS", "ieMAS (Low Density) + GS", "ieMAS (High Density) + GS", "ieMAS (Perfect Markers) + GS"), ylabel=FALSE)
  gain20_MAS <- plotGeneticGain(geneticGain.df, "20", c("GS", "PS", "ieMAS (Low Density) + GS", "ieMAS (High Density) + GS", "ieMAS (Perfect Markers) + GS"), ylabel=FALSE)

  (gain5_MAS | gain10_MAS | gain20_MAS) + plot_layout(guides = "collect", axes = "collect")
  ggplot2::ggsave(filename = paste0("GainAll_Boxplot_", GS_MODEL, ".jpg"),
                  path=output_dir,
                  device = "jpg",
                  width=6.5,
                  height=2.5,
                  dpi=600)
  ggplot2::ggsave(filename = paste0("GainAll_Boxplot_", GS_MODEL, ".pdf"),
                  path=output_dir,
                  device = "pdf",
                  width=6.5,
                  height=2.5)
  
  minMaxGain <- geneticGain.df %>%
    dplyr::summarize(min=min(gain),
                     max=max(gain))
  
  minGain <- minMaxGain$min
  maxGain <- minMaxGain$max
  
  # Plot isoeliteness against genetic gain after varying numbers of cycles
  plotGeneticGainVsIe <- function(df, selType, title, ylabel=TRUE) {
    df %>%
      dplyr::filter(type == "Admixed", sel==selType) %>%
      ggplot(aes(x = meanIe, y = gain, color = cycles)) +
      geom_point(size=0.5) +
      geom_smooth(method="lm", se=FALSE, aes(color=cycles), linewidth=0.4) +
      sig_cor +
      labs(
        title  = title,
        x = "Mean Isoeliteness",
        y = "Genetic Gain"
      ) + 
      scale_color_manual(
        name = "Number of Cycles",
        values = c("black", "grey40", "grey70")
      ) +
      scale_y_continuous(limits=c(minGain, maxGain)) +
      theme +
      theme(
        axis.title.y = if (!ylabel) element_blank() else element_text(),
        axis.text.y = if (!ylabel) element_blank() else element_text()
      )
  }
  
  
  (plotGeneticGainVsIe(geneticGain.df, "GS", "GS") |
      plotGeneticGainVsIe(geneticGain.df, "ieMAS (High Density) + GS", "ieMAS + GS", FALSE)) +
    plot_layout(guides = "collect", axes = "collect")
  
  ggplot2::ggsave(filename = paste0("Gain_vs_Ie_", GS_MODEL, ".jpg"),
                  path=output_dir,
                  device = "jpg",
                  width=6.5,
                  height=3,
                  dpi=600)
  ggplot2::ggsave(filename = paste0("Gain_vs_Ie_", GS_MODEL, ".pdf"),
                  path=output_dir,
                  device = "pdf",
                  width=6.5,
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
  
  ggplot2::ggsave(filename = paste0("Heterozygosity_", GS_MODEL, ".jpg"),
                  path=output_dir,
                  device = "jpg",
                  width=3.5,
                  height=3,
                  dpi=600)
  ggplot2::ggsave(filename = paste0("Heterozygosity_", GS_MODEL, ".pdf"),
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
    dplyr::filter(sel %in% c("GS", "PS", "highResMAS")) %>%
    ggplot(aes(x = c, y = ie, group = pop)) +
    geom_line() +
    geom_point(aes(fill = pt_fill, shape = pop), stroke = 0.5, size = 1) +
    scale_shape_manual(
      name = "Population",
      values = c(
        "Admixed GS"   = 16,
        "Admixed PS"    = 21,
        "Admixed highResMAS" = 21,
        "Unadmixed GS" = 17,
        "Unadmixed PS"  = 24,
        "Unadmixed highResMAS" = 24
      ),
      guide = guide_legend(override.aes = list(
        fill  = c("black", "white", "grey", "black", "white", "grey"),
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
    write.csv(file.path(output_dir, "gwp_accuracy.csv"), quote=FALSE)
  
  geneticGain.df %>%
    dplyr::filter(type=="Admixed") %>%
    dplyr::group_by(qtl, sel, cycles, ie_cat) %>%
    dplyr::summarize(meanGain=mean(gain),
                     varGain=var(gain))%>%
    write.csv(file.path(output_dir, "genetic_gain.csv"), quote=FALSE)
  
  
  # ANOVA
  sig.df <- geneticGain.df %>%
    dplyr::filter(type=="Admixed",
                  qtl==10,
                  sel %in% c("GS", "highResMAS"))
  
  
  # Test whether there is more variance in GS than MAS
  var.test(gain ~ sel, data=sig.df %>% filter(cycles==5,
                                              ie_cat=="Low"))
  var.test(gain ~ sel, data=sig.df %>% filter(cycles==5,
                                              ie_cat=="Moderate"))
  var.test(gain ~ sel, data=sig.df %>% filter(cycles==5,
                                              ie_cat=="High"))
  
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
  anova(lm(gain ~ sel, data=sig.df %>% filter(cycles==5,
                                              ie_cat == "High")))
  anova(lm(gain ~ sel, data=sig.df %>% filter(cycles==10,
                                              ie_cat == "High")))
  anova(lm(gain ~ sel, data=sig.df %>% filter(cycles==20,
                                              ie_cat == "High")))
}




