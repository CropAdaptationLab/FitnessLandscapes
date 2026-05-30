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

# setwd("~/Documents/CSU/FitnessLandscapes/output/GWP/Sim_2026-05-28_16_22")
# output_dir <- getwd()
# RIL.df <- rbind(read.csv("RRBLUP/ril_results.csv"))
#                read.csv("GBLUP/ril_results.csv"))
# RS.df <- rbind(read.csv("RRBLUP/rs_results.csv"))
#               read.csv("GBLUP/rs_results.csv"))

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

scale_color_gwp <- scale_color_manual(name = "Test Family",
                                values = c("CV (Both Subpopulations)" = "navy",
                                           "CV (Single Subpopulation)" = "steelblue2",
                                           "Unadmixed" = "#7EBD4F",
                                           "Admixed" = "goldenrod1",
                                           "Cross-population" = "#A66A28"),
                                breaks = c("CV (Both Subpopulations)",
                                           "CV (Single Subpopulation)",
                                           "Unadmixed",
                                           "Admixed", 
                                           "Cross-population"))

color_scheme <- c(
  "PS" = "#8E44C6",
  "ieMAS (High Density) + PS" = "#eec7ff",
  "GS" = "#2EBD65",
  "GS (No Update)" = "#b7edba",
  "ieMAS (Low Density) + GS" = "#FFD4A8",
  "ieMAS (High Density) + GS" = "#F06A00",
  "ieMAS (High Density) + GS (No Update)" = "#C14B00",
  "ieMAS (Perfect Markers) + GS" = "#7A2900")

scale_fill_sel <- scale_fill_manual(name = "Selection",
                                      values = color_scheme)


scale_color_sel <- scale_color_manual(name = "Selection",
                                    values = color_scheme)

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
    type %in% c("rP1", "rP2") ~ "Unadmixed",
    type %in% c("rP1P2", "rP2P1") ~ "Cross-population",
    type %in% c("rCvP1", "rCvP2") ~ "CV (Single Subpopulation)",
    type == "rCvP1P2" ~ "CV (Both Subpopulations)",
    type == "rAdmixed" ~ "Admixed",
    TRUE ~ type
  ),
  pt_fill = dplyr::case_when(
    r > 0 ~ "pos",
    r <= 0 ~ "neg"
  ))

gwp.df$type <- factor(gwp.df$type,
                      levels=c(
                        "CV (Both Subpopulations)",
                        "CV (Single Subpopulation)",
                        "Unadmixed",
                        "Admixed",
                        "Cross-population"))

# PLOT TABLES
gwp.df %>%
  dplyr::group_by(qtl, type, model) %>%
  dplyr::summarize(meanR=mean(r)) %>%
  write.csv(file.path(output_dir, "gwp_accuracy.csv"), quote=FALSE)

for (GS_MODEL in c("RRBLUP")) {
  # Plot the GWP results by test scenario type
  gwp.df %>% dplyr::filter(model == GS_MODEL) %>%
    ggplot(aes(x=type, y=r, color=type)) +
    # 95% CI
    stat_summary(
      fun.data = mean_cl_normal, 
      geom = "errorbar", 
      position = position_dodge(width = 0.75),
      width = 0.5,
      linewidth = 1
    ) +
    # Point estimate
    stat_summary(
      fun = mean, 
      geom = "point", 
      position = position_dodge(width = 0.75),
      size = 1.5
    ) +
    #stat_compare_means(
    #  aes(group=model),
    #  label = "p.signif",
    #  hide.ns = FALSE,
    #) +
    scale_color_gwp +
    scale_y_continuous(expand=c(0,0.1)) +
    labs(x="Scenario", y=expression("GWP accuracy for breeding fitness (r"[MG]*")")) +
    theme +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggplot2::ggsave(filename = paste0("gwp_accuracy_", GS_MODEL, ".jpg"),
                  path=output_dir,
                  device = "jpg",
                  width=6,
                  height=5,
                  dpi=600)
  ggplot2::ggsave(filename = paste0("gwp_accuracy_", GS_MODEL, ".pdf"),
                  path=output_dir,
                  device = "pdf",
                  width=6,
                  height=5)
  
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
      y = expression("GWP accuracy for breeding fitness (r"[MG]*")"),
    ) +
    theme +
    theme(legend.position = "none")
  
  
  ggplot2::ggsave(filename = paste0("gwp_ie_scatterheat_", GS_MODEL, ".jpg"),
                  path=output_dir,
                  device = "jpg",
                  width=5,
                  height=4,
                  dpi=600)
  ggplot2::ggsave(filename = paste0("gwp_ie_scatterheat_", GS_MODEL, ".pdf"),
                  path=output_dir,
                  device = "pdf",
                  width=5,
                  height=4)
  
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
    dplyr::group_by(founder, rep, type, sel) %>%
    dplyr::mutate(gain = wGV - wGV[1]) %>%
    dplyr::mutate(sel = case_when(
      sel == "highResMAS" ~ "ieMAS (High Density) + GS",
      sel == "lowResMAS" ~ "ieMAS (Low Density) + GS",
      sel == "ieMAS" ~ "ieMAS (Perfect Markers) + GS",
      sel == "GSnoUpdate" ~ "GS (No Update)",
      sel == "PSieMAS" ~ "ieMAS (High Density) + PS",
      sel == "highResMASnoUpdate" ~ "ieMAS (High Density) + GS (No Update)",
      TRUE ~ sel
    )) %>%
    dplyr::mutate(
      pt_fill = color_scheme[sel]
    ) %>%
    dplyr::mutate(pop=paste0(type, " ", sel)) %>%
    dplyr::mutate(
      type=factor(type, levels=c("Admixed",
                                 "Unadmixed")),
      sel=factor(sel, levels=c("PS",
                               "ieMAS (High Density) + PS",
                               "GS (No Update)",
                               "GS",
                               "ieMAS (Low Density) + GS",
                               "ieMAS (High Density) + GS (No Update)",
                               "ieMAS (High Density) + GS",
                               "ieMAS (Perfect Markers) + GS")),
      pop=factor(pop, levels=c("Admixed PS",
                               "Admixed ieMAS (High Density) + PS",
                               "Admixed GS (No Update)",
                               "Admixed GS",
                               "Admixed ieMAS (Low Density) + GS",
                               "Admixed ieMAS (High Density) + GS (No Update)",
                               "Admixed ieMAS (High Density) + GS",
                               "Admixed ieMAS (Perfect Markers) + GS",
                               "Unadmixed PS",
                               "Unadmixed ieMAS (High Density) + PS",
                               "Unadmixed GS (No Update)",
                               "Unadmixed GS",
                               "Unadmixed ieMAS (Low Density) + GS",
                               "Unadmixed ieMAS (High Density) + GS (No Update)",
                               "Unadmixed ieMAS (High Density) + GS",
                               "Unadmixed ieMAS (Perfect Markers) + GS")))
  

  # Calculate average breeding fitness per pop per cycle
  cycleMean.df <- gs.df %>%
    dplyr::group_by(pop, type, sel, c, pt_fill) %>%
    dplyr::summarize(gain_CI  = 1.96*(sd(gain) / sqrt(n())),
                     gain = mean(gain),
                     genHt = mean(genome_het),
                     attHt = mean(attained_het),
                     desHt = mean(desired_het),
                     r = mean(r),
                     gvar = mean(gvar))
  
  
  # Plot for QTL == 10 only
  cycleMean.df %>%
    dplyr::filter(sel %in% c("GS", "PS")) %>%
    ggplot(aes(x = c, y = gain, group = pop, color=pt_fill)) +
    geom_line() +
    geom_errorbar(aes(ymin = gain - gain_CI, ymax = gain + gain_CI),
                  width = 0.2,
                  linewidth = 0.3,
                  alpha = 0.6) +
    geom_point(aes(fill = ifelse(type == "Unadmixed", "white", pt_fill),
                   shape = pop),
               stroke = 0.5,
               size = 1) +
    scale_shape_manual(
      name = "Population",
      values = c(
        "Admixed PS"    = 21,
        "Admixed GS"   = 21,
        "Unadmixed PS"  = 21,
        "Unadmixed GS" = 21
      ),
      guide = guide_legend(override.aes = list(
        fill  = c(color_scheme[["PS"]], color_scheme[["GS"]], "white", "white"),
        color = c(color_scheme[["PS"]], color_scheme[["GS"]], color_scheme[["PS"]], color_scheme[["GS"]]),
        size  = 3
      ))
    ) +
    scale_fill_identity() +
    scale_color_identity() +
    scale_x_continuous(breaks=seq(from=0, to=20, by=2)) +
    #scale_y_continuous(limits = c(minCycleW, maxCycleW)) +
    labs(
      x = "Year",
      y = "Genetic Gain\n(Realized Yield units)"
    ) +
    theme

  ggplot2::ggsave(filename = paste0("GS_vs_PS_", GS_MODEL, ".jpg"),
                  path=output_dir,
                  device = "jpg",
                  width=6,
                  height=5,
                  dpi=600)
  ggplot2::ggsave(filename = paste0("GS_vs_PS_", GS_MODEL, ".pdf"),
                  path=output_dir,
                  device = "pdf",
                  width=6,
                  height=5)
  
  
  upper75 <- quantile(gs.df$isoElite[gs.df$type == "Admixed"], 0.75)
  lower25 <- quantile(gs.df$isoElite[gs.df$type == "Admixed"], 0.25)
  gs.df <- gs.df %>%
    dplyr::mutate(
      ie_cat = dplyr::case_when(
        isoElite >= upper75 ~ "High",
        isoElite > lower25 & isoElite < upper75 ~ "Moderate",
        isoElite <=  lower25 ~ "Low"
      ),
      ie_cat=factor(ie_cat, levels=c("Low", "Moderate", "High")))
  
  # Show the relationship between mean isoeliteness and breeding fitness in RIL family
  gs.df %>%
    dplyr::filter(type == "Admixed") %>%
    dplyr::filter(c == 0) %>%
    dplyr::group_by(founder, rep) %>%
    dplyr::summarize(ie=mean(isoElite),
                     w=mean(wGV)) %>%
    dplyr::ungroup() %>%
    ggplot(aes(x=ie, y=w)) +
    geom_point(size=0.5) +
    geom_smooth(method="lm", se=FALSE, linewidth=0.4, color= "black") +
    guides(color = guide_legend(override.aes = list(shape = 16, linetype = 1, size=2))) +
    sig_cor +
    labs(x="Mean Isoeliteness",
         y="Breeding Fitness of RIL Family") +
    theme +
    theme(
      legend.position = "none",
    )
  
  ggplot2::ggsave(filename = paste0("w_vs_ie_", GS_MODEL, ".jpg"),
                  path=output_dir,
                  device = "jpg",
                  width=5,
                  height=4,
                  dpi=600)
  ggplot2::ggsave(filename = paste0("w_vs_ie_", GS_MODEL, ".pdf"),
                  path=output_dir,
                  device = "pdf",
                  width=5,
                  height=4)
  
  # PLOT EACH REPLICATE CURVE
  # Group by replicate
  ie_cats.df <- gs.df %>%
    dplyr::filter(type == "Admixed") %>%
    dplyr::group_by(ie_cat, sel, c) %>%
    dplyr::summarize(gain_CI  = 1.96*(sd(gain) / sqrt(n())),
                     gain = mean(gain),
                     w_CI  = 1.96*(sd(wGV) / sqrt(n())),
                     w = mean(wGV)) %>%
    dplyr::mutate(ie_cat = factor(ie_cat, levels = c("Low", "Moderate", "High"))) %>%
    dplyr::mutate(pop = paste0(sel, " ", ie_cat))
  
  minCycleGain <- min(ie_cats.df$gain-ie_cats.df$gain_CI)
  maxCycleGain <- max(ie_cats.df$gain+ie_cats.df$gain_CI)

  minCycleW <- min(ie_cats.df$w-ie_cats.df$w_CI)
  maxCycleW <- max(ie_cats.df$w+ie_cats.df$w_CI)
  
  plotAllCurvesGain <- function(selType, ylabel = TRUE) {
    ie_cats.df %>%
      dplyr::filter(sel == selType) %>%
      ggplot(aes(x = c, y = gain, group = ie_cat, color = ie_cat)) +
      geom_line(linewidth = 0.8) +
      geom_point(size=1) +
      geom_errorbar(aes(ymin = gain - gain_CI, ymax = gain + gain_CI),
                    width = 0.2,
                    linewidth = 0.3,
                    alpha = 0.6) +
      scale_color_manual(
        name = "Mean Isoeliteness",
        values =c(
          "Low" = "#cc0000",
          "Moderate" = "grey40",
          "High" = "#0957d6"
          )) +
      labs(
        x = "Year",
        y = "Genetic Gain\n(Realized Yield units)",
        title = selType
      ) +
      scale_x_continuous(breaks=seq(from=0, to=20, by=2)) +
      scale_y_continuous(limits = c(minCycleGain, maxCycleGain)) +
      theme +
      theme(
        axis.title.y = if (!ylabel) element_blank() else element_text(),
        axis.text.y = if (!ylabel) element_blank() else element_text()
      )
  }
  
  ((plotAllCurvesGain("GS") | plotAllCurvesGain("ieMAS (High Density) + GS", FALSE)) /
      (plotAllCurvesGain("PS")  | plotAllCurvesGain("ieMAS (High Density) + PS", FALSE))) +
    plot_layout(guides = "collect", axes = "collect")
  
  ggplot2::ggsave(filename = paste0("GS_vs_MAS_Curves_Gain_", GS_MODEL, ".jpg"),
                  path=output_dir,
                  device = "jpg",
                  width=7,
                  height=5,
                  dpi=600)
  ggplot2::ggsave(filename = paste0("GS_vs_MAS_Curves_Gain_", GS_MODEL, ".pdf"),
                  path=output_dir,
                  device = "pdf",
                  width=7,
                  height=5)
  
  plotAllCurvesW <- function(selType, ylabel = TRUE) {
    ie_cats.df %>%
      dplyr::filter(sel == selType) %>%
      ggplot(aes(x = c, y = w, group = ie_cat, color = ie_cat)) +
      geom_line(linewidth = 0.8) +
      geom_point(size=1) +
      geom_errorbar(aes(ymin = w - w_CI, ymax = w + w_CI),
                    width = 0.2,
                    linewidth = 0.3,
                    alpha = 0.6) +
      scale_color_manual(
        name = "Mean Isoeliteness",
        values =c(
          "Low" = "#cc0000",
          "Moderate" = "grey40",
          "High" = "#0957d6"
        )) +
      labs(
        x = "Year",
        y = "Breeding Fitness\n(Realized Yield)",
        title = selType
      ) +
      scale_x_continuous(breaks=seq(from=0, to=20, by=2)) +
      scale_y_continuous(limits = c(minCycleW, maxCycleW)) +
      theme +
      theme(
        axis.title.y = if (!ylabel) element_blank() else element_text(),
        axis.text.y = if (!ylabel) element_blank() else element_text()
      )
  }
  
  ((plotAllCurvesW("GS") | plotAllCurvesW("ieMAS (High Density) + GS", FALSE)) /
      (plotAllCurvesW("PS")  | plotAllCurvesW("ieMAS (High Density) + PS", FALSE))) +
    plot_layout(guides = "collect", axes = "collect")
  
  ggplot2::ggsave(filename = paste0("GS_vs_MAS_Curves_W_", GS_MODEL, ".jpg"),
                  path=output_dir,
                  device = "jpg",
                  width=7,
                  height=5,
                  dpi=600)
  ggplot2::ggsave(filename = paste0("GS_vs_MAS_Curves_W_", GS_MODEL, ".pdf"),
                  path=output_dir,
                  device = "pdf",
                  width=7,
                  height=5)

  plotGainByIe <- function(ie, ylabel=TRUE) {
    ie_cats.df %>%
      dplyr::filter(ie_cat == ie) %>%
      ggplot(aes(x = c, y = gain, group = sel, color=sel)) +
      geom_line(linewidth = 0.3) +
      geom_errorbar(aes(ymin = gain - gain_CI,
                        ymax = gain + gain_CI),
                    width = 0.2,
                    linewidth = 0.2,
                    alpha = 0.6) +
      labs(
        x = "Year",
        y = "Genetic Gain\n(Realized Yield units)",
        title = ie
      ) +
      scale_x_continuous(breaks=seq(from=0, to=20, by=2)) +
      scale_y_continuous(limits = c(minCycleGain, maxCycleGain)) +
      scale_color_sel +
      theme +
      theme(
        axis.title.y = if (!ylabel) element_blank() else element_text(),
        axis.text.y = if (!ylabel) element_blank() else element_text()
      )
  }
  
  (plotGainByIe("Low") | plotGainByIe("Moderate", FALSE) | plotGainByIe("High", FALSE)) + plot_layout(guides = "collect", axes = "collect")
  
  ggplot2::ggsave(filename = paste0("Methods_Curves_Gain_", GS_MODEL, ".jpg"),
                  path=output_dir,
                  device = "jpg",
                  width=10,
                  height=4,
                  dpi=600)
  ggplot2::ggsave(filename = paste0("Methods_Curves_Gain_", GS_MODEL, ".pdf"),
                  path=output_dir,
                  device = "pdf",
                  width=10,
                  height=4)
  
  plotWByIe <- function(ie, ylabel=TRUE) {
    ie_cats.df %>%
      dplyr::filter(ie_cat == ie) %>%
      ggplot(aes(x = c, y = w, group = sel, color=sel)) +
      geom_line(linewidth = 0.3) +
      geom_errorbar(aes(ymin = w - w_CI,
                        ymax = w + w_CI),
                    width = 0.2,
                    linewidth = 0.2,
                    alpha = 0.6) +
      labs(
        x = "Year",
        y = "Breeding Fitness\n(Realized Yield)",
        title = ie
      ) +
      scale_x_continuous(breaks=seq(from=0, to=20, by=2)) +
      scale_y_continuous(limits = c(minCycleW, maxCycleW)) +
      scale_color_sel +
      theme +
      theme(
        axis.title.y = if (!ylabel) element_blank() else element_text(),
        axis.text.y = if (!ylabel) element_blank() else element_text()
      )
  }
  
  (plotWByIe("Low") | plotWByIe("Moderate", FALSE) | plotWByIe("High", FALSE)) + plot_layout(guides = "collect", axes = "collect")
  
  ggplot2::ggsave(filename = paste0("Methods_Curves_W_", GS_MODEL, ".jpg"),
                  path=output_dir,
                  device = "jpg",
                  width=10,
                  height=4,
                  dpi=600)
  ggplot2::ggsave(filename = paste0("Methods_Curves_W_", GS_MODEL, ".pdf"),
                  path=output_dir,
                  device = "pdf",
                  width=10,
                  height=4)
  
  
  geneticGain.df <- gs.df %>%
    dplyr::filter(type == "Admixed") %>%
    dplyr::group_by(founder, rep, type, sel, ie_cat) %>%
    dplyr::summarize(
      "gain_5"  = gain[c == 5],
      "gain_10" = gain[c == 10],
      "gain_20" = gain[c == 20],
      "w_5"  = w[c == 5],
      "w_10" = w[c == 10],
      "w_20" = w[c == 20],
      ie = mean(isoElite),
      .groups = "drop"
    ) %>%
    drop_na() %>%
    tidyr::pivot_longer(cols=c("gain_5", "gain_10", "gain_20", "w_5", "w_10", "w_20"),
                        names_to=c(".value", "cycles"),
                        names_sep = "_") %>%
    dplyr::mutate(cycles=factor(cycles, levels=c("5", "10", "20")))
  
  # calculate confidence intervals and means of each category to determine axis limits
  gain_stats <- geneticGain.df %>%
    dplyr::group_by(ie_cat, cycles, sel) %>%
    dplyr::summarize(
      mean_gain = mean(gain, na.rm = TRUE),
      mean_w = mean(w, na.rm = TRUE),
      ci_gain = qnorm(0.975) * sd(gain, na.rm = TRUE) / sqrt(n()),
      ci_w = qnorm(0.975) * sd(w, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
  
  minDiscreteGain <- min(gain_stats$mean_gain - gain_stats$ci_gain)
  maxDiscreteGain <- max(gain_stats$mean_gain + gain_stats$ci_gain)
  minDiscreteW <- min(gain_stats$mean_w - gain_stats$ci_w)
  maxDiscreteW <- max(gain_stats$mean_w + gain_stats$ci_w)
  
  # Plot the genetic gain per replicate as a boxplot
  plotGainByCycle <- function(df, cycle, ylabel=TRUE) {
    df %>%
      dplyr::filter(cycles==cycle) %>%
      ggplot(aes(x = ie_cat, y = gain, color=sel, group=sel)) +
      # 95% CI
      stat_summary(
        fun.data = mean_cl_normal, 
        geom = "errorbar", 
        position = position_dodge(width = 0.75),
        width = 0.5,
        linewidth = 1
      ) +
      # Point estimate
      stat_summary(
        fun = mean, 
        geom = "point", 
        position = position_dodge(width = 0.75),
        size = 1.5
      ) +
      scale_y_continuous(limits = c(minDiscreteGain, maxDiscreteGain)) +
      labs(
        title  = paste0("Cycles: ", cycle),
        x = "Mean Isoeliteness",
        y = "Genetic Gain\n(Realized Yield units)"
      ) + 
      scale_color_sel +
      theme +
      theme(
        axis.title.y = if (!ylabel) element_blank() else element_text(),
        axis.text.y = if (!ylabel) element_blank() else element_text()
      )
  }
  
  gain5 <- plotGainByCycle(geneticGain.df, "5")
  gain10 <- plotGainByCycle(geneticGain.df, "10", ylabel=FALSE)
  gain20 <- plotGainByCycle(geneticGain.df, "20", ylabel=FALSE)
  
  (gain5 | gain10 | gain20) + plot_layout(guides = "collect", axes = "collect")
  ggplot2::ggsave(filename = paste0("Methods_Discrete_Cycles_Gain_", GS_MODEL, ".jpg"),
                  path=output_dir,
                  device = "jpg",
                  width=10,
                  height=4,
                  dpi=600)
  ggplot2::ggsave(filename = paste0("Methods_Discrete_Cycles_Gain_", GS_MODEL, ".pdf"),
                  path=output_dir,
                  device = "pdf",
                  width=10,
                  height=4)
  
  # Plot the genetic gain per replicate as a boxplot
  plotWByCycle <- function(df, cycle, ylabel=TRUE) {
    df %>%
      dplyr::filter(cycles==cycle) %>%
      ggplot(aes(x = ie_cat, y = w, color=sel, group=sel)) +
      # 95% CI
      stat_summary(
        fun.data = mean_cl_normal, 
        geom = "errorbar", 
        position = position_dodge(width = 0.75),
        width = 0.5,
        linewidth = 1
      ) +
      # Point estimate
      stat_summary(
        fun = mean, 
        geom = "point", 
        position = position_dodge(width = 0.75),
        size = 1.5
      ) +
      scale_y_continuous(limits = c(minDiscreteW, maxDiscreteW)) +
      labs(
        title  = paste0("Cycles: ", cycle),
        x = "Mean Isoeliteness",
        y = "Breeding Fitness\n(Realized Yield)"
      ) + 
      scale_color_sel +
      theme +
      theme(
        axis.title.y = if (!ylabel) element_blank() else element_text(),
        axis.text.y = if (!ylabel) element_blank() else element_text()
      )
  }
  
  w5 <- plotWByCycle(geneticGain.df, "5")
  w10 <- plotWByCycle(geneticGain.df, "10", ylabel=FALSE)
  w20 <- plotWByCycle(geneticGain.df, "20", ylabel=FALSE)
  
  (w5 | w10 | w20) + plot_layout(guides = "collect", axes = "collect")
  ggplot2::ggsave(filename = paste0("Methods_Discrete_Cycles_W_", GS_MODEL, ".jpg"),
                  path=output_dir,
                  device = "jpg",
                  width=10,
                  height=4,
                  dpi=600)
  ggplot2::ggsave(filename = paste0("Methods_Discrete_Cycles_W_", GS_MODEL, ".pdf"),
                  path=output_dir,
                  device = "pdf",
                  width=10,
                  height=4)
  
  # Plot the genetic gain per replicate as a boxplot
  plotGainByIe <- function(df, cycle, ylabel=TRUE) {
    df %>%
      dplyr::filter(cycles==cycle) %>%
      ggplot(aes(x = sel, y = gain, color=ie_cat, group=ie_cat)) +
      # 95% CI
      stat_summary(
        fun.data = mean_cl_normal, 
        geom = "errorbar", 
        position = position_dodge(width = 0.75),
        width = 0.5,
        linewidth = 1
      ) +
      stat_compare_means(
        aes(group=ie_cat),
        label = "p.signif",
        method = "anova",
        hide.ns = TRUE,
        label.y = maxCycleGain*0.9,
        show.legend = FALSE
      ) +
      # Point estimate
      stat_summary(
        fun = mean, 
        geom = "point", 
        position = position_dodge(width = 0.75),
        size = 1.5
      ) +
      scale_y_continuous(limits = c(minDiscreteGain, maxDiscreteGain)) +
      labs(
        title  = paste0("Cycles: ", cycle),
        x = "Mean Isoeliteness",
        y = "Genetic Gain"
      ) + 
      scale_color_manual(
        name = "Mean Isoeliteness",
        values =c(
          "Low" = "#cc0000",
          "Moderate" = "grey40",
          "High" = "#0957d6"
        )) +
      theme +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = if (!ylabel) element_blank() else element_text(),
        axis.text.y = if (!ylabel) element_blank() else element_text()
      )
  }
  
  gain5_MAS <- plotGainByIe(geneticGain.df, "5")
  gain10_MAS <- plotGainByIe(geneticGain.df, "10", ylabel=FALSE)
  gain20_MAS <- plotGainByIe(geneticGain.df, "20", ylabel=FALSE)

  (gain5_MAS | gain10_MAS | gain20_MAS) + plot_layout(guides = "collect", axes = "collect")
  ggplot2::ggsave(filename = paste0("DiscreteCycles_Gain_", GS_MODEL, ".jpg"),
                  path=output_dir,
                  device = "jpg",
                  width=10,
                  height=4,
                  dpi=600)
  ggplot2::ggsave(filename = paste0("DiscreteCycles_Gain_", GS_MODEL, ".pdf"),
                  path=output_dir,
                  device = "pdf",
                  width=6.5,
                  height=2.5)
  
  # Plot the breeding fitness per replicate as a boxplot
  plotWByIe <- function(df, cycle, ylabel=TRUE) {
    df %>%
      dplyr::filter(cycles==cycle) %>%
      ggplot(aes(x = sel, y = w, color=ie_cat, group=ie_cat)) +
      # 95% CI
      stat_summary(
        fun.data = mean_cl_normal, 
        geom = "errorbar", 
        position = position_dodge(width = 0.75),
        width = 0.5,
        linewidth = 1
      ) +
      stat_compare_means(
        aes(group=ie_cat),
        label = "p.signif",
        method = "anova",
        hide.ns = TRUE,
        label.y = maxDiscreteW*0.9,
        show.legend = FALSE
      ) +
      # Point estimate
      stat_summary(
        fun = mean, 
        geom = "point", 
        position = position_dodge(width = 0.75),
        size = 1.5
      ) +
      scale_y_continuous(limits = c(minDiscreteW, maxDiscreteW)) +
      labs(
        title  = paste0("Cycles: ", cycle),
        x = "Mean Isoeliteness",
        y = "Breeding Fitness"
      ) + 
      scale_color_manual(
        name = "Mean Isoeliteness",
        values =c(
          "Low" = "#cc0000",
          "Moderate" = "grey40",
          "High" = "#0957d6"
        )) +
      theme +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = if (!ylabel) element_blank() else element_text(),
        axis.text.y = if (!ylabel) element_blank() else element_text()
      )
  }
  
  w5_MAS <- plotWByIe(geneticGain.df, "5")
  w10_MAS <- plotWByIe(geneticGain.df, "10", ylabel=FALSE)
  w20_MAS <- plotWByIe(geneticGain.df, "20", ylabel=FALSE)
  
  (w5_MAS | w10_MAS | w20_MAS) + plot_layout(guides = "collect", axes = "collect")
  ggplot2::ggsave(filename = paste0("DiscreteCycles_W_", GS_MODEL, ".jpg"),
                  path=output_dir,
                  device = "jpg",
                  width=10,
                  height=4,
                  dpi=600)
  ggplot2::ggsave(filename = paste0("DiscreteCycles_W_", GS_MODEL, ".pdf"),
                  path=output_dir,
                  device = "pdf",
                  width=6.5,
                  height=2.5)
  
  # DETERMINE CORRELATION BETWEEN ISOELITENESS AND GAIN AT EACH CYCLE
  ie_cor.df <- gs.df %>%
    dplyr::filter(type == "Admixed") %>%
    # Grouping only by population selection and chromosome/condition
    dplyr::filter(c > 0) %>%
    dplyr::group_by(sel, c) %>%
    dplyr::summarize(
      rGain = cor(gain, isoElite, method = "pearson", use = "complete.obs"),
      rW = cor(w, isoElite, method = "pearson", use = "complete.obs"),
      sample_size = n(),
      .groups     = "drop"
    )
  
  ie_cor.df %>%
    ggplot(aes(x = c, y = rGain, group = sel, color=sel)) +
    geom_line(linewidth = 0.3) +
    geom_point(size=1) +
    labs(
      x = "Year",
      y = "Correlation (r) between\nIsoeliteness and Genetic Gain",
    ) +
    scale_x_continuous(breaks=seq(from=0, to=20, by=2)) +
    scale_color_sel +
    theme
  
  ggplot2::ggsave(filename = paste0("IeCorrelation_Gain_", GS_MODEL, ".jpg"),
                  path=output_dir,
                  device = "jpg",
                  width=5,
                  height=4,
                  dpi=600)
  ggplot2::ggsave(filename = paste0("IeCorrelation_Gain_", GS_MODEL, ".pdf"),
                  path=output_dir,
                  device = "pdf",
                  width=5,
                  height=4)
  
  ie_cor.df %>%
    ggplot(aes(x = c, y = rW, group = sel, color=sel)) +
    geom_line(linewidth = 0.3) +
    geom_point(size=1) +
    labs(
      x = "Year",
      y = "Correlation (r) between\nIsoeliteness and Breeding Fitness",
    ) +
    scale_x_continuous(breaks=seq(from=0, to=20, by=2)) +
    scale_color_sel +
    theme
  
  ggplot2::ggsave(filename = paste0("IeCorrelation_W_", GS_MODEL, ".jpg"),
                  path=output_dir,
                  device = "jpg",
                  width=5,
                  height=4,
                  dpi=600)
  ggplot2::ggsave(filename = paste0("IeCorrelation_W_", GS_MODEL, ".pdf"),
                  path=output_dir,
                  device = "pdf",
                  width=5,
                  height=4)
  
  # Plot isoeliteness against genetic gain after varying numbers of cycles
  plotGeneticGainVsIe <- function(df, selType, title, ylabel=TRUE) {
    df %>%
      dplyr::filter(sel==selType) %>%
      ggplot(aes(x = ie, y = gain, color = cycles)) +
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
      scale_y_continuous(limits=c(minDiscreteGain, maxDiscreteGain)) +
      theme +
      theme(
        axis.title.y = if (!ylabel) element_blank() else element_text(),
        axis.text.y = if (!ylabel) element_blank() else element_text()
      )
  }

  ((plotGeneticGainVsIe(geneticGain.df, "GS", "GS") |
      plotGeneticGainVsIe(geneticGain.df, "ieMAS (High Density) + GS", "ieMAS + GS", FALSE)) / 
  (plotGeneticGainVsIe(geneticGain.df, "PS", "PS") |
       plotGeneticGainVsIe(geneticGain.df, "ieMAS (High Density) + PS", "ieMAS + PS", FALSE))) +
    plot_layout(guides = "collect", axes = "collect")
  
  ggplot2::ggsave(filename = paste0("Ie_vs_Gain_", GS_MODEL, ".jpg"),
                  path=output_dir,
                  device = "jpg",
                  width=6.5,
                  height=3,
                  dpi=600)
  ggplot2::ggsave(filename = paste0("Ie_vs_Gain_", GS_MODEL, ".pdf"),
                  path=output_dir,
                  device = "pdf",
                  width=6.5,
                  height=3)
  
  # Plot isoeliteness against genetic gain after varying numbers of cycles
  plotWVsIe <- function(df, selType, title, ylabel=TRUE) {
    df %>%
      dplyr::filter(sel==selType) %>%
      ggplot(aes(x = ie, y = w, color = cycles)) +
      geom_point(size=0.5) +
      geom_smooth(method="lm", se=FALSE, aes(color=cycles), linewidth=0.4) +
      sig_cor +
      labs(
        title  = title,
        x = "Mean Isoeliteness",
        y = "Breeding Fitness"
      ) + 
      scale_color_manual(
        name = "Number of Cycles",
        values = c("black", "grey40", "grey70")
      ) +
      scale_y_continuous(limits=c(minDiscreteW, maxDiscreteW)) +
      theme +
      theme(
        axis.title.y = if (!ylabel) element_blank() else element_text(),
        axis.text.y = if (!ylabel) element_blank() else element_text()
      )
  }
  
  ((plotWVsIe(geneticGain.df, "GS", "GS") |
      plotWVsIe(geneticGain.df, "ieMAS (High Density) + GS", "ieMAS + GS", FALSE)) / 
      (plotWVsIe(geneticGain.df, "PS", "PS") |
         plotWVsIe(geneticGain.df, "ieMAS (High Density) + PS", "ieMAS + PS", FALSE))) +
    plot_layout(guides = "collect", axes = "collect")
  
  ggplot2::ggsave(filename = paste0("Ie_vs_W_", GS_MODEL, ".jpg"),
                  path=output_dir,
                  device = "jpg",
                  width=6.5,
                  height=3,
                  dpi=600)
  ggplot2::ggsave(filename = paste0("Ie_vs_W_", GS_MODEL, ".pdf"),
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
    dplyr::filter(c > 1) %>%
    ggplot(aes(x = c, y = het, color = qtlType)) +
    geom_line() +
    geom_point() +
    labs(
      x = "Year",
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
  
  # Plot R over cycles
  cycleMean.df %>%
    dplyr::filter(pop=="Admixed GS") %>%
    ggplot(aes(x = c, y = r)) +
    geom_line() +
    geom_point() +
    labs(
      x = "Year",
      y = "GWP accuracy for breeding fitness (r)"
    ) +
    theme
  ggplot2::ggsave(filename = paste0("r_per_cycle_", GS_MODEL, ".jpg"),
                  path=output_dir,
                  device = "jpg",
                  width=3.5,
                  height=3,
                  dpi=600)
  ggplot2::ggsave(filename = paste0("r_per_cycle_", GS_MODEL, ".pdf"),
                  path=output_dir,
                  device = "pdf",
                  width=3.5,
                  height=3)
  
  geneticGain.df %>%
    dplyr::filter(type=="Admixed") %>%
    dplyr::group_by(sel, cycles, ie_cat) %>%
    dplyr::summarize(meanGain=mean(gain),
                     varGain=var(gain))%>%
    write.csv(file.path(output_dir, "genetic_gain.csv"), quote=FALSE)
  
}

# ANOVA
sig.df <- geneticGain.df %>%
  dplyr::filter(type=="Admixed",
                sel %in% c("GS", "ieMAS (High Density) + GS"))


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
