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
library(purrr)
library(tibble)
library(tidyr)


setwd("~/Documents/CSU/FitnessLandscapes/output/GWP/cycle1_noupdate_merged_6_22_1237reps")
#setwd("~/Documents/CSU/FitnessLandscapes/output/GWP/cycle1_update_6_17_400reps")
output_dir <- getwd()

#write.table(RS.df, file.path(output_dir, "RRBLUP/rs_results.csv"), col.names=TRUE, quote=FALSE, sep=",")
#write.table(RIL.df, file.path(output_dir, "RRBLUP/ril_results.csv"), col.names=TRUE, quote=FALSE, sep=",")

RIL.df <- rbind(read.csv("RRBLUP/ril_results.csv"))
#                read.csv("../cycle1_noupdate_6_21_465reps/RRBLUP/ril_results.csv"))

RS.df <- rbind(read.csv("RRBLUP/rs_results.csv"))
#               read.csv("../cycle1_noupdate_6_21_465reps/RRBLUP/rs_results.csv"))

read_costs <- function(fname) {
  costs.df <- read.csv(fname, check.names=FALSE) %>%
    pivot_longer(
      cols = -c(season, year),
      names_to = "sel",
      values_to = "cost"
    )
  return (costs.df)
}

costs_high_labor_high_seq <- read_costs("../../../data/costs_high_labor_high_seq.csv") %>%
  dplyr::rename(cost_high_high = cost)
costs_high_labor_low_seq <- read_costs("../../../data/costs_high_labor_low_seq.csv") %>%
  dplyr::rename(cost_high_low = cost)
costs_low_labor_high_seq <- read_costs("../../../data/costs_low_labor_high_seq.csv") %>%
  dplyr::rename(cost_low_high = cost)
costs_low_labor_low_seq <- read_costs("../../../data/costs_low_labor_low_seq.csv") %>%
  dplyr::rename(cost_low_low = cost)

# The first columns should be retained
join_cols <- names(costs_high_labor_high_seq)[1:3]

costs.df <- list(
  costs_high_labor_high_seq,
  costs_high_labor_low_seq,
  costs_low_labor_high_seq,
  costs_low_labor_low_seq
) %>%
  reduce(full_join, by = join_cols)


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

scale_color_gwp <- scale_color_manual(name = "Scenario",
                                values = c("CV (Single Subpopulation)" = "steelblue2",
                                           "CV (Both Subpopulations)" = "navy",
                                           "Unadmixed" = "#7EBD4F",
                                           "Admixed" = "goldenrod1",
                                           "Cross-population" = "#A66A28"),
                                breaks = c("CV (Single Subpopulation)",
                                           "CV (Both Subpopulations)",
                                           "Unadmixed",
                                           "Admixed", 
                                           "Cross-population"))

color_scheme <- c(
  "PS" = "#8E44C6",
  "isoMAS-PS" = "#eec7ff",
  "GS" = "#2EBD65",
  "GS (No Update)" = "#b7edba",
  "isoMAS-GS (Low Density)" = "#FFD4A8",
  "isoMAS-GS" = "#F06A00",
  "isoMAS-GS (No Update)" = "#C14B00",
  "isoMAS-GS (Perfect Markers)" = "#7A2900")

scale_fill_sel <- scale_fill_manual(name = "Scheme",
                                      values = color_scheme)

scale_color_sel <- scale_color_manual(name = "Scheme",
                                      values = color_scheme,
                                      limits = c("PS",
                                                 "isoMAS-PS",
                                                 "GS",
                                                 "isoMAS-GS"))
                                                 #"isoMAS-GS (Low Density)",
                                                 #"isoMAS-GS (Perfect Markers)",
                                                 #"GS (No Update)",
                                                 #"isoMAS-GS (No Update)"))

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
                        "CV (Single Subpopulation)",
                        "CV (Both Subpopulations)",
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
    labs(x="Scenario", y=expression("GWP accuracy for realized yield (r"[MG]*")")) +
    theme +
    theme(
      axis.text.x = element_blank())
  
  ggplot2::ggsave(filename = paste0("gwp_accuracy_", GS_MODEL, ".jpg"),
                  path=output_dir,
                  device = "jpg",
                  width=5,
                  height=4,
                  dpi=600)
  ggplot2::ggsave(filename = paste0("gwp_accuracy_", GS_MODEL, ".pdf"),
                  path=output_dir,
                  device = "pdf",
                  width=5,
                  height=4)
  
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
      y = expression("GWP accuracy for realized yield (r"[MG]*")"),
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
         y = expression("GWP accuracy for realized yield (r"[MG]*")")) +
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
  
  (ie_des | att_des) + plot_annotation(tag_levels='a',
                                     theme = theme(plot.tag = element_text(family = "Helvetica", size = 6)))
  
  ggplot2::ggsave(filename = paste0("gwp_ieDes_", GS_MODEL, ".jpg"),
                  path=output_dir,
                  device = "jpg",
                  width=8,
                  height=4,
                  dpi=600)
  ggplot2::ggsave(filename = paste0("gwp_ieDes_", GS_MODEL, ".pdf"),
                  path=output_dir,
                  device = "pdf",
                  width=8,
                  height=4)
  highIe <- 0.9
  #quantile(gs.df$isoElite[gs.df$type == "Admixed"], 0.8)
  lowIe <- 0.8 #
  #quantile(gs.df$isoElite[gs.df$type == "Admixed"], 0.2)
  
  startingVals <- RS.df %>%
    dplyr::filter(model == GS_MODEL) %>%
    dplyr::filter(type == "Admixed") %>%
    dplyr::filter(c == 0)
  
  mean(startingVals$wGV)
  
  # Create a unique identifier for the RILTYPE_SELECTION TYPE
  gs.df <- RS.df %>%
    dplyr::filter(model == GS_MODEL) %>%
    dplyr::filter(
      case_when(
        sel %in% c("PS", "PS_ieMAS") ~ season %% 6 == 0,
        TRUE ~ season %% 2 == 0
      )
    ) %>%
    dplyr::mutate(year = season/2) %>%
    dplyr::filter(year <= 18) %>%
    dplyr::group_by(founder, rep, type, sel) %>%
    dplyr::filter(wGV[1] >= n.minW & wGV[1] <= n.maxW) %>%
    dplyr::mutate(gain = wGV - wGV[1]) %>%
    dplyr::mutate(sel = case_when(
      sel == "ieMAS" ~ "isoMAS-GS",
      sel == "ieMAS_low" ~ "isoMAS-GS (Low Density)",
      sel == "ieMAS_perfect" ~ "isoMAS-GS (Perfect Markers)",
      sel == "GS_noUpdate" ~ "GS (No Update)",
      sel == "PS_ieMAS" ~ "isoMAS-PS",
      sel == "ieMAS_noUpdate" ~ "isoMAS-GS (No Update)",
      TRUE ~ sel
    )) %>%
    dplyr::mutate(
      ie_cat = dplyr::case_when(
        isoElite >= highIe ~ "High",
        isoElite > lowIe & isoElite < highIe ~ "Moderate",
        isoElite <=  lowIe ~ "Low"
      ),
      ie_cat=factor(ie_cat, levels=c("Low", "Moderate", "High"))) %>%
    dplyr::mutate(
      pt_fill = color_scheme[sel]
    ) %>%
    dplyr::mutate(pop=paste0(type, " ", sel)) %>%
    dplyr::mutate(
      type=factor(type, levels=c("Admixed",
                                 "Unadmixed")),
      sel=factor(sel, levels=c("PS",
                               "isoMAS-PS",
                               "GS",
                               "isoMAS-GS (Low Density)",
                               "isoMAS-GS",
                               "isoMAS-GS (Perfect Markers)",
                               "GS (No Update)",
                               "isoMAS-GS (No Update)")),
      pop=factor(pop, levels=c("Admixed PS",
                               "Admixed isoMAS-PS",
                               "Admixed GS",
                               "Admixed isoMAS-GS (Low Density)",
                               "Admixed isoMAS-GS",
                               "Admixed isoMAS-GS (Perfect Markers)",
                               "Admixed GS (No Update)",
                               "Admixed isoMAS-GS (No Update)",
                               "Unadmixed PS",
                               "Unadmixed isoMAS-PS",
                               "Unadmixed GS",
                               "Unadmixed isoMAS-GS (Low Density)",
                               "Unadmixed isoMAS-GS",
                               "Unadmixed isoMAS-GS (Perfect Markers)",
                               "Unadmixed GS (No Update)",
                               "Unadmixed isoMAS-GS (No Update)")))
  
  gs.df %>% 
    filter(type == "Admixed") %>% 
    group_by(founder) %>%
    summarize(unique_founders = n_distinct(founder)) %>%
    nrow()
  
  gs.df %>%
    filter(type == "Admixed") %>%
    group_by(ie_cat) %>%
    summarize(unique_founders = n_distinct(founder))
  
  # Show the relationship between mean isoeliteness and realized yield in RIL family
  gs.df %>%
    dplyr::filter(type == "Admixed") %>%
    dplyr::filter(c == 0) %>%
    dplyr::group_by(founder, rep) %>%
    dplyr::summarize(ie=mean(isoElite),
                     wGV=mean(wGV)) %>%
    dplyr::ungroup() %>%
    ggplot(aes(x=ie, y=wGV)) +
    geom_point(size=0.5) +
    geom_smooth(method="lm", se=FALSE, linewidth=0.4, color= "black") +
    guides(color = guide_legend(override.aes = list(shape = 16, linetype = 1, size=2))) +
    sig_cor +
    labs(x="Mean Isoeliteness",
         y="Mean Realized Yield in RIL Family") +
    theme +
    theme(
      legend.position = "none",
    )
  
  ggplot2::ggsave(filename = paste0("w_vs_ie_filtered_", GS_MODEL, ".jpg"),
                  path=output_dir,
                  device = "jpg",
                  width=5,
                  height=4,
                  dpi=600)
  ggplot2::ggsave(filename = paste0("w_vs_ie_filtered_", GS_MODEL, ".pdf"),
                  path=output_dir,
                  device = "pdf",
                  width=5,
                  height=4)

  # Calculate average realized yield per pop per cycle
  yearMean.df <- gs.df %>%
    dplyr::group_by(pop, type, sel, year, ie_cat, pt_fill) %>%
    dplyr::summarize(gain_CI  = 1.96*(sd(gain) / sqrt(n())),
                     gain = mean(gain),
                     wGV_CI  = 1.96*(sd(wGV) / sqrt(n())),
                     wGV = mean(wGV),
                     genHt = mean(genome_het),
                     attHt = mean(attained_het),
                     desHt = mean(desired_het),
                     r_CI = 1.96*(sd(r) / sqrt(n())),
                     r = mean(r, na.rm=TRUE),
                     gvar = mean(gvar, na.rm=TRUE)) %>%
    dplyr::left_join(costs.df, by = c("year", "sel")) %>%
    dplyr::mutate(unit_cost_high_high = cost_high_high / gain / 1000,
                  unit_cost_high_low = cost_high_low / gain / 1000,
                  unit_cost_low_high = cost_low_high / gain / 1000,
                  unit_cost_low_low = cost_low_low / gain / 1000)
    
  
  minYearGain <- min(yearMean.df$gain-yearMean.df$gain_CI, na.rm=TRUE)
  maxYearGain <- max(yearMean.df$gain+yearMean.df$gain_CI, na.rm=TRUE)
  
  minYearW <- min(yearMean.df$wGV-yearMean.df$wGV_CI, na.rm=TRUE)
  maxYearW <- max(yearMean.df$wGV+yearMean.df$wGV_CI, na.rm=TRUE)
  
  # Plot for QTL == 10 only
  yearMeanByIe <- function(ie) {
    yearMean.df %>%
      dplyr::filter(sel %in% c("GS", "PS")) %>%
      dplyr::filter(ie_cat == ie) %>%
      ggplot(aes(x = year, y = gain, group = pop, color=pt_fill)) +
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
          "Admixed PS" = 21,
          #"Admixed ieMAS + PS" = 21,
          "Admixed GS" = 21,
          #"Admixed ieMAS + GS" = 21,
          "Unadmixed PS" = 21,
          #"Unadmixed ieMAS + PS" = 21,
          "Unadmixed GS" = 21
          #"Unadmixed ieMAS + GS" = 21
        ),
        guide = guide_legend(override.aes = list(
          fill  = c(color_scheme[["PS"]],
                    #color_scheme[["ieMAS + PS"]],
                    color_scheme[["GS"]],
                    #color_scheme[["ieMAS + GS"]],
                    "white",
                    #"white",
                    #"white",
                    "white"),
          color = c(color_scheme[["PS"]],
                    #color_scheme[["ieMAS + PS"]],
                    color_scheme[["GS"]],
                    #color_scheme[["ieMAS + GS"]],
                    color_scheme[["PS"]],
                    #color_scheme[["ieMAS + PS"]],
                    color_scheme[["GS"]]),
          #color_scheme[["ieMAS + GS"]]),
          size  = 3
        ))
      ) +
      scale_fill_identity() +
      scale_color_identity() +
      scale_x_continuous(breaks=seq(from=0, to=n.Y, by=3)) +
      #scale_y_continuous(limits = c(minYearW, maxYearW)) +
      labs(
        title = paste0(ie, " Founder isoMAS"),
        x = "Year",
        y = "Genetic Gain\n(Realized Yield Units)"
      ) +
      theme
  }
  yearMeanByIe("High")

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

  
  plotAllCurvesGain <- function(selType, ylabel = TRUE) {
    yearMean.df %>%
      dplyr::filter(type == "Admixed") %>%
      dplyr::filter(sel == selType) %>%
      ggplot(aes(x = year, y = gain, group = ie_cat, color = ie_cat)) +
      geom_line(linewidth = 0.8) +
      geom_point(size=1) +
      geom_errorbar(aes(ymin = gain - gain_CI, ymax = gain + gain_CI),
                    width = 0.2,
                    linewidth = 0.3,
                    alpha = 0.6) +
      scale_color_manual(
        name = "Founder isoMAS",
        values =c(
          "Low" = "#cc0000",
          "Moderate" = "grey40",
          "High" = "#0957d6"
        )) +
      labs(
        x = "Year",
        y = "Genetic Gain\n(Realized Yield Units)",
        title = selType
      ) +
      scale_x_continuous(breaks=seq(from=0, to=n.Y, by=3)) +
      scale_y_continuous(limits = c(-2, maxYearGain)) +
      theme +
      theme(
        axis.title.y = if (!ylabel) element_blank() else element_text(),
        axis.text.y = if (!ylabel) element_blank() else element_text()
      )
  }
  
  (plotAllCurvesGain("GS") | (plotAllCurvesGain("PS", FALSE))) +
    plot_layout(guides = "collect", axes = "collect") +
    plot_annotation(tag_levels='a',
                    theme = theme(plot.tag = element_text(family = "Helvetica", size = 6)))
  
  
  ggplot2::ggsave(filename = paste0("Founder_isoMAS_", GS_MODEL, ".jpg"),
                  path=output_dir,
                  device = "jpg",
                  width=8,
                  height=4,
                  dpi=600)
  ggplot2::ggsave(filename = paste0("Founder_isoMAS_", GS_MODEL, ".pdf"),
                  path=output_dir,
                  device = "pdf",
                  width=8,
                  height=4)
  
  (plotAllCurvesGain("isoMAS-GS") | (plotAllCurvesGain("isoMAS-PS", FALSE))) +
    plot_layout(guides = "collect", axes = "collect") +
    plot_annotation(tag_levels='a',
                    theme = theme(plot.tag = element_text(family = "Helvetica", size = 6)))
  
  
  ggplot2::ggsave(filename = paste0("Founder_isoMAS_Supp_", GS_MODEL, ".jpg"),
                  path=output_dir,
                  device = "jpg",
                  width=8,
                  height=4,
                  dpi=600)
  ggplot2::ggsave(filename = paste0("Founder_isoMAS_Supp_", GS_MODEL, ".pdf"),
                  path=output_dir,
                  device = "pdf",
                  width=8,
                  height=4)
  
  plotGainByIe <- function(df, ie, withControls=FALSE, ylabel=TRUE) {
    if (!withControls) {
      df <- df %>% dplyr::filter(sel %in% c("GS", "PS", "isoMAS-PS", "isoMAS-GS"))
    }
    df %>%
      dplyr::filter(type == "Admixed") %>%
      dplyr::filter(ie_cat == ie) %>%
      ggplot(aes(x = year, y = gain, group = sel, color=sel)) +
      geom_line(linewidth = 0.5) +
      geom_errorbar(aes(ymin = gain - gain_CI,
                        ymax = gain + gain_CI),
                    width = 0.2,
                    linewidth = 0.2,
                    alpha = 0.6) +
      labs(
        x = "Year",
        y = "Genetic Gain\n(Realized Yield Units)",
        title = paste0(ie, " Founder isoMAS"),
      ) +
      scale_x_continuous(breaks=seq(from=0, to=n.Y, by=3)) +
      scale_y_continuous(limits = c(-2, maxYearGain)) +
      scale_color_sel +
      theme +
      theme(
        axis.title.y = if (!ylabel) element_blank() else element_text(),
        axis.text.y = if (!ylabel) element_blank() else element_text()
      )
  }

  (plotGainByIe(yearMean.df, "Low", withControls=TRUE) |
      plotGainByIe(yearMean.df, "Moderate", withControls=TRUE, ylabel=FALSE) |
      plotGainByIe(yearMean.df, "High", withControls=TRUE, ylabel=FALSE)) +
    plot_layout(guides = "collect", axes = "collect") +
    plot_annotation(tag_levels='a',
                    theme = theme(plot.tag = element_text(family = "Helvetica", size = 6)))
  
  ggplot2::ggsave(filename = paste0("Methods_Curves_Gain_Controls_", GS_MODEL, ".jpg"),
                  path=output_dir,
                  device = "jpg",
                  width=10,
                  height=4,
                  dpi=600)
  ggplot2::ggsave(filename = paste0("Methods_Curves_Gain_Controls_", GS_MODEL, ".pdf"),
                  path=output_dir,
                  device = "pdf",
                  width=10,
                  height=4)

  plotCostByIe <- function(df, ie, withControls=FALSE, cost_var="unit_cost_low_high", ylabel=TRUE) {
    if (cost_var == "unit_cost_low_high") {
      yrange <- c(0,16)
    } else if (cost_var == "unit_cost_low_low") {
      yrange <- c(0, 16)
    } else if (cost_var == "unit_cost_high_high") {
      yrange <- c(7, 31)
    } else if (cost_var == "unit_cost_high_low") {
      yrange <- c(7, 31)
    }
    
    df %>% dplyr::filter(sel %in% c("GS", "PS", "isoMAS-PS", "isoMAS-GS")) %>%
      dplyr::filter(year > 0) %>%
      dplyr::filter(type == "Admixed") %>%
      dplyr::filter(ie_cat == ie) %>%
      ggplot(aes(x = year, y = .data[[cost_var]], group = sel, color=sel)) +
      geom_line(linewidth = 0.5) +
      labs(
        x = "Year",
        y = "Cost / Units Gain\n(In $1000s of USD)",
        title = paste0(ie, " Founder isoMAS"),
      ) +
      scale_x_continuous(breaks=seq(from=0, to=n.Y, by=3)) +
      scale_y_continuous(limits = yrange) +
      scale_color_sel +
      theme +
      theme(
        axis.title.y = if (!ylabel) element_blank() else element_text(),
        #axis.text.y = if (!ylabel) element_blank() else element_text()
      )
  }
  top_row <- (plotGainByIe(yearMean.df, "Low") |
                plotGainByIe(yearMean.df, "Moderate", ylabel=FALSE) |
                plotGainByIe(yearMean.df, "High", ylabel=FALSE))
  
  bottom_row <- (plotCostByIe(yearMean.df, "Low") |
                   plotCostByIe(yearMean.df, "Moderate", ylabel=FALSE) |
                   plotCostByIe(yearMean.df, "High", ylabel=FALSE)) & 
    theme(legend.position = "none")
  (top_row / bottom_row) +
    plot_layout(guides = "collect", axes = "collect") +
    plot_annotation(tag_levels='a',
                    theme = theme(plot.tag = element_text(family = "Helvetica", size = 6)))
  
  ggplot2::ggsave(filename = paste0("Methods_Curves_Gain_", GS_MODEL, ".jpg"),
                  path=output_dir,
                  device = "jpg",
                  width=9,
                  height=7,
                  dpi=600)
  ggplot2::ggsave(filename = paste0("Methods_Curves_Gain_", GS_MODEL, ".pdf"),
                  path=output_dir,
                  device = "pdf",
                  width=9,
                  height=7)
  
  cost_low_high <- (plotCostByIe(yearMean.df, "Low", cost_var="unit_cost_low_high") |
                   plotCostByIe(yearMean.df, "Moderate", cost_var="unit_cost_low_high", ylabel=FALSE) |
                   plotCostByIe(yearMean.df, "High", cost_var="unit_cost_low_high", ylabel=FALSE)) & 
    theme(legend.position = "none")
  
  cost_low_low <- (plotCostByIe(yearMean.df, "Low", cost_var="unit_cost_low_low") |
                       plotCostByIe(yearMean.df, "Moderate", cost_var="unit_cost_low_low", ylabel=FALSE) |
                       plotCostByIe(yearMean.df, "High", cost_var="unit_cost_low_low", ylabel=FALSE)) & 
    theme(legend.position = "none")
  cost_high_high <- (plotCostByIe(yearMean.df, "Low", cost_var="unit_cost_high_high") |
                       plotCostByIe(yearMean.df, "Moderate", cost_var="unit_cost_high_high", ylabel=FALSE) |
                       plotCostByIe(yearMean.df, "High", cost_var="unit_cost_high_high", ylabel=FALSE)) & 
    theme(legend.position = "none")
  cost_high_low <- (plotCostByIe(yearMean.df, "Low", cost_var="unit_cost_high_low") |
                       plotCostByIe(yearMean.df, "Moderate", cost_var="unit_cost_high_low", ylabel=FALSE) |
                       plotCostByIe(yearMean.df, "High", cost_var="unit_cost_high_low", ylabel=FALSE)) & 
    theme(legend.position = "none")
  (cost_low_high / cost_low_low / cost_high_high / cost_high_low) +
    plot_layout(guides = "collect", axes = "collect") +
    plot_annotation(tag_levels='a',
                    theme = theme(plot.tag = element_text(family = "Helvetica", size = 6)))
  
  ggplot2::ggsave(filename = paste0("Costs_", GS_MODEL, ".jpg"),
                  path=output_dir,
                  device = "jpg",
                  width=9,
                  height=12,
                  dpi=600)
  ggplot2::ggsave(filename = paste0("Costs_", GS_MODEL, ".pdf"),
                  path=output_dir,
                  device = "pdf",
                  width=9,
                  height=12)
  
  geneticGainAll.df <- gs.df %>%
    dplyr::filter(year %in% c(3, 9, 18)) %>% 
    dplyr::group_by(founder, rep, type, sel, ie_cat, year) %>%
    dplyr::summarize(
      gain = mean(gain, na.rm = TRUE),
      wGV  = mean(wGV, na.rm = TRUE),
      ie   = mean(isoElite, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(year = factor(year, levels = c("3", "9", "18")))
  
  geneticGain.df <- geneticGainAll.df %>%
    dplyr::filter(type == "Admixed")
  
  global_bounds <- geneticGain.df %>%
    dplyr::group_by(ie_cat, year, sel) %>%
    dplyr::summarize(
      margin = 1.96 * (sd(gain, na.rm = TRUE) / sqrt(n())),
      low    = mean(gain, na.rm = TRUE) - margin,
      high   = mean(gain, na.rm = TRUE) + margin,
      .groups = "drop"
    )

  minDiscreteGain <- min(global_bounds$low, na.rm=TRUE)
  maxDiscreteGain <- max(global_bounds$high, na.rm=TRUE)

  # Plot the genetic gain per replicate as a boxplot
  plotGainByYear <- function(df, years, ylabel=TRUE) {
    df %>%
      dplyr::filter(year == years) %>%
      ggplot(aes(x = ie_cat, y = gain, color=sel, group=sel)) +
      # 95% CI
      stat_summary(
        fun.data = mean_cl_normal, 
        geom = "errorbar", 
        position = position_dodge(width = 0.75),
        width = 0.5,
        linewidth = 0.5
      ) +
      # Point estimate
      stat_summary(
        fun = mean, 
        geom = "point", 
        position = position_dodge(width = 0.75),
        size = 1
      ) +
      coord_cartesian(ylim = c(minDiscreteGain, maxDiscreteGain)) +
      labs(
        title  = paste0("Years: ", years),
        x = "Founder isoMAS",
        y = "Genetic Gain\n(Realized Yield Units)"
      ) + 
      scale_color_sel +
      theme +
      theme(
        axis.title.y = if (!ylabel) element_blank() else element_text(),
        axis.text.y = if (!ylabel) element_blank() else element_text()
      )
  }
  
  w3 <- plotGainByYear(geneticGain.df, "3")
  w9 <- plotGainByYear(geneticGain.df, "9", FALSE)
  w18 <- plotGainByYear(geneticGain.df, "18", FALSE)
  
  (w3 | w9 | w18) + plot_layout(guides = "collect", axes = "collect") +
    plot_annotation(tag_levels='a',
                    theme = theme(plot.tag = element_text(family = "Helvetica", size = 6)))
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
  
  # Plot the realized yield per replicate as a boxplot
  plotGainByIe <- function(df, years, ylabel=TRUE) {
    df %>%
      dplyr::filter(year==years) %>%
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
        label.y = maxDiscreteGain*0.97,
        show.legend = FALSE
      ) +
      # Point estimate
      stat_summary(
        fun = mean, 
        geom = "point", 
        position = position_dodge(width = 0.75),
        size = 1.5
      ) +
      coord_cartesian(ylim = c(minDiscreteGain, maxDiscreteGain)) +
      labs(
        title  = paste0("Years: ", years),
        x = "Mean Isoeliteness",
        y = "Genetic Gain\n(Realized Yield Units)"
      ) + 
      scale_color_manual(
        name = "Founder isoMAS",
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
  
  w3_MAS <- plotGainByIe(geneticGain.df, "3")
  w9_MAS <- plotGainByIe(geneticGain.df, "9", ylabel=FALSE)
  w18_MAS <- plotGainByIe(geneticGain.df, "18", ylabel=FALSE)
  
  (w3_MAS | w9_MAS | w18_MAS) + plot_layout(guides = "collect", axes = "collect")
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
  
  # DETERMINE CORRELATION BETWEEN ISOELITENESS AND GAIN AT EACH CYCLE
  ie_cor.df <- gs.df %>%
    dplyr::filter(type == "Admixed") %>%
    # Grouping only by population selection and chromosome/condition
    dplyr::filter(year >= 3) %>%
    dplyr::group_by(sel, year) %>%
    dplyr::summarize(
      rGain = cor(gain, isoElite, method = "pearson", use = "complete.obs"),
      rW = cor(w, isoElite, method = "pearson", use = "complete.obs"),
      sample_size = n(),
      .groups     = "drop"
    )
  
  ie_cor.df %>%
    ggplot(aes(x = year, y = rGain, group = sel, color=sel)) +
    geom_line(linewidth = 0.3) +
    geom_point(size=1) +
    labs(
      x = "Year",
      y = "Correlation (r) between\nIsoeliteness and Genetic Gain",
    ) +
    scale_x_continuous(breaks=seq(from=0, to=n.Y, by=3)) +
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
  
  minHt <- min(c(yearMean.df$genHt, yearMean.df$attHt, yearMean.df$desHt))
  maxHt <- max(c(yearMean.df$genHt, yearMean.df$attHt, yearMean.df$desHt))
  
  ht.df <- yearMean.df %>%
    pivot_longer(cols=c(genHt, attHt, desHt),
                 names_to="qtlType",
                 values_to="het")
  ht.df$qtlType <- factor(ht.df$qtlType,
                          levels=c("attHt", "desHt", "genHt"))
  
  hetPlot <- function(targetSel, targetIe, ylabel=TRUE) {
    ht.df %>%
      dplyr::filter(sel == targetSel,
                    type == "Admixed",
                    ie_cat == targetIe) %>%
      dplyr::filter(year > 0) %>%
      ggplot(aes(x = year, y = het, color = qtlType)) +
      geom_line() +
      geom_point() +
      labs(
        title = paste0(targetIe),
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
      coord_cartesian(ylim = c(0, 0.08)) +
      theme +
      theme(
        axis.title.y = if (!ylabel) element_blank() else element_text(),
        axis.text.y = if (!ylabel) element_blank() else element_text()
      )
  }
  ((hetPlot("GS", "Low") | hetPlot("GS", "Moderate", FALSE) | hetPlot("GS", "High", FALSE)) /
    (hetPlot("isoMAS-GS", "Low") | hetPlot("isoMAS-GS", "Moderate", FALSE) | hetPlot("isoMAS-GS", "High", FALSE))) +
    plot_layout(guides = "collect", axes = "collect")
  
  ggplot2::ggsave(filename = paste0("Heterozygosity_", GS_MODEL, ".jpg"),
                  path=output_dir,
                  device = "jpg",
                  width=8,
                  height=6,
                  dpi=600)
  ggplot2::ggsave(filename = paste0("Heterozygosity_", GS_MODEL, ".pdf"),
                  path=output_dir,
                  device = "pdf",
                  width=8,
                  height=6)
  # Plot R over cycles
  
  plotRPerCycle <- function(ie, ylabel=TRUE) {
    yearMean.df %>%
      dplyr::filter(ie_cat==ie) %>%
      dplyr::filter(type == "Admixed") %>%
      dplyr::filter(sel %in% c("isoMAS-GS", "GS")) %>%
      ggplot(aes(x = year, y = r, color = sel, group = sel)) +
      geom_line(na.rm = TRUE) +
      geom_point(na.rm = TRUE) +
      geom_errorbar(aes(ymin = r - r_CI,
                        ymax = r + r_CI),
                    width = 0.2,
                    linewidth = 0.2,
                    alpha = 0.6) +
      labs(
        title = paste0(ie, " Founder isoMAS"),
        x = "Year",
        y = expression("GWP accuracy for realized yield (r"[MG]*")"),
      ) +
      coord_cartesian(ylim = c(-0.15, 0.55)) +
      scale_color_sel +
      theme +
      theme(
        axis.title.y = if (!ylabel) element_blank() else element_text(),
        axis.text.y = if (!ylabel) element_blank() else element_text()
      )
  }
  
  (plotRPerCycle("Low") | plotRPerCycle("Moderate", FALSE) | plotRPerCycle("High", FALSE)) +
  plot_layout(guides = "collect", axes = "collect") +
    plot_annotation(tag_levels='a',
                    theme = theme(plot.tag = element_text(family = "Helvetica", size = 6)))
  ggplot2::ggsave(filename = paste0("r_per_cycle_", GS_MODEL, ".jpg"),
                  path=output_dir,
                  device = "jpg",
                  width=10,
                  height=4,
                  dpi=600)
  ggplot2::ggsave(filename = paste0("r_per_cycle_", GS_MODEL, ".pdf"),
                  path=output_dir,
                  device = "pdf",
                  width=10,
                  height=4)
  
  geneticGainAll.df %>%
    dplyr::group_by(type, sel, year, ie_cat) %>%
    dplyr::summarize(meanGain=mean(gain),
                     varGain=var(gain))%>%
    write.csv(file.path(output_dir, "genetic_gain.csv"), quote=FALSE)
  
  
  
}

# ANOVA
sig.df <- geneticGain.df %>%
  dplyr::filter(type=="Admixed",
                sel %in% c("GS", "ieMAS (High Density) + GS"))

# Test whether there is more variance in GS than MAS
var.test(gain ~ sel, data=sig.df %>% filter(year==3,
                                            ie_cat=="Low"))
var.test(gain ~ sel, data=sig.df %>% filter(cycles==5,
                                            ie_cat=="Moderate"))
var.test(gain ~ sel, data=sig.df %>% filter(cycles==5,
                                            ie_cat=="High"))

# SIGNIFICANT
anova(lm(gain ~ sel, data=sig.df %>% filter(year==3,
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
