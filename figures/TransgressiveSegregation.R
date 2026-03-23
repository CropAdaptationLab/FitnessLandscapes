# Title: TRANSGRESSIVE SEGREGATION
# Author: Ted Monyak
# Description: Plots the excessive transgressive segregation for each of the traits,
# comparing the purelines to the admixed RIL families
# Assumes that CreateIndependentPops has been run, and purelines have been created
# from each subpopulation

if (saveTraitPlots) {
  theme <- theme_minimal(base_size = 8,
                base_family="Helvetica") +
    theme(
    plot.margin = margin(0,0,0,5, "pt"),
    axis.text.x = element_text(angle = 0, hjust=1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.text = element_text(),
    legend.title = element_text(),
    legend.key = element_rect(linewidth=0.05),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = "black"),
    aspect.ratio = 1)
  
  # Join all the phenotypes together
  all_pheno <- rbind(RIL_pheno, pureline_pheno)
  
  # Create density plots for each of the attained traits
  ev1 <- all_pheno %>%
    ggplot(aes(Trait1, fill=Population, color=Population)) +
    scale_color_manual(values=c("#CC0000", "#3C78D8", "gold"),
                       aesthetics=c("color", "fill")) +
    geom_density(size=0.5, alpha=0.3) +
    xlim(-n.initTraitVal,n.initTraitVal) +
    labs(x="Attained Trait 1", y="Density") +
    annotate("text",
             label=paste0("EV = ", round(ev_T1,2)),
             color="black",
             x = -Inf,
             y = Inf,
             hjust=-0.1,
             vjust=1.5,
             family="Helvetica",
             size=2) +
    theme
  
  ev2 <- all_pheno %>%
    ggplot(aes(Trait2, fill=Population, color=Population)) +
    scale_color_manual(values=c("#CC0000", "#3C78D8", "gold"),
                       aesthetics=c("color", "fill")) +
    geom_density(size=0.5, alpha=0.3) +
    xlim(-n.initTraitVal,n.initTraitVal) +
    labs(x="Attained Trait 2", y="Density") +
    annotate("text",
             label=paste0("EV = ", round(ev_T2,2)),
             color="black",
             x = -Inf,
             y = Inf,
             hjust=-0.1,
             vjust=1.5,
             family="Helvetica",
             size=2) +
    theme
  
  # Suitability density plot
  evSuit <- all_pheno %>%
    ggplot(aes(Suitability, fill=Population, color=Population)) +
    scale_color_manual(values=c("#CC0000", "#3C78D8", "gold"),
                       aesthetics=c("color", "fill")) +
    geom_density(size=0.5, alpha=0.3) +
    labs(x="Suitability", y="Density") +
    annotate("text",
             label=paste0("EV = ", round(ev_Suit,2)),
             color="black",
             x = -Inf,
             y = Inf,
             hjust=-0.1,
             vjust=1.5,
             family="Helvetica",
             size=2) +
    theme
  
  # Yield potential density plot
  ev3 <- all_pheno %>%
    ggplot(aes(Trait3, fill=Population, color=Population)) +
    scale_color_manual(values=c("#CC0000", "#3C78D8", "gold"),
                       aesthetics=c("color", "fill")) +
    geom_density(size=0.5, alpha=0.3) +
    labs(x="Desired Trait", y="Density") +
    annotate("text",
             label=paste0("EV = ", round(ev_T3,2)),
             color="black",
             x = -Inf,
             y = Inf,
             hjust=-0.1,
             vjust=1.5,
             family="Helvetica",
             size=2) +
    theme
  
  # Breeding fitness density plot
  evW <- all_pheno %>%
    ggplot(aes(W, fill=Population, color=Population)) +
    scale_color_manual(values=c("#CC0000", "#3C78D8", "gold"),
                       aesthetics=c("color", "fill")) +
    geom_density(size=0.5, alpha=0.3) +
    labs(x="Breeding Fitness", y="Density") +
    annotate("text",
             label=paste0("EV = ", round(ev_W,2)),
             color="black",
             x = -Inf,
             y = Inf,
             hjust=-0.1,
             vjust=1.5,
             family="Helvetica",
             size=2) +
    theme

  # Join all the plots together in a multipanel
  (ev1|ev2|evSuit|ev3|evW) +
    plot_layout(guides='collect', axes='collect') +
    plot_annotation(tag_levels='a')
  ggplot2::ggsave(filename = "trans_seg.jpg",
                  path=ril_dir,
                  device = "jpg",
                  width=6.5,
                  height=2,
                  dpi=600)
  ggplot2::ggsave(filename = "trans_seg.pdf",
                  path=ril_dir,
                  device = "pdf",
                  width=6.5,
                  height=2)
  
  plotTraitArchitecture(pop=RIL, trait=1, popName="RIL Trait 1")
  ggplot2::ggsave(filename = "RIL_traitarchitecture_1.jpg",
                  path=ril_dir,
                  device = "jpg",
                  width=10,
                  height=7)
  
  
  plotTraitArchitecture(pop=RIL, trait=2, popName="RIL Trait 2")
  ggplot2::ggsave(filename = "RIL_traitarchitecture_2.jpg",
                  path=ril_dir,
                  device = "jpg",
                  width=10,
                  height=7)
  
  
  
}
