# Title: CROSS POPULATION
# Author: Ted Monyak
# Description: This script produces summary plots for cross population prediction
# with data generated in GwpPipeline.R

# Pull out the admixed RIL family results
gwp_admixed <- res.df %>%
  dplyr::filter(type=="Admixed") %>%
  dplyr::select(isoElite, gwpR, qtl, type)

# Ungroup the data from both unadmixed families and designate them by the
# subpopulation from which they were derived
gwp_p1 <- res.df %>%
  dplyr::filter(type=="Unadmixed") %>%
  dplyr::select(isoElite, gwpR, qtl) %>%
  dplyr::mutate(type="Subpopulation 1")

gwp_p2 <- res.df %>%
  dplyr::filter(type=="Unadmixed") %>%
  dplyr::select(isoElite, gwpR_P2, qtl) %>%
  dplyr::rename(
    "gwpR"=gwpR_P2
  ) %>%
  dplyr::mutate(type="Subpopulation 2")

# Merge data from all RIL families
gwp.df <- dplyr::bind_rows(gwp_admixed, gwp_p1, gwp_p2) %>%
  drop_na()

gwp.df$type <- factor(gwp.df$type,
                      levels=c("Subpopulation 1", "Subpopulation 2", "Admixed"))
gwp.df$qtl <- as.factor(gwp.df$qtl)

# Test population should only have 2 categories
gwp.df <- gwp.df %>%
  dplyr::mutate(test_pop = ifelse(type=="Admixed", "Admixed", "Unadmixed"))

theme <- theme_minimal(base_size = 8,
                       base_family="Helvetica") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(),
    panel.background = element_rect(fill = "white", color = "black"),
    plot.margin= unit(c(10,10,10,10), unit="pt"),
    legend.position = "bottom",
    legend.title.position = "left",
    legend.direction="vertical")

scale_fill_gwp <- scale_fill_manual(name = "Training Population",
                                values = c("Admixed" = "gold2",
                                           "Subpopulation 1" = "#CC0000",
                                           "Subpopulation 2" = "#3C78D8"),
                                breaks = c("Subpopulation 1", "Subpopulation 2", "Admixed"))

scale_color <- scale_color_manual(name = "QTL per\nAttained Trait",
                                  values = c("10" = "#4A1A6B",
                                             "20" = "#9B59B6",
                                             "50" = "#D7B8F3"))

# Plot the GWP results by RIL family type
g <- gwp.df %>%
  ggplot(aes(x=qtl, y=gwpR, fill=type)) +
  geom_boxplot(
    outlier.shape=NA,
    position = position_dodge(width = 1),
    linewidth = 0.1,
  ) +
  stat_compare_means(
    aes(group = test_pop),
    label = "p.signif",
    hide.ns = FALSE,
  ) +
  scale_fill_gwp +
  scale_y_continuous(expand=c(0,0.1)) +
  labs(x="QTL per Attained Trait", y="GWP accuracy for breeding fitness (r)", title="Breeding Fitness") +
  theme

# Plot isoeliteness against GWP accuracy
gIE <- gwp.df %>%
  dplyr::filter(type=="Admixed") %>%
  ggplot(aes(x=isoElite, y=gwpR)) +
  geom_point(aes(color=qtl), alpha=0.6) +
  geom_smooth(method="lm", se=FALSE, color="black", linewidth=0.4) +
  stat_cor(method="pearson", 
           aes(label=paste(
             sub("R", "r", after_stat(r.label)), 
             after_stat(p.label), 
             sep="~`,`~"
           ))) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  scale_color +
  labs(x="Mean Isoeliteness", y="GWP accuracy for breeding fitness (r)", title="Breeding Fitness") +
  theme
gIE
(g|gIE) + plot_annotation(tag_levels=list(c('a', 'b')))
  
ggplot2::ggsave(filename = "gwp_ie.jpg",
                path=output_dir,
                device = "jpg",
                width=6,
                height=4,
                dpi=600)
ggplot2::ggsave(filename = "gwp_ie.pdf",
                path=output_dir,
                device = "pdf",
                width=6,
                height=4)
