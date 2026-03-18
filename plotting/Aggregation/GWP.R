

gwp_1.df <- res.df %>%
  dplyr::select(isoElite_T1, gwpR_T1, gwpR_T1_Pop2, type, qtl) %>%
  dplyr::rename(
    "isoElite" = isoElite_T1,
    "gwpR" = gwpR_T1,
    "gwpR_Pop2" = gwpR_T1_Pop2
  )

gwp_2.df <- res.df %>%
  dplyr::select(isoElite_T2, gwpR_T2, gwpR_T2_Pop2, type, qtl) %>%
  dplyr::rename(
    "isoElite" = isoElite_T2,
    "gwpR" = gwpR_T2,
    "gwpR_Pop2" = gwpR_T2_Pop2
  )

gwp_att.df <- dplyr::bind_rows(gwp_1.df, gwp_2.df)

gwp_att_admixed <- gwp_att.df %>%
  dplyr::filter(type=="Admixed") %>%
  dplyr::select(-gwpR_Pop2)

gwp_att_pop1 <- gwp_att.df %>%
  dplyr::filter(type=="Unadmixed") %>%
  dplyr::select(isoElite, gwpR, qtl) %>%
  dplyr::mutate(type="Subpopulation 1")

gwp_att_pop2 <- gwp_att.df %>%
  dplyr::filter(type=="Unadmixed") %>%
  dplyr::select(isoElite, gwpR_Pop2, qtl) %>%
  dplyr::rename(
    "gwpR"=gwpR_Pop2
  ) %>%
  dplyr::mutate(type="Subpopulation 2")

gwp_att.df <- dplyr::bind_rows(gwp_att_admixed, gwp_att_pop1, gwp_att_pop2) %>%
  drop_na()
gwp_att.df$type <- as.factor(gwp_att.df$type)
gwp_att.df$qtl <- as.factor(gwp_att.df$qtl)


gwp_w.df <- res.df %>%
  dplyr::select(isoElite_Att, gwpR_W, gwpR_W_Pop2, type, qtl) %>%
  dplyr::rename("isoElite" = isoElite_Att)

gwp_w_admixed <- gwp_w.df %>%
  dplyr::filter(type=="Admixed") %>%
  dplyr::select(-gwpR_W_Pop2)

gwp_w_pop1 <- gwp_w.df %>%
  dplyr::filter(type=="Unadmixed") %>%
  dplyr::select(isoElite, gwpR_W, qtl) %>%
  dplyr::mutate(type="Subpopulation 1")

gwp_w_pop2 <- gwp_w.df %>%
  dplyr::filter(type=="Unadmixed") %>%
  dplyr::select(isoElite, gwpR_W_Pop2, qtl) %>%
  dplyr::rename(
    "gwpR_W"=gwpR_W_Pop2
  ) %>%
  dplyr::mutate(type="Subpopulation 2")

gwp_w.df <- dplyr::bind_rows(gwp_w_admixed, gwp_w_pop1, gwp_w_pop2) %>%
  drop_na()

gwp_w.df$type <- factor(gwp_w.df$type,
                        levels=c("Subpopulation 1", "Subpopulation 2", "Admixed"))
gwp_w.df$qtl <- as.factor(gwp_w.df$qtl)
gwp_w.df <- gwp_w.df %>%
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
#gwp_att.df %>%
#  ggplot(aes(x=qtl, y=gwpR, fill=type)) +
#  geom_boxplot(
#    outlier.shape=NA,
#    position = position_dodge(width = 1),
#    linewidth = 0.1,
#  ) +
#  stat_compare_means(
#    aes(group = test_pop),
#    label = "p.signif",
#    hide.ns = FALSE,
#  ) +
#  scale_fill_gwp +
#  scale_y_continuous(expand=c(0,0.1)) +
#  labs(x="QTL per Attained Trait", y="GWP Accuracy (r)", title="Attained Traits") +
#  theme

#ggplot2::ggsave(filename = paste0("gwp_att.jpg"),
#                path=output_dir,
#                device = "jpg",
#                width=3.5,
#                height=3,
#                dpi=600)

g <- gwp_w.df %>%
  ggplot(aes(x=qtl, y=gwpR_W, fill=type)) +
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
g
gIE <- gwp_w.df %>%
  dplyr::filter(type=="Admixed") %>%
  ggplot(aes(x=isoElite, y=gwpR_W)) +
  geom_point(aes(color=qtl), alpha=a, size=s) +
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
  
ggplot2::ggsave(filename = paste0("gwp_ie.jpg"),
                path=output_dir,
                device = "jpg",
                width=6,
                height=4,
                dpi=600)
