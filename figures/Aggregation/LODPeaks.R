theme <- theme_minimal(base_size = 8,
                         base_family="Helvetica") +
  theme(
    legend.position = "bottom",
    legend.title.position = "left",
    legend.direction="vertical",
    legend.text=element_text(size=5),
    legend.title=element_text(size=8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = "black"),
    axis.text = element_text(),
    plot.margin=unit(c(0,0,0,0), unit="pt"))

s <- 0.2

sig_cor <- stat_cor(method="pearson",
                    aes(label=paste(
                      sub("R", "r", after_stat(r.label)),
                      ifelse(after_stat(p) < 0.001, '"***"',
                             ifelse(after_stat(p) < 0.01,  '"**"',
                                    ifelse(after_stat(p) < 0.05,  '"*"', '"ns"'))),
                      sep="~"
                    )),
                    parse=TRUE)

qMax <- max(c(res.df$nLod_T1, res.df$nLod_T2, res.df$nLod_T3, res.df$nLod_W, res.df$nLod_Suit), na.rm=TRUE, inf.rm=TRUE)

q1 <- res.df %>%
  ggplot(aes(x=qtl, y=nLod_T1, fill=type)) +
  geom_boxplot(
    outlier.shape=NA,
    position = position_dodge(width = 1),
    linewidth = 0.1,
  ) +
  stat_compare_means(aes(group=type),
                     method="t.test",
                     label="p.signif",
                     hide.ns=FALSE,
                     label.y = c(10.5,10.5)) +
  scale_y_continuous(
    limits=c(0,qMax),
    expand=c(0,0.2)) +
  scale_fill +
  labs(x="QTL per Attained Trait", y="Significant Peaks") +
  theme

q2 <- res.df %>%
  ggplot(aes(x=qtl, y=nLod_T2, fill=type)) +
  geom_boxplot(
    outlier.shape=NA,
    position = position_dodge(width = 1),
    linewidth = 0.1,
  ) +
  stat_compare_means(aes(group=type),
                     method="t.test",
                     label="p.signif",
                     hide.ns=FALSE,
                     label.y = c(10.5,10.5)) +
  scale_y_continuous(
    limits=c(0,qMax),
    expand=c(0,0.2)) +
  scale_fill +
  labs(x="QTL per Attained Trait", y="Significant Peaks") +
  theme

qSuit <- res.df %>%
  ggplot(aes(x=qtl, y=nLod_Suit, fill=type)) +
  geom_boxplot(
    outlier.shape=NA,
    position = position_dodge(width = 1),
    linewidth = 0.1,
  ) +
  stat_compare_means(aes(group=type),
                     method="t.test",
                     label="p.signif",
                     hide.ns=FALSE,
                     label.y = c(10.5,10.5)) +
  scale_y_continuous(
    limits=c(0,qMax),
    expand=c(0,0.2)) +
  scale_fill +
  labs(x="QTL per Attained Trait", y="Significant Peaks") +
  theme

q3 <- res.df %>%
  ggplot(aes(x=qtl, y=nLod_T3, fill=type)) +
  geom_boxplot(
    outlier.shape=NA,
    position = position_dodge(width = 1),
    linewidth = 0.1,
  ) +
  stat_compare_means(aes(group=type),
                     method="t.test",
                     label="p.signif",
                     hide.ns=FALSE,
                     label.y = c(10.5,10.5)) +
  scale_y_continuous(
    limits=c(0,qMax),
    expand=c(0,0.2)) +
  scale_fill +
  labs(x="QTL per Attained Trait", y="Significant Peaks") +
  theme

qW <- res.df %>%
  ggplot(aes(x=qtl, y=nLod_W, fill=type)) +
  geom_boxplot(
    outlier.shape=NA,
    position = position_dodge(width = 1),
    linewidth = 0.1,
  ) +
  stat_compare_means(aes(group=type),
                     method="t.test",
                     label="p.signif",
                     hide.ns=FALSE,
                     label.y = c(10.5,10.5)) +
  scale_y_continuous(
    limits=c(0,qMax),
    expand=c(0,0.2)) +
  scale_fill +
  labs(x="QTL per Attained Trait", y="Significant Peaks") +
  theme

q <- (q1/q2/qSuit/q3/qW) +
  plot_layout(guides='collect') &
  theme(legend.position='bottom')

iMin <- 0.4
iMax <- 1

# Create ev_T plots
isoLod1 <- res.df %>%
  dplyr::filter(type=="Admixed") %>%
  ggplot(aes(x=isoElite_T1, y=nLod_T1)) +
  geom_point(aes(color=qtl), alpha=a, size=s) +
  geom_smooth(method="lm", se=FALSE, color="black", linewidth=0.4) +
  sig_cor +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  scale_x_continuous(
    limits=c(iMin,iMax),
    breaks=seq(0.4,1, by=0.2)) +
  scale_y_continuous(
    limits=c(0,qMax),
    expand=c(0,0.2)) +
  scale_color +
  labs(x="Attained Trait 1 Isoeliteness", y="Significant Peaks") +
  theme

isoLod2 <- res.df %>%
  dplyr::filter(type=="Admixed") %>%
  ggplot(aes(x=isoElite_T2, y=nLod_T2)) +
  geom_point(aes(color=qtl), alpha=a, size=s) +
  geom_smooth(method="lm", se=FALSE, color="black", linewidth=0.4) +
  sig_cor +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  scale_x_continuous(
    limits=c(iMin,iMax),
    breaks=seq(0.4,1, by=0.2)) +
  scale_y_continuous(
    limits=c(0,qMax),
    expand=c(0,0.2)) +
  scale_color +
  labs(x="Attained Trait 2 Isoeliteness", y="Significant Peaks") +
  theme

isoLodSuit <- res.df %>%
  dplyr::filter(type=="Admixed") %>%
  ggplot(aes(x=isoElite_Att, y=nLod_Suit)) +
  geom_point(aes(color=qtl), alpha=a, size=s) +
  geom_smooth(method="lm", se=FALSE, color="black", linewidth=0.4) +
  sig_cor +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  scale_x_continuous(
    limits=c(iMin,iMax),
    breaks=seq(0.4,1, by=0.2)) +
  scale_y_continuous(
    limits=c(0,qMax),
    expand=c(0,0.2)) +
  scale_color +
  labs(x="Mean Isoeliteness", y="Significant Peaks") +
  theme

isoLod3 <- res.df %>%
  dplyr::filter(type=="Admixed") %>%
  ggplot(aes(x=isoElite_Att, y=nLod_T3)) +
  geom_point(aes(color=qtl), alpha=a, size=s) +
  geom_smooth(method="lm", se=FALSE, color="black", linewidth=0.4) +
  sig_cor +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  scale_x_continuous(
    limits=c(iMin,iMax),
    breaks=seq(0.4,1, by=0.2)) +
  scale_y_continuous(
    limits=c(0,qMax),
    expand=c(0,0.2)) +
  scale_color +
  labs(x="Mean Isoeliteness", y="Significant Peaks") +
  theme

isoLodW <- res.df %>%
  dplyr::filter(type=="Admixed") %>%
  ggplot(aes(x=isoElite_Att, y=nLod_W)) +
  geom_point(aes(color=qtl), alpha=a, size=s) +
  geom_smooth(method="lm", se=FALSE, color="black", linewidth=0.4) +
  sig_cor +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  scale_x_continuous(
    limits=c(iMin,iMax),
    breaks=seq(0.4,1, by=0.2)) +
  scale_y_continuous(
    limits=c(0,qMax),
    expand=c(0,0.2)) +
  scale_color +
  labs(x="Mean Isoeliteness", y="Significant Peaks") +
  theme

li <- (isoLod1/isoLod2/isoLodSuit/isoLod3/isoLodW) +
  plot_layout(guides='collect') &
  theme(legend.position='bottom')

(q | li) + plot_annotation(tag_levels=list(c('f', 'g', 'h', 'i\ ', 'j', 'k', 'l\ ', 'm', 'n', 'o')),
                           theme = theme(plot.tag = element_text(family = "Helvetica", size = 6)))

ggplot2::ggsave(filename = paste0("lod_peaks.jpg"),
                path=output_dir,
                device = "jpg",
                width=3.5, # use 6.5 for 4 panels
                height=8,
                dpi=600)

ggplot2::ggsave(filename = paste0("lod_peaks.pdf"),
                path=output_dir,
                device = "pdf",
                width=3.5, # use 6.5 for 4 panels
                height=8)
theme <- theme + theme(
    plot.margin=unit(c(0,0,0,10), unit="pt"),
    legend.position = "right",
    legend.title.position = "top",
    legend.direction="vertical")

# Create emergent lod peaks plot
intW <- res.df %>%
  ggplot(aes(x=qtl, y=nLod_Int, fill=type)) +
  geom_boxplot(
    outlier.shape=NA,
    position = position_dodge(width = 1),
    linewidth = 0.1,
  ) +
  stat_compare_means(aes(group=type),
                     method="t.test",
                     label="p.signif",
                     hide.ns=FALSE,
                     label.y=c(7,7)) +
  scale_y_continuous(
    limits=c(0,9),
    breaks=seq(0,8, by=2)) +
  scale_fill +
  labs(x="QTL per Attained Trait", y="Breeding Fitness Interactions") +
  theme

ggplot2::ggsave(filename = file.path(output_dir, paste0("interaction.jpg")),
                device = "jpg",
                height=2,
                width=3,
                dpi=600)

ggplot2::ggsave(filename = file.path(output_dir, paste0("interaction.pdf")),
                device = "pdf",
                height=2,
                width=3)

intW_ie <- res.df %>%
  dplyr::filter(type=="Admixed") %>%
  ggplot(aes(x=isoElite_Att, y=nLod_Int)) +
  geom_point(aes(color=qtl), alpha=a, size=s) +
  geom_smooth(method="lm", se=FALSE, color="black", linewidth=0.4) +
  sig_cor +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  scale_color +
  labs(x="Mean Isoeliteness", y="Breeding Fitness Interactions") +
  theme

intW_ie
fname <- file.path(output_dir, paste0("interaction_ie.jpg"))
ggplot2::ggsave(filename = fname,
                device = "jpg",
                height=2,
                width=3,
                dpi=600)
fname <- file.path(output_dir, paste0("interaction_ie.pdf"))
ggplot2::ggsave(filename = fname,
                device = "pdf",
                height=2,
                width=3)
