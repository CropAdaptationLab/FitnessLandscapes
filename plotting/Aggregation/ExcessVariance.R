theme <- theme_minimal(base_size = 8,
                         base_family="Helvetica") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = "black"),
    plot.margin = margin(2,0,0,2, "pt"),
    axis.title.x = element_text(size=7),
    axis.title.y = element_text(size=7),
    plot.title = element_text(hjust = 0.5, size=7),
    aspect.ratio=1)

s <- 0.2

vMin <- min(c(res.df$ev_T1, res.df$ev_T2, res.df$ev_T3, res.df$ev_Suit, res.df$ev_W), na.rm=TRUE, inf.rm=TRUE)
vMax <- max(c(res.df$ev_T1, res.df$ev_T2, res.df$ev_T3, res.df$ev_Suit, res.df$ev_W), na.rm=TRUE, inf.rm=TRUE)

# Create ev_T plots
v1 <- res.df %>%
  dplyr::filter(!is.na(ev_T1) & !is.infinite(ev_T1)) %>%
  ggplot(aes(x=qtl, y=ev_T1, fill=type)) +
  geom_boxplot(
    outlier.shape=NA,
    position = position_dodge(width = 1),
    linewidth = 0.1,
  ) +
  scale_y_continuous(
    limits=c(vMin,vMax),
    expand=c(0,0.1)) +
  scale_fill +
  labs(x="QTL per\nAttained Trait", y="EV", title="Attained Trait 1") +
  theme
v2 <- res.df %>%
  dplyr::filter(!is.na(ev_T2) & !is.infinite(ev_T2)) %>%
  ggplot(aes(x=qtl, y=ev_T2, fill=type)) +
  geom_boxplot(
    outlier.shape=NA,
    position = position_dodge(width = 1),
    linewidth = 0.1,
  ) +
  scale_y_continuous(
    limits=c(vMin,vMax),
    expand=c(0,0.1)) +
  scale_fill +
  labs(x="QTL per\nAttained Trait", y="EV", title="Attained Trait 2") +
  theme

vSuit <- res.df %>%
  dplyr::filter(!is.na(ev_Suit) & !is.infinite(ev_Suit)) %>%
  ggplot(aes(x=qtl, y=ev_Suit, fill=type)) +
  geom_boxplot(
    outlier.shape=NA,
    position = position_dodge(width = 1),
    linewidth = 0.1,
  ) +
  scale_y_continuous(
    limits=c(vMin,vMax),
    expand=c(0,0.1)) +
  scale_fill +
  labs(x="QTL per\nAttained Trait", y="EV", title="Suitability") +
  theme

v3 <- res.df %>%
  dplyr::filter(!is.na(ev_T3) & !is.infinite(ev_T3)) %>%
  ggplot(aes(x=qtl, y=ev_T3, fill=type)) +
  geom_boxplot(
    outlier.shape=NA,
    position = position_dodge(width = 1),
    linewidth = 0.1,
  ) +
  scale_y_continuous(
    limits=c(vMin,vMax),
    expand=c(0,0.1)) +
  scale_fill +
  labs(x="QTL per\nAttained Trait", y="EV", title="Desired Trait") +
  theme

vW <- res.df %>%
  dplyr::filter(!is.na(ev_W) & !is.infinite(ev_W)) %>%
  ggplot(aes(x=qtl, y=ev_W, fill=type)) +
  geom_boxplot(
    outlier.shape=NA,
    position = position_dodge(width = 1),
    linewidth = 0.1,
  ) +
  scale_y_continuous(
    limits=c(vMin,vMax),
    expand=c(0,0.1)) +
  scale_fill +
  labs(x="QTL per\nAttained Trait", y="EV", title="Breeding Fitness") +
  theme

v <- (v1|v2|vSuit|v3|vW) +
  plot_layout(guides='collect', axes='collect')

# SCATTERPLOTS

iMin <- 0.4
iMax <- 1

# Create ev_T plots
isoVar1 <- res.df %>%
  dplyr::filter(type=="Admixed") %>%
  dplyr::filter(!is.na(ev_T1) & !is.infinite(ev_T1)) %>%
  ggplot(aes(x=isoElite_T1, y=ev_T1, color=qtl)) +
  geom_point(alpha=a, size=s) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  scale_x_continuous(
    limits=c(iMin,iMax),
    breaks=seq(0.4,1, by=0.2)) +
  scale_y_continuous(
    limits=c(vMin,vMax*1.05),
    expand=c(0,0.1)) +
  scale_color +
  labs(x="Attained Trait 1\nIsoeliteness", y="EV", title="Attained Trait 1") +
  theme
isoVar2 <- res.df %>%
  dplyr::filter(type=="Admixed") %>%
  dplyr::filter(!is.na(ev_T2) & !is.infinite(ev_T2)) %>%
  ggplot(aes(x=isoElite_T2, y=ev_T2, color=qtl)) +
  geom_point(alpha=a, size=s) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  scale_x_continuous(
    limits=c(iMin,iMax),
    breaks=seq(0.4,1, by=0.2)) +
  scale_y_continuous(
    limits=c(vMin,vMax*1.05),
    expand=c(0,0.1)) +
  scale_color +
  labs(x="Attained Trait 2\nIsoeliteness", y="EV", title="Attained Trait 2") +
  theme

isoVarSuit <- res.df %>%
  dplyr::filter(type=="Admixed") %>%
  dplyr::filter(!is.na(ev_Suit) & !is.infinite(ev_Suit)) %>%
  ggplot(aes(x=isoElite_Att, y=ev_Suit, color=qtl)) +
  geom_point(alpha=a, size=s) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  scale_x_continuous(
    limits=c(iMin,iMax),
    breaks=seq(0.4,1, by=0.2)) +
  scale_y_continuous(
    limits=c(vMin,vMax*1.05),
    expand=c(0,0.1)) +
  scale_color +
  labs(x="Mean\nIsoeliteness", y="EV", title="Suitability") +
  theme

isoVar3 <- res.df %>%
  dplyr::filter(type=="Admixed") %>%
  dplyr::filter(!is.na(ev_T3) & !is.infinite(ev_T3)) %>%
  ggplot(aes(x=isoElite_Att, y=ev_T3, color=qtl)) +
  geom_point(alpha=a, size=s) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  scale_x_continuous(
    limits=c(iMin,iMax),
    breaks=seq(0.4,1, by=0.2)) +
  scale_y_continuous(
    limits=c(vMin,vMax*1.05),
    expand=c(0,0.1)) +
  scale_color +
  labs(x="Mean\nIsoeliteness", y="EV", title="Desired Trait") +
  theme

isoVarW <- res.df %>%
  dplyr::filter(type=="Admixed") %>%
  dplyr::filter(!is.na(ev_W) & !is.infinite(ev_W)) %>%
  ggplot(aes(x=isoElite_Att, y=ev_W, color=qtl)) +
  geom_point(alpha=a, size=s) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  scale_x_continuous(
    limits=c(iMin,iMax),
    breaks=seq(0.4,1, by=0.2)) +
  scale_y_continuous(
    limits=c(vMin,vMax*1.05),
    expand=c(0,0.1)) +
  scale_color +
  labs(x="Mean\nIsoeliteness", y="EV", title="Breeding Fitness") +
  theme


vi <- (isoVar1|isoVar2|isoVarSuit|isoVar3|isoVarW) +
  plot_layout(guides='collect', axes='collect')

(v / vi) + plot_annotation(tag_levels=list(c('f', 'g', 'h', 'i\ ', 'j', 'k', 'l\ ', 'm', 'n', 'o')),
                           theme = theme(plot.tag = element_text(family = "Helvetica", size = 6)))
ggplot2::ggsave(filename = paste0("excess_variance.jpg"),
                path=output_dir,
                device = "jpg",
                width=6.5, # use 6.5 for 4 panels
                height=3,
                dpi=600)

ggplot2::ggsave(filename = paste0("excess_variance.pdf"),
                path=output_dir,
                device = "pdf",
                width=6.5, # use 6.5 for 4 panels
                height=3)
