# Title: ALLOGENIC RANK
# Author: Ted Monyak
# Description: Determine, for a particular set of parameters, the rank of the largest
# segregating QTL from the initial founder trait architecture series that is segregating in
# an admixed RIL family
# Allogenic = differentially fixed (i.e. not isogenic)
# This script should be called only from AggregateFigures.R

theme <- theme_minimal(base_size = 8,
                       base_family="Helvetica") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(),
    panel.background = element_rect(fill = "white", color = "black"),
    plot.margin= unit(c(10,10,10,10), unit="pt"))

# Aggregate attained traits 1 and 2
allo_1.df <- res.df %>%
  dplyr::filter(type=="Admixed") %>%
  dplyr::select(qtl, relRank_T1, isoElite_T1) %>%
  dplyr::rename("relRank"=relRank_T1,
                "isoElite"=isoElite_T1)
allo_2.df <- res.df %>%
  dplyr::filter(type=="Admixed") %>%
  dplyr::select(qtl, relRank_T2, isoElite_T2) %>%
  dplyr::rename("relRank"=relRank_T2,
                "isoElite"=isoElite_T2)
allo.df <- rbind(
  allo_1.df[is.finite(allo_1.df$relRank), ],
  allo_2.df[is.finite(allo_2.df$relRank), ]
)

# Make a boxplot of the mean relative ranks, sorted by number of initial QTL
allo.df %>%
  ggplot(aes(x=qtl, y=relRank, fill=qtl)) +
  geom_boxplot(
    outlier.shape = NA,
    position = position_dodge(width = 1)
  ) +
  scale_fill_qtl +
  scale_y_continuous(
    limits=c(0,10),
    expand=c(0,1)) +
  stat_compare_means(method="t.test",
                     comparisons = list(c("10", "20"), c("20", "50")),
                     label="p.signif",
                     hide.ns=FALSE,
                     label.y = c(9,9),
                     bracket.size = 0) +
  labs(x="QTL per Attained Trait ", y="Mean Rank\nFirst Allogenic QTL") +
  theme

ggplot2::ggsave(filename = file.path(output_dir, "relRank.jpg"),
                device = "jpg",
                height=3,
                width=3,
                dpi=600)
ggplot2::ggsave(filename = file.path(output_dir, "relRank.pdf"),
                device = "pdf",
                height=3,
                width=3)
