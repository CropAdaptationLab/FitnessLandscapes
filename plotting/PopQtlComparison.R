two_small <- read.csv("Ne_500_qtl_2_h2_0.1_aa_0_2026-01-28_15_56/sim_results.csv")
five_small <- read.csv("Ne_500_qtl_5_h2_0.1_aa_0_2026-01-28_17_13/sim_results.csv")
two_large <- read.csv("Ne_1000_qtl_2_h2_0.1_aa_0_2026-01-28_18_32/sim_results.csv")
five_large <- read.csv("Ne_1000_qtl_5_h2_0.1_aa_0_2026-01-28_20_11/sim_results.csv")

res.df <- rbind(two_small, five_small, two_large, five_large)
res.df$popSize <- as.factor(res.df$popSize)
res.df$qtl <- as.factor(res.df$qtl)
res.df$aa <- as.factor(res.df$aa)
res.df <- res.df %>%
  dplyr::mutate(diff=nSigQtlInt-nSigQtlFit)
res.df <- res.df %>%
  filter(popSize==1000)


theme <- theme(
  axis.title.x = element_blank(),
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title.y = element_text(family="Helvetica", size=12),
  axis.text.y = element_text(angle = 0, hjust=1, size=10),
  plot.title = element_text(family="Helvetica", size=14, hjust = 0.5),
  legend.text = element_text(family="Helvetica", size=10),
  legend.title = element_text(family="Helvetica", size=12),
  legend.key = element_rect(linewidth=0.05),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_rect(fill = "white", color = "black"),
  plot.margin= unit(c(10,10,10,10), unit="pt"),
  aspect.ratio = 1)

qtl_ie <- res.df %>%
  ggplot(aes(x=qtl, y=ie, fill=qtl)) +
  geom_violin(trim=TRUE) +
  geom_boxplot(width=0.1, outlier.shape=NA, show.legend=FALSE) +
  stat_compare_means(method="t.test",
                     label="p.signif",
                     hide.ns=FALSE) +
  #scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_fill_manual(name = "QTL",
                    values = c("2" = "yellow2",
                               "5" = "lightcoral")) +
  labs(title="#QTL vs Iso-eliteness Score Pop Size 1000",
       y="Iso-eliteness score") +
  theme

qtl_ie
fname <- file.path(output_dir, "qtl_ie_pop1000.jpg")
ggplot2::ggsave(filename = fname,
                device = "jpg")

qtl_var <- res.df %>%
  ggplot(aes(x=qtl, y=varRat, fill=qtl)) +
  geom_violin(trim=TRUE) +
  geom_boxplot(width=0.1, outlier.shape=NA, show.legend=FALSE) +
  stat_compare_means(method="t.test",
                     label="p.signif",
                     hide.ns=FALSE) +
  #scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_fill_manual(name = "QTL",
                    values = c("2" = "yellow2",
                               "5" = "lightcoral")) +
  labs(title="#QTL vs Variance Ratio Pop Size 1000",
       y="Iso-eliteness score") +
  theme

qtl_var
fname <- file.path(output_dir, "qtl_var_pop1000.jpg")
ggplot2::ggsave(filename = fname,
                device = "jpg")

qtl_sigQtl <- res.df %>%
  ggplot(aes(x=qtl, y=nSigQtl, fill=qtl)) +
  geom_violin(trim=TRUE) +
  geom_boxplot(width=0.1, outlier.shape=NA, show.legend=FALSE) +
  stat_compare_means(method="t.test",
                     label="p.signif",
                     hide.ns=FALSE) +
  #scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_fill_manual(name = "QTL",
                    values = c("2" = "yellow2",
                               "5" = "lightcoral")) +
  labs(title="#QTL vs # Sig Qtl Pop Size 1000",
       y="# LOD peaks") +
  theme

qtl_sigQtl
fname <- file.path(output_dir, "qtl_lod_pop1000.jpg")
ggplot2::ggsave(filename = fname,
                device = "jpg")

qtl_diff <- res.df %>%
  ggplot(aes(x=qtl, y=diff, fill=qtl)) +
  geom_violin(trim=TRUE) +
  geom_boxplot(width=0.1, outlier.shape=NA, show.legend=FALSE) +
  stat_compare_means(method="t.test",
                     label="p.signif",
                     hide.ns=FALSE) +
  #scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_fill_manual(name = "QTL",
                    values = c("2" = "yellow2",
                               "5" = "lightcoral")) +
  labs(title="#QTL vs # Sig Qtl Pop Size 1000",
       y="Difference in Interaction vs Fitness Peaks") +
  theme

qtl_diff
fname <- file.path(output_dir, "qtl_diff_pop1000.jpg")
ggplot2::ggsave(filename = fname,
                device = "jpg")
