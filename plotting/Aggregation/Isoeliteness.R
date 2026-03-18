theme <- theme_minimal(base_size = 8,
                         base_family="Helvetica") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(),
    panel.background = element_rect(fill = "white", color = "black"),
    plot.margin= unit(c(10,10,10,10), unit="pt"))

# Create iso-eliteness plot
res.df %>%
  dplyr::filter(!is.na(isoElite_Att) & !is.infinite(isoElite_Att)) %>%
  ggplot(aes(x=qtl, y=isoElite_Att, fill=type)) +
  geom_boxplot(
    outlier.shape=NA,
    position = position_dodge(width = 1),
    linewidth = 0.1,
  ) +
  stat_compare_means(aes(group=type),
                     method="t.test",
                     label="p.signif",
                     hide.ns=FALSE) +
  scale_fill +
  scale_y_continuous(expand=c(0,0.1)) +
  labs(x="QTL per Attained Trait", y="Attained Trait Isoeliteness") +
  theme
output_dir

ggplot2::ggsave(filename = file.path(output_dir, paste0("isoeliteness.jpg")),
                device = "jpg",
                height=2,
                width=3.5,
                dpi=600)
ggplot2::ggsave(filename = file.path(output_dir, paste0("isoeliteness.pdf")),
                device = "pdf",
                height=2,
                width=3.5)