# Title: LOD PEAKS
# Author: Ted Monyak
# Description: This script, which should only be run from AggregateFigures.R,
# creates agregated figures of the number of LOD peaks, from single and epistatic
# linkage mapping, detected for each combination of parameters

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

# The size of the geom_points
s <- 0.2

# Determine a string representation of the correlation and significance
sig_cor <- stat_cor(method="pearson",
                    aes(label=paste(
                      sub("R", "r", after_stat(r.label)),
                      ifelse(after_stat(p) < 0.001, '"***"',
                             ifelse(after_stat(p) < 0.01,  '"**"',
                                    ifelse(after_stat(p) < 0.05,  '"*"', '"ns"'))),
                      sep="~"
                    )),
                    parse=TRUE)

# Find the maximum number of LOD peaks to have a consistent y axis
qMax <- max(c(res.df$nLod_T1,
              res.df$nLod_T2,
              res.df$nLod_T3,
              res.df$nLod_W,
              res.df$nLod_Suit), na.rm=TRUE, inf.rm=TRUE)

# Make a boxplot for each kind of LOD peak
# y_var: The column in res.df starting with nLod, for each phenotype
make_boxplot <- function(y_var) {
  res.df %>%
    ggplot(aes(x = qtl, y = .data[[y_var]], fill = type)) +
    geom_boxplot(outlier.shape = NA, position = position_dodge(width = 1), linewidth = 0.1) +
    stat_compare_means(aes(group = type), method = "t.test", label = "p.signif",
                       hide.ns = FALSE, label.y = c(10.5, 10.5)) +
    scale_y_continuous(limits = c(0, qMax), expand = c(0, 0.2)) +
    scale_fill +
    labs(x = "QTL per Attained Trait", y = "Significant Peaks") +
    theme
}

# Make a multipanel of each kind of LOD peak
y_vars <- c("nLod_T1", "nLod_T2", "nLod_Suit", "nLod_T3", "nLod_W")
plots  <- lapply(y_vars, make_boxplot)
q <- wrap_plots(plots) +
  plot_layout(ncol = 1, guides = "collect") &
  theme(legend.position = "bottom")

# x-axis limits for isoeliteness
iMin <- 0.4
iMax <- 1

# Make a scatter plot of isoeliteness against 
# x_var: The specific isoeliteness metric
# y_var: The nLod peak to extract, one for each phenotype
# x_label: The name of the isoeliteness metric
make_isolod <- function(x_var, y_var, x_label) {
  res.df %>%
    dplyr::filter(type == "Admixed") %>%
    ggplot(aes(x = .data[[x_var]], y = .data[[y_var]])) +
    geom_point(aes(color = qtl), alpha = a, size = s) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.4) +
    sig_cor +
    guides(color = guide_legend(override.aes = list(size = 2))) +
    scale_x_continuous(limits = c(iMin, iMax), breaks = seq(0.4, 1, by = 0.2)) +
    scale_y_continuous(limits = c(0, qMax), expand = c(0, 0.2)) +
    scale_color +
    labs(x = x_label, y = "Significant Peaks") +
    theme
}

# Each of the plots to generate
plot_params <- list(
  list("isoElite_T1", "nLod_T1", "Attained Trait 1 Isoeliteness"),
  list("isoElite_T2", "nLod_T2", "Attained Trait 2 Isoeliteness"),
  list("isoElite_Att", "nLod_Suit", "Mean Isoeliteness"),
  list("isoElite_Att", "nLod_T3", "Mean Isoeliteness"),
  list("isoElite_Att", "nLod_W", "Mean Isoeliteness")
)

iso_plots <- lapply(plot_params, function(p) make_isolod(p[[1]], p[[2]], p[[3]]))

li <- wrap_plots(iso_plots) +
  plot_layout(ncol = 1, guides = "collect") &
  theme(legend.position = "bottom")

(q | li) + plot_annotation(tag_levels=list(c('f', 'g', 'h', 'i\ ', 'j', 'k', 'l\ ', 'm', 'n', 'o')),
                           theme = theme(plot.tag = element_text(family = "Helvetica", size = 6)))

ggplot2::ggsave(filename = paste0("lod_peaks.jpg"),
                path=output_dir,
                device = "jpg",
                width=3.5,
                height=8,
                dpi=600)

ggplot2::ggsave(filename = paste0("lod_peaks.pdf"),
                path=output_dir,
                device = "pdf",
                width=3.5,
                height=8)
theme <- theme + theme(
    plot.margin=unit(c(0,0,0,10), unit="pt"),
    legend.position = "right",
    legend.title.position = "top",
    legend.direction="vertical")

# Create emergent interaction LOD peaks plot
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

# Plot isoeliteness against the number of interaction LOD peaks
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
