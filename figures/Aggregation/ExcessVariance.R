# Title: EXCESS VARIANCE
# Author: Ted Monyak
# Description: This function aggregates all of the transgressive segregation
# data from QtlMonteCarlo, quantified with 'excess variance' (EV)
# Run this script only from AggregateFigures.R

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

# Function for plotting the excess variance by number of QTL
# y_var: The column in res.df for the 'ev' of each phenotype
# title: The phenotype
make_ev_box <- function(y_var, title) {
  res.df %>%
    dplyr::filter(!is.na(.data[[y_var]]) & !is.infinite(.data[[y_var]])) %>%
    ggplot(aes(x = qtl, y = .data[[y_var]], fill = type)) +
    geom_boxplot(outlier.shape = NA, position = position_dodge(width = 1), linewidth = 0.1) +
    scale_y_continuous(limits = c(vMin, vMax), expand = c(0, 0.1)) +
    scale_fill +
    labs(x = "QTL per\nAttained Trait", y = "EV", title = title) +
    theme
}

# x-axis limits for isoeliteness
iMin <- 0.4
iMax <- 1

# Plot isoeliteness against excess variance
# x_var: The isoeliteness metric
# y_var: The 'ev' metric
# x_label: The name of the isoeliteness metric
# title: The title for the plot
make_ev_scatter <- function(x_var, y_var, x_label, title) {
  res.df %>%
    dplyr::filter(type == "Admixed") %>%
    dplyr::filter(!is.na(.data[[y_var]]) & !is.infinite(.data[[y_var]])) %>%
    ggplot(aes(x = .data[[x_var]], y = .data[[y_var]], color = qtl)) +
    geom_point(alpha = a, size = s) +
    guides(color = guide_legend(override.aes = list(size = 2))) +
    scale_x_continuous(limits = c(iMin, iMax), breaks = seq(0.4, 1, by = 0.2)) +
    scale_y_continuous(limits = c(vMin, vMax * 1.05), expand = c(0, 0.1)) +
    scale_color +
    labs(x = x_label, y = "EV", title = title) +
    theme
}

plot_params <- list(
  list(y = "ev_T1", x = "isoElite_T1", x_lab = "Attained Trait 1\nIsoeliteness", title = "Attained Trait 1"),
  list(y = "ev_T2", x = "isoElite_T2", x_lab = "Attained Trait 2\nIsoeliteness", title = "Attained Trait 2"),
  list(y = "ev_Suit", x = "isoElite_Att", x_lab = "Mean\nIsoeliteness", title = "Suitability"),
  list(y = "ev_T3", x = "isoElite_Att", x_lab = "Mean\nIsoeliteness", title = "Desired Trait"),
  list(y = "ev_W", x = "isoElite_Att", x_lab = "Mean\nIsoeliteness", title = "Breeding Fitness")
)

box_plots <- lapply(plot_params, function(p) make_ev_box(p$y, p$title))
scatter_plots <- lapply(plot_params, function(p) make_ev_scatter(p$x, p$y, p$x_lab, p$title))

v  <- wrap_plots(box_plots) + plot_layout(ncol = 5, guides = "collect", axes = "collect")
vi <- wrap_plots(scatter_plots) + plot_layout(ncol = 5, guides = "collect", axes = "collect")

(v / vi) + plot_annotation(tag_levels=list(c('f', 'g', 'h', 'i\ ', 'j', 'k', 'l\ ', 'm', 'n', 'o')),
                           theme = theme(plot.tag = element_text(family = "Helvetica", size = 6)))
ggplot2::ggsave(filename = paste0("excess_variance.jpg"),
                path=output_dir,
                device = "jpg",
                width=6.5,
                height=3,
                dpi=600)

ggplot2::ggsave(filename = paste0("excess_variance.pdf"),
                path=output_dir,
                device = "pdf",
                width=6.5,
                height=3)
