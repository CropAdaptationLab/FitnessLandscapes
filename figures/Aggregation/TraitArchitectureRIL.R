# Title: TRAIT ARCHITECTURE RIL
# Author: Ted Monyak
# Description: This creates a multipanel where the trait architecture of admixed
# RIL families (by number of QTL per attained trait) is overlaied on the founder trait architecture
# Assumes that SimulationPipeline and TraitArchitectureFounder have already been run
library(dplyr)
library(ggplot2)
library(ggpubr)
library(grid)
library(patchwork)

# Designate the type as 'ril'
effectSizes.df <- effectSizes.df %>%
  dplyr::mutate(type="ril") %>%
  dplyr::select(rank, eff_size, qtl, type)

init_df <- read.csv("~/Documents/CSU/FitnessLandscapes/output/TraitArchitecture/initial_architecture.csv")

# Designate the type as 'initial'
init_df <- init_df %>% 
  dplyr::mutate(type="initial") %>%
  dplyr::select(rank, eff_size, qtl, type)

# Join the RIL and Iniital dataframes
merged <- rbind(effectSizes.df, init_df)
merged$type <- factor(merged$type, levels=c("initial", "ril"))

# Get the mean effect size at each rank
# There are only 20 samples per initial architecture, so filter above 19
df_summary <- merged %>%
  dplyr::group_by(qtl, type, rank) %>%
  dplyr::summarize(mean_eff_size = mean(eff_size), n = n(), .groups = "drop") %>%
  dplyr::filter(n>1) # Increase filter with more data?

df_summary$rank <- as.numeric(df_summary$rank)

geometric_theme <- theme_minimal(base_size=7,
                                 base_family="Helvetica") +
  theme(
    plot.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = "black"),
    legend.text = element_text(size=8),
    plot.margin= unit(c(0,0,0,0), unit="pt"),
    legend.position = "bottom",
    legend.direction="horizontal")

# This is the line of best fit for the initial trait architecture
init_equation <- stat_fit_tidy(method="nls",
                               method.args=list(formula=y ~ (a * rho^(x-1)),
                                                start=c(a=1,rho=0.9),
                                                algorithm="port"), # Use port algorithm so it will converge
                               data = ~filter(.x, type == "initial"),
                               label.x=0.95,
                               label.y=0.95,
                               aes(label=sprintf("alpha[k]~`=`~%.2g~`*`~%.2g^{k-1}",
                                                 after_stat(a_estimate),
                                                 after_stat(rho_estimate))),
                               parse=TRUE,
                               size=2.5,
                               family="Helvetica")

# This is the line of best fit for the RIL trait architecture
# A 'b' intercept is added to account for the fact that the initial trait architecture
# does not converge to zero
ril_equation <- stat_fit_tidy(method="nls",
                                   method.args=list(formula=y ~ (a * rho^(x-1)) + b,
                                                    start=c(a=1,rho=0.9, b=0),
                                                    algorithm="port"), # Use port algorithm so it will converge
                                   data = ~filter(.x, type == "ril"),
                                   label.x=0.95,
                                   label.y=0.8,
                                   aes(label=sprintf("alpha[k]~`=`~%.2g~`*`~%.2g^{k-1} ~`+`~%.2g",
                                                     after_stat(a_estimate),
                                                     after_stat(rho_estimate),
                                                     after_stat(b_estimate))),
                                   parse=TRUE,
                                   size=2.5,
                                   family="Helvetica")

# Get the axis ranges
max_x <- max(df_summary$rank)
max_y <- max(df_summary$mean_eff_size)

# Generates an overlaid plot of the two geometric series
# Filters by nQtl
plot_series <- function(nQtl) {
  
  # Fit the NLS model to series data
  fit <- nls(mean_eff_size ~ a * rho^(rank-1) + b,
             data = filter(df_summary, type == "ril", qtl==nQtl),
             start = c(a=1, rho=0.9, b=0),
             algorithm = "port")
  
  a_hat   <- coef(fit)["a"]
  rho_hat <- coef(fit)["rho"]
  b_hat <- coef(fit)["b"]
  
  # Calculate the split line between stochastic and deterministic
  ril_max <- df_summary %>%
    filter(type == "ril", qtl==nQtl) %>%
    pull(mean_eff_size) %>%
    max()

  df_summary %>%
    dplyr::filter(qtl==nQtl) %>%
    ggplot(aes(x = rank, y = mean_eff_size, color = type)) +
    annotate("rect", xmin=-Inf, xmax=Inf, ymin=ril_max, ymax=Inf,
             fill="lightblue", alpha=0.2) +
    annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=ril_max,
             fill="lightyellow", alpha=0.2) +
    geom_hline(yintercept=ril_max, linetype="dashed", linewidth=0.3, color="black") +
    annotate("text", x=Inf, y=ril_max - (max_y * 0.03), 
             label="stochastic", hjust=1.1, vjust=1, 
             size=2, color="black", family="Helvetica") +
    annotate("text", x=Inf, y=ril_max + (max_y * 0.03), 
             label="deterministic", hjust=1.1, vjust=0, 
             size=2, color="black", family="Helvetica") +
    stat_function(fun = function(x) (a_hat * rho_hat^(x-1)) + b_hat,
                  color = "black", linewidth = 0.3, linetype = "dotted") +
    geom_point() +
    init_equation +
    ril_equation +
    scale_color_manual(values = c("initial" = "#808080", "ril" = "black"),
                       labels = c("initial" = "Initial Architecture", "ril" = "RIL Architecture"),
                       name=NULL) +
    guides(color = guide_legend(override.aes = list(shape = 16, linetype = 0, size=2))) +
    labs(x="Rank", y="Mean Allelic\nSubstitution Effect") +
    xlim(0,max_x) +
    scale_y_continuous(limits = c(0, max_y), expand = c(0, 0.05)) +
    geometric_theme
}

qtl1_plot <- plot_series(1)
qtl2_plot <- plot_series(2)
qtl5_plot <- plot_series(5)

labelfont <- gpar(fontsize=8,
                  fontfamily="Helvetica")

qtl1_label <- wrap_elements(panel = textGrob('10 QTL per Attained Trait',
                                             rot=0,
                                             gp=labelfont),
                            ignore_tag = TRUE)
qtl2_label <- wrap_elements(panel = textGrob('20 QTL per Attained Trait',
                                             rot=0,
                                             gp=labelfont),
                            ignore_tag = TRUE)
qtl5_label <- wrap_elements(panel = textGrob('50 QTL per Attained Trait',
                                             rot=0,
                                             gp=labelfont),
                            ignore_tag = TRUE)

p <- qtl1_label + qtl2_label + qtl5_label +
  qtl1_plot + qtl2_plot + qtl5_plot +
  plot_layout(guides='collect',
              nrow=2,
              ncol=3,
              heights=c(0.5,3),
              widths=c(3,3,3)) +
  plot_annotation(tag_levels=list(c('a', 'b', 'c'))) & theme(plot.tag = element_text(size = 10), legend.position="bottom")

p

ggplot2::ggsave(filename = file.path(output_dir, "emergent_architecture.jpg"),
                device = "jpg",
                height=2.5,
                width=6.5,
                dpi=600)
ggplot2::ggsave(filename = file.path(output_dir, "emergent_architecture.pdf"),
                device = "pdf",
                height=2.5,
                width=6.5)