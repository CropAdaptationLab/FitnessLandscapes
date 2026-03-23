# Title: FIXATION ORDER
# Author: Ted Monyak
# Description: This script overlays the mean fixation order of alleles along
# an adaptive walk over the founder trait architecture
# Assumes that SimulationPipeline and TraitArchitectureFounder have already been run
library(dplyr)
library(ggplot2)
library(ggpubr)
library(grid)
library(patchwork)

# Designate the type as 'fixation'
fixedAlleles.df <- fixedAlleles.df %>%
  dplyr::mutate(type="fixation") %>%
  dplyr::select(rank, eff_size, qtl, type)

# Read in the founder trait architecture
init_df <- read.csv("~/Documents/CSU/FitnessLandscapes/output/TraitArchitecture/initial_architecture.csv")

# Designate the type as 'initial'
init_df <- init_df %>% 
  dplyr::mutate(type="initial") %>%
  dplyr::select(rank, eff_size, qtl, type)

merged <- rbind(fixedAlleles.df, init_df)
merged$type <- factor(merged$type, levels=c("initial", "fixation"))

# Get the mean effect size at each rank
df_summary <- merged %>%
  dplyr::group_by(qtl, type, rank) %>%
  dplyr::summarize(mean_eff_size = mean(eff_size), n = n(), .groups = "drop") %>%
  dplyr::filter((type=="fixation" & n > 100) | type=="initial")

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

# Geometric series line of best fit for the initial trait architecture
init_equation <- stat_fit_tidy(method="nls",
                          method.args=list(formula=y ~ (a * rho^(x-1)),
                                           start=c(a=1,rho=0.9),
                                           algorithm="port"),
                          data = ~filter(.x, type == "initial"),
                          label.x=0.95,
                          label.y=0.95,
                          aes(label=sprintf("alpha[k]~`=`~%.2g~`*`~%.2g^{k-1}",
                                            after_stat(a_estimate),
                                            after_stat(rho_estimate))),
                          parse=TRUE,
                          size=2.5,
                          family="Helvetica")

# Geometric series line of best fit for the fixation order of alleles
fixation_equation <- stat_fit_tidy(method="nls",
                               method.args=list(formula=y ~ (a * rho^(x-1)),
                                                start=c(a=1,rho=0.9),
                                                algorithm="port"),
                               data = ~filter(.x, type == "fixation"),
                               label.x=0.95,
                               label.y=0.8,
                               aes(label=sprintf("alpha[k]~`=`~%.2g~`*`~%.2g^{k-1}",
                                                 after_stat(a_estimate),
                                                 after_stat(rho_estimate))),
                               parse=TRUE,
                               size=2.5,
                               family="Helvetica")
max_x <- max(df_summary$rank)
max_y <- max(df_summary$mean_eff_size)

# Generates an overlaid plot of the two geometric series
# Filters by nQtl
plot_series <- function(nQtl) {
  # Fit the NLS model to fixation data
  fit <- nls(mean_eff_size ~ a * rho^(rank-1),
             data = filter(df_summary, type == "fixation", qtl==nQtl),
             start = c(a=1, rho=0.9),
             algorithm = "port")
  
  a_hat   <- coef(fit)["a"]
  rho_hat <- coef(fit)["rho"]

  df_summary %>%
    dplyr::filter(qtl==nQtl) %>%
    ggplot(aes(x = rank, y = mean_eff_size, color = type)) +
    stat_function(fun = function(x) a_hat * rho_hat^(x-1),
                  color = "black", linewidth = 0.3, linetype = "dotted") +
    geom_point() +
    init_equation +
    fixation_equation +
    scale_color_manual(values = c("initial" = "#808080", "fixation" = "black"),
                       labels = c("initial" = "Initial Architecture", "fixation" = "Fixation Order"),
                       name=NULL) +
    guides(color = guide_legend(override.aes = list(shape = 16, linetype = 0, size=2))) +
    labs(x="Rank", y="Mean Allele\nSubstitution Effect") +
    xlim(0,max_x) +
    ylim(0,max_y) +
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
  plot_annotation(tag_levels=list(c('c', 'd', 'e'))) & theme(plot.tag = element_text(size = 14), legend.position="bottom")

p

ggplot2::ggsave(filename = file.path(output_dir, "fixation_order.jpg"),
                device = "jpg",
                height=2.5,
                width=6.5,
                dpi=600)
ggplot2::ggsave(filename = file.path(output_dir, "fixation_order.pdf"),
                device = "pdf",
                height=2.5,
                width=6.5)