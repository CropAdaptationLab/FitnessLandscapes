library(dplyr)
library(ggplot2)
library(ggpubr)
library(grid)
library(patchwork)

setwd("~/Documents/CSU/Lab-Notebooks/members/Ted/breeding_sims/output/QtlMonteCarlo/effect_sizes")
output_dir <- file.path(getwd(), "../../Geometric/3_15")
if (!dir.exists(output_dir)) dir.create(output_dir)

# QTL: 1 V: 0.4
# 50 pop resets x 2 sims (100)
qtl1_v4_100 <- read.csv("qtl_1_Ne_1000_var_0.4_2026-02-26_21_53/fixed_alleles.csv")
qtl1_v4_100$qtl <- 1
qtl1_v4_100$var <- 0.4

# 200 pop resets x 2 sims (400)
qtl1_v4_400 <- read.csv("qtl_1_Ne_1000_var_0.4_2026-02-27_22_27/fixed_alleles.csv")
qtl1_v4_400$qtl <- 1
qtl1_v4_400$var <- 0.4


# QTL: 2 V: 0.4
# 50 pop resets x 2 sims (100)
qtl2_v4_100 <- read.csv("qtl_2_Ne_1000_var_0.4_2026-02-27_00_24/fixed_alleles.csv")
qtl2_v4_100$qtl <- 2
qtl2_v4_100$var <- 0.4

# 100 pop resets x 2 sims (200)
qtl2_v4_200_a <- read.csv("qtl_2_Ne_1000_var_0.4_2026-02-28_08_42/fixed_alleles.csv")
qtl2_v4_200_a$qtl <- 2
qtl2_v4_200_a$var <- 0.4

# 100 pop resets x 2 sims (200)
qtl2_v4_200_b <- read.csv("qtl_2_Ne_1000_var_0.4_2026-03-01_21_25/fixed_alleles.csv")
qtl2_v4_200_b$qtl <- 2
qtl2_v4_200_b$var <- 0.4


# QTL: 5 V: 0.4
# 50 pop resets x 2 sims (100)
qtl5_v4_100 <- read.csv("qtl_5_Ne_1000_var_0.4_2026-02-27_02_55/fixed_alleles.csv")
qtl5_v4_100$qtl <- 5
qtl5_v4_100$var <- 0.4

# 100 pop resets x 2 sims (200)
qtl5_v4_200_a <- read.csv("qtl_5_Ne_1000_var_0.4_2026-02-28_13_52/fixed_alleles.csv")
qtl5_v4_200_a$qtl <- 5
qtl5_v4_200_a$var <- 0.4

# 100 pop resets x 2 sims (200)
qtl5_v4_200_b <- read.csv("qtl_5_Ne_1000_var_0.4_2026-03-02_02_39/fixed_alleles.csv")
qtl5_v4_200_b$qtl <- 5
qtl5_v4_200_b$var <- 0.4


fixed_df <- rbind(qtl1_v4_100,
                  qtl1_v4_400,
                  qtl2_v4_100,
                  qtl2_v4_200_a,
                  qtl2_v4_200_b,
                  qtl5_v4_100,
                  qtl5_v4_200_a,
                  qtl5_v4_200_b)

#source("Plotting/AggregateTraitArchitecture.R")

fixed_df <- fixed_df %>%
  dplyr::mutate(type="fixation") %>%
  dplyr::rename("rank"="order_fixed") %>%
  dplyr::select(rank, eff_size, qtl, type)

init_df <- read.csv(file.path(output_dir, "initial_architecture.csv"))
  
init_df <- init_df %>% 
  dplyr::mutate(type="initial") %>%
  dplyr::select(rank, eff_size, qtl, type)

merged <- rbind(fixed_df, init_df)
merged$type <- factor(merged$type, levels=c("initial", "fixation"))

df_summary <- merged %>%
  dplyr::group_by(qtl, type, rank) %>%
  dplyr::summarize(mean_eff_size = exp(mean(log(eff_size))), n = n(), .groups = "drop") %>% #mean(eff_size)
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

qtl1_v4_plot <- plot_series(1)
qtl2_v4_plot <- plot_series(2)
qtl5_v4_plot <- plot_series(5)

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

v2_label <- wrap_elements(panel = textGrob('Initial Additive\nVariance: 0.2',
                                           rot=90,
                                           gp=labelfont),
                          ignore_tag = TRUE)

v4_label <- wrap_elements(panel = textGrob('Initial Additive\nVariance: 0.4',
                                           rot=90,
                                           gp=labelfont),
                          ignore_tag = TRUE)

p <- qtl1_label + qtl2_label + qtl5_label +
  qtl1_v4_plot + qtl2_v4_plot + qtl5_v4_plot +
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

write.table(merged, file.path(output_dir, "fixationorder.csv"), col.names=TRUE, quote=FALSE, sep=",")
