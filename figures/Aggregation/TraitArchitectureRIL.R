library(dplyr)
library(ggplot2)
library(ggpubr)
library(grid)
library(patchwork)

setwd("~/Documents/CSU/Lab-Notebooks/members/Ted/breeding_sims/output/QtlMonteCarlo/effect_sizes")
output_dir <- file.path(getwd(), "../../Geometric/3_15")
if (!dir.exists(output_dir)) dir.create(output_dir)

# Still needed to get to 500:
# qtl1_v2: 460
# qtl2_v2: 460
# qtl5_v2: 460

# QTL: 1 V: 0.2
# 20 pop resets x 2 sims (40)
qtl1_v2_40 <- read.csv("qtl_1_Ne_1000_var_0.2_2026-02-27_06_23/effect_size.csv")
qtl1_v2_40$qtl <- 1
qtl1_v2_40$var <- 0.2


# QTL: 2 V: 0.2
# 20 pop resets x 2 sims (40)
qtl2_v2_40 <- read.csv("qtl_2_Ne_1000_var_0.2_2026-02-27_07_21/effect_size.csv")
qtl2_v2_40$qtl <- 2
qtl2_v2_40$var <- 0.2


# QTL: 5 V: 0.2
# 20 pop resets x 2 sims (40)
qtl5_v2_40 <- read.csv("qtl_5_Ne_1000_var_0.2_2026-02-27_08_21/effect_size.csv")
qtl5_v2_40$qtl <- 5
qtl5_v2_40$var <- 0.2


# QTL: 1 V: 0.4
# 50 pop resets x 2 sims (100)
qtl1_v4_100 <- read.csv("qtl_1_Ne_1000_var_0.4_2026-02-26_21_53/effect_size.csv")
qtl1_v4_100$qtl <- 1
qtl1_v4_100$var <- 0.4

# 200 pop resets x 2 sims (400)
qtl1_v4_400 <- read.csv("qtl_1_Ne_1000_var_0.4_2026-02-27_22_27/effect_size.csv")
qtl1_v4_400$qtl <- 1
qtl1_v4_400$var <- 0.4


# QTL: 2 V: 0.4
# 50 pop resets x 2 sims (100)
qtl2_v4_100 <- read.csv("qtl_2_Ne_1000_var_0.4_2026-02-27_00_24/effect_size.csv")
qtl2_v4_100$qtl <- 2
qtl2_v4_100$var <- 0.4

# 100 pop resets x 2 sims (200)
qtl2_v4_200_a <- read.csv("qtl_2_Ne_1000_var_0.4_2026-02-28_08_42/effect_size.csv")
qtl2_v4_200_a$qtl <- 2
qtl2_v4_200_a$var <- 0.4

# 100 pop resets x 2 sims (200)
qtl2_v4_200_b <- read.csv("qtl_2_Ne_1000_var_0.4_2026-03-01_21_25/effect_size.csv")
qtl2_v4_200_b$qtl <- 2
qtl2_v4_200_b$var <- 0.4


# QTL: 5 V: 0.4
# 50 pop resets x 2 sims (100)
qtl5_v4_100 <- read.csv("qtl_5_Ne_1000_var_0.4_2026-02-27_02_55/effect_size.csv")
qtl5_v4_100$qtl <- 5
qtl5_v4_100$var <- 0.4

# 100 pop resets x 2 sims (200)
qtl5_v4_200_a <- read.csv("qtl_5_Ne_1000_var_0.4_2026-02-28_13_52/effect_size.csv")
qtl5_v4_200_a$qtl <- 5
qtl5_v4_200_a$var <- 0.4

# 100 pop resets x 2 sims (200)
qtl5_v4_200_b <- read.csv("qtl_5_Ne_1000_var_0.4_2026-03-02_02_39/effect_size.csv")
qtl5_v4_200_b$qtl <- 5
qtl5_v4_200_b$var <- 0.4

ril_df <- rbind(qtl1_v4_100,
                  qtl1_v4_400,
                  qtl2_v4_100,
                  qtl2_v4_200_a,
                  qtl2_v4_200_b,
                  qtl5_v4_100,
                  qtl5_v4_200_a,
                  qtl5_v4_200_b,
                  qtl1_v2_40,
                  qtl2_v2_40,
                  qtl5_v2_40)

#source("Plotting/AggregateTraitArchitecture.R")

ril_df <- ril_df %>%
  dplyr::mutate(type="ril") %>%
  dplyr::select(rank, eff_size, qtl, var, type)

init_df <- read.csv(file.path(output_dir, "initial_architecture.csv"))

init_df <- init_df %>% 
  dplyr::mutate(type="initial") %>%
  dplyr::select(rank, eff_size, qtl, var, type)

merged <- rbind(ril_df, init_df)
merged$type <- factor(merged$type, levels=c("initial", "ril"))

# There are only 20 samples per init architecture, so filter above 19
df_summary <- merged %>%
  group_by(qtl, var, type, rank) %>%
  summarise(mean_eff_size = mean(eff_size), n = n(), .groups = "drop") %>%
  filter(n>1) # Increase filter with more data?

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
                                   method.args=list(formula=y ~ (a * rho^(x-1)) + b,
                                                    start=c(a=1,rho=0.9, b=0),
                                                    algorithm="port"),
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

max_x <- max(df_summary$rank)
max_y <- max(df_summary$mean_eff_size)
plot_series <- function(nQtl, nVar) {
  
  # Fit the NLS model to fixation data
  fit <- nls(mean_eff_size ~ a * rho^(rank-1) + b,
             data = filter(df_summary, type == "ril", qtl==nQtl, var==nVar),
             start = c(a=1, rho=0.9, b=0),
             algorithm = "port")
  
  a_hat   <- coef(fit)["a"]
  rho_hat <- coef(fit)["rho"]
  b_hat <- coef(fit)["b"]
  
  # Calculate the split line
  ril_max <- df_summary %>%
    filter(type == "ril", qtl==nQtl, var==nVar) %>%
    pull(mean_eff_size) %>%
    max()

  df_summary %>%
    dplyr::filter(qtl==nQtl, var==nVar) %>%
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
    fixation_equation +
    scale_color_manual(values = c("initial" = "#808080", "ril" = "black"),
                       labels = c("initial" = "Initial Architecture", "ril" = "RIL Architecture"),
                       name=NULL) +
    guides(color = guide_legend(override.aes = list(shape = 16, linetype = 0, size=2))) +
    labs(x="Rank", y="Mean Allelic\nSubstitution Effect") +
    xlim(0,max_x) +
    scale_y_continuous(limits = c(0, max_y), expand = c(0, 0.05)) +
    geometric_theme
}
#qtl1_v2_plot <- plot_series(1, 0.2)
#qtl2_v2_plot <- plot_series(2, 0.2)
#qtl5_v2_plot <- plot_series(5, 0.2)
qtl1_v4_plot <- plot_series(1, 0.4)
qtl2_v4_plot <- plot_series(2, 0.4)
qtl5_v4_plot <- plot_series(5, 0.4)

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

write.table(merged, file.path(output_dir, "emergentarchitecture.csv"), col.names=TRUE, quote=FALSE, sep=",")
