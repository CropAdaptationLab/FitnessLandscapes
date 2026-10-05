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

fixedAlleles.df <- read.csv("fixationorder.csv")

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

merged <- merged %>%
  dplyr::mutate(qtl = qtl*10) %>%
  dplyr::mutate(qtl = as.factor(qtl))

# Get the mean effect size at each rank
df_summary <- merged %>%
  dplyr::filter(type == "fixation") %>%
  dplyr::group_by(qtl, type, rank) %>%
  dplyr::summarize(mean_eff_size = mean(eff_size), n = n(), .groups = "drop") %>%
  dplyr::filter(n > 500)

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

# Geometric series line of best fit for the fixation order of alleles
fixation_equation <- stat_fit_tidy(method="nls",
                               method.args=list(formula=y ~ (a * rho^(x-1)) + b,
                                                start=c(a=1,rho=0.9, b=0),
                                                algorithm="port"),
                               data = ~filter(.x, type == "fixation"),
                               label.x=0.95,
                               aes(label=sprintf("alpha[k]~`=`~%.2g~`*`~%.2g^{k-1} ~`+`~%.2g",
                                                 after_stat(a_estimate),
                                                 after_stat(rho_estimate),
                                                 after_stat(b_estimate))),
                               parse=TRUE,
                               size=2.5,
                               family="Helvetica")
max_x <- max(df_summary$rank)
max_y <- max(df_summary$mean_eff_size)

# Generates an overlaid plot of the two geometric series
# Fit the NLS model to fixation data
fit <- nls(mean_eff_size ~ a * rho^(rank-1) + b,
           data = filter(df_summary, type == "fixation", qtl==10),
           start = c(a=1, rho=0.9, b=0),
           algorithm = "port")

a_hat_10   <- coef(fit)["a"]
rho_hat_10 <- coef(fit)["rho"]
b_hat_10 <- coef(fit)["b"]

fit <- nls(mean_eff_size ~ a * rho^(rank-1) + b,
           data = filter(df_summary, type == "fixation", qtl==20),
           start = c(a=1, rho=0.9, b=0),
           algorithm = "port")

a_hat_20   <- coef(fit)["a"]
rho_hat_20 <- coef(fit)["rho"]
b_hat_20 <- coef(fit)["b"]

fit <- nls(mean_eff_size ~ a * rho^(rank-1) + b,
           data = filter(df_summary, type == "fixation", qtl==50),
           start = c(a=1, rho=0.9, b=0),
           algorithm = "port")

a_hat_50   <- coef(fit)["a"]
rho_hat_50 <- coef(fit)["rho"]
b_hat_50 <- coef(fit)["b"]

df_summary %>%
  dplyr::filter(type=="fixation") %>%
  ggplot(aes(x = rank, y = mean_eff_size, color = qtl)) +
  stat_function(fun = function(x) (a_hat_10 * rho_hat_10^(x-1)) + b_hat_10,
                color = "black", linewidth = 0.3, linetype = "dotted") +
  stat_function(fun = function(x) (a_hat_20 * rho_hat_20^(x-1)) + b_hat_20,
                color = "grey40", linewidth = 0.3, linetype = "dotted") +
  stat_function(fun = function(x) (a_hat_50 * rho_hat_50^(x-1)) + b_hat_50,
                color = "grey70", linewidth = 0.3, linetype = "dotted") +
  geom_point() +
  fixation_equation +
  scale_color +
  guides(color = guide_legend(override.aes = list(shape = 16, linetype = 0, size=2))) +
  labs(x="Order of Fixation", y="Mean Allele\nSubstitution Effect") +
  xlim(0,max_x) +
  ylim(0.1,max_y) +
  geometric_theme

ggplot2::ggsave(filename = file.path(output_dir, "fixation_order.jpg"),
                device = "jpg",
                height=2.5,
                width=3,
                dpi=600)
ggplot2::ggsave(filename = file.path(output_dir, "fixation_order.pdf"),
                device = "pdf",
                height=2.5,
                width=3)