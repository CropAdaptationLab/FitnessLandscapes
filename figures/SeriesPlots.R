# Title: SERIES PLOTS
# Author: Ted Monyak
# Description: Plot the series of allelic substitution effect sizes for the
# fixation order of alleles, and the initial RIL architecture
# Run this from QtlMonteCarlo.R

library(ggpubr)
library(ggplot2)
library(ggpmisc)

equation <- stat_fit_tidy(method="nls",
                          method.args=list(formula=y ~ (a * rho^(x-1)),start=c(a=1,rho=0.9), algorithm="port"),
                          label.x="right",
                          label.y="top",
                          aes(label=sprintf("a[n]~`=`~%.2g~`*`~%.2g^{n-1}",
                                            after_stat(a_estimate),
                                            after_stat(rho_estimate))),
                          parse=TRUE,
                          size=7,
                          family="Helvetica")

# Save a plot of the average effect sizes, sorted by rank
if (saveEffectSizes) {
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

  effectSizes.df %>%
    group_by(rank) %>%
    summarize(meanEffectSize = mean(eff_size)) %>%
    ggplot(aes(x=rank, y=meanEffectSize)) +
    geom_point() +
    theme +
    labs(title="",
         x="Rank", y="Mean Effect Size") +
    ril_equation
  
  ggplot2::ggsave(filename = file.path(output_dir, "average_effect_size_RIL.jpg"),
                  device = "jpg",
                  height=6,
                  width=6,
                  units="in",
                  dpi=600)
  
  ggplot2::ggsave(filename = file.path(output_dir, "average_effect_size_RIL.pdf"),
                  device = "pdf",
                  height=6,
                  width=6,
                  units="in")
}

# Save a plot of fixation order of alleles, sorted by rank
if (saveFixationOrder) {
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

  
  # Determine the average additive effect size at each 'step'
  fixedAlleles.df %>%
    group_by(order_fixed) %>%
    summarize(meanEffectSize = mean(eff_size), n=n()) %>%
    ggplot(aes(x=order_fixed, y=meanEffectSize)) +
    geom_point() +
    theme +
    labs(x="Order of Fixation",
         y="Allele Substitution Effect") +
    fixation_equation
  
  ggplot2::ggsave(filename = file.path(output_dir, "average_effect_size_fixed.jpg"),
                  device = "jpg",
                  height=6,
                  width=6,
                  units="in",
                  dpi=600)
  ggplot2::ggsave(filename = file.path(output_dir, "average_effect_size_fixed.pdf"),
                  device = "pdf",
                  height=6,
                  width=6,
                  units="in")
}