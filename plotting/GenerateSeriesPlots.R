# Title: Generate Plots
# Author Ted Monyak
# Description:
# This script generates plots for the dataframes generated in QtlMonteCarlo.R

theme <- theme(
  axis.title.x = element_text(family="Helvetica", size=18),
  axis.text.x = element_text(angle = 0, hjust=1, size=16),
  axis.title.y = element_text(family="Helvetica", size=18),
  axis.text.y = element_text(angle = 0, hjust=1, size=16),
  plot.title = element_text(family="Helvetica", size=18, hjust = 0.5),
  legend.text = element_text(family="Helvetica", size=14),
  legend.title = element_text(family="Helvetica", size=16),
  legend.key = element_rect(linewidth=0.05),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_rect(fill = "white", color = "black"),
  plot.margin= unit(c(10,10,10,10), unit="pt"),
  aspect.ratio = 1)

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

# Save a plot of the average effect sizes, sorted by effect size
if (saveEffectSizes) {
  effectSizes.df %>%
    group_by(rank) %>%
    dplyr::summarize(meanEffectSize = mean(eff_size)) %>%
    ggplot(aes(x=rank, y=meanEffectSize)) +
    geom_point(size=3) +
    theme +
    labs(title="",
         x="Rank", y="Mean Effect Size") +
    equation
  
  fname <- file.path(base_dir, paste0(base_fname, "average_effect_size_RIL.jpg"))
  ggplot2::ggsave(filename = fname,
                  device = "jpg",
                  height=6,
                  width=6,
                  units="in")
  write.table(effectSizes.df, file.path(base_dir, "effect_size.csv"), col.names=TRUE, quote=FALSE, sep=",")
}

if (saveFixationOrder) {
  # Determine the average additive effect size at each 'step'
  fixedAlleles.df %>%
    group_by(order_fixed) %>%
    dplyr::summarize(meanEffectSize = mean(eff_size), n=n()) %>%
    filter(n>100) %>%
    ggplot(aes(x=order_fixed, y=meanEffectSize)) +
      geom_point(size=3) +
      theme +
      labs(x="Order of Fixation",
           y="Mean Additive Effect") +
    equation
  
  fname <- file.path(base_dir, paste0(base_fname, "average_effect_size_fixed.jpg"))
  ggplot2::ggsave(filename = fname,
                  device = "jpg",
                  height=6,
                  width=6,
                  units="in")
  write.table(fixedAlleles.df, file.path(base_dir, "fixed_alleles.csv"), col.names=TRUE, quote=FALSE, sep=",")
}