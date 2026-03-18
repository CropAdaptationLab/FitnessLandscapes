# Title: ALLELE ORDERING
# Author: Ted Monyak
# Description: This script aggregates the results of several allele fixation runs

library("ggpubr")
library("grid")
library("psych")
setwd("~/Documents/CSU/R/BreedingSims")

#subpop50_qtl20 <- read.csv(file.path(getwd(), "Output/AverageEffectSize/subpop50_qtl20/effSize.csv"))
#subpop50_qtl2 <- read.csv(file.path(getwd(), "Output/AverageEffectSize/subpop50_qtl2/effSize.csv"))
#subpop500_qtl2 <- read.csv(file.path(getwd(), "Output/AverageEffectSize/subpop500_qtl2/effSize.csv"))
#dfs <- list(subpop50_qtl20, subpop50_qtl2, subpop500_qtl2)
#output_dir1 <- "~/Documents/CSU/R/BreedingSims/Output/AverageEffectSize/subpop50_qtl20"
#output_dir2 <- "~/Documents/CSU/R/BreedingSims/Output/AverageEffectSize/subpop50_qtl2"
#output_dir3 <- "~/Documents/CSU/R/BreedingSims/Output/AverageEffectSize/subpop500_qtl2"

#subpop50_qtl20_h30 <- read.csv(file.path(getwd(), "Output/AverageEffectSize/subpop50_qtl20_h2_30/effSize.csv"))
#subpop50_qtl2_h30 <- read.csv(file.path(getwd(), "Output/AverageEffectSize/subpop50_qtl2_h2_30/effSize.csv"))
#subpop500_qtl2_h30 <- read.csv(file.path(getwd(), "Output/AverageEffectSize/subpop500_qtl2_h2_30/effSize.csv"))
#dfs <- list(subpop50_qtl20_h30, subpop50_qtl2_h30, subpop500_qtl2_h30)
#output_dir1 <- file.path(getwd(), "Output/AverageEffectSize/subpop50_qtl20_h2_30")
#output_dir2 <- file.path(getwd(), "Output/AverageEffectSize/subpop50_qtl2_h2_30")
#output_dir3 <- file.path(getwd(), "Output/AverageEffectSize/subpop500_qtl2_h2_30")

output_dir1 <- file.path(getwd(), "Output/AverageEffectSize/2025-02-05_21_11_36")
subpop50_qtl2 <- read.csv(file.path(output_dir1, "effSize.csv"))
output_dir2 <- file.path(getwd(), "Output/AverageEffectSize/2025-02-05_23_22_28")
subpop500_qtl2 <- read.csv(file.path(output_dir2, "effSize.csv"))
dfs <- list(subpop50_qtl2, subpop500_qtl2)

output_dirs <- c(output_dir1, output_dir2)

theme <- theme(
  axis.text.x = element_text(angle = 0, hjust=1, size=16),
  axis.text.y = element_text(angle = 0, hjust=1, size=16),
  axis.title.x = element_text(family="Helvetica", size=30),
  axis.title.y = element_text(family="Helvetica", size=30),
  plot.title = element_text(family="Helvetica", size=32, hjust = 0.5),
  legend.text = element_text(family="Helvetica", size=10),
  legend.title = element_text(family="Helvetica", size=12),
  legend.key = element_rect(linewidth=0.05),
  plot.caption = element_text(family="Helvetica", size=10, hjust = 0.5),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_rect(fill = "white", color = "black"),
  aspect.ratio = 1)

removeOutliers <- "false"
for (i in 1:length(output_dirs)) {
  df <- dfs[[i]]
  for (N in c(100)) {
    # ORDER
    order.df <- df %>%
      group_by(orderFixed) %>%
      #filter(!(abs(effectSize - median(effectSize)) > 2*sd(effectSize))) %>%
      summarize(meanEffectSize = mean(effectSize), n=n()) %>%
      filter(n > N)
    
    gO <- ggplot(order.df, aes(x=orderFixed, y=meanEffectSize)) +
      geom_point() +
      theme +
      ylim(0.065, 0.1) +
      labs(title="Subpopulation Size = 50", x="Order Fixed", y="Mean Effect Size") +
      geom_smooth(formula = y ~ x, method="loess", color="black")
    
    ggplot2::ggsave(filename = paste0("effect_size_order.pdf"),
                    path=output_dirs[i],
                    device = "pdf",
                    width=8, height=8, units="in")

    # FITNESS
    fit.df <- df %>%
      mutate(smoothedFit=round(fitness,digits=1)/-4) %>%
      group_by(smoothedFit) %>%
      #filter(!(abs(effectSize - median(effectSize)) > 2*sd(effectSize))) %>%
      summarize(meanEffectSize = mean(effectSize), n=n()) %>%
      filter(n > N)
    
    gF <- ggplot(fit.df, aes(x=smoothedFit, y=meanEffectSize)) +
      geom_point() +
      scale_x_reverse() +
      theme +
      labs(title="", x="Normalized Distance from Fitness Optimum", y="Mean Effect Size") +
      geom_smooth(formula = y ~ x, method="loess", color="black")
    
    ggplot2::ggsave(filename = paste0("effect_size_fitness.jpg"),
                    path=output_dirs[i],
                    device = "jpg",
                    width=3, height=3, units="in")
    
    # GENERATION
    gen.df <- df %>%
      group_by(gen) %>%
      #filter(!(abs(effectSize - median(effectSize)) > 2*sd(effectSize))) %>%
      summarize(meanEffectSize = mean(effectSize), n=n()) %>%
      filter(n > N)
    
    gG <- ggplot(gen.df, aes(x=gen, y=meanEffectSize)) +
      geom_point() +
      theme +
      labs(title="", x="Generation", y="Mean Effect Size") +
      geom_smooth(formula = y ~ x, method="loess", color="black")
    
    ggplot2::ggsave(filename = paste0("effect_size_generation.jpg"),
                    path=output_dirs[i],
                    device = "jpg",
                    width=3, height=3, units="in")
  }
}