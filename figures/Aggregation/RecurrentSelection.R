# Title: RECURRENT SELECTION
# Author: Ted Monyak
# Description: Aggregate results from the recurrent selection simulations

library(ggplot2)

gainOverCycles <- function(df, nQtl) {
  df %>%
    dplyr::filter(qtl==nQtl) %>%
    dplyr::mutate(group=paste0("RIL=", ril, ", sel=", sel)) %>%
    ggplot(aes(x = c, y = meanW, color = group, group = group)) +
      geom_line() +
      geom_point() +
      scale_color_manual(
        values = c(
          "RIL=Admixed, sel=GARS"     = "#E41A1C",
          "RIL=Admixed, sel=PRS"      = "#377EB8",
          "RIL=Unadmixed, sel=GARS"   = "#4DAF4A",
          "RIL=Unadmixed, sel=PRS"    = "#FF7F00"
        )
      ) +
      labs(
        title  = paste("QTL =", nQtl),
        x      = "c",
        y      = "meanW",
        color  = "Group"
      ) + 
      theme_minimal()
}

gainOverCycles.df <-  gs.df %>%
  dplyr::group_by(qtl, ril, c, sel) %>%
  dplyr::summarize(meanW = mean(w))

gainOverCycles(gainOverCycles.df, 10)
gainOverCycles(gainOverCycles.df, 20)
gainOverCycles(gainOverCycles.df, 50)

plotGeneticGain <- function(df, nQtl) {
  df %>%
    dplyr::filter(qtl==nQtl) %>%
    dplyr::mutate(group=paste0("RIL=", ril, ", sel=", sel)) %>%
    ggplot(aes(x = group, y = geneticGain)) +
    geom_boxplot() +
    labs(
      title  = paste("QTL =", nQtl),
      x      = "Group",
      y      = "Genetic Gain"
    ) + 
    theme_minimal()
}

geneticGain.df <- gs.df %>%
  dplyr::group_by(qtl, pop, rep, ril, sel) %>%
  dplyr::summarize(
    geneticGain = w[c == 10] - w[c == 1],
    .groups = "drop"
  )


plotGeneticGain(geneticGain.df, 10)
plotGeneticGain(geneticGain.df, 20)
plotGeneticGain(geneticGain.df, 50)