selectHaplo <- function(pop, qtl1, qtl2, parent1, parent2, fitCalc) {
  pop_0_0 <- vector(mode = "numeric", length = 0)
  pop_0_1 <- vector(mode = "numeric", length = 0)
  pop_0_2 <- vector(mode = "numeric", length = 0)
  pop_1_0 <- vector(mode = "numeric", length = 0)
  pop_1_1 <- vector(mode = "numeric", length = 0)
  pop_1_2 <- vector(mode = "numeric", length = 0)
  pop_2_0 <- vector(mode = "numeric", length = 0)
  pop_2_1 <- vector(mode = "numeric", length = 0)
  pop_2_2 <- vector(mode = "numeric", length = 0)

  qtls <- getUniqueQtl(pop)
  if (!segLocus(qtls[,qtl1]) || !segLocus(qtls[,qtl2])) {
    return (NA)
  }
  pheno <- pheno(pop)
  for (i in 1:nInd(pop)){
    geno_1 <- qtls[i, qtl1]
    geno_2 <- qtls[i, qtl2]
    fitness <- fitCalc(pheno[i,1], pheno[i,2], pheno[i,3])
    if (geno_1 == 0 & geno_2 == 0) {
      pop_0_0 <- append(pop_0_0, fitness)
    }
    if (geno_1 == 0 & geno_2 == 1) {
      pop_0_1 <- append(pop_0_1, fitness)
    }
    if (geno_1 == 0 & geno_2 == 2) {
      pop_0_2 <- append(pop_0_2, fitness)
    }
    if (geno_1 == 1 & geno_2 == 0) {
      pop_1_0 <- append(pop_1_0, fitness)
    }
    if (geno_1 == 1 & geno_2 == 1) {
      pop_1_1 <- append(pop_1_1, fitness)
    }
    if (geno_1 == 1 & geno_2 == 2) {
      pop_1_2 <- append(pop_1_2, fitness)
    }
    if (geno_1 == 2 & geno_2 == 0) {
      pop_2_0 <- append(pop_2_0, fitness)
    }
    if (geno_1 == 2 & geno_2 == 1) {
      pop_2_1 <- append(pop_2_1, fitness)
    }
    if (geno_1 == 2 & geno_2 == 2) {
      pop_2_2 <- append(pop_2_2, fitness)
    }
  }
  if (length(pop_0_0) == 0 ||
      length(pop_0_2) == 0 ||
      length(pop_2_0) == 0 ||
      length(pop_2_2) == 0) {
    return (NA)
  }
  
  hap_mat <- matrix(c(mean(pop_0_0), mean(pop_0_2), # QTL1: 00
                      mean(pop_2_0), mean(pop_2_2)), # QTL2: 11
                    nrow = 2,
                    byrow = TRUE)
  
  haplo1 <- getUniqueQtl(parent1) %>%
    dplyr::select(all_of(c(qtl1, qtl2)))
  
  # This assumes that parents 1 and 2 are fixed for the opposite alleles at both QTLs
  # Parent 1's allele at QTL1
  p1_qtl1 <- haplo1[1,1]
  
  # Parent 1's allele at QTL2
  p1_qtl2 <- haplo1[1,2]
  
  if (p1_qtl1 == 0) {
    rownames(hap_mat) <- c("P1", "P2")
  } else {
    rownames(hap_mat) <- c("P2", "P1")
  }
  
  if (p1_qtl2 == 0) {
    colnames(hap_mat) <- c("P1", "P2")
  } else {
    colnames(hap_mat) <- c("P2", "P1")
  }

  return (hap_mat)
}

plot_reaction_norm <- function(pop, qtl1, qtl2, parent1, parent2, fitCalc, save_dir) {
  hap_mat <- selectHaplo(pop, qtl1, qtl2, parent1, parent2, fitCalc)
  
  qtl1_labels <- colnames(hap_mat) # Rows - will be different lines
  qtl2_labels <- rownames(hap_mat)  # Columns - x-axis
  
  df <- data.frame(
    qtl1 = rep(qtl1_labels, each = 2),
    qtl2 = rep(qtl2_labels, times = 2),
    fitness = as.vector(t(hap_mat))
  )
  
  rn <- ggplot(df, aes(x = qtl2, y = fitness, group = qtl1, color = qtl1)) +
    geom_line(size = 1) +
    geom_point() +
    scale_color_manual(values = c("P1" = "#CC0000", "P2"="#3C78D8"),
                       labels = c("Parental Type 1", "Parental Type 2")) +
    scale_x_discrete(labels = c("Parental\nType 1", "Parental\nType 2"),
                     expand = c(0.1, 0.1)) +
    labs(
      x = paste0(qtl2, " Genotype"),
      y = "Breeding Fitness",
      color = paste0(qtl1, " Genotype")
    ) +
    theme_minimal(base_size = 10,
                  base_family="Helvetica") +
    theme(
      plot.title = element_blank(),
      axis.title.y = element_text(margin=margin(t=0, r=10, b=0, l=10, unit="pt")),
      legend.position = "right",
      legend.key.spacing.y = unit(5, "pt"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      #axis.line = element_line(color = "black"),
      panel.background = element_rect(fill = "white", color = "black")
    )  +
    ylim(min(df$fitness)*0.97, max(df$fitness)*1.03)

  ggplot2::ggsave(filename = file.path(save_dir, "rxn_norm.jpg"),
                  device = "jpg",
                  height=2,
                  width=3.5,
                  units="in",
                  dpi=600)
  ggplot2::ggsave(filename = file.path(save_dir, "rxn_norm.pdf"),
                  device = "pdf",
                  height=2,
                  width=3.5)
}

