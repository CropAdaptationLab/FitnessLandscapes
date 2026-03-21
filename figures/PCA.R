# Title: PCA
# Author: Ted Monyak
# Description: this file contains functions for plotting the principal component analysis of populations
library(ade4)
library(akima)
library(factoextra)
library(ggdensity)
library(stats)

#sampled_inds <- unlist(sampled_inds)
#sampled_inds_geno <- rbind(
#  sampled_inds[[1]],
#  sampled_inds[[2]],
#  sampled_inds[[3]],
#  sampled_inds[[4]]
#)

#plotPCA_pop(pca_pops=pops, project_pop=createNAM())

# Plots a PCA of a mapping population
plotPCA_pop <- function(pca_pops, project_pop) {
  pca_geno <- pullSnpGeno(pops[[1]])
  for (p in 2:length(pops)) {
    pca_geno <- rbind(pca_geno,
                      pullSnpGeno(pops[[p]]))
  }
  pca_geno <- rbind(pca_geno, pullSnpGeno(founderPop))
  
  project_geno <- pullSnpGeno(RIL)
  pca <- prcomp(pca_geno)
  VAF <- summary(pca)$importance[2,1:2] * 100
  colnames(project_geno) <- colnames(pca_geno)
  pc_pred <- predict(pca, project_geno)

  pheno.df <- data.frame(Suit=gaussianFitFunc(pheno(RIL)))
  
  r2_per_pc <- sapply(1:ncol(PCs), function(i) {
    m <- lm(pheno.df$Suit ~ pc_pred[, i])
    summary(m)$r.squared
  })
  
  names(r2_per_pc) <- colnames(PCs)
  r2_sorted <- sort(r2_per_pc, decreasing = TRUE)
  
  head(r2_sorted, 20)
  
  df.PCA = data.frame(
    Suit = pheno.df$Suit,
    "PC1" = pc_pred[, 829],
    "PC2" = pc_pred[, 1514]
  )
  
  

  ggplot(pheno.df, aes(x=Suit)) +
    geom_density()
  
  # For sampled pops
  #df.PCA$pop <- as.factor(df.PCA$pop)
  #ggplot(df.PCA, aes(x = PC1, y = PC2, color=pop)) +
  #  geom_line()

  #plot_NAM_family(df.PCA)
  
  plot_NAM_suit(df.PCA)
  
  smoothed <- generate_landscape(df.PCA)

  render_2d_landscape(smoothed)
  
  df.PCA$Suit_smoothed <- fields::interp.surface(
    smoothed, 
    cbind(df.PCA$PC1, df.PCA$PC2)
  )
  df.PCA$Suit_smoothed <- df.PCA$Suit_smoothed + 0.01
  
  p <- render_3d_landscape(smoothed, df.PCA, families)
  p
}

# Plots a PCA of all of the subpopulations
# pops: a list of AlphaSim populations
# save_dir: where to write the plot
# gen: the generation of the populations (for naming)
plotPCA_RIL <- function(pops, save_dir, gen) {
  # Collect the genotype and phenotype data of each population
  subpops_geno <- pullSnpGeno(pops[[1]])
  if (length(pops) > 1) {
    for (p in 2:length(pops)){
      subpops_geno <- rbind(subpops_geno, pullSnpGeno(pops[[p]]))
    }
  }
  
  RIL_geno <- pullSnpGeno(RIL)
  
  pca <- prcomp(geno)
  PCs <- pca$x
  
  pheno.df <- data.frame(fitness=gaussianLandraceFitFunc(rbind(pheno(pops[[1]]), pheno(pops[[2]]))))
  
  r2_per_pc <- sapply(1:ncol(PCs), function(i) {
    m <- lm(pheno.df$fitness ~ PCs[, i])
    summary(m)$r.squared
  })
  
  names(r2_per_pc) <- colnames(PCs)
  r2_sorted <- sort(r2_per_pc, decreasing = TRUE)
  
  head(r2_sorted, 20)
  
  
  RIL_geno <- getUniqueQtl(RIL, traits=c(1,2))
  subpops_geno <- rbind(getUniqueQtl(pop1, traits=c(1,2)),
                        getUniqueQtl(pop2, traits=c(1,2)))
  
  
  pca <- prcomp(geno)
  PCs <- pca$x
  
  
  RIL_pheno.df <- data.frame(fitness=gaussianLandraceFitFunc(pheno(RIL)),
                             haplo=haplo$cat)
  
  pheno.df <- data.frame(fitness=gaussianLandraceFitFunc(rbind(pheno(pops[[1]]), pheno(pops[[2]]))))
  
  r2_per_pc <- sapply(1:ncol(PCs), function(i) {
    m <- lm(pheno.df$fitness ~ PCs[, i])
    summary(m)$r.squared
  })
  
  names(r2_per_pc) <- colnames(PCs)
  r2_sorted <- sort(r2_per_pc, decreasing = TRUE)
  
  head(r2_sorted, 20)
  
  # Calculate the fitness of each individual
  fitness.df <- as.data.frame(pheno) %>%
    dplyr::mutate(fitness = fitCalc(Trait1, Trait2)) %>%
    dplyr::select(fitness)
  
  # Plot the PCA
  VAF <- summary(pca)$importance[2,1:2] * 100
  
  df.PCA = data.frame(
    "Fitness" = fitness.df$fitness,
    "PC1" = pca$x[,1],
    "PC2" = pca$x[,2]
  )

  ggplot(df.PCA, aes(x = PC1, y = PC2)) +
    ggtitle("Population structure") +
    geom_point(aes(colour = Fitness)) +
    scale_color_viridis(alpha=0.9) +
    xlab(paste("PC1: ", round(VAF[1]), "%", sep = "")) +
    ylab(paste("PC2: ", round(VAF[2]), "%", sep = ""))
  
  fname <- file.path(save_dir, paste0("PCA_RIL", gen, ".jpg"))
  ggplot2::ggsave(filename = fname,
                  device = "jpg",
                  width=8,
                  height=8)
}

plot_1D_PCA <- function(RIL, pop1, pop2, parent1, parent2, qtl1, qtl2, save_dir) {
  # Make a larger RIL to reduce noise
  RIL <- self(RIL, nProgeny=5)
  RIL_geno <- getUniqueQtl(RIL, traits=c(1,2))
  subpops_geno <- rbind(getUniqueQtl(pop1, traits=c(1,2)),
                        getUniqueQtl(pop2, traits=c(1,2)))

  haplo1 <- getUniqueQtl(parent1) %>%
    dplyr::select(all_of(c(qtl1, qtl2)))
  
  haplo2 <- getUniqueQtl(parent2) %>%
    dplyr::select(all_of(c(qtl1, qtl2)))
  
  categorize <- function(allele1, allele2) {
    if (all(c(allele1, allele2) == haplo1)) {
      return("P1")
    } else if (all(c(allele1, allele2) == haplo2)) {
      return("P2")
    }
    return("R")
  }

  haplo <- RIL_geno %>%
    dplyr::select(all_of(c(qtl1,qtl2))) %>%
    setNames(c("qtl1", "qtl2")) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(cat=categorize(qtl1, qtl2)) %>%
    dplyr::ungroup()
  
  RIL_pheno.df <- data.frame(fitness=breedingFitness(pheno(RIL)),
                             haplo=haplo$cat)
  
  pca <- prcomp(subpops_geno)
  VAF <- summary(pca)$importance[2,1] * 100
  pc_RIL <- predict(pca, RIL_geno)
  
  df.PCA = data.frame(
    "Fitness" = RIL_pheno.df$fitness,
    "PC1" = pc_RIL[, 1],
    "PC2" = pc_RIL[, 2],
    "Haplotype" = RIL_pheno.df$haplo)
  
  ggplot(df.PCA, aes(x = PC1, y = Fitness)) +
    geom_hdr(aes(fill = Haplotype), alpha = 0.3, color = NA, probs = seq(0.1, 0.9, by = 0.1)) +
    geom_point(aes(color = Haplotype), alpha = 0.1, size = 0.8) +
    geom_smooth(
      #aes(group = Haplotype, color = Haplotype),
      color="black", # comment this and uncomment above to get haplo-specific lines
      method = "loess",
      se = FALSE,
      linewidth = 0.8
    ) +
    scale_color_manual(values = c("#CC0000", "#3C78D8", "gold"),,
                       labels = c("Parental Type 1", "Parental Type 2", "Recombinant"),
                       name="Haplotype at\nEpistatic Loci") +
    scale_fill_manual(values = c("#CC0000", "#3C78D8", "gold"),
                      labels = c("Parental Type 1", "Parental Type 2", "Recombinant"),
                      name="Haplotype at\nEpistatic Loci") +
    guides(alpha="none") +
    xlab("PC1 of Oligogenic Component") +
    ylab("Breeding Fitness") +
    theme_minimal(base_size = 14,
                  base_family="Helvetica") +
    theme(
      plot.title = element_blank(),
      axis.title.y = element_text(margin=margin(t=0, r=10, b=0, l=10, unit="pt")),
      axis.text.x = element_blank(),
      legend.position = "right",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white")
    )
  ggplot2::ggsave(filename = file.path(save_dir, "PCA_Bridge.jpg"),
                  device = "jpg",
                  width=5.5,
                  height=3.5,
                  dpi=600)
  
  ggplot2::ggsave(filename = file.path(save_dir, "PCA_Bridge.pdf"),
                  device = "pdf",
                  width=5.5,
                  height=3.5)
}

plot_NAM_family <- function(pca_plot_df, VAF, pcx=1, pcy=2) {
  # Plot PCA with colors by group
  xcol <- colnames(pca_plot_df)[pcx]
  ycol <- colnames(pca_plot_df)[pcy]
  
  ggplot(pca_plot_df, aes(x = .data[[xcol]], y = .data[[ycol]], color = family_color)) +
    geom_point(data = subset(pca_plot_df, grepl("_RIL$", Family)),
              size = 0.5, shape = 16) +
    geom_point(data = subset(pca_plot_df, !grepl("_RIL$", Family)),
               size = 3, shape = 18) +
    scale_color_identity(guide="legend",
                         breaks = pca_plot_df$family_color,
                         labels = pca_plot_df$Family,
                         name="Family") +
    guides(color = guide_legend(override.aes = list(size = 2))) +
    labs(x = paste0("PC", pcx),
         y = paste0("PC", pcy)) +
    #labs(x = paste0("PC", pcx, ": ", round(VAF[pcx]), "%"),
    #     y = paste0("PC", pcy, ": ", round(VAF[pcy]), "%")) +
    theme_minimal(base_family="Helvetica", base_size=8) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 1),
          axis.title.x = element_text(hjust = 0.5),
          axis.title.y = element_text(hjust = 0.5),
          axis.text = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "white", color = "black"),
          aspect.ratio=1)
}

plot_NAM_suit <- function(pca_plot_df, VAF, pcx=1, pcy=2) {
  # Plot PCA with colors by group
  xcol <- colnames(pca_plot_df)[pcx]
  ycol <- colnames(pca_plot_df)[pcy]
  
  ggplot(pca_plot_df, aes(x = .data[[xcol]], y = .data[[ycol]])) +
    geom_point(data = subset(pca_plot_df, grepl("_RIL$", Family)),
               size = 0.5, shape = 16,
               aes(colour = Suit, alpha=0.9), show.legend = c(color = TRUE, size = FALSE, alpha = FALSE)) +
    geom_point(data = subset(pca_plot_df, !grepl("_RIL$", Family)),
               size = 3, shape = 18,
               aes(colour = Suit, alpha=0.9), show.legend = c(color = TRUE, size = FALSE, alpha = FALSE)) +
    #geom_point(size=0.2) +
    scale_color_viridis(name="Suitability") +
    #scale_color_distiller(palette = "Greys", name = "Suitability", guide="legend") +
    #labs(x = paste0("PC", pcx, ": ", round(VAF[pcx]), "%"),
    #     y = paste0("PC", pcy, ": ", round(VAF[pcy]), "%")) +
    labs(x = paste0("PC", pcx),
         y = paste0("PC", pcy)) +
    theme_minimal(base_family="Helvetica", base_size=8) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 1),
          axis.title.x = element_text(hjust = 0.5),
          axis.title.y = element_text(hjust = 0.5),
          axis.text = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = "white", color = "black"),
          aspect.ratio=1)
}



