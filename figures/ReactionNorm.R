# TITLE: REACTION NORM
# Author: Ted Monyak
# Description: This file generates a reaction norm plot

# Create a reaction norm plot based on the genotypes at 2 specified markers
# Assumes that parents 1 and 2 are fixed for the opposite alleles at both markers
# pop: The admixed RIL family
# m1: The first marker in a putative interaction
# m2: The second marker in a putative interaction
# parent1: One of the parents from the biparental cross
# parent2: The other parent from the biparental cross
# suitFunc: See Fitness.R
# ril_dir: Where to save the resulting plot
plotReactionNorm <- function(pop, m1, m2, parent1, parent2, suitFunc, snpChip, ril_dir) {
  
  # Determine the mean fitness of each haplotype w.r.t. two markers
  # Returns: a 2x2 matrix storing the mean fitness for each haplotype,
  # With the colnames and rownames designating which parental haplotype
  # The alleles correspond to
  getHaploMeans <- function(pop, m1, m2, parent1, parent2, snpChip, suitFunc) {
    # Vectors to store the phenotypes of each haplotype
    pop_0_0 <- vector(mode = "numeric", length = 0)
    pop_0_2 <- vector(mode = "numeric", length = 0)
    pop_2_0 <- vector(mode = "numeric", length = 0)
    pop_2_2 <- vector(mode = "numeric", length = 0)
    
    geno <- as.data.frame(pullSnpGeno(pop, snpChip))
    # Check that both QTL are segregating
    if (!segLocus(geno[,m1]) || !segLocus(geno[,m2])) {
      return (NA)
    }
    pheno <- pheno(pop)
    # Iterate through each individual to determine its haplotype
    for (i in 1:nInd(pop)){
      # Determine the genotype for that individual for each marker
      geno_1 <- geno[i, m1]
      geno_2 <- geno[i, m2]
      # Calculate W
      fitness <- calculateBreedingFitness(t1=pheno[i,1], t2=pheno[i,2], t3=pheno[i,3], suitFunc=suitFunc)
      if (geno_1 == 0 & geno_2 == 0) {
        pop_0_0 <- append(pop_0_0, fitness)
      }
      if (geno_1 == 0 & geno_2 == 2) {
        pop_0_2 <- append(pop_0_2, fitness)
      }
      if (geno_1 == 2 & geno_2 == 0) {
        pop_2_0 <- append(pop_2_0, fitness)
      }
      if (geno_1 == 2 & geno_2 == 2) {
        pop_2_2 <- append(pop_2_2, fitness)
      }
    }
    # There must be some individuals with each haplotype to get a rxn norm
    if (length(pop_0_0) == 0 ||
        length(pop_0_2) == 0 ||
        length(pop_2_0) == 0 ||
        length(pop_2_2) == 0) {
      return (NA)
    }
    
    # Create a matrix with the mean values of each haplotype
    hap_mat <- matrix(c(mean(pop_0_0), mean(pop_0_2), # M1: 00
                        mean(pop_2_0), mean(pop_2_2)), # M2: 11
                      nrow = 2,
                      byrow = TRUE)
    
    # Get parent 1's haplotype
    haplo1 <- pullSnpGeno(parent1, snpChip) %>%
      dplyr::select(all_of(c(m1, m2)))
    
    # Parent 1's allele at M1
    p1_m1 <- haplo1[1,1]
    
    # Parent 1's allele at M2
    p1_m2 <- haplo1[1,2]
    
    # Assign the rownames based on parent 1s haplotype at the first QTL
    if (p1_m1 == 0) {
      rownames(hap_mat) <- c("P1", "P2")
    } else {
      rownames(hap_mat) <- c("P2", "P1")
    }
    
    # Assign the column names based on parent 1s haplotype at the second QTL
    if (p1_m2 == 0) {
      colnames(hap_mat) <- c("P1", "P2")
    } else {
      colnames(hap_mat) <- c("P2", "P1")
    }
    
    return (hap_mat)
  }
  hap_mat <- getHaploMeans(pop, m1, m2, parent1, parent2, snpChip, suitFunc)
  if (is.na(hap_mat)) {
    return (NA)
  }
  
  m1_labels <- colnames(hap_mat) # Y-axis
  m2_labels <- rownames(hap_mat) # X-axis
  
  # Store the mean fitness as a long dataframe
  df <- data.frame(
    m1 = rep(m1_labels, each = 2),
    m2 = rep(m2_labels, times = 2),
    fitness = as.vector(t(hap_mat))
  )
  
  # Create a rxn norm plot
  rn <- ggplot(df, aes(x = m2, y = fitness, group = m1, color = m1)) +
    geom_line(size = 1) +
    geom_point() +
    scale_color_manual(values = c("P1" = "#CC0000", "P2"="#3C78D8"),
                       labels = c("Parental Type 1", "Parental Type 2")) +
    scale_x_discrete(labels = c("Parental\nType 1", "Parental\nType 2"),
                     expand = c(0.1, 0.1)) +
    labs(
      x = paste0(m2, " Genotype"),
      y = "Breeding Fitness",
      color = paste0(m1, " Genotype")
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
      panel.background = element_rect(fill = "white", color = "black")
    )  +
    ylim(min(df$fitness)*0.97, max(df$fitness)*1.03)

  ggplot2::ggsave(filename = file.path(ril_dir, "rxn_norm.jpg"),
                  device = "jpg",
                  height=2,
                  width=3.5,
                  units="in",
                  dpi=600)
  ggplot2::ggsave(filename = file.path(ril_dir, "rxn_norm.pdf"),
                  device = "pdf",
                  height=2,
                  width=3.5)
}

