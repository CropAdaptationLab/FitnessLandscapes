# Title: QUANT GEN
# Author: Ted Monyak
# Description: contains functions for calculating various quant gen metrics


# Calculates 'isoeliteness' of the specified trait as the weighted sum of
# the allelic differences between individuals multiplied by QTL effect sizes
# ind1: one of the individuals
# ind2: the other individual
# founderPop: the original founderPop from which both individuals came
# trait: the index of the AlphaSim phenotype to measure
# Returns: a float between 0 (allo-elite) and 1 (iso-elite)
isoEliteness <- function(ind1, ind2, founderPop, trait) {
  # Get the effect sizes of all original QTL from the population
  eff_sizes <- singleTraitArchitecture(founderPop, trait) %>%
    tibble::column_to_rownames(var="id")
  # Get the sum of all the QTL effect sizes, to normalize effect sizes
  tot_eff_size <- sum(eff_sizes$eff_size)
  
  # Join the qtl genotypes of both individuals
  qtlGeno <- as.data.frame(t(rbind(pullQtlGeno(ind1, trait),
                                   pullQtlGeno(ind2, trait))))
  colnames(qtlGeno) <- c("ind1", "ind2")
  qtlGeno$ind1 <- as.numeric(qtlGeno$ind1)
  qtlGeno$ind2 <- as.numeric(qtlGeno$ind2)
  
  # Determine the allelic differences between the individuals, and a weighted score
  # of each difference
  iso_elite.df <- merge(eff_sizes, qtlGeno, by="row.names", all=FALSE) %>%
    dplyr::mutate(diff=1-(abs(ind1-ind2) /2)) %>%
    dplyr::mutate(val=diff*eff_size/tot_eff_size)
  
  # Return the weighted sum
  return (sum(iso_elite.df$val))
}

# Calculates Hamming Distance - the number of genotypic differences
# only considering QTL underlying the specified trait
# ind1: one of the individuals
# ind2: the other individual
# founderPop: the original founderPop from which both individuals came
# trait: The AlphaSimR index of the phenotype
# Returns: an integer hamming distance
hammingDistance <- function(ind1, ind2, trait) {
  # Join the qtl genotypes of both individuals and calculate the number of different
  # alleles
  qtlGeno <- as.data.frame(t(rbind(pullQtlGeno(ind1,trait),
                                   pullQtlGeno(ind2,trait))))
  colnames(qtlGeno) <- c("ind1", "ind2")
  qtlGeno$ind1 <- as.numeric(qtlGeno$ind1)
  qtlGeno$ind2 <- as.numeric(qtlGeno$ind2)
  
  qtlGeno <- qtlGeno %>%
    dplyr::mutate(diff=(abs(ind1-ind2)))
  # Return the weighted sum
  return (sum(qtlGeno$diff))
}

# Calculates the EV in the admixed RIL family for a particular trait
# RIL_pheno: a dataframe with the phenotypes in the RIL family
# pureline_pheno: a dataframe with the combined phenotypes of the two purelines
# parents used to create a RIL family
# Returns: a float, EV, of the magnitude of increase of the variance for that phenotype
# in the RIL family compared to the purelines
excessVariance <- function(RIL_pheno, pureline_pheno) {
  return ((var(RIL_pheno)-var(pureline_pheno)) / var(pureline_pheno))
}


# Calculate mean heterozygosity across all loci
# Returns: a float between 0 and 1
meanHetLocus <- function(geno) {
  nInd <- nrow(geno)
  # Count the number of heterozygotes per locus
  het <- apply(geno, MARGIN=2, FUN= function(x) sum(x==1)/nInd)
  return (mean(het))
}

# Calculate the mean isoeliteness across the entire population
# TODO ADD DESCRIPTION
popIsoeliteness <- function(pop) {
  nInd <- nInd(pop)
  # Get the effect sizes of all original QTL from the population
  eff_sizes <- getQtlEffectSizes(founderPop)

  # Get the sum of all the QTL effect sizes, to normalize effect sizes
  tot_eff_size <- sum(eff_sizes$eff_size)
  
  # Get allele frequencies
  p <- as.data.frame(apply(getUniqueQtl(pop),
                           MARGIN=2,
                           FUN=function(x)
                             (sum(x==n.allele)/nInd) + ((sum(x==1)/nInd)/2))) %>%
    tibble::rownames_to_column(var="id")
  
  colnames(p) <- c("id", "freq")
  
  # Determine the degree of segregation
  # Multiply the effect size by degree of segregation
  iso_elite.df <- merge(eff_sizes, p, by="id", all=FALSE) %>%
    dplyr::mutate(seg=abs(0.5-freq)/0.5) %>%
    dplyr::mutate(val=seg*eff_size/tot_eff_size)
  
  # Return the weighted sum
  return (sum(iso_elite.df$val))
}

# Get the haplotypes of each individual in the population based on 2 markers
# and the parental haplotypes found in parent1 and parent2
# If useQtls is set to true, use QTL markers only, otherwise use the markers
# stored in the snpChip
# Returns a list of size nInd(rilGeno) of 'P1', 'P2', and 'R'
getHaplos <- function(rilGeno, m1, m2, parent1, parent2, snpChip=2, useQtls=TRUE) {
  # Determine the haplotypes of parent1 and parent 2 w.r.t. the two QTL markers
  if (m1 == m2) warning("m1 and m2 are the same marker: ", m1)

  if (useQtls) {
    haplo1 <- getUniqueQtl(parent1) %>%
      dplyr::select(all_of(c(m1, m2)))
    haplo2 <- getUniqueQtl(parent2) %>%
      dplyr::select(all_of(c(m1, m2)))
  } else {
    haplo1 <- as.data.frame(pullSnpGeno(parent1, snpChip)) %>%
      dplyr::select(all_of(c(m1, m2)))
    haplo2 <- as.data.frame(pullSnpGeno(parent2, snpChip)) %>%
      dplyr::select(all_of(c(m1, m2)))
  }

  
  # Determine whether a haplotype corresponds to parent 1, parent 2, or neither (recombinant)
  # allele1: The allele at the first locus
  # allele2: The allele at the second locus
  categorize <- function(allele1, allele2) {
    if (all(c(allele1, allele2) == haplo1)) {
      return("P1")
    } else if (all(c(allele1, allele2) == haplo2)) {
      return("P2")
    }
    return("R")
  }
  
  # Determine the parental haplotypes of each individual in the RIL family  
  haplo <- tibble(m1 = rilGeno[[m1]],
                  m2 = rilGeno[[m2]]) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(cat = categorize(m1, m2)) %>%
    dplyr::ungroup()
  
  return (haplo$cat)
}
