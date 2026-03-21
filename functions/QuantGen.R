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

# Train an RRBLUP model to predict one of the traits
# trainPop: the training population
# testPop: the test population
# trait: the AlphaSimR phenotype index
# Return correlation (r) between the EBVs and the actual genetic values in the test pop
evaluateGWP <- function(trainPop, testPop, trait) {
  # Update phenotype to have heritability associated with breeding programs
  trainPop <- setPheno(trainPop, h2=c(n.h2Breeding, n.h2Breeding, n.yieldH2Breeding))
  # Train the model
  model <- fastRRBLUP(trainPop, traits=trait, use="pheno")
  # Set the estimated breeding values
  testPop <- setEBV(testPop, model)
  # Determine the correlation between genetic values and estimated breeding values
  # in the test population
  r <- cor(gv(testPop), ebv(testPop))[trait]
  return (r)
}

# Train an RRBLUP model to predict breeding fitness
# trainPop: the training population
# testPop: the test population
# trait: the AlphaSimR phenotype index
# Return correlation (r) between the EBVs and the actual genetic values in the test pop
evaluateGWP_W <- function(trainPop, testPop) {
  # Update phenotype to have heritabilities associated with breeding programs
  trainPop <- setPheno(trainPop, h2=c(n.h2Breeding, n.h2Breeding, n.yieldH2Breeding))
  model <- fastRRBLUP(trainPop, traits=calculateW_GWP, use="gv")
  testPop <- setEBV(testPop, model)
  # Determine the correlation between genetic values and estimated breeding values
  # in the test population
  r <- cor(calculateW_GWP(gv(testPop)), ebv(testPop))[1]
  return (r)
}
