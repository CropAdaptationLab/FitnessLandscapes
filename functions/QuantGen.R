library(genetics)

# TODO remove this funciton
# Calculates 'isogenicity' as a weighted sum of the genotypic differences between
# two individuals multiplied by the effect size of each QTL
# ind1: one of the individuals
# ind2: the other individual
# founderPop: the original founderPop from
# Returns: a float between 0 (allo-elite) and 1 (iso-elite)
isogenicity <- function(ind1, ind2, founderPop) {
  # Get the effect sizes of all original QTL from the population
  eff_sizes <- getQtlEffectSizes(founderPop)
  # Get the sum of all the QTL effect sizes, to normalize effect sizes
  tot_eff_size <- sum(eff_sizes$eff_size)
  
  # Join the qtl genotypes of both individuals
  qtlGeno <- as.data.frame(t(rbind(getUniqueQtl(ind1), getUniqueQtl(ind2))))
  
  # Determine the allelic differences between the individuals, and a weighted score
  # of each difference
  iso_elite.df <- merge(eff_sizes, qtlGeno, by="row.names", all=FALSE) %>%
    dplyr::mutate(diff=1-(abs(V2-V1)/2)) %>%
    dplyr::mutate(val=diff*eff_size/tot_eff_size)
  
  # Return the weighted sum
  return (sum(iso_elite.df$val))
}

# Calculates 'isogenicity' of the desired trait
isoEliteness <- function(ind1, ind2, founderPop, trait) {
  # Get the effect sizes of all original QTL from the population
  eff_sizes <- getSingleTraitEffectSizes(founderPop, trait) %>%
    dplyr::mutate(eff_size=abs(eff_size))
  # Get the sum of all the QTL effect sizes, to normalize effect sizes
  tot_eff_size <- sum(eff_sizes$eff_size)
  
  # Join the qtl genotypes of both individuals
  qtlGeno <- as.data.frame(t(rbind(pullQtlGeno(ind1, trait), pullQtlGeno(ind2, trait))))
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
# only considering QTL
# ind1: one of the individuals
# ind2: the other individual
# founderPop: the original founderPop from
# Returns: an integer hamming distance
hammingDistance <- function(ind1, ind2, trait) {
  # Join the qtl genotypes of both individuals and calculate the number of different
  # alleles
  qtlGeno <- as.data.frame(t(rbind(pullQtlGeno(ind1,trait), pullQtlGeno(ind2,trait))))
  colnames(qtlGeno) <- c("ind1", "ind2")
  qtlGeno$ind1 <- as.numeric(qtlGeno$ind1)
  qtlGeno$ind2 <- as.numeric(qtlGeno$ind2)
  
  qtlGeno <- qtlGeno %>%
    dplyr::mutate(diff=(abs(ind1-ind2)))
  # Return the weighted sum
  return (sum(qtlGeno$diff))
}

excessVariance <- function(RIL_pheno, pureline_pheno) {
  return ((var(RIL_pheno)-var(pureline_pheno)) / var(pureline_pheno))
}

# Train an RRBLUP model on trainPop
# Evaluate it on testPop
# Return correlation between the EBVs and the actual genetic values
evaluateGWP <- function(trainPop, testPop, trait) {
  # Update phenotype to have heritabilities associated with breeding programs
  trainPop <- setPheno(trainPop, h2=c(n.h2Breeding, n.h2Breeding, n.yieldH2Breeding))
  model <- fastRRBLUP(trainPop, traits=trait, use="pheno")
  testPop <- setEBV(testPop, model)
  r <- cor(gv(testPop), ebv(testPop))[trait]
  return (r)
}

# Train an RRBLUP model on trainPop
# Evaluate it on testPop
# Return correlation between the EBVs and the actual genetic values
evaluateGWP_W <- function(trainPop, testPop) {
  # Update phenotype to have heritabilities associated with breeding programs
  trainPop <- setPheno(trainPop, h2=c(n.h2Breeding, n.h2Breeding, n.yieldH2Breeding))
  model <- fastRRBLUP(trainPop, traits=calculateW_GWP, use="gv")
  testPop <- setEBV(testPop, model)
  r <- cor(calculateW_GWP(gv(testPop)), ebv(testPop))[1]
  return (r)
}

evaluateGPRR <- function(trainPop, testPop, trait) {
  library(rrBLUP)
  trainPop@pheno[1,]
  model <- mixed.solve(trainPop@pheno[,1], Z=pullSnpGeno(trainPop))
  cor(bv(trainPop)[,1], model$u)
  model$u
  # Update phenotype to have heritabilities associated with breeding programs
  trainPop <- setPheno(trainPop, h2=c(n.h2Breeding, n.h2Breeding, n.yieldH2Breeding))
  model <- fastRRBLUP(trainPop, traits=trait, use="pheno")
  testPop <- setEBV(testPop, model)
  r <- cor(gv(testPop), ebv(testPop))[trait]
  return (r)
  
  M <- matrix(rep(0,200*1000),200,1000)
  for (i in 1:200) {
    M[i,] <- ifelse(runif(1000)<0.5,-1,1)
  }
  u <- rnorm(1000)
  g <- as.vector(crossprod(t(M),u))
  h2 <- 0.5 #heritability
  y <- g + rnorm(200,mean=0,sd=sqrt((1-h2)/h2*var(g)))
  K <- A.mat(M)
}
