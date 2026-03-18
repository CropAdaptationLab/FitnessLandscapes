# Title: EPISTASIS
# Author: Ted Monyak
# Description: This file contains functions for calculating breeding values and
# genetic variance, and particularly epistatic variance w.r.t. fitness

# Calculates the coded genotypic values
# g_11: the genotypic value of the '11' homozygote
# g_10: the genotypic value of the heterozygote
# g_00: the genotypic value of the '00' homozygote
# Returns: a list of 2 items: 'a' and 'd'
coded_gvs <- function(g_11, g_10, g_00) {
  # Calculate the mid-parent value
  P <- (g_11 + g_00) / 2
  a <- g_11 - P # coded genotypic value of the '11' homozygote
  d <- g_10 - P # coded genotypic value of the heterozygote
  return (c(a, d))
}

# Calculate alpha
# p: allele frequency of '1'
# q: allele frequency of '0'
# a: coded genotypic value of '11' genotype
# d: coded genotypic value of the heterozygote
# Returns: a list of 3 items: alpha_1, alpha_0, and alpha
allele_substitution <- function(p, q, a, d) {
  alpha_1 <- q * (a+(d*(q-p))) # alpha for allele 1
  alpha_0 <- -p * (a+(d*(q-p))) # alpha for allele 0
  alpha <- alpha_1 - alpha_0 # substitution value
  return (c(alpha_1, alpha_0, alpha))
}

# Calculate additive variance and dominance deviation variance
# p: allele frequency of '1'
# q: allele frequency of '0'
# a: coded genotypic value of '11' genotype
# d: coded gentoypic value of the heterozygote
# Returns: a list of two items: additive variance and dominance deviation
genetic_variance <- function(p, q, a, d) {
  Va <- 2*p*q*((a+(d*(q-p)))^2) # Va = 2pq(a+d(q-p))^2
  Vd <- 4*(p^2)*(q^2)*(d^2) # Vd = 4p^2q^2d^2
  return(c(Va, Vd))
}

# Alternative method for calculating breeding values using allele substitution effect
# Not used currently
# p: allele frequency of '1'
# q: allele frequency of '0'
# d: coded gentoypic value of the heterozygote
# Returns a list of three breeding values
breeding_values <- function(alpha_1, alpha_0, p, q, d) {
  b_11 <- alpha_1 + alpha_1 - 2*(q^2)*d 
  b_10 <- alpha_1 + alpha_0 + 2*p*q*d
  b_00 <- alpha_0 + alpha_0 - 2*(p^2)*d
  return (c(b_11,b_10,b_00))
}

# Calculate the genetic variance between 2 loci
# merged.df: a dataframe with nInd rows and nQtl + nTraits columns
# qtlA: locus 1
# qtlB: locus 2
# trait: one of 'fitness', 'Trait1', or 'Trait2'
# Returns: a list of four values: Va, Vd, Vi, alpha (Va, Vd, and alpha are just for qtlA)
calc_gvar <- function(merged.df, qtlA, qtlB, trait) {
  # Calculate the mean phenotype
  mu <- mean(merged.df[[trait]])
  
  # Aggregate the dataframe to calculate the genotypic values based on deviation from the mean
  means_AB <- merged.df %>%
    dplyr::group_by(!!sym(qtlA), !!sym(qtlB)) %>%
    summarize(g=mean(!!sym(trait))-mu,
              n=n(),
              .groups = 'drop') %>%
    dplyr::mutate(freq=n/sum(n)) %>%
    complete(
      !!sym(qtlA) := 0:2,
      !!sym(qtlB) := 0:2,
      fill = list()
    ) %>%
    dplyr::mutate(across(everything(), ~replace_na(.x, 0)))
  
  # For qtlA, calculate the weighted mean of each genotypic value (11, 10, 00)
  # independent of qtlB
  means_A <- means_AB %>%
    dplyr::group_by(.data[[qtlA]]) %>%
    summarise(g = weighted.mean(g, w = n),
              n=sum(n),
              freq = sum(freq)) %>%
    column_to_rownames(var = {{qtlA}}) %>%
    dplyr::mutate(across(everything(), ~replace_na(.x, 0)))
  
  # Calculate allelic frequencies
  p <- means_A["2", "freq"] + means_A["1", "freq"]/2 # p is the '1' allele
  q <- means_A["0", "freq"] + means_A["1", "freq"]/2 # q is the '0' allele
  
  # For qtlB, calculate the weighted mean of each genotypic value (11, 10, 00)
  # independent of qtlA
  means_B <- means_AB %>%
    dplyr::group_by(.data[[qtlB]]) %>%
    summarise(g = weighted.mean(g, w = n),
              n= sum(n),
              freq = sum(freq)) %>%
    column_to_rownames(var = {{qtlB}}) %>%
    dplyr::mutate(across(everything(), ~replace_na(.x, 0)))
  
  # Calculate coded genotpic values (a, d)
  a_d <- coded_gvs(means_A["2","g"], means_A["1","g"], means_A["0","g"])
  a <- a_d[1]
  d <- a_d[2]
  
  # Calculate allele substitution effect
  alphas <- allele_substitution(p, q, a, d)
  alpha <- alphas[3]
  #alpha <- 0
  
  # Equivalent way of getting alpha
  #lm(fitness ~ `5_941`, merged.df)$coefficients[2]
  
  # Calculate breeding values
  # bvs_alt <- calc_bvs_alt(alpha_1, alpha_0, p, q, d)
  
  # Calculate genetic variances Va and Vd
  gvars <- genetic_variance(p, q, a, d)
  Va <- gvars[1]
  Vd <- gvars[2]
  
  # Calculate the epistatic effect for each combination of loci as:
  # The difference between the genotypic value and the sum of the genotypic values
  # of the two qtls independently
  means_AB <- means_AB %>%
    dplyr::mutate(I=g - means_A[as.character(.data[[qtlA]]), "g"] - means_B[as.character(.data[[qtlB]]), "g"]) %>%
    dplyr::mutate(VarI=(I^2)*freq)
  
  # Calculate a weighted average based on allele frequencies
  Vi <- sum(means_AB$VarI)
  
  return(c(Va, Vd, Vi, alpha))
}

# Calculate the total genetic variance for a population
# pop: The population
# calc_fitness: whether to calculate for the 'fitness' trait
# calc_traits: whether to calculate for the standard traits
# remove_hets: For mapping populations, throw away all individuals with heterozygosity
# Returns: a 3x3 dataframe with the following columns: Va, Vd, Vi, and the following
# rows: Fitness, Trait1, Trait2
population_gvar <- function(pop, calc_fitness=TRUE, calc_traits=FALSE, remove_hets=FALSE) {
  
  # Default to return
  default_result <- data.frame(Va=c(0,0,0),
                       Vd=c(0,0,0),
                       Vi=c(0,0,0),
                       row.names=c("Fitness", "Trait1", "Trait2"))

  # Return default result if not calculating fitness or traits
  if(!(calc_fitness | calc_traits)) {
    return (default_result)
  }
  
  # Randomly mate the population for 2 generations to put it at HWE
  hwe_pop <- randCross(pop, nCrosses=nInd(pop))
  hwe_pop <- randCross(hwe_pop, nCrosses=nInd(hwe_pop))
  
  # Obtain the genotype data
  qtl.df <- as.data.frame(getUniqueQtl(hwe_pop))
  
  # Remove all non-segregating loci
  qtl.df <- qtl.df[, sapply(qtl.df, segLocus)]
  
  # If there are no segregating loci, there is no genetic variance, so return the default
  if(ncol(qtl.df) == 0) {
    return (default_result)
  }
  
  # Remove all het individuals in a RIL
  if (remove_hets) {
    qtl.df <- qtl.df %>%
      rowwise() %>%
      mutate(het = any(c_across(everything()) == 1)) %>%
      ungroup() %>%
      filter(het==FALSE) %>%
      select(-het)
  }

  # Get the names of the qtls
  qtls <- colnames(qtl.df)
  # Nuber of qtl
  nQtl <- length(qtls)
  # Get the phenotypes for each of the traits
  pheno.df <- as.data.frame(pheno(hwe_pop))
  
  # Perform an inner join on the qtl (genotype) data and the phenotype data
  merged.df <- qtl.df %>%
    tibble::rownames_to_column(var = "row_name") %>%
    inner_join(
      pheno.df %>% tibble::rownames_to_column(var = "row_name"),
      by = "row_name"
    ) %>%
    tibble::column_to_rownames(var = "row_name") %>%
    arrange(as.integer(row.names(.))) %>%
    dplyr::mutate(fitness = calculateFitnessTwoTrait(Trait1, Trait2))
  
  # Create an nQtl x nQtl matrix to store the epistatic fitness variance at each combination
  # of qtl. Only the bottom half of the matrix will be populated
  Vi_mat_fit <- matrix(data=NA, nrow=nQtl, ncol=nQtl)
  rownames(Vi_mat_fit) <- qtls
  colnames(Vi_mat_fit) <- qtls

  # Create duplicate versions of this matrix to store the epistatic variance for
  # the two other traits
  Vi_mat_trait1 <- Vi_mat_fit
  Vi_mat_trait2 <- Vi_mat_fit

  # Create a vector to store the additive fitness variance for each qtl
  Va_vec_fit <- numeric(nQtl)
  # Create duplicate versions for each trait
  Va_vec_trait1 <- Va_vec_fit
  Va_vec_trait2 <- Va_vec_fit
  
  # Create a vector to store the fitness dominance deviation for each qtl
  Vd_vec_fit <- numeric(nQtl)
  # Create duplicate versions for each trait
  Vd_vec_trait1 <- Vd_vec_fit
  Vd_vec_trait2 <- Vd_vec_fit

  # Create a vector to store the allele substitution fitess effect for each qtl
  alpha_vec_fit <- numeric(nQtl)
  # Create duplicate versions for each trait
  alpha_vec_trait1 <- alpha_vec_fit
  alpha_vec_trait2 <- alpha_vec_fit
  
  # A nested for-loop for assessing each qtl x qtl combination
  for(a_idx in 1:nQtl) {
    qtlA <- qtls[a_idx]
    for(b_idx in 1:nQtl) {
      qtlB <- qtls[b_idx]
      # Only create the lower diagonal of the matrix, to avoid redundancy
      if (qtlA == qtlB) break
      
      # Only calculate fitness variance if specified
      if (calc_fitness) {
        vars_fit <- calc_gvar(merged.df, qtlA, qtlB, "fitness")
        Va_vec_fit[a_idx] <- vars_fit[1]
        Vd_vec_fit[a_idx] <- vars_fit[2]
        Vi_mat_fit[qtlA, qtlB] = vars_fit[3]
        alpha_vec_fit[a_idx] <- vars_fit[4]
      } else { # Just fill in zeroes
        Va_vec_fit[a_idx] = 0
        Vd_vec_fit[a_idx] = 0
        Vi_mat_fit[qtlA, qtlB] = 0
        alpha_vec_fit[a_idx] = 0
      }
      
      # Only calculate trait variance if specified
      if (calc_traits) {
        vars_trait1 <- calc_gvar(merged.df, qtlA, qtlB, "Trait1")
        Va_vec_trait1[a_idx] <- vars_trait1[1]
        Vd_vec_trait1[a_idx] <- vars_trait1[2]
        Vi_mat_trait1[qtlA, qtlB] = vars_trait1[3]
        
        vars_trait2 <- calc_gvar(merged.df, qtlA, qtlB, "Trait2")
        Va_vec_trait2[a_idx] <- vars_trait2[1]
        Vd_vec_trait2[a_idx] <- vars_trait2[2]
        Vi_mat_trait2[qtlA, qtlB] = vars_trait2[3]
      } else { # Just fill in zeroes
        Va_vec_trait1[a_idx] = 0
        Vd_vec_trait1[a_idx] = 0
        Vi_mat_trait1[qtlA, qtlB] = 0
        alpha_vec_trait1[a_idx] = 0
        Va_vec_trait2[a_idx] = 0
        Vd_vec_trait2[a_idx] = 0
        Vi_mat_trait2[qtlA, qtlB] = 0
        alpha_vec_trait2[a_idx] = 0
      }
    }
  }
  
  # Calculate total variance across all loci by summing
  Va_fit <- sum(Va_vec_fit)
  Vd_fit <- sum(Vd_vec_fit)
  Vi_fit <- sum(Vi_mat_fit[lower.tri(Vi_mat_fit, diag = FALSE)])
  
  Va_trait1 <- sum(Va_vec_trait1)
  Vd_trait1 <- sum(Vd_vec_trait1)
  Vi_trait1 <- sum(Vi_mat_trait1[lower.tri(Vi_mat_trait1, diag = FALSE)])
  
  Va_trait2 <- sum(Va_vec_trait2)
  Vd_trait2 <- sum(Vd_vec_trait2)
  Vi_trait2 <- sum(Vi_mat_trait2[lower.tri(Vi_mat_trait2, diag = FALSE)])
  
  # Create a 3 x 3 dataframe to return all 9 values
  return (data.frame(Va=c(Va_fit, Va_trait1, Va_trait2),
                     Vd=c(Vd_fit, Vd_trait1, Vd_trait2),
                     Vi=c(Vi_fit, Vi_trait1, Vi_trait2),
                     row.names=c("Fitness", "Trait1", "Trait2")))
}
