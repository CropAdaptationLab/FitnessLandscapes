# Title: TRAIT ARCHITECTURE
# Author: Ted Monyak
# Description: Contains functions for examining and plotting trait architectures

# Determine whether a population is segregating at a locus
# locus: a list of genotypes from a population
# Returns: true if a locus is segregating within a population
segLocus <- function(locus) {
  return (length(unique(locus)) > 1)
}
# Get the QTL of all traits in a population
# Since QTL can be overlapping between traits, we must filter out the
# duplicated QTL
# pop: an AlphaSimR population
# traits: a vector specifying the indices of traits to acquire
# Returns: a dataframe of the QTL (where columns are the QTL and rows are the individuals)
getUniqueQtl <- function(pop, traits=c(1,2)) {
  # Get the qtl from trait 1
  qtlGeno <- pullQtlGeno(pop, traits[1])
  # Get the qtl from the other specified traits
  if (length(traits) > 1) {
    for (t in 2:length(traits)) {
      qtlGeno <- cbind(qtlGeno, pullQtlGeno(pop, traits[t]))
    }
  }
  qtlGeno <- as.data.frame(qtlGeno)
  # Remove duplicate QTL and order by QTL name
  qtlGeno <- qtlGeno %>%
    dplyr::select(which(!duplicated(names(.)))) %>%
    dplyr::select(sort(names(.)))

  # Set the rownames to be 1:nInd
  rownames(qtlGeno) <- c(1:nrow(qtlGeno))
  return (qtlGeno)
}


# Gets the absolute value of the effect size for each qtl
# Assumes there are 2 traits with QTL of additive effects
# Calculates effect size as sqrt(e1^2 + e2^2) (where eN is effect size of trait N)
# pop: The population
# Returns a dataframe with the qtl indices as the rows, and one column: 'eff_size'
getQtlEffectSizes <- function(pop) {
  # For each trait:
  # Merge the QTL dataframes with the dataframe of the effect sizes, and filter
  # out rows 1:pop@nInd, which are the genotypes, which we don't care about
  e1 <- data.frame(rbind(pullQtlGeno(pop, 1),
                         SP$traits[[1]]@addEff)[pop@nInd+1,])
  e2 <- data.frame(rbind(pullQtlGeno(pop, 2),
                         SP$traits[[2]]@addEff)[pop@nInd+1,])

  # Outer join the two dataframes
  eff_sizes <- merge(e1, e2, by="row.names", all=TRUE)
  # Change all NA values to zero
  eff_sizes[is.na(eff_sizes)] <- 0
  colnames(eff_sizes) <- c("snp", "eff1", "eff2")
  # Create a 3rd column called eff_size which is sqrt(e1^2 + e2^2) and move 1st
  # column back to the index spot
  eff_sizes <- eff_sizes %>%
    mutate(eff_sizes, eff_size=sqrt(eff1^2 + eff2^2)) %>%
    column_to_rownames(., var="snp")
  # Only store the generated effect sizes from both traits
  eff_sizes <- eff_sizes[3]
  return(eff_sizes)
}

# Get the QTL effect sizes for the desired trait
getSingleTraitEffectSizes <- function(pop, trait) {
  eff_sizes <- data.frame(rbind(pullQtlGeno(pop, trait),
                         SP$traits[[trait]]@addEff)[pop@nInd+1,])
  # Change all NA values to zero
  eff_sizes[is.na(eff_sizes)] <- 0
  colnames(eff_sizes) <- c("eff_size")
  
  return(eff_sizes)
}

# Calculates the fitness effect size of the QTL
# qtl: the name of the QTL in the form chr_loc
# merged.df: a dataframe with nInd rows and nQtl + nTraits columns (including a 'fitness' column)
# Returns: an effect size > 0, or NA if the locus is not segregating in the
# population, or it isn't a QTL
getFitnessEffectSize <- function(merged.df, qtl) {
  # Calculate the mean phenotype
  mu <- mean(merged.df$fitness)
  
  # Aggregate the dataframe to calculate the genotypic values based on deviation from the mean
  means <- merged.df %>%
    dplyr::group_by(!!sym(qtl)) %>%
    summarize(g=mean(fitness)-mu,
              n=n(),
              .groups = 'drop') %>%
    dplyr::mutate(freq=n/sum(n)) %>%
    complete(
      !!sym(qtl) := 0:2,
      fill = list()
    ) %>%
    dplyr::mutate(across(everything(), ~replace_na(.x, 0))) %>%
    column_to_rownames(var = {{qtl}})
  
  # Calculate the coded genotypic values
  a_d <- coded_gvs(means["2","g"], means["1","g"], means["0","g"])
  return (a_d[[1]])
}

# Calculates the trait architecture based on the genotypic values of individuals
# in the population, i.e. the mean difference in fitness between individuals
# homozygous for each genotype of each QTL
# pop: An AlphaSimR population
# Returns a dataframe with the following columns: 'id', and 'eff_size'
getFitnessTraitArchitecture <- function (pop) {
  eff_sizes <- data.frame(id = character(nLoci),
                          eff_size=numeric(nLoci))
  
  # Obtain the genotype data
  qtl.df <- as.data.frame(getUniqueQtl(pop))
  
  # Remove all non-segregating loci
  qtl.df <- qtl.df[, sapply(qtl.df, segLocus)]
  
  # If there are no segregating loci, there is no genetic variance, so return the default
  if(ncol(qtl.df) == 0) {
    return (eff_sizes)
  }
  
  # Get names of qtl
  qtls <- colnames(geno)
  nQtl <- length(qtls)
  popSize <- nrow(qtl.df)
  
  # Get the phenotypes for each of the traits
  pheno.df <- as.data.frame(pheno(pop))
  
  # Perform an inner join on the qtl (genotype) data and the phenotype data
  merged.df <- qtl.df %>%
    tibble::rownames_to_column(var = "row_name") %>%
    inner_join(
      pheno.df %>% tibble::rownames_to_column(var = "row_name"),
      by = "row_name"
    ) %>%
    tibble::column_to_rownames(var = "row_name") %>%
    arrange(as.integer(row.names(.))) %>%
    dplyr::mutate(fitness = fitCalc(Trait1, Trait2))
  
  # Calculate the effect size for each allele
  for (l in 1:nQtl) {
    id <- qtls[l]
    eff_sizes$id[l] <- id
    # locus is all of the genotypes in the population at a given locus
    locus = geno[,l]
    eff_sizes$eff_size[l] <- getFitnessEffectSize(merged.df, id)
  }
  # Filter out all effect sizes of zero
  eff_sizes <- eff_sizes[apply(eff_sizes!=0, 1, all),]
  # Filter out all NA effect sizes
  eff_sizes <- na.omit(eff_sizes)
  return (eff_sizes)
}

# Determines the genetic architecture of a trait in a population
# pop: the population of interest
# trait: the index of the trait to get an effect size for (only matters if method
# type is 'Additive')
# Returns a dataframe with the following columns: 'id', and 'eff_size'
singleTraitArchitecture <- function(pop, trait) {
  # Determine the segregating alleles
  qtlGeno <- pullQtlGeno(pop, trait)
  segLoci <- as.data.frame(apply(qtlGeno, MARGIN=2, FUN=segLocus))
  colnames(segLoci) <- "seg"
  # Join the seg table with the effect size table
  
  eff_sizes <- data.frame(rbind(qtlGeno,
                                SP$traits[[trait]]@addEff)[pop@nInd+1,])
  
  eff_sizes <- merge(segLoci, eff_sizes, by="row.names")
  colnames(eff_sizes) <- c("id", "seg", "eff_size")
  
  # Filter out any non-segregating alleles
  eff_sizes <- eff_sizes %>%
    filter(seg) %>%
    dplyr::select(-seg) %>%
    dplyr::mutate(eff_size=abs(eff_size)) %>%
    dplyr::arrange(desc(eff_size))
  return (eff_sizes)
}


# Determines the genetic architecture of a trait in a population
# pop: the population of interest
# methodType: one of 'Additive' (for determining the additive effect) or
# 'Fitness' (for determining the effect of the allele on fitness)
# trait: the index of the trait to get an effect size for (only matters if method
# type is 'Additive')
# Returns a dataframe with the following columns: 'id', and 'eff_size'
twoTraitArchitecture <- function(pop, methodType="Additive") {
  # Determine the segregating alleles
  segLoci <- as.data.frame(apply(getUniqueQtl(pop, traits=c(1,2)), MARGIN=2, FUN=segLocus))
  colnames(segLoci) <- "seg"
  # Join the seg table with the effect size table
  eff_sizes <- merge(segLoci, getQtlEffectSizes(pop),by="row.names")
  colnames(eff_sizes)[1] <- "id"
  # Filter out any non-segregating alleles
  eff_sizes <- eff_sizes %>% filter(seg)
  eff_sizes <- eff_sizes[c(1,3)]
  if (methodType == "Additive") {
    return (eff_sizes)
  } else {
    eff_sizes <- eff_sizes %>%
      dplyr::mutate(eff_size=eff_size^2)
    return (eff_sizes)
  }
}

# This function returns a dataframe of the effect sizes of segregating QTL 
# in the population, sorted in descending order of effect size,
# with the columns: id, eff_size, and rank
sortedEffectSizes <- function(pop, traits=c(1,2), methodType="Additive") {
  # Get the effect sizes
  if (length(traits) == 2) {
    effSizes <- twoTraitArchitecture(pop, methodType)
  } else if (length(traits) == 1) {
    effSizes <- singleTraitArchitecture(pop, trait=traits[[1]])
  }
  # If there are no segregating QTL, return the empty dataframe
  if (nrow(effSizes) == 0) {
    return(effSizes)
  }
  # Put in descending order
  effSizes <- effSizes %>%
    arrange(desc(eff_size))
  # Add a new column called rank which goes from 1 to the number of rows
  effSizes$rank <- seq.int(nrow(effSizes))
  return (effSizes)
}

# This function will return the additive effect sizes for QTL and the genetic location in morgans
# pop: The population in question
# traits: a vector specifying the indices of traits to acquire (plus fitness)
# Returns: a dataframe with columns:
#   "snp" (the qtl id)
#   "trait" (the trait - Trait 1, Trait 2, Trait 3, Fitness)
#   "eff_size" (the effect size)
#   "pos" (the genetic location in morgans)
#   "chr", the chromosome
getPerTraitQtlEffectSizesAndLocations <- function(pop, traits=c(1,2)) {
  # Get all segregating loci
  segLoci <- as.data.frame(apply(getUniqueQtl(pop, traits), MARGIN=2, FUN=segLocus)) %>%
    rename_with(~ "seg", .cols=1) %>%
    rownames_to_column("snp")
  
  # Get the additive effect sizes for all the QTL
  e1 <- data.frame(eff_trait1=SP$traits[[1]]@addEff)
  rownames(e1) <- colnames(pullQtlGeno(pop,1))
  
  e2 <- data.frame(eff_trait2=SP$traits[[2]]@addEff)
  rownames(e2) <- colnames(pullQtlGeno(pop,2))
  
  if (length(traits) > 2) {
    e3 <- data.frame(eff_trait3=SP$traits[[3]]@addEff)
    rownames(e3) <- colnames(pullQtlGeno(pop,3))
  }
  
  # Outer join the two dataframes
  eff_sizes <- merge(e1, e2, by="row.names", all=TRUE) %>%
    column_to_rownames(var="Row.names")
  
  if (length(traits) > 2) {
    eff_sizes <- merge(eff_sizes, e3, by="row.names", all=TRUE)
  }
  
  # Change all NA values to zero
  eff_sizes[is.na(eff_sizes)] <- 0
  
  eff_sizes <- rownames_to_column(eff_sizes, "snp")
  
  # TODO figure out how to add fitness QTL
  # Create a 3rd column called Fitness which is e1^2 + e2^2
  eff_sizes <- eff_sizes %>%
    #mutate(eff_sizes, eff_fitness=fitCalc(eff_trait1, eff_trait2)) %>%
    pivot_longer(
      cols = starts_with("eff"),
      names_to = "trait",
      values_to = "eff_size",
    ) %>% mutate(trait = case_when(
      trait == "eff_trait1" ~ "Trait 1",
      trait == "eff_trait2" ~ "Trait 2",
      trait == "eff_trait3" ~ "Trait 3",
      trait == "eff_fitness" ~ "Fitness",
      TRUE ~ trait
    ))
  
  # Drop all non-segregating markers and make effect sizes positive
  eff_sizes <- eff_sizes %>%
    left_join(segLoci, select(marker, seg), by="snp") %>%
    dplyr::filter(seg==TRUE) %>%
    dplyr::select(-seg) %>%
    dplyr::mutate(eff_size=abs(eff_size)) %>%
    dplyr::filter(eff_size>0)
  
  # Get the genetic map, and flatten it so it isn't grouped by chromosome
  genMap <- as.data.frame(unlist(SP$genMap)) %>%
    rename_with(~ "pos", .cols=1) %>%
    rownames_to_column("snp") %>%
    mutate(snp = sub(".*\\.", "", snp)) # Get ids in the form chr_loc
  
  # Merge effect sizes and locations on the snp id
  eff_sizes <- eff_sizes %>%
    left_join(genMap, select(snp, pos), by="snp") %>%
    mutate(chr = as.numeric(sub("_.*", "", snp))) # Determine the chromosome of each qtl
  return(eff_sizes)
}

# This function will create a ggplot of the trait architecture
# pop: The population to calculate trait architecture for
# methodType: to use for the calculation of effect size. Either 'Additive', for
# the underlying additive effect size, or 'Fitness', to determine the effect of
# each allele on fitness
# trait: the index of the trait under question
# Returns: a ggplot
plotTraitArchitecture <- function(pop, traits=c(1,2), popName="", methodType="Additive") {
  if (length(traits) == 1) {
    eff_sizes <- singleTraitArchitecture(pop, traits[1])
  } else if (length(traits) == 2) {
    eff_sizes <- twoTraitArchitecture(pop)
  }
  
  # Create a plot with the effect sizes ranked in descending order
  g <- ggplot(data=eff_sizes, aes(x=reorder(id, -eff_size), y=eff_size)) +
    geom_bar(stat="identity") +
    labs(x = "Variant Id", y = "Effect Size", title=popName) +
    theme(axis.text.x = element_text(angle = 45,
                                     vjust = 0.5,
                                     hjust=1,
                                     size = 6,
                                     margin = margin(b = 10)),
          axis.text.y = element_text(margin = margin(l=10, r=10)))
  return (g)
}


