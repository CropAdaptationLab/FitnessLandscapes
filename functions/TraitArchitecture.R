# Title: TRAIT ARCHITECTURE
# Author: Ted Monyak
# Description: Contains functions for examining and plotting trait architectures

# Determine whether a population is segregating at a locus
# locus: a list of genotypes from a population
# Returns: true if a locus is segregating within a population
segLocus <- function(locus) {
  return (length(unique(locus)) > 1)
}
# Get the QTL of all specified traits
# Since QTL can be overlapping between traits, we must filter out the
# duplicated QTL
# pop: an AlphaSimR population
# traits: a vector specifying the indices traits (defaults to both attained traits)
# Returns: a dataframe of the QTL (where columns are the QTL and rows are the individuals)
getUniqueQtl <- function(pop, traits=c(1,2)) {
  # Get the qtl from trait 1
  qtlGeno <- pullQtlGeno(pop, traits[1])
  # Get the qtl from the other traits
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
    dplyr::arrange(desc(eff_size)) %>%
    rownames_to_column("rank")
  return (eff_sizes)
}

# This function will return the additive effect sizes for QTL and the genetic location in centimorgans
# pop: The population in question
# traits: a vector specifying the indices of traits to acquire (plus fitness)
# Returns: a dataframe with columns:
#   "snp" (the qtl id)
#   "trait" (the trait - Trait 1, Trait 2, Trait 3, Fitness)
#   "eff_size" (the effect size)
#   "pos" (the genetic location in morgans)
#   "chr", the chromosome
getQtlEffectSizes <- function(pop, traits=c(1,2)) {
  # Get all segregating loci
  segLoci <- as.data.frame(apply(getUniqueQtl(pop, traits), MARGIN=2, FUN=segLocus)) %>%
    rename_with(~ "seg", .cols=1) %>%
    rownames_to_column("id")
  
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
  
  eff_sizes <- rownames_to_column(eff_sizes, "id")
  eff_sizes <- eff_sizes %>%
    pivot_longer(
      cols = starts_with("eff"),
      names_to = "trait",
      values_to = "eff_size",
    ) %>% mutate(trait = case_when(
      trait == "eff_trait1" ~ "Trait 1",
      trait == "eff_trait2" ~ "Trait 2",
      trait == "eff_trait3" ~ "Trait 3",
      TRUE ~ trait
    ))
  
  # Drop all non-segregating markers and make effect sizes positive
  eff_sizes <- eff_sizes %>%
    left_join(segLoci, select(marker, id), by="id") %>%
    dplyr::filter(seg==TRUE) %>%
    dplyr::select(-seg) %>%
    dplyr::mutate(eff_size=abs(eff_size)) %>%
    dplyr::filter(eff_size>0)
  
  # Get the genetic map, and flatten it so it isn't grouped by chromosome
  genMap <- as.data.frame(unlist(SP$genMap)) %>%
    rename_with(~ "pos", .cols=1) %>%
    rownames_to_column("id") %>%
    mutate(id = sub(".*\\.", "", id)) # Get ids in the form chr_loc
  
  # Merge effect sizes and locations on the snp id
  eff_sizes <- eff_sizes %>%
    left_join(genMap, select(id, pos), by="id") %>%
    mutate(chr = as.numeric(sub("_.*", "", id))) %>% # Determine the chromosome of each qtl
    dplyr::mutate(pos = pos*n.genMapLen)
  return(eff_sizes)
}

# This function will create a ggplot of the trait architecture
# pop: The population to calculate trait architecture for
# trait: the index of the trait to plot
# Returns: a barplot of the effect sizes in descending order
plotTraitArchitecture <- function(pop, trait=1, popName="") {
  eff_sizes <- singleTraitArchitecture(pop, trait)
  
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


