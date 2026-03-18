# Title: POP GEN
# Author: Ted Monyak
# Description: This file contains functions for calculating population genetics
# statistics, such as Fst

# This function calculates the mean FST (fixation index) for a metapopulation
# FST is calculated as (Ht-mean(Hs))/Ht, where
# Ht is the expected heterozygosity in the metapopulation (at each locus)
# Hs is the expected heterozygosity in a subpopulation
# Pops is a list of the populations in the metapopulation
# Returns the average FST value over all of the markers
FST <- function(pops) {
  # This function calcluates the expected heterozygosity at a locus with the 
  # equation 2pq, where p and q are the allele frequencies of the '2' and '0'
  # alleles, respectively
  calcExpHet <- function(locus, pop) {
    popSize <- nInd(pop)
    # The frequency of the p allele is the frequency of the '22' homozygotes plus
    # half of the frequency of the heterozygotes
    p <- (sum(locus==2)/popSize) + ((sum(locus==1)/popSize)/2)
    q <- 1 - p
    return (2 * p * q)
  }
  
  # Calculate the expected heterozygosity in the whole metapopulation
  # Returns a vector of size n.markers*1
  calcHt <- function(pops) {
    # Merge all the subpopulations into one metapopulation
    metaPop <- pops[[1]]
    for (p in 2:length(pops)) {
      metaPop <- c(metaPop, pops[[p]])
    }
    # Get the marker genotype data
    metaPopGeno <- pullSnpGeno(metaPop)
    # Calculate per-locus expected heterozygosity by applying calcExpHet to each
    # column in the dataframe
    ht <- apply(metaPopGeno, MARGIN=2, FUN=calcExpHet,pop=metaPop)
    return (ht)
  }
  
  # Calculate the mean expected heterozygosity across each subpopulation
  # Returns a vector of size n.markers*1
  calcMeanHs <- function(pops) {
    # Create a result dataframe
    all_hs <- data.frame(matrix(ncol=length(pops), nrow=ncol(pullSnpGeno(pops[[1]]))))
    # Set the column names as 1, 2, etc.
    colnames(all_hs) <- c(1:length(pops))
    # Iterate through each subpopulation
    for (p in 1:length(pops)) {
      pop <- pops[[p]]
      # Calculate per-locus heterozygosity for each column in the genotype dataframe
      hs <- apply(pullSnpGeno(pop), MARGIN=2, FUN=calcExpHet, pop=pop)
      # Update the result dataframe
      all_hs[,p] <- hs
    }
    # Return a vector of the mean of each row (locus) in the result dataframe
    return (rowMeans(all_hs))
  }
  
  # Calculate Ht
  ht <- as.data.frame(calcHt(pops))
  # Calculate mean Hs
  meanHs <- as.data.frame(calcMeanHs(pops))
  # Merge the dataframes
  combined <- cbind(ht, meanHs)
  colnames(combined) <- c("HT", "HS")
  # Calculate the Fst at each 
  combined <- combined %>% mutate(FST=((HT-HS)/HT)) %>% replace(is.na(.), 0)
  return (mean(combined$FST))
}