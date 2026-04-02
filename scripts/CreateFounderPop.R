# Title: CREATE FOUNDER POPULATION
# Author: Ted Monyak
# Description: This script creates a founder population to use in subsequent simulations

# For reproducibility
set.seed(123)
# Create initial population
founders = runMacs(
  nInd=n.popSize,
  nChr=n.chr,
  segSites=n.segSites,
  nThreads=n.cores
)

SP <- SimParam$new(founders)

# Allow markers and segmentation sites to overlap
SP$restrSegSites(overlap = T)

# Assign the number of QTL per chromosome

# If n.L is divisible by n.chr, assign an equal number of QTL per chromosome
if (n.L %% n.chr == 0) {
  qtlPerChr <- rep(n.L/n.chr, times=n.chr)
} else if (n.L >= n.chr) {
  # Use 'stars and bars' to assign a random number of QTL to each chr
  cuts <- sort(sample(1:(n.L-1), n.chr-1, replace = FALSE))
  qtlPerChr <- diff(c(0, cuts, n.L))
} else {
  # There are fewer QTL than chromosomes
  qtlPerChr <- sample(c(rep(1,n.L), rep(0, n.chr-n.L)))
}

# Adds an additive trait with an explicitly-defined genetic architecture
# Effect sizes will follow a geometric distribution
# Effect sizes are scaled to achieve a desired genetic variance
addGeometricAdditiveTrait <- function() {

  # Randomly select markers to be QTL
  # Selected QTL are randomly chosen from the loci with a MAF above n.minMAF
  # QTL are sorted by position
  # Returns: a list of selected markers
  selectQtl <- function() {
    qtls <- c()
    # Select qtlPerChr markers per chromosome
    for (chr in 1:n.chr) {
      # The number of QTL on this chromosome
      nQtl <- qtlPerChr[chr]
      if (nQtl == 0) next
  
      # The initial genotypes for all markers
      geno <- pullSegSiteGeno(founders, chr)
      
      # Get the minor allele frequency of a locus
      getMAF <- function(locus) {
        # The frequency of the '1' allele
        allele_1_freq <- (sum(locus==n.allele)/nInd(founders)) +
          ((sum(locus==1)/nInd(founders))/2)
        # Return the MAF
        return (min(allele_1_freq, 1-allele_1_freq))
      }
      
      # Get MAFs for each locus
      mafs <- apply(geno, MARGIN=2, FUN=getMAF)
      # Retain only markers with a MAF above the threshold
      eligible_qtls <- names(mafs[mafs >= n.minMAF])
      # Select the random QTLs from these eligible markers
      rand_qtls <- sample(eligible_qtls, nQtl)
      
      # Sort by location
      locs <- sort(sapply(rand_qtls, function(x) as.numeric(sub(".*_", "", x))))
      # Re-add the chromosome to the marker name
      chr_qtls <- sapply(locs, function(x) paste0(as.character(chr), "_", as.character(x)))
      # Add the qtls from this chromsoome to the list
      qtls <- c(qtls, chr_qtls)
    }
    return(unname(qtls))
  }
  
  # Create a geometric series
  # a: Scaling term
  # r: Decay term
  # n: Length of the series
  # Returns: a vector of size n
  geom_series <- function(a, r, n) {
    res <- a * r ** c(0:(n-1))
    return (res)
  }
  sample(c(0,1), size=1) == 1
  
  # Select random QTL
  markerNames <- selectQtl()
  # Create the geometric series of effect sizes
  addEff <- geom_series(a=1,r=n.a,n=n.L)
  # Randomly assign some effects as negative
  addEff <- sapply(addEff, function(x) if(sample(c(0,1), size=1) == 1) -x else x)
  
  # Randomize the order of effect sizes
  addEff <- sample(addEff)
  SP$importTrait(markerNames, addEff)
  
}

# Calculate the decay parameter based on a formula from Lande and Thompson 1990
n.a <- (n.L-1) / (n.L+1)

# Additive G > T
if (n.relAA == 0) { 
  if (SAMPLING=="normal") {
    SP$addTraitA(mean=n.initTraitVal, var=n.var, nQtlPerChr=qtlPerChr)
    SP$addTraitA(mean=n.initTraitVal, var=n.var, nQtlPerChr=qtlPerChr)
  } else if (SAMPLING=="gamma") {
    SP$addTraitA(mean=n.initTraitVal, var=n.var, nQtlPerChr=qtlPerChr, gamma=TRUE, shape=n.shape)
    SP$addTraitA(mean=n.initTraitVal, var=n.var, nQtlPerChr=qtlPerChr, gamma=TRUE, shape=n.shape)
  } else if (SAMPLING=="geometric") {
    addGeometricAdditiveTrait()
    addGeometricAdditiveTrait()
    SP$rescaleTraits(mean=c(n.initTraitVal,n.initTraitVal),var=c(n.var,n.var),varEnv = c(0,0),varGxE = c(0,0))
  }
  
  # Add a yield potential polygenic trait
  SP$addTraitA(mean=n.initYieldVal, var=n.yieldVar, nQtlPerChr=n.yieldQtlPerChr)
} else { # Epistatic G > T
  SP$addTraitAE(mean=n.initTraitVal, var=n.var, nQtlPerChr=qtlPerChr, relAA=n.relAA, useVarA=FALSE)
  SP$addTraitAE(mean=n.initTraitVal, var=n.var, nQtlPerChr=qtlPerChr, relAA=n.relAA, useVarA=FALSE)
  
  # Add a yield potential polygenic trait
  SP$addTraitAE(mean=n.initYieldVal, var=n.yieldVar, nQtlPerChr=n.yieldQtlPerChr, relAA=n.relAA, useVarA=FALSE)
}

# Set heritability for each of the traits
SP$setVarE(h2=c(n.h2, n.h2, n.yieldH2))

# Add a SNP chip with n.markers*n.chr markers
SP$addSnpChip(nSnpPerChr=n.markers)
SP$addSnpChip(nSnpPerChr=n.GSmarkers)

# Create base population
founderPop <- newPop(founders, simParam = SP)

# Save the initial trait architecture
if (saveFitnessPlots) {
  plotTraitArchitecture(founderPop, trait=1, popName=paste0("QTL: ", n.L))
  ggplot2::ggsave(filename = paste0("founderarchitecture.jpg"),
                  path=pop_dir,
                  device = "jpg",
                  width=10,
                  height=7)
}
