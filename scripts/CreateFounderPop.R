# Title: CREATE FOUNDER POPULATION
# Author: Ted Monyak
# Description: This script creates a founder population to use in subsequent simulations

# For reproducibility
set.seed(123)
# Create initial population
if (basicPop) {
  # default founder parameters
  founders = runMacs(
    nInd=n.popSize,
    nChr=n.chr,
    segSites=n.segSites,
    nThreads=n.cores
  )
} else {
  # NOT BEING USED
  # This is from a Brian Rice script
  founders = runMacs(nInd=n.ne,
                     nChr=10,
                     segSites=n.segSites,
                     inbred=TRUE,
                     manualCommand = paste(
                       "1000000000 -t", #Physical length 1e8 base pairs
                       2.5/1E8*(4*n.ne), #Mutation rate adjusted for Ne
                       "-r",1/1E8*(4*n.ne), #Recombination rate adjusted for Ne
                       "-eN",10/(4*n.ne),100/n.ne), #Modeling Ne=100 at 10 generations ago
                     manualGenLen=1, nThreads=n.cores) #Genetic length 1 Morgan
}

SP <- SimParam$new(founders)

# Allow markers and segmentation sites to overlap
SP$restrSegSites(overlap = T)

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
    # Select n.qtlPerChr markers per chromosome
    for (chr in 1:n.chr) {
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
      rand_qtls <- sample(eligible_qtls, n.qtlPerChr)
      
      # Sort by location
      locs <- sort(sapply(rand_qtls, function(x) as.numeric(sub(".*_", "", x))))
      # Re-add the chromosome to the marker name
      chr_qtls <- sapply(locs, function(x) paste0(as.character(chr), "_", as.character(x)))
      # Add the qtls from this chromsoome to the list
      qtls <- c(qtls, chr_qtls)
    }
    return(unname(qtls))
  }
  
  geom_series <- function(a, r, n) {
    res <- a * r ** c(0:(n-1))
    return (res)
  }
  sample(c(0,1), size=1) == 1
  
  markerNames <- selectQtl()
  addEff <- geom_series(a=1,r=n.a,n=n.L)
  # Randomly assign some effects as negative
  addEff <- sapply(addEff, function(x) if(sample(c(0,1), size=1) == 1) -x else x)
  
  #varScale <- n.var#*(1-n.R)
  #addEff <- addEff * varScale
  addEff <- sample(addEff)
  SP$importTrait(markerNames, addEff)
  
}

n.L <- n.qtlPerChr * n.chr
n.a <- (n.L-1) / (n.L+1)

if (n.relAA == 0) { # Additive G > T
  if (SAMPLING=="normal") {
    SP$addTraitA(mean=n.initTraitVal, var=n.var, nQtlPerChr=n.qtlPerChr)
    SP$addTraitA(mean=n.initTraitVal, var=n.var, nQtlPerChr=n.qtlPerChr)
  } else if (SAMPLING=="gamma") {
    # Higher shape parameter means more oligogenic
    SP$addTraitA(mean=n.initTraitVal, var=n.var, nQtlPerChr=n.qtlPerChr, gamma=TRUE, shape=n.shape)
    SP$addTraitA(mean=n.initTraitVal, var=n.var, nQtlPerChr=n.qtlPerChr, gamma=TRUE, shape=n.shape)
  } else if (SAMPLING=="geometric") {
    addGeometricAdditiveTrait()
    addGeometricAdditiveTrait()
    SP$rescaleTraits(mean=c(n.initTraitVal,n.initTraitVal),var=c(n.var,n.var),varEnv = c(0,0),varGxE = c(0,0))
  }
  
  # Add a yield potential polygenic trait
  SP$addTraitA(mean=n.initYieldVal, var=n.yieldVar, nQtlPerChr=n.yieldQtlPerChr)
} else { # Epistatic G > T
  SP$addTraitAE(mean=n.initTraitVal, var=n.var, nQtlPerChr=n.qtlPerChr, relAA=n.relAA, useVarA=FALSE)
  SP$addTraitAE(mean=n.initTraitVal, var=n.var, nQtlPerChr=n.qtlPerChr, relAA=n.relAA, useVarA=FALSE)
  
  # Add a yield potential polygenic trait
  SP$addTraitAE(mean=n.initYieldVal, var=n.yieldVar, nQtlPerChr=n.yieldQtlPerChr, relAA=n.relAA, useVarA=FALSE)
}

# Set heritability for each of the traits
SP$setVarE(h2=c(n.h2, n.h2, n.yieldH2))

# Add a SNP chip with n.markers*n.chr markers
SP$addSnpChip(nSnpPerChr=n.markers)

# Create base population
founderPop <- newPop(founders, simParam = SP)

if (saveFitnessPlots) {
  plotTraitArchitecture(founderPop, traits=c(1,2), popName=paste0("QTL: ", n.L))
  ggplot2::ggsave(filename = paste0("founderarchitecture.jpg"),
                  path=pop_dir,
                  device = "jpg",
                  width=10,
                  height=7)
}
#alpha_1 <- singleTraitArchitecture(founderPop, 2)[1,1]

#singleTraitArchitecture(RIL, 1) %>%
#  dplyr::mutate(scaled=(eff_size/0.015)*(1-0.9)) %>%
#  summarize(sum=sum(scaled))

#varG(RIL)
#sum(eff_sizes$scaled)
