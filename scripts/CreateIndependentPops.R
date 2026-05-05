# Title: CREATE INDEPENDENT POPULATIONS
# Author: Ted Monyak
# Description: This script creates n.nPops independent sub-populations from an initial founder
# population and has them follow independent adaptive walks to a fitness optimum

# Assumes that CreateFounderPop.R has been run already
# Each replication should create a new subpop_dir, where this data is stored

# founderPop is created in CreateFounderPop.R
pop <- founderPop
# Check if # independent pops * # individuals per pop does not exceed the
# number of individuals in the founder pop
if (n.nPops*n.subPopSize > nInd(pop)) {
  stop(paste0("Population of size ", nInd(pop), " not large enough to create ",
              n.nPops, " subpopulations of size ", n.subPopSize, "."))
}

# Burn-in
for (gen in 1:n.burnInGens) {
  pop <- selectCross(pop, use="rand", nInd=nInd(pop), nCrosses=nInd(pop))
}

# Create a random vector of size n.pops, with a random order of sub-population ids
randVec <- sample(rep(c(1:n.nPops), times=n.popSize/n.nPops))

# Create n.nPops sub populations
pops <- vector(mode="list", length=n.nPops)
# Lists for storing results
fit_dfs <- list()
subpop_dirs <- list()

# Sample individuals for the genotype-to-fitness landscape
#sampled_inds <- matrix(data=pullSnpGeno(founderPop, snpChip=2),
#                       nrow=nInd(founderPop),
#                       ncol=(n.GSmarkers*n.chr))

sampled_inds <- founderPop

for (p in 1:n.nPops) {
  # Assign subpop 'p'
  pops[[p]] <- selectInd(pop, trait=selectSubPop, selectTop=TRUE, nInd=n.subPopSize, idx=p, randVec=randVec)
  
  # Get the names of all the QTLs
  qtl <- colnames(getUniqueQtl(pops[[p]]))
  
  # Create a dataframe of all zeros where the columns are the QTL ids, and the # rows is the # of generations
  # Add six to account for the number of other datapoints being stored (colnames)
  fit.df <- data.frame(matrix(ncol=length(qtl)+6, nrow=0))
  colnames(fit.df) <- c("gen",
                        "traitVal1",
                        "traitVal2",
                        "suit",
                        "meanSuit",
                        "yieldPotential",
                        qtl)
  fit_dfs[[p]] <- fit.df
  
  # Create subpop dirs
  subpop_dir <- file.path(rep_dir, paste0("Subpop_", p))
  if (!dir.exists(subpop_dir)) dir.create(subpop_dir)
  subpop_dirs[[p]] <- subpop_dir
}

# Stores the order in which an allele is fixed along an adaptive walk
# Row 1 stores trait 1 and row 2 stores trait 2
# Each column is a subpopulation
fixation_idx <- matrix(1, nrow=2, ncol=n.nPops)

# Get the effect sizes of each qtl
qtlEff.df <- getQtlEffectSizes(pop)

# Each population follows an adaptive walk for a maximum of n.gens generations
# Each will terminate once it is within n.margin of the fitness optimum
for (gen in 1:n.gens) {
  for (p in 1:length(pops)) {
    pop <- pops[[p]]
    
    # Get mean phenotypes for this generation
    pheno <- as.data.frame(pheno(pop)) %>%
      dplyr::mutate(suit=suitFunc(Trait1, Trait2)) %>%
      dplyr::mutate(id=pop@id)
    meanSuit <- mean(pheno$suit)
    traitVal1 <- mean(pheno$Trait1)
    traitVal2 <- mean(pheno$Trait2)
    yieldPotential <- mean(pheno$Trait3)
    suit <- suitFunc(traitVal1, traitVal2)
    
    # Get the qtl genotype data
    qtlGeno <- getUniqueQtl(pop)
    
    # The new data to add
    newRow <- data.frame(
      gen=gen,
      traitVal1=traitVal1,
      traitVal2=traitVal2,
      suit=suit,
      meanSuit=meanSuit,
      yieldPotential=yieldPotential)
    
    # Sample individuals closest to the mean for each of traits 1 and 2
    if (sampleInds) {
      #sample <- selectInd(pop, use="rand", nInd=n.selInds)
      #sampled_inds <- rbind(sampled_inds,pullSnpGeno(sample, snpChip=2))
      
      # The midpoint of the number of individuals in the population
      mid <- nInd(pop)/2
      
      # Sort phenotypes by attained trait 1
      trait1 <- pheno %>%
        arrange(Trait1)
      # Get the middle individuals w.r.t. trait 1
      middle_t1 <- trait1[(mid-(window/2)):(mid+(window/2)),]  
      trait2 <- pheno %>%
        arrange(Trait2)
      # Get the middle individuals w.r.t. trait 2
      middle_t2 <- trait2[(mid-(window/2)):(mid+(window/2)),]
      
      # Find the individuals with middle phenotypes for both of the traits
      common_ids <- intersect(middle_t1$id,
                              middle_t2$id)
      
      # Get the phenotypes of the selected individuals
      sampleNum <- min(length(common_ids), n.selInds)
      sampleIds <- sample(pheno$id, size=sampleNum, replace=FALSE)
      sampled_inds <- c(sampled_inds, pop[sampleIds])
    }
    
    # Select n.m percentage of individuals to migrate in a stepping stone model
    if (n.m > 0) {
      # Number of individuals that are migrating
      migrationSize <- n.subPopSize * n.m
      
      # Randomly designate individuals with a 1 (staying in the population), 2 (moving left), or 3 (moving right)
      randVec <- sample(c(rep(1, nInd(pop) - migrationSize), rep(2, migrationSize/2), rep(3, migrationSize/2)))
      migrateLeft <- selectInd(pop, trait=selectSubPop, selectTop=TRUE, nInd=migrationSize/2, idx=2, randVec=randVec)
      migrateRight <- selectInd(pop, trait=selectSubPop, selectTop=TRUE, nInd=migrationSize/2, idx=3, randVec=randVec)
      pop <- selectInd(pop, trait=selectSubPop, selectTop=TRUE, nInd=nInd(pop) - migrationSize, idx=1, randVec=randVec)
      
      
      if (p > 1) { # If this is not the first population, migrate to the population to the left
        pops[[p-1]] <- c(pops[[p-1]], migrateLeft)
      } else { # If this is the first population, migrate to the last population
        pops[[n.nPops]] <- c(pops[[n.nPops]], migrateLeft)
      }
      
      if (p < length(pops)) { # If this is not the last population, migrate to the population to the right
        pops[[p+1]] <- c(pops[[p+1]], migrateRight)
      } else { # If this is the last population, migrate to the first population
        pops[[1]] <- c(pops[[1]], migrateRight)
      }
    }
    
    # At each stage, select the top individuals according to how close each 
    # is from the fitness optimum, according to selectionPheno
    pop <- selectCross(pop, trait=selectionPheno, nInd=nInd(pop)*n.selProp, nCrosses=nInd(pop))
  
    if (saveAllelePlots | saveFixationOrder) {
      # The allele frequencies of the previous generation
      prevAlleleFreq <- as.data.frame(apply(qtlGeno,
                                            MARGIN=2,
                                            FUN=function(x)
                                              (sum(x==n.allele)/n.subPopSize) + ((sum(x==1)/n.subPopSize)/2)))
      newRow <- cbind(newRow, t(prevAlleleFreq))
      
      # Join the allele frequencies with the effect sizes
      prevAlleleFreq <- as.data.frame(prevAlleleFreq) %>%
        dplyr::rename("prevFreq"=1) %>%
        tibble::rownames_to_column(var="id") %>%
        dplyr::left_join(qtlEff.df, by="id")
      
      # Determine which of the new alleles were fixed
      # newQtlGeno: The QTL marker data in this generation
      # trait: The AlphaSim phenotype index
      # Returns: a dataframe of alleles that were fixed this generation
      getNewFixedAlleles <- function(newQtlGeno, trait) {
        # Get the number of alleles that have been fixed
        fixation_order <- fixation_idx[trait, p]
        fixedAlleles <- as.data.frame(apply(newQtlGeno,
                                            MARGIN=2,
                                            FUN=function(x)
                                              (sum(x==n.allele)/n.subPopSize) + ((sum(x==1)/n.subPopSize)/2))) %>%
          dplyr::rename("newFreq"=1) %>%
          tibble::rownames_to_column(var="id") %>%
          dplyr::inner_join(prevAlleleFreq, by="id") %>% # Join new allele frequencies with previous
          dplyr::mutate(fixed=((newFreq==1 | newFreq==0) & newFreq != prevFreq)) %>% # Determine which were fixed
          dplyr::filter(fixed==TRUE) %>%
          dplyr::mutate(rank=fixation_order,
                        subpop=p,
                        qtl=n.L) %>%
          dplyr::select(c(id,eff_size,rank,subpop, qtl))
        
        return(fixedAlleles)
      }
        
      
      # Determine which trait1 qtl were fixed this generation
      trait <- 1
      fixed_t1 <- getNewFixedAlleles(pullQtlGeno(pop, trait), trait)
      # The fixation index is incremented by the number of QTL fixed for that trait this generation
      fixation_idx[trait,p] <- fixation_idx[trait,p] + nrow(fixed_t1)
      
      # Determine which trait2 qtl were fixed this generation
      trait <- 2
      fixed_t2 <- getNewFixedAlleles(pullQtlGeno(pop, trait), trait)
      fixation_idx[trait,p] <- fixation_idx[trait,p] + nrow(fixed_t2)
      
      fixedAlleles.df <- rbind(fixedAlleles.df, fixed_t1, fixed_t2)
    }
    # Add the new data
    fit_dfs[[p]] <- rbind(fit_dfs[[p]], newRow)
    
    pops[[p]] <- pop

  } # end pops
} # end n.gens
