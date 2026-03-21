# Title: GENO CONVERSIONS
# Author: Ted Monyak
# Description: Contains functions for converting AlphaSimR genotype encodings
# into encodings for rqtl and rrblup

# Create a genotype encoded using the ABH encoding for use in the rqtl package
# Uses the following documentation as reference:
# ?read.cross
# https://bitbucket.org/tasseladmin/tassel-5-source/wiki/UserManual/GenosToABH/GenosToABHPlugin
# Assumes a RIL family has been created
# RIL: the RIL family to get the genotype encoding
# parentA: one parent that was crossed to create the biparental RIL
# parentB: the other parent that was crossed to create the biparental RIL
# genoEnc: a list of encodings where the first element is missing (NA), the second element
# is homozygous for parent A (AA), the third element is heterozygous (AB), and the
# fourth element is homozygous for parent B (BB)
# Returns: a dataframe of i individuals x m snps
getABHGeno <- function(RIL, parentA, parentB, genoEnc=c(NA,1,NA,2), snpChip=1) {
  # Get the marker genotype for parent A and B
  aSnps <- pullSnpGeno(parentA, snpChip) 
  bSnps <- pullSnpGeno(parentB, snpChip) 
  # Transpose the snp dataframes and merge them, so the dataframe has 1 column
  # for each parent, and a row for each snp
  allSnps <- as.data.frame(t(rbind(aSnps,bSnps)))
  colnames(allSnps) <- c('parentA', 'parentB')

  # Find the monomorphic snps (where parentA geno matches parentB geno)
  monoSnps <- allSnps %>% filter(parentA==parentB)
  monoMarkers <- rownames(monoSnps)

  # Find the snps where parentA or parentB is heterozygous
  hetSnps <- allSnps %>% filter(parentA==1|parentB==1)
  hetMarkers <- rownames(hetSnps)

  # Get the marker genotype for the RIL family
  rilSnps <- pullSnpGeno(RIL, snpChip)

  # Create the output dataframe
  output <- matrix(0, nrow=nrow(rilSnps), ncol=ncol(rilSnps))
  rownames(output) <- rownames(rilSnps)
  colnames(output) <- colnames(rilSnps)
  # Iterate through all snps
  for (m in 1:ncol(rilSnps)) {
    # If the snp is monomorphic in the parents, set it as 'missing'
    if (colnames(rilSnps)[m] %in% monoMarkers) {
      output[,m] <- genoEnc[1]
      next
    }
    # If the snp is heterozygous in one of the parents, set it as 'missing'
    if (colnames(rilSnps)[m] %in% hetMarkers) {
      output[,m] <- genoEnc[1]
      next
    }
    # Iterate through all individuals
    for (i in 1:nrow(rilSnps)) {
      rilGeno <- rilSnps[i,m]
      aGeno <- aSnps[1,m]
      bGeno <- bSnps[1,m]
      if (rilGeno == 1) {
        # if the RIL is heterozygous, set it as AB
        output[i,m] <- genoEnc[3]
      } else if (rilGeno == aGeno) {
        # if the RIL matches parentA, set it as AA
        output[i,m] <- genoEnc[2]
      } else if (rilGeno == bGeno) {
        # if the RIL matches parentB, set it as BB
        output[i,m] <- genoEnc[4]
      }
    }
  }
  return (output)
}

# Creates a Cross object to be used in rqtl
# pop: the population (assumed to be a biparental RIL)
# parentA: one parent crossed to create the RIL
# parentB: the other parent crossed to create the RIL
# popType: 'riself' (can also be 'f2', 'bc'). See ?read.cross and rqtl documentation
# genoEnc: a list of encodings where the first element is missing (NA), the second element
# is homozygous for parent A (AA), the third element is heterozygous (AB), and the
# fourth element is homozygous for parent B (BB)
# Returns an rqtl Cross object. See more at https://rqtl.org/tutorials/rqtltour.pdf
getCross <- function(pop, parentA, parentB, popType="riself", genoEnc=c(NA,1,NA,2), snpChip=1) {
  # Get the genetic distance map from AlphaSim
  genMap <- as.data.frame(unlist(SP$genMap))
  # Remove the chromosome prefix from the rownames, leaving them in the format chr_loc
  rownames(genMap) <- sub(".*\\.", "", rownames(genMap))
  colnames(genMap) <- "pos"

  # Get the ABH encoded genotypes of the RIL population (NA,1,2,3)
  snpGeno <- as.data.frame(t(getABHGeno(pop, parentA, parentB, genoEnc, snpChip)))

  # Join the genotype data with the genetic distance map
  snpGeno <- merge(snpGeno, genMap, by="row.names")
  
  # Add chromosome and site number to order by physical location
  snpGeno <- snpGeno %>%
    dplyr::mutate(chr = strtoi(sub("_.*", "", Row.names))) %>%
    dplyr::mutate(site = strtoi(sub(".*_", "", Row.names))) %>%
    arrange(chr, site) %>%
    dplyr::select(-one_of(c("chr", "site")))
  rownames(snpGeno) <- snpGeno$Row.names
  
  # Remove the first column of the dataframe (which is a duplicate of the row names)
  snpGeno <- snpGeno[c(-1)]
  # Transpose, to get the columns as the snps, and the rows as the individuals
  snpGeno <- as.data.frame(t(snpGeno))

  # Change this if using a different number of chromosomes
  chrs <- as.character(c(1:n.chr))

  # Pull the phenotype data from AlphaSim
  pheno <- data.frame(pheno=pop@pheno) %>%
    dplyr::mutate(pheno.Suitability=suitFunc(pheno.Trait1, pheno.Trait2),
                  pheno.W=calculateBreedingFitness(pheno.Trait1, pheno.Trait2, pheno.Trait3))
  
  geno <- list()
  # Iterate through each chromosome
  for (c in chrs) {
    # Get only the snps from chromosome 'c'
    snps <- snpGeno[,grepl(paste0(c,"_"), colnames(snpGeno))]
    # Remove the last row in the dataframe containing the genetic locations
    chr_geno <- data.matrix(snps[c(-length(rownames(snps))),])
    # Get the last row in the snps data frame with the genetic locations
    chr_snps <- snps["pos",]
    # Since the genetic locations are stored in morgans, convert them to
    # cM, as expected by rqtl
    chr_snps <- sapply(chr_snps, function(x) as.numeric(x)*100)
    geno[[c]] <- list()
    # A is an autosome
    attributes(geno[[c]])$class<-"A"
    # Conform to the rqtl cross formatting
    geno[[c]]$"data" <- chr_geno
    geno[[c]]$"map" <- chr_snps
  }
  cross <- list(geno=geno,pheno=pheno)
  class(cross) <- c(popType,"cross")
  return (cross)
}

# Function for creating a genotype encoding that conforms to the rrBLUP standard
# pop: the population to convert
# Returns: a dataframe with the following columns: 'snp' (the marker id), 'chr'
# (the chromosome), 'pos' (the genetic location), and a column for every individual
rrBlupGeno <- function(pop) {
  # Get marker data from the population
  snpGeno <- as.data.frame(pullSnpGeno(pop))
  # Filter out monomorphic snps
  snpGeno <- snpGeno[sapply(snpGeno, function(x) length(unique(x)) > 1)]

  # Obtain the genetic map and flatten it so it isn't grouped by chromosome
  genMap <- as.data.frame(unlist(SP$genMap))

  # Remove the prefix from each line, leaving the format chr_loc
  rownames(genMap) <- sub(".*\\.", "", rownames(genMap))

  # set the column names to be 'snp' and 'pos' (genetic location)
  colnames(genMap) <- "pos"
  genMap <- rownames_to_column(genMap, "snp")

  # Create the output dataframe with one row for each individual
  geno = data.frame(
    snp = colnames(snpGeno),
    chr = numeric(length(colnames(snpGeno))),
    t(snpGeno - 1)
  )
  # Set the column names of each individual to be a simple index from 1:popSize
  colnames(geno)[-c(1:2)] = 1:pop@nInd
  # Set the chromosome value to be the 'snp' id that precedes the '_'
  geno <- geno %>%
    mutate(chr = sub("_.*", "", snp))

  # Join the genotype data with the genetic map data
  geno <- merge(x=geno, y=genMap, by="snp")
  # Re-order the columns so that the 'pos' column becomes third
  geno <- geno[,c("snp", "chr", "pos", c(1:pop@nInd))]
  geno$chr <- as.integer(geno$chr)
  return (geno)
}


