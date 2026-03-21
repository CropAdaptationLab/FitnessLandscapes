# Title: QTL Mapping
# Author: Ted Monyak
# Description: Contains functions for doing linkage mapping

# Determines the number of significant LOD peaks in a linkage map
# Also plots the linkage map and a PCA plot of the RIL
# RIL: a recombinant inbred line family
# parent1: one of the parents crossed to create the RIL
# parent2: the other parent crossed to create the RIL
# save_dir: a directory in which to save the plots
# is homozygous for parent A (AA), the third element is heterozygous (AB), and the
# fourth element is homozygous for parent B (BB)
# Returns: a list of 3, with the number of significant LOD peaks per trait (Trait1, Trait2, Fitness)
getLodPeaks <- function(RIL, parent1, parent2, save_dir) {
  # Create a cross object used in rqtl
  cross <- getCross(RIL, parent1, parent2, "riself", snpChip = 1)
  cross <- drop.nullmarkers(cross)

  # Check to see if there are no markers once null markers have been dropped
  if (length(cross$geno) == 0) {
    return (0)
  }
  
  # TODO: Remove this block?
  # Create a vector of size n.chr to check whether each chromosome has the minimum
  # number of markers required to do linkage mapping. If not, return '0'
  markersPerChr <- rep(n.minMarkers, times=n.chr)
  if (any(nmar(cross) < n.minMarkers)) {
    return (0)
  }
  
  # R/qtl2
  cross2 <- convert2cross2(cross)
  map <- insert_pseudomarkers(cross2$gmap, step=n.step)
  probs <- calc_genoprob(cross2, map, error_prob=n.errorProb)
  s1.phes <- cross2$pheno[,c(1,2,3,4,5)]
  s1.output  <- scan1(genoprobs=probs, pheno=s1.phes, cores=n.cores)
  s1.perm <- scan1perm(genoprobs=probs, pheno=s1.phes, n_perm=1000, cores=n.cores)
  thresholds <- summary(s1.perm, alpha=0.05)
  
  # Use a LOD support interval of 5
  sigQtl <- qtl2::find_peaks(scan1_output=s1.output, map=map, threshold=thresholds, peakdrop = 5)
  if (saveQtlPlots) {
    plotLinkageMap(RIL, s1.output=qtl2toqtl1(s1.output, map), s1.perm=s1.perm, save_dir)
  }
  
  # Create a list storing the number of peaks per trait
  res <- c(sum(sigQtl$lodindex==1),
           sum(sigQtl$lodindex==2),
           sum(sigQtl$lodindex==3),
           sum(sigQtl$lodindex==4),
           sum(sigQtl$lodindex==5))
  return (res)
}

# Calculate the number of significant interactions between QTL
# Runs scantwo() and processes the result, calculating the number of significant
# interaction effects, by evaluating whether the LOD score between the peak
# LOD values on each chromosome is above a pre-determined threshold
# Known inaccuracies:
# - uses predefined [conservative] significance thresholds to avoid expensive permutation tests
# - only considers the top LOD peak on each chromosome, which may undercount
#   LOD peaks if there are multiple QTL per chromomsome
# RIL: a recombinant inbred line family
# parent1: one of the parents crossed to create the RIL
# parent2: the other parent crossed to create the RIL
# trait: an integer (1: Trait 1, 2: Trait 2: 3: Yield, 4: Suitability, 5: Breeding Fitness)
# save_dir: a directory in which to save the plots
# Returns: the number of significant LOD peaks
epistaticLodPeaks <- function(RIL, parent1, parent2, trait, save_dir) {
  cross <- getCross(RIL, parent1, parent2, "riself", snpChip=1)
  cross <- drop.nullmarkers(cross)
  # Check to see if there are no markers once null markers have been dropped
  if (length(cross$geno) == 0) {
    return (0)
  }
  
  s2.phes <- phenames(cross)[trait]
  # Run QTL mapping5
  cross <- calc.genoprob(cross, step=2, error.prob=n.errorProb)
  
  s2.output <- scantwo(cross, pheno.col=s2.phes, method=n.mappingMethod,
                   n.cluster=n.cores, verbose=FALSE)
  #s2.perm <- scantwo(cross, pheno.col=s2.phes, method=n.mappingMethod,
  #                  n.perm=1000, n.cluster=n.cores)
  #summary(s2.perm)

  # These predefined thresholds were recommended by AlphaSimR:
  # https://rqtl.org/tutorials/new_summary_scantwo.pdf
  # Use these to avoid the significant computational time required to do 1000
  # permutations
  # These were obtained by looking at 10 independendent 1000 perm tests
  thresholds <- c(6.6,5.2,4.2,5.0,3.1)
  # This function only calculates based on the pairs of peak values on each chromosome
  # It does not account for multiple QTL per chromosome
  peaks <- summary(object=s2.output, thresholds=thresholds, what="int")
  
  findQtl <- function(chrId, pos) {
    qtls <- getQtlEffectSizes(RIL) %>%
      dplyr::filter(chr==chrId)
    
    if (nrow(qtls) == 0) {
      return(NA)
    }
    
    chrMap <- as.data.frame(SP$genMap[[chrId]])
    colnames(chrMap) <- c("pos")
    
    chrMap <- chrMap %>%
      dplyr::mutate(pos=pos*n.genMapLen)
    
    snp <- rownames(chrMap)[which.min(abs(chrMap$pos - pos))]
    snp_pos <- as.numeric(sub(".*_", "", snp))

    qtl_id <- qtls$id[which.min(abs(qtls$pos - snp_pos))]
    
    return(qtl_id)
  }
  
  if (nrow(peaks) > 0) {
    peaks <- peaks %>%
      dplyr::mutate(chr1=as.numeric(chr1),
                    chr2=as.numeric(chr2),
                    pos1=as.numeric(pos1),
                    pos2=as.numeric(pos2)) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(qtl1 = findQtl(chr1, pos1),
                    qtl2 = findQtl(chr2, pos2)) %>%
      dplyr::ungroup()
  }

  if (saveQtlPlots) {
    plot2DLinkageMap(RIL=RIL, s2.output=s2.output, save_dir=save_dir)
  }
  return(peaks)
}


