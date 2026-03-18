

## NOT USED BELOW HERE

alleleEffectSizeEstimate <- function(chr, pos, marker, probs, map, trait) {
  g <- maxmarg(probs, map, chr=chr, pos=pos, return_char=TRUE)
  #par(mar=c(4.1, 4.1, 0.6, 0.6))
  #plot_pxg(g,cross2$pheno[,1], ylab="Trait 1 Phenotype")
  phenos <- data.frame(geno=g,
                       pheno=cross2$pheno[,trait]) %>%
    filter(geno =="AA" | geno == "BB")
  phenos$geno <- as.factor(phenos$geno)
  
  ggplot(phenos, aes(x=geno, y=pheno, fill=geno)) +
    geom_boxplot(na.rm=TRUE, outliers=FALSE) +
    stat_compare_means(method="t.test",
                       label="p.signif",
                       hide.ns=FALSE) +
    labs(title=paste0("Marker ", marker), y="Trait 1") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank())
}


LDplot <- function(RIL) {
  geno <- as.data.frame(pullSegSiteGeno(RIL))
  maf <- geno %>%
    summarise(across(everything(), ~{
      n_2 <- sum(. == "2")
      n_1 <- sum(. == "1")
      n_0 <- sum(. == "0")
      total <- n_2 + n_1 + n_0
      
      freq_1 <- (2 * n_2 + n_1) / (2 * total)
      freq_0 <- (2 * n_0 + n_1) / (2 * total)
      
      min(freq_1, freq_0)
    }))
  keep <- as.numeric(maf[1, ]) >= 0.05
  geno <- geno[, keep]
  
  nsnp <- ncol(geno) # Number of snp loci
  
  # # Empty matrix to hold LD values
  c.mat <- matrix(NA, nrow = nsnp, ncol = nsnp)
  rownames(c.mat) <- colnames(c.mat) <- colnames(geno) # Add row and column names
  
  diag(c.mat) <- 1 # Make diagonal values = 1
  
  # For loop to compute pairwise r^2 between all loci
  for (i in seq_len(nsnp - 1)) {
    
    for (j in (i + 1):nsnp) {
      
      c.mat[i, j] <- c.mat[j, i] <- round(cor(geno[, i], geno[, j]) ** 2, 2)
      
    }
  }
  
  # Check for symmetry in computed LD matrix
  isSymmetric(c.mat)
  
  n.genMapLen <- 100
  # Get the genetic map, and flatten it so it isn't grouped by chromosome
  genMap <- as.data.frame(unlist(SP$genMap)) %>%
    rename_with(~ "pos", .cols=1) %>%
    rownames_to_column("snp") %>%
    dplyr::mutate(locus = sub(".*\\.", "", snp), # Get ids in the form chr_loc
                  chr = as.numeric(sub("_.*", "", locus)), # Determine the chromosome of each qtl
                  pos = pos*n.genMapLen) %>%
    dplyr::select(-snp)
  
  genMap <- genMap[keep, ]
  
  # Get unique chromosomes and calculate cumulative positions
  chr_info <- data.frame(chr = seq(1:n.chr)) %>%
    mutate(
      chr_length = n.genMapLen,
      cumul_start = cumsum(c(0, head(chr_length, -1))),
      cumul_mid = cumul_start + n.genMapLen/ 2
    )
  
  # Add cumulative positions to effect sizes
  genMap <- genMap %>%
    left_join(chr_info %>% dplyr::select(chr, cumul_start), by = "chr") %>%
    dplyr::mutate(cumul_pos = pos + cumul_start)
  
  # Plot a tiny part of the LD matrix
  #LD.plot(c.mat, snp.positions = genMap$cumul_pos)
  
  #' Color palette for plotting LD heatmap
  rgb.palette <- colorRampPalette(rev(c("blue", "yellow", "red")), 
                                  space = "rgb")
  
  # It accepts raw geno data or numerical format or already computed LD values
  myLDmap <- LDheatmap::LDheatmap(c.mat, 
                                  genetic.distances = genMap$pos,
                                  color = rgb.palette(20),
                                  distances = 'genetic',
                                  geneMapLabelX = 0.4,
                                  geneMapLabelY = 0.75,
                                  geneMapLocation = 0.1,
                                  flip = TRUE,
                                  text = FALSE,
                                  SNP.name = colnames(geno))
  myLDmap
  
  #geno_str <- geno %>%
  #  mutate(across(everything(), ~case_when(
  #    . == 2 ~ "D/D",
  #    . == 1 ~ "D/I",
  #    . == 0 ~ "I/I"
  #  )))
  #g <- makeGenotypes(geno_str)
  #ld <- genetics::LD(g)
}

plotAlleleEffect <- function() {
  alleleEffectSizeEstimate(3, 93.58591, "3_926", probs, map, trait=1)
  alleleEffectSizeEstimate(1, 18.786996, "1_186", probs, map, trait=2)
  View(SP$traits)
  View(pullQtlGeno(RIL, 1))
  #hap1 <- pullQtlHaplo(pop, 1)
  #SP$qtlIndex
  #SP$traits
  #matrix(sample.int(sum(n.qtlPerChr),sum(n.qtlPerChr)),ncol=2)
  
  #View(SP$traits)
  #SP$traits[[1]]@epiEff
}

gwas <- function() {
  pop <- founderPop
  # Get genotype in rrBLUP encoded format
  geno <- getGwasGeno(pop)
  # Get phenotype data
  # dataframe has 3 columns: id, pheno.Trait1, pheno.Trait2
  pheno = data.frame(id = 1:pop@nInd,
                     pheno  = pop@pheno)
  model.1 <- GWAS(pheno[,c(1,2)], geno, plot = T, n.core=n.cores)
  model.2 <- GWAS(pheno[,c(1,3)], geno, plot = T, n.core=n.cores)
  model.3 <- GWAS(pheno[,c(1,4)], geno, plot = T, n.core=n.cores)
  
  geno
  
  cor(pop@pheno[,2],pop@pheno[,3])
  
  # Run GWAS for 2 traits
  save_dir <- file.path(output_dir, "GWAS")
  if (!dir.exists(save_dir)) dir.create(save_dir)
  save_dir <- file.path(save_dir, format(Sys.time(), "%F_%H_%M_%S"))
  if (!dir.exists(save_dir)) dir.create(save_dir)
  #fname <- file.path(save_dir, "gwas.pdf")
  #pdf(fname)
  #dev.new()
  
  library(qqman)
  # select the 2nd column for trait 1
  
  #colnames(model.1) <- c("SNP", "CHR", "BP", "P")
  
  
  str(model.1)
  manhattan(model.1,  ylim=c(0,10))
  #dev.off()
  #plot(model.1)
  
  fname <- file.path(save_dir, "gwas_trait2_NAM.pdf")
  pdf(fname)
  # select the 3rd column for trait 2
  dev.off()
  
  
  model.pc = GWAS(pheno[,c(1,trait+1)], geno, n.PC = 3, plot = T)
  
  model1$Trait1 = 10 ^ (-model1$pheno.Trait1)
  model2$Trait1 = 10 ^ (-model2$pheno.Trait1)
  
  model1$Trait2 = 10 ^ (-model1$pheno.Trait2)
  model2$Trait2 = 10 ^ (-model2$pheno.Trait2)
  
  qtl1 = as.vector(getQtlMap(trait = 1)$id)
  qtl2 = as.vector(getQtlMap(trait = 2)$id)
  par(mfrow = c(2, 1))
  manhattan(model1, chr = "chr", bp = "pos", p = "pheno.Trait2", snp = "snp", highlight = qtl2,
            main = "Marker", ylim = c(0,10))
  manhattan(model2, chr = "chr", bp = "pos", p = "pheno.Trait2", snp = "snp", highlight = qtl2,
            main = "Marker + principal components", ylim = c(0,10))
  # Check QQ-plot
  par(mfrow = c(2, 1))
  qq(model1$Trait1, main = "Marker")
  qq(model2$Trait1, main = "Marker + principal components")
  
  
  geno <- as.data.frame(pullQtlGeno(founderPop, 1))
  genoVar <- geno %>%
    summarize(across(where(is.numeric), var, na.rm = TRUE))
  
  x <- c(1:length(geno))
  x <- c(1:nInd(founderPop))
  
  
  y_orig <- SP$traits[[2]]@epiEff
  y <- matrix(nrow=20, ncol=2)
  y[1:10, ] <- y_orig[, c(1, 3)]
  y[11:20, ] <- y_orig[, c(2, 3)]
  
  x <- y[,1]
  eff <- y[,2]
  plot(c(1:n.popSize),c(aa(pop)[,2]))
}
