# Title: GLOBAL VARIABLES
# Author: Ted Monyak
# Description: This file contains default parameters for use in all breeding simulations

# POPULATION
n.popSize <- 2000 # Num. individuals in the founder population
n.subPopSize <- 1000 # Num individuals in the independent sub-populations
n.m <- 0 # Migration rate. This percentage of individuals will migrate to adjacent populations each generation in a stepping stone model
n.ne <- n.popSize # Effective population size
n.segSites <- 1000 # Initial # segregating alleles per chromosome in the population
n.markers <- 1000 # Num. loci per chromosome to include on the SNP chip 
n.chr <- 10 # Num. chromosomes (sorghum has 10)
n.minMAF <- 0.3 # Minimum minor allele frequency to use when assigning random QTL for the founder population
n.genMapLen <- 100 # Centimorgans per chromosome

# TRAITS
n.qtlPerChr <- 2 # Number of qtl per chromosome, per trait
n.h2 <- 0.1 # Narrow-sense heritability for each acquired trait for landrace adaptation
n.h2Breeding <- 0.8 # Narrow-sense heritability for each acquired trait for breeding adaptation
n.initTraitVal <- 2 # Starting value for each of the two traits
n.var <- 0.4 # Initial variance for each trait. This is used as the 'rate' parameter in a gamma distribution
n.shape <- 1 # Initial shape for the gamma distribution for each trait
n.allele <- 2 # This is the allele for which to track frequency over time
n.relAA <- 0 # The relative value of additive-by-additive variance compared to additive variance w/ allele freq. 0.5
SAMPLING <- "geometric" # Draw from either a "gamma", "normal" of effect sizes, or use a "geometric" distribution

# YIELD
n.yieldQtlPerChr <- 50 # No. QTL per chromosome for yield
n.yieldVar <- 50 # Initial variance for yield
n.yieldH2 <- 0.05 # h2 for yield
n.yieldH2Breeding <- 0.3 # h2 for yield at the breeding stage
n.initYieldVal <- 50 # Starting mean value for yield

# ADAPTIVE WALKS
n.margin <- 0 # Populations will terminate their adaptive walks once they reach this fitness value
n.gens <- 200 # maximum number of generations for an adaptive walk
n.burnInGens <- 5 # number of burn-in generations for founder population
n.nPops <- 2 # number of independent subpopulations to create
n.selProp <- 0.1 # % of the population to advance during main adaptive walk
n.selR <- 1 # the r value to use in the geometric series for a decaying selection intensity. Set to 1 for a non-decaying intensity, and decrease this to increase the rate of decay

# SIMULATIONS
n.popResets <- 1 # number of times to reset the founder population in a simulation
n.reps <- 1 # number of replications to run per pop reset

# MAPPING POPULATIONS
n.RILs <- 250 # number of RILs to create
n.indPerRIL <- 4 # number of replicates in each RIL

# PLOTTING
saveQtlPlots <- FALSE # Save linkage map plots
saveTraitPlots <- FALSE # Save trait distribution plots
saveAllelePlots <- FALSE # Save allele frequency plots
saveFitnessPlots <- FALSE # Save adaptive walk and fitness landscape plots
saveFixationOrder <- TRUE # Save the fixation order of alleles
saveEffectSizes <- TRUE # Save the effect sizes in RIL families

# BREEDING
n.F1 <- 1000 # Initial number of individuals to select from each landraces to purify
n.F2 <- 2000
n.F3 <- 1000
n.F4 <- 500
n.F5 <- 250
n.F6 <- 125
n.F7 <- 60
n.F8 <- 30
n.F9 <- 15
n.F10 <- 10

# QTL Mapping Parameters
n.mappingMethod <- "hk" # One of 'hk' (Haley-Knott), 'em', 'imp' (imputation), or others. See ?scanone
n.errorProb <- 0 # assumed genotyping error rate. See ?calc.genoprob
n.step <- 1 # max distance b/t calculated genotype probabilities. See ?calc.genoprob
n.cores <- 8 # Number of cores to use for QTL mapping and population generation
n.minMarkers <- 5 # Minimum number of markers per chromosome required to do QTL mapping after all het and monomorphic markers have been removed

# Function for returning all parameters (including ones that have been updated)
# for writing to the params.txt file along with graphs
getParams <- function() {
  n.df <- data.frame(
    popSize=n.popSize,
    subPopSize=n.subPopSize,
    ne=n.ne,
    segSites=n.segSites,
    markers=n.markers,
    chr=n.chr,
    minMAF=n.minMAF,
    qtlPerChr=n.qtlPerChr,
    h2=n.h2,
    initTraitVal=n.initTraitVal,
    var=n.var,
    shape=n.shape,
    a=n.a,
    L=n.L,
    relAA=n.relAA,
    sampling=SAMPLING,
    allele=n.allele,
    margin=n.margin,
    selProp=n.selProp,
    r=n.selR,
    gens=n.gens,
    burnInGens=n.burnInGens,
    reps=n.reps,
    popResets=n.popResets,
    RILs=n.RILs,
    indPerRIL=n.indPerRIL,
    nPops=n.nPops,
    mappingMethod=n.mappingMethod,
    errorProb=n.errorProb,
    step=n.step,
    cores=n.cores,
    F1=n.F1,
    F2=n.F2,
    F3=n.F3,
    F4=n.F4,
    F5=n.F5,
    F6=n.F6,
    F7=n.F7,
    F8=n.F8,
    yieldQtlPerChr=n.yieldQtlPerChr,
    yieldVar=n.yieldVar,
    yieldH2=n.yieldH2,
    yieldH2Breeding=n.yieldH2Breeding,
    initYieldVal=n.initYieldVal
  )
  return (t(n.df))
}
