
source("Scripts/CreateFounderPop.R")

n.tests <- 10

results <- data.frame(
  bulkMean=c(),
  bulkVar=c(),
  pedMean=c(),
  pedVar=c(),
  hybridMean=c(),
  hybridVar=c()
)
for (i in 1:n.tests) {
  fitFunc <- gaussianLandraceFitFunc
  fitCalc <- calculateFitnessGaussian

  SP$setVarE(h2=c(n.h2, n.h2, n.yieldH2))
  source("Scripts/CreateIndependentPops.R")
  SP$setVarE(h2=c(n.h2Breeding, n.h2Breeding, n.yieldH2Breeding))
  fitFunc <- gaussianFitFunc
  
  
  bulkPure1 <- makePurelinesBulk(pops[[1]])
  bulkPure2 <- makePurelinesBulk(pops[[2]])
  bulk1 <- bulkPure1[[1]]
  bulk2 <- bulkPure1[[2]]
  bulk3 <- bulkPure2[[1]]
  bulk4 <- bulkPure2[[2]]
  
  pedPure1 <- makePurelinesPed(pops[[1]])
  pedPure2 <- makePurelinesPed(pops[[2]])
  ped1 <- pedPure1[[1]]
  ped2 <- pedPure1[[2]]
  ped3 <- pedPure2[[1]]
  ped4 <- pedPure2[[2]]
  
  hybPure1 <- makePurelinesHybrid(pops[[1]])
  hybPure2 <- makePurelinesHybrid(pops[[2]])
  hyb1 <- hybPure1[[1]]
  hyb2 <- hybPure1[[2]]
  hyb3 <- hybPure2[[1]]
  hyb4 <- hybPure2[[2]]
  
  newRow <- data.frame(
    bulkMean=c(meanP(bulk1)[1], meanP(bulk1)[2], meanP(bulk2)[1], meanP(bulk2)[2], meanP(bulk3)[1], meanP(bulk3)[2], meanP(bulk4)[1], meanP(bulk4)[2]),
    bulkVar=c(varP(bulk1)[1,1], varP(bulk1)[2,2], varP(bulk2)[1,1], varP(bulk2)[2,2], varP(bulk3)[1,1], varP(bulk3)[2,2], varP(bulk4)[1,1], varP(bulk4)[2,2]),
    pedMean=c(meanP(ped1)[1], meanP(ped1)[2], meanP(ped2)[1], meanP(ped2)[2], meanP(ped3)[1], meanP(ped3)[2], meanP(ped4)[1], meanP(ped4)[2]),
    pedVar=c(varP(ped1)[1,1], varP(ped1)[2,2], varP(ped2)[1,1], varP(ped2)[2,2], varP(ped3)[1,1], varP(ped3)[2,2], varP(ped4)[1,1], varP(ped4)[2,2]),
    hybridMean=c(meanP(hyb1)[1], meanP(hyb1)[2], meanP(hyb2)[1], meanP(hyb2)[2], meanP(hyb3)[1], meanP(hyb3)[2], meanP(hyb4)[1], meanP(hyb4)[2]),
    hybridVar=c(varP(hyb1)[1,1], varP(hyb1)[2,2], varP(hyb2)[1,1], varP(hyb2)[2,2], varP(hyb3)[1,1], varP(hyb3)[2,2], varP(hyb4)[1,1], varP(hyb4)[2,2])
  )
  
  results <- rbind(results, newRow)
}

results %>%
  summarise(across(where(is.numeric), 
                   list(var = \(x) var(x, na.rm = TRUE))))
