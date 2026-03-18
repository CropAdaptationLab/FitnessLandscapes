rm(list = ls())

# Set to false if running on Alpine
runLocal = TRUE

if (runLocal) {
  setwd("~/Documents/CSU/R/BreedingSims")
  output_dir <- file.path(getwd(), "Output")
} else {
  setwd("/pl/active/Morris_CSU/Ted_Monyak/BreedingSims")
  output_dir <- ("/scratch/alpine/c837220672@colostate.edu/Output")
}

library(AlphaSimR)
library(ade4)
library(devtools)
library(dplyr)
library(ggpubr)
library(ggplot2)
library(factoextra)
library(patchwork)
library(plotly)
library(purrr)
library(reshape2)
library(snow)
library(tibble)
library(tidyr)
library(viridis)
library(zoo)

source("Functions/Fitness.R")
source("Functions/TraitArchitecture.R")
source("Scripts/GlobalParameters.R")
source("Scripts/CreateFounderPop.R")

pop <- founderPop

geno <- getUniqueQtl(pop)
pheno <- pheno(pop)

n.exploreGens <- 200

pop <- founderPop
n.nPops <- 10
n.selProp <- 0.8
# Create a random vector of size n.pops, with a random order of sub-population ids
randVec <- sample(rep(c(1:n.nPops), times=n.popSize/n.nPops))

# Create n.nPops sub populations
pops <- vector(mode="list", length=n.nPops)
for (p in 1:n.nPops) {
  pops[[p]] <- selectInd(pop, trait=selectSubPop, selectTop=TRUE, nInd=n.subPopSize, idx=p, randVec=randVec)
}

#fig <- plot_ly()
#fig <- fig %>% layout(legend=list(title=list(text='Population Size')),
#                      scene = list(xaxis = list(title = "Trait A"),
#                                   yaxis = list(title = "Trait B"),
#                                   zaxis = list(title = "Fitness"),
#                                   aspectmode='cube')) %>% hide_colorbar()

for (p in 1:length(pops)) {
  print(p)
  pop <- pops[[p]]
  fit.df <- data.frame(gen=1:n.exploreGens,
                       fitness=numeric(n.burnInGens),
                       traitValA=numeric(n.burnInGens),
                       traitValB=numeric(n.burnInGens))
  for (gen in 1:n.exploreGens) {
    randPop <- selectInd(pop, nInd=5, use="rand")
    fit.df$fitness[gen] <- mean(twoTraitFitFunc(pheno(randPop)))
    fit.df$traitValA[gen] <- meanP(randPop)[1]
    fit.df$traitValB[gen] <- meanP(randPop)[2]
    #pop <- selectCross(pop, nInd = n.subPopSize, nCrosses = n.subPopSize, use = "rand")
    
    pop <- selectCross(pop, trait=twoTraitFitFunc, nInd=nInd(pop)*n.selProp, nCrosses=nInd(pop))
    geno <- rbind(geno, getUniqueQtl(randPop))
    pheno <- rbind(pheno, pheno(randPop))
  }
  #fig <- add_trace(
  #  fig,
  #  fit.df,
  #  name = p,
  #  x = fit.df$traitValA,
  #  y = fit.df$traitValB,
  #  z = fit.df$fitness,
  #  type = 'scatter3d',
  #  mode = 'lines',
  #  opacity = 1,
  #  color = p,
  #  line = list(width = 5))
}
#fig

PCA  = dudi.pca(df = geno, center = T, scale = F, scannf = F, nf = 2)
(VAF = 100 * PCA$eig[1:2] / sum(PCA$eig)) # variance explained
df.PCA = data.frame(
  "Pop" = c(rep("base",times=1000), (rep(c(1:n.nPops),each=5*n.exploreGens))),
  "PC1" = PCA$l1$RS1,
  "PC2" = PCA$l1$RS2)

ggplot(df.PCA, aes(x = PC1, y = PC2)) +
  geom_point(aes(color=factor(Pop))) +
  ggtitle("Population structure") +
  xlab(paste("Pcomp1: ", round(VAF[1], 2), "%", sep = "")) +
  ylab(paste("Pcomp2: ", round(VAF[2], 2), "%", sep = ""))


head(sort(df.PCA$PC1))
head(sort(PCA$l1$RS1))
length(sort(PCA$l1$RS1))
rows <- sort(unique(round(PCA$l1$RS2, digits=2)))
cols <- sort(unique(round(PCA$l1$RS1, digits=2)))

res_mat <- matrix(NA, nrow=length(rows), ncol=length(cols))
rownames(res_mat) <- rows
colnames(res_mat) <- cols

for (i in 1:nrow(geno)) {
  #geno[i,] %*% PCA$c1$CS1
  x <- round(PCA$l1$RS1[i],digits=2)
  y <- round(PCA$l1$RS2[i],digits=2)
  z <- calculateFitnessTwoTrait(pheno[i,1], pheno[i,2])
  res_mat[as.character(y),as.character(x)] <- z
}
rownames(res_mat) <- rows
colnames(res_mat) <- cols
res_mat.df <- data.frame(res_mat)
a <- na.approx(res_mat.df)

p <- plot_ly(x=rows,
             y=cols,
             z=a,
             type='surface',
             colors = "PuBuGn",
             opacity=1) %>%
  layout(scene = list(xaxis = list(title = "PC 1"),
                      yaxis = list(title = "PC 2"),
                      zaxis = list(title = "Fitness"),
                      aspectmode='cube'))
p
return (p) 