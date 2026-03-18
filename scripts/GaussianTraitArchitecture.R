# LOOKS LIKE WE WANT VAR = 0.4, SHAPE = 1
n.var <- 0.4
n.initTraitVal <- 2
n.initYieldVal <- 50
n.yieldVar <- 50
n.yieldQtlPerChr <- 20
n.gens <- 200 
n.popSize <- 500
normalDist <- FALSE
n.popResets <- 5
shapes <- c(1,2,3,4,5)
vars <- c(0.2,0.4,0.6,1)
#allQtl <- data.frame(id=c(),
#                     eff_size=c(),
#                     shape=c())

for (v in 1:length(vars)) {
  print(paste0("VAR: ", v))
  n.var <- vars[v]
  for (s in 1:length(shapes)) {
    print(paste0("Shape: ", s))
    n.shape <- shapes[s]
    for (r in 1:n.popResets) {
      print(paste0("Pop: ", r))
      source("Scripts/CreateFounderPop.R")
      traitArch <- twoTraitArchitecture(founderPop)
      traitArch['shape'] <- n.shape
      traitArch['var'] <- n.var
      allQtl <- rbind(allQtl, traitArch)
    }
  }
}

allQtl$shape <- as.numeric(allQtl$shape)
allQtl$var <- as.factor(allQtl$var)

allQtl %>%
  dplyr::filter(shape==1) %>%
  ggplot(aes(x=eff_size, fill=var, color=var)) +
  geom_density(size=1, alpha=0.1) +
  labs(title="Shape 1", x="Effect Size", y="Density")

allQtl %>%
  dplyr::filter(shape==2) %>%
  ggplot(aes(x=eff_size, fill=var, color=var)) +
  geom_density(size=1, alpha=0.1) +
  labs(title="Shape 2", x="Effect Size", y="Density")

allQtl %>%
  dplyr::filter(shape==3) %>%
  ggplot(aes(x=eff_size, fill=var, color=var)) +
  geom_density(size=1, alpha=0.1) +
  labs(title="Shape 3", x="Effect Size", y="Density")

allQtl %>%
  dplyr::filter(shape==4) %>%
  ggplot(aes(x=eff_size, fill=var, color=var)) +
  geom_density(size=1, alpha=0.1) +
  labs(title="Shape 4", x="Effect Size", y="Density")

allQtl %>%
  dplyr::filter(shape==5) %>%
  ggplot(aes(x=eff_size, fill=var, color=var)) +
  geom_density(size=1, alpha=0.1) +
  labs(title="Shape 5", x="Effect Size", y="Density")


allQtl$var <- as.factor(allQtl$var)
allQtl$shape <- as.factor(allQtl$shape)

allQtl %>%
  dplyr::filter(var==0.2) %>%
  ggplot(aes(x=eff_size, fill=shape, color=shape)) +
  geom_density(size=1, alpha=0.1) +
  labs(title="Var 0.2", x="Effect Size", y="Density")

allQtl %>%
  dplyr::filter(var==0.4) %>%
  ggplot(aes(x=eff_size, fill=shape, color=shape)) +
  geom_density(size=1, alpha=0.1) +
  labs(title="Var 0.4", x="Effect Size", y="Density")

allQtl %>%
  dplyr::filter(var==0.6) %>%
  ggplot(aes(x=eff_size, fill=shape, color=shape)) +
  geom_density(size=1, alpha=0.1) +
  labs(title="Var 0.6", x="Effect Size", y="Density")

allQtl %>%
  dplyr::filter(var==1) %>%
  ggplot(aes(x=eff_size, fill=shape, color=shape)) +
  geom_density(size=1, alpha=0.1) +
  labs(title="Var 1", x="Effect Size", y="Density")

write.table(allQtl, file.path("~/Documents/CSU/Lab-Notebooks/members/Ted/breeding_sims/output/TraitArchitecture/allQtl.csv"), col.names=NA, quote=FALSE, sep=",")
