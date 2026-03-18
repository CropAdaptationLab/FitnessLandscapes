library(dplyr)
library(stringr)
library(tibble)
library(qtl)
library(qtl2)

#setwd("/pl/active/Morris_CSU/Ted_Monyak/Lab-Notebooks/members/Ted/breeding_sims")
setwd("~/Documents/CSU/Lab-Notebooks/members/Ted/breeding_sims")
#source("functions/Fitness.R")

save_dir <- file.path(getwd(), "output/NAM")

codes.df <- read.csv("data/NAM_population_codes.csv") %>%
  dplyr::rename(Taxa="Plot_ID")

codes.df <- rbind(codes.df,
                  data.frame(
                    Taxa=c("P898012", "SC265", "SC283", "Macia", "SC35", "SC1103", "Segeolane", "SC1345", "Ajabsido", "SC971", "RTx430"),
                    Pedigree=c("P898012", "SC265", "SC283", "Macia", "SC35", "SC1103", "Segeolane", "SC1345", "Ajabsido", "SC971", "RTx430"),
                    LIB=c("SAMEA12302610", "ISIM", "IWNY", "ISGP", "IZGS", "SAMEA12302682", "ISHA", "IUIN", "ISGR", "IXYF", "IFMC")
                  ),
                  data.frame(
                    Taxa=c("PI656057", "PI533766", "PI533869", "PI565121", "PI534133", "PI576434", "PI656023", "PI597980", "PI656015", "PI656111", "PI655996"),
                    Pedigree=c("P898012", "SC265", "SC283", "Macia", "SC35", "SC1103", "Segeolane", "SC1345", "Ajabsido", "SC971", "RTx430"),
                    LIB=c("SAMEA12302610", "ISIM", "IWNY", "ISGP", "IZGS", "SAMEA12302682", "ISHA", "IUIN", "ISGR", "IXYF", "IFMC")
                  ))

height.df <- read.csv("data/Plant_height_NAM.csv") %>%
  dplyr::filter(Location=="Manhattan") %>%
  dplyr::rename("PH"=Plant_height) %>%
  dplyr::select(Taxa,PH) %>%
  dplyr::group_by(Taxa) %>%
  dplyr::summarize(PH=mean(PH, na.rm=TRUE),
                   .groups='drop') %>%
  inner_join(codes.df) %>%
  dplyr::select(-c(Taxa, Pedigree))


ft.df <- read.csv("data/Flowering_time_NAM_2015.csv") %>%
  dplyr::filter(Year=="2015") %>%
  dplyr::rename("FT"=Flowering_time) %>%
  dplyr::select(Taxa,FT) %>%
  dplyr::group_by(Taxa) %>%
  dplyr::summarize(FT=mean(FT, na.rm=TRUE),
                   .groups='drop') %>%
  replace(., . == "SC35-14E", "SC35") %>%
  replace(., . == "Segaolane", "Segeolane") %>%
  inner_join(codes.df) %>%
  dplyr::select(-c(Taxa, Pedigree))
pheno.df <- inner_join(height.df, ft.df)

#setwd("/pl/active/Morris_CSU/Ted_Monyak/NAM")
saveRDS(pheno.df, file="data/pheno.rds")

#pheno.df <- pheno.df %>% dplyr::mutate(
#  Pedigree=sub(".*/", "", Pedigree),
#  Pedigree=sub(" .*$", "", Pedigree)
#)

means.df <- pheno.df %>%
  dplyr::group_by(Pedigree) %>%
  dplyr::summarize(ft=mean(FT),
                   ht=mean(PH))


FT_opt <- 71
PH_opt <- 150

FT_s <- sd(pheno.df$FT)
PH_s <- sd(pheno.df$PH)

pheno.df <- pheno.df %>%
  dplyr::mutate(FT_norm=(FT_opt-FT)/(FT_s),
                PH_norm=(PH_opt-PH)/(PH_s))

pheno.df <- pheno.df %>%
  dplyr::filter(Pedigree=="P898012") %>%
  dplyr::mutate(LIB=str_replace_all(LIB, ":", ".")) %>%
  dplyr::mutate(Suitability=exp(-(FT_norm^2 + PH_norm^2) / (16)))

ggplot(pheno.df, aes(x=Suitability)) +
  geom_density()


p898012 <- read.cross(format="csv",
                      dir="data",
                      file="P898012.cross_new.csv",
                      genotypes=c("AltAlt","TxTx"))
ids <- data.frame(p898012[["pheno"]][["id"]]) %>%
  rownames_to_column(var="idx")

colnames(ids)[2] <- "LIB"
ids$LIB <- as.character(ids$LIB)
ids$idx <- as.numeric(ids$idx)

merged.df <- inner_join(pheno.df, ids, by="LIB") %>%
  dplyr::arrange(idx)

# 187th index is missing

remove_genotype <- function(cross, chr, idx=187) {
  geno <- cross[["geno"]][[as.character(chr)]][["data"]]
  geno <- geno[-idx,]
  cross[["geno"]][[as.character(chr)]][["data"]] <- geno
  return (cross)
}

for (chr in 1:10) {
  p898012 <- remove_genotype(p898012, chr)  
}

p898012[["pheno"]] <- merged.df %>% dplyr::select(FT_norm, PH_norm, Suitability)
cross <- p898012
cross <- drop.nullmarkers(cross)

# Check to see if there are no markers once null markers have been dropped
if (length(cross$geno) == 0) {
  return (0)
}
phenames(cross)
n.step <- 1
n.cores <- 8
#n.errorProb <- 
# R/qtl1:
# Get the phenotype names
## Run QTL mapping
phes <- phenames(cross)[c(1,2,3)]
cross <- calc.genoprob(cross, step=n.step)
out.hk <- scanone(cross, pheno.col=phes, method="hk",
                  n.cluster=4)
operm.hk <- scanone(cross, pheno.col=phes, method="hk",
                    n.perm=1000,n.cluster=n.cores)
plot(out.hk, lodcolumn = c(1,2,3))

s2.phes <- phenames(cross)[3]
s2.output <- scantwo(cross, pheno.col=s2.phes, method="hk",
                     n.cluster=n.cores, verbose=FALSE)

summary(s2.output, what="int")
plot(s2.output)
#cross2 <- convert2cross2(cross)
#map <- insert_pseudomarkers(cross2$gmap, step=n.step)
#probs <- calc_genoprob(cross2, map, error_prob=n.errorProb)
#s1.phes <- cross2$pheno[,c(1,2,3,4,5)]
#s1.output  <- scan1(genoprobs=probs, pheno=s1.phes, cores=n.cores)
#s1.perm <- scan1perm(genoprobs=probs, pheno=s1.phes, n_perm=1000, cores=n.cores)
#thresholds <- summary(s1.perm, alpha=0.05)

#plotLinkageMap(RIL, s1.output=out.hk, s1.perm=operm.hk, save_dir)
lod_matrix <- s2.output$lod

# Mask the lower triangle (set to NA so it's ignored)
lod_upper <- t(lod_matrix)
lod_upper[upper.tri(lod_upper, diag = TRUE)] <- NA

best <- which(lod_upper == max(lod_upper, na.rm=TRUE), arr.ind = TRUE)

mar1 <- find.marker(cross, chr=6, pos=55)
mar2 <- find.marker(cross, chr=10, pos=51)
g1 <- pull.geno(cross)[, mar1]
g2 <- pull.geno(cross)[, mar2]
pheno <- cross$pheno[, 3]

df <- data.frame(
  pheno = pheno,
  geno1 = factor(g1, levels = c(1, 2), labels = c("AA", "AB")),
  geno2 = factor(g2, levels = c(1, 2), labels = c("AA", "AB"))
)
df <- df[complete.cases(df), ]

means_df <- df %>%
  group_by(geno1, geno2) %>%
  summarise(mean_pheno = mean(pheno), se = sd(pheno)/sqrt(n()), .groups = "drop")

ggplot(means_df, aes(x = geno2, y = mean_pheno, color = geno1, group = geno1)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  labs(
    title = "Reaction Norm: Chr6 @ 55cM × Chr10 @ 51cM",
    x = paste("Genotype at", mar2, "(Chr10)"),
    y = "Mean Phenotype",
    color = paste("Genotype at\n", mar1, "(Chr6)")
  ) +
  theme_classic(base_size = 14) +
  scale_color_manual(values = c("#E41A1C", "#377EB8"))

qtl1 <- "Chr 6"
qtl2 <- "Chr 10"
ggplot(means_df, aes(x = geno2, y = mean_pheno, group = geno1, color = geno1)) +
  geom_line(size = 1) +
  geom_point() +
  scale_color_manual(values = c("1" = "#CC0000", "2"="#3C78D8")) +
  scale_x_discrete(labels = c("1", "2"),
                   expand = c(0.1, 0.1)) +
  labs(
    x = paste0(qtl2, " Genotype"),
    y = "Suitability",
    color = paste0(qtl1, " Genotype")
  ) +
  theme_minimal(base_size = 10,
                base_family="Helvetica") +
  theme(
    plot.title = element_blank(),
    axis.title.y = element_text(margin=margin(t=0, r=10, b=0, l=10, unit="pt")),
    legend.position = "right",
    legend.key.spacing.y = unit(5, "pt"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    #axis.line = element_line(color = "black"),
    panel.background = element_rect(fill = "white", color = "black")
  )
save_dir
  fname <- file.path(save_dir, "rxn_norm.jpg")
  ggplot2::ggsave(filename = fname,
                  device = "jpg",
                  height=2,
                  width=3.5,
                  units="in",
                  dpi=600)
  