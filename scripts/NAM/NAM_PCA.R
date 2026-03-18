library(akima)
library(ggplot2)
library(dplyr)
library(fields)
library(RColorBrewer)
library(parallel)
library(plotly)
library(tibble)
library(tidyr)
library(viridis)

setwd("~/Documents/CSU/Lab-Notebooks/members/Ted/breeding_sims/data/pca_sap")
output_dir <- "~/Documents/CSU/Lab-Notebooks/members/Ted/breeding_sims/output/NAM/pca_sap"

source("../../plotting/PCA.R")
source("../../plotting/FitnessLandscapes.R")

pca_df <- readRDS("pca_df.rds")
VAF <- readRDS('vaf.rds')
full_pca <- readRDS("full_pca_projection.rds")
pca_df <- rownames_to_column(as.data.frame(full_pca), var="LIB")

metadata <- read.csv("../metadata_NAM.txt", header = T,sep = "\t")
metadata <- metadata %>% dplyr::select(LIB, Family)

founders <- read.csv("../founders.csv", header=T)
metadata <- rbind(metadata, founders)

# filtering in case there are some missing
metadata <- metadata[metadata$LIB %in% pca_df$LIB, ]
metadata <- metadata[match(pca_df$LIB, metadata$LIB), ]

# Merge PCA results with the coloring column
full_pca_plot_df<- left_join(pca_df, metadata)

# DROP IZPY and IZPZ (duplicates)
full_pca_plot_df <- full_pca_plot_df[!(full_pca_plot_df$LIB %in% c("IZPY", "IZPZ")),]

family_colors <- data.frame(
  Family=c("Macia", "P898012", "Segeolane", "SC35", "SC283", "SC265", "SC1345", "Ajabsido", "SC971", "SC1103", "RTx430"),
  family_color=c("#da0001", "#70cd30", "#aacc4b", "#902d27", "#e9a425", "#9dcce9", "#8c2eea", "#1b1bfb", "#df7d71", "#c5b192", "black")
)
RIL_family_names <- paste0(family_colors$Family,"_RIL")
family_colors <- rbind(family_colors,
                       data.frame(Family=RIL_family_names,
                                  family_color=family_colors$family_color))
full_pca_plot_df <- left_join(full_pca_plot_df, family_colors, by="Family")

founder_genos <- tail(full_pca_plot_df, n=11)
pheno.df <- readRDS("../pheno.rds")
pheno.df <- pheno.df[pheno.df$LIB %in% full_pca_plot_df$LIB,]
#pheno.df <- pheno.df[match(pca_plot_df$LIB, pheno.df$LIB),]

full_pca_plot_df <- inner_join(full_pca_plot_df, pheno.df)

means.df <- full_pca_plot_df %>%
  dplyr::group_by(Family) %>%
  dplyr::summarize(ft=mean(FT),
                   ft_sd=sd(FT),
                   ht=mean(PH),
                   ht_sd=sd(PH))

FT_mean <- mean(full_pca_plot_df$FT)
PH_mean <- mean(full_pca_plot_df$PH)

FT_s <- sd(full_pca_plot_df$FT)
PH_s <- sd(full_pca_plot_df$PH)

opt.mat <- matrix(data=c(
  73, 134,
  68, 112,
  69, 115,
  70, 129,
  72, 113,
  72.5, 106.5),
  nrow=6,
  ncol=2,
  byrow=TRUE
)

rownames(opt.mat) <- c("P898012_RIL", "SC265_RIL", "SC283_RIL", "Macia_RIL", "SC35_RIL", "RTx430")
colnames(opt.mat) <- c("FT", "PH")

# Just use RTx430 as optimal values
m <- 2
full_pca_plot_df <- full_pca_plot_df %>%
  dplyr::mutate(
    FT_norm = (opt.mat["RTx430", "FT"] - FT)/(FT_s*m),
    PH_norm = (opt.mat["RTx430", "PH"] - PH)/(PH_s*m),
    Suit=exp(-(FT_norm^2 + PH_norm^2) /2)
  )
#plot_families <- unique(full_pca_plot_df$Family)
plot_families <- c("SC35_RIL", "SC283_RIL", "Macia_RIL", "SC265_RIL", "Ajabsido_RIL",
                   "SC35", "SC283", "Macia", "SC265", "Ajabsido", "RTx430")

pca_plot_df <- full_pca_plot_df %>%
  dplyr::filter(Family %in% plot_families)

plot_NAM_trait_distribution(pca_plot_df)
ggplot2::ggsave(filename = file.path(output_dir, "trait_distr.jpg"),
                device = "jpg",
                width=5.5,
                height=2.5,
                dpi=600)
ggplot2::ggsave(filename = file.path(output_dir, "trait_distr.pdf"),
                device = "pdf",
                width=5.5,
                height=2.5)

pca_plot_df <- pca_plot_df %>%
  dplyr::select(-c(LIB, PH, FT))

plot_NAM_suit_distribution(pca_plot_df)
ggplot2::ggsave(filename = file.path(output_dir, "suit_distr.jpg"),
                device = "jpg",
                width=3.5,
                height=2.5,
                dpi=600)
ggplot2::ggsave(filename = file.path(output_dir, "suit_distr.pdf"),
                device = "pdf",
                width=3.5,
                height=2.5)



#r2_interactions <- mclapply(pc_pairs, function(pair) {
#  m <- lm(Suit ~ pca_plot_df[[pair[1]]] * pca_plot_df[[pair[2]]], 
#          data = pca_plot_df)
#  summary(m)$r.squared
#}, mc.cores = n_cores)

#r2_interactions <- unlist(r2_interactions)
#names(r2_interactions) <- sapply(pc_pairs, paste, collapse = ":")

# Results
#best_interaction <- names(which.max(r2_interactions))
#best_r2          <- max(r2_interactions)

#cat("Best PC interaction:", best_interaction, "\n")
#cat("R-squared:", round(best_r2, 4), "\n")

#print(head(sort(r2_interactions, decreasing = TRUE), 10))

r2_per_pc <- sapply(1:(ncol(pca_plot_df)-6), function(i) {
  m <- lm(pca_plot_df$FT_norm ~ pca_plot_df[, i])
  summary(m)$r.squared
})

names(r2_per_pc) <- colnames(pca_plot_df)[1:length(r2_per_pc)]
r2_sorted <- sort(r2_per_pc, decreasing = TRUE)
head(r2_sorted, 20)

# FT: 3
# HT: 6

# Suitability: 260, 182
PCX <- 1
PCY <- 2
for (PCX in c(1, 3)) {
  for (PCY in c(2, 6)) {
    plot_NAM_family(pca_plot_df, VAF, pcx=PCX, pcy=PCY)
    ggplot2::ggsave(filename = file.path(output_dir, paste0("PCA_fam_PC", PCX, "_", PCY, ".jpg")),
                    device = "jpg",
                    width=3.5,
                    height=2.5,
                    dpi=600)
    ggplot2::ggsave(filename = file.path(output_dir, paste0("PCA_fam_PC", PCX, "_", PCY, ".pdf")),
                    device = "pdf",
                    width=3.5,
                    height=2.5)
    
    ggplot2::ggsave(filename = file.path(output_dir, paste0("PCA_fam_PC", PCX, "_", PCY, ".svg")),
                    device = "svg",
                    width=3.5,
                    height=2.5)
    
    plot_NAM_suit(pca_plot_df, VAF, pcx=PCX, pcy=PCY)
    ggplot2::ggsave(filename = file.path(output_dir, paste0("noisy_landscape_", PCX, "_", PCY, ".jpg")),
                    device = "jpg",
                    width=3.5,
                    height=2.5,
                    dpi=600)
    
    ggplot2::ggsave(filename = file.path(output_dir, paste0("noisy_landscape_", PCX, "_", PCY, ".pdf")),
                    device = "pdf",
                    width=3.5,
                    height=2.5)
    ggplot2::ggsave(filename = file.path(output_dir, paste0("noisy_landscape_", PCX, "_", PCY, ".svg")),
                    device = "svg",
                    width=3.5,
                    height=2.5)
    
    
    smoothed <- generate_landscape(pca_plot_df, PCX, PCY, thetas=c(1:15))
    render_2d_landscape(smoothed, VAF, PCX, PCY)
    ggplot2::ggsave(filename = file.path(output_dir,  paste0("2d_smoothed_landscape_", PCX, "_", PCY, ".jpg")),
                    device = "jpg",
                    width=3.5,
                    height=2.5,
                    dpi=600)
    ggplot2::ggsave(filename = file.path(output_dir,  paste0("2d_smoothed_landscape_", PCX, "_", PCY, ".pdf")),
                    device = "pdf",
                    width=3.5,
                    height=2.5)
    ggplot2::ggsave(filename = file.path(output_dir,  paste0("2d_smoothed_landscape_", PCX, "_", PCY, ".svg")),
                    device = "svg",
                    width=3.5,
                    height=2.5)
    
    pca_plot_df$Suitability <- fields::interp.surface(
      smoothed, 
      cbind(pca_plot_df[,PCX], pca_plot_df[,PCY])
    )
    pca_plot_df$Suitability <- pca_plot_df$Suitability + 0.005
    
    wright <- wright_landscape(smoothed, PCX, PCY, window=3)
    wright
    htmlwidgets::saveWidget(as_widget(wright), file.path(output_dir, "wright.html"))
    
    p <- render_3d_landscape(smoothed, pca_plot_df, plot_families, founders, PCX, PCY)
    p
    htmlwidgets::saveWidget(as_widget(p), file.path(output_dir, paste0("landscape_", PCX, "_", PCY, ".html")))
    
    pi <- render_individuals(smoothed, pca_plot_df, plot_families, founders, PCX, PCY)
    pi
    htmlwidgets::saveWidget(as_widget(pi), file.path(output_dir, paste0("individuals_", PCX, "_", PCY, ".html")))
  }
}

