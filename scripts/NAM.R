# Title: NAM_PCA
# Author: Ted Monyak
# Description: This script uses the sorghum NAM resource to generate an
# empirical genotype-to-fitness landscape

library(ggplot2)
library(dplyr)
library(fields)
library(RColorBrewer)
library(parallel)
library(plotly)
library(stringr)
library(tibble)
library(tidyr)
library(viridis)

setwd("~/Documents/CSU/FitnessLandscapes")
output_dir <- "~/Documents/CSU/FitnessLandscapes/output/NAM/pca_sap"
save_dir <- file.path(getwd(), "output/NAM/pca_sap")

source("figures/NAM_Plots.R")
source("figures/G2FLandscapes.R")

# Load the phenotype data, if not already loaded
phenotypeData <- function() {
  # Check if the file already exists
  if (file.exists("NAM/pheno.rds")) {
    return (readRDS("NAM/pheno.rds"))
  }
  
  # Contains the LIB information for each NAM sample
  codes.df <- read.csv("NAM/NAM_population_codes.csv") %>%
    dplyr::rename(Taxa="Plot_ID")
  
  # Append the LIBs of the founder lines
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
  
  # Read in the plant height data from Manhattan 2015
  height.df <- read.csv("NAM/Plant_height_NAM.csv") %>%
    dplyr::filter(Location=="Manhattan") %>%
    dplyr::rename("PH"=Plant_height) %>%
    dplyr::select(Taxa,PH) %>%
    dplyr::group_by(Taxa) %>%
    dplyr::summarize(PH=mean(PH, na.rm=TRUE),
                     .groups='drop') %>%
    inner_join(codes.df) %>%
    dplyr::select(-c(Taxa, Pedigree))
  
  # Read in the flowering time data from Manhattan 2015
  ft.df <- read.csv("NAM/Flowering_time_NAM_2015.csv") %>%
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
  
  # Join the height and flowering time dataframes
  pheno.df <- inner_join(height.df, ft.df)
  
  saveRDS(pheno.df, "NAM/pheno.df")
  return (pheno.df)
}

# Phenotype data
pheno.df <- phenotypeData()

# PCA data
# pca_df <- readRDS("data/pca_sap/pca_df.rds")
VAF <- readRDS("NAM/pca_sap/vaf.rds")
full_pca <- readRDS("NAM/pca_sap/full_pca_projection.rds")
pca_df <- rownames_to_column(as.data.frame(full_pca), var="LIB")

# Metadata
metadata <- read.csv("NAM/metadata_NAM.txt", header = T,sep = "\t")
metadata <- metadata %>% dplyr::select(LIB, Family)

# Founder line metadata
founders <- read.csv("NAM/founders.csv", header=T)
metadata <- rbind(metadata, founders)

# Filter in case there are missing samples
metadata <- metadata[metadata$LIB %in% pca_df$LIB, ]
metadata <- metadata[match(pca_df$LIB, metadata$LIB), ]

# Merge PCA results with metadata
full_pca_plot_df<- left_join(pca_df, metadata)

# DROP IZPY and IZPZ (duplicates)
full_pca_plot_df <- full_pca_plot_df[!(full_pca_plot_df$LIB %in% c("IZPY", "IZPZ")),]

# Use the same color scheme from Bouchet et al. 2017
family_colors <- data.frame(
  Family=c("Macia", "P898012", "Segeolane", "SC35", "SC283", "SC265", "SC1345", "Ajabsido", "SC971", "SC1103", "RTx430"),
  family_color=c("#da0001", "#70cd30", "#aacc4b", "#902d27", "#e9a425", "#9dcce9", "#8c2eea", "#1b1bfb", "#df7d71", "#c5b192", "black")
)
# Duplicate the colors for the RIL family samples
RIL_family_names <- paste0(family_colors$Family,"_RIL")
family_colors <- rbind(family_colors,
                       data.frame(Family=RIL_family_names,
                                  family_color=family_colors$family_color))

# Join the color code with metadata
full_pca_plot_df <- left_join(full_pca_plot_df, family_colors, by="Family")

# Only retain the phenotypes for which there are genotypes
pheno.df <- pheno.df[pheno.df$LIB %in% full_pca_plot_df$LIB,]

# Join PCA data with phenotype data
full_pca_plot_df <- inner_join(full_pca_plot_df, pheno.df)

# Get the mean phenotypes for RTx430
RTx430.df <- full_pca_plot_df %>%
  dplyr::group_by(Family) %>%
  dplyr::summarize(ft=mean(FT),
                   ph=mean(PH)) %>%
  dplyr::filter(Family == "RTx430")

# Get the standard deviation of each trait across all families
FT_s <- sd(full_pca_plot_df$FT)
PH_s <- sd(full_pca_plot_df$PH)

# Divide by this many standard deviations
m <- 2
# Normalize the data by subtracting the phenotypes of RTx430 and
# dividing by 2*standard deviation of the original trait value
full_pca_plot_df <- full_pca_plot_df %>%
  dplyr::mutate(
    FT_norm = (RTx430.df$ft[1] - FT)/(FT_s*m),
    PH_norm = (RTx430.df$ph[1] - PH)/(PH_s*m),
    Suit=exp(-(FT_norm^2 + PH_norm^2) /2)
  )
# Uncomment to plot all NAM families
# plot_families <- unique(full_pca_plot_df$Family)

# Select just these families to plot
plot_families <- c("SC35_RIL", "SC283_RIL", "Macia_RIL", "SC265_RIL", "Ajabsido_RIL",
                   "SC35", "SC283", "Macia", "SC265", "Ajabsido", "RTx430")

# Filter data by the selected families
pca_plot_df <- full_pca_plot_df %>%
  dplyr::filter(Family %in% plot_families)

# Plot the distribution of phenotypes
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

# Remove the phenotype data from the dataframe
pca_plot_df <- pca_plot_df %>%
  dplyr::select(-c(LIB, PH, FT))

# Plot the distribution of suitability values
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


# Determine the PCs that explain the most variation for each trait
# phenotype: one of 'FT_norm', 'PH_norm', or 'Suit'
# Returns: a named vector of PCs and R-squared
getTopPC <- function(phenotype) {
  # The last 6 columns are not PCs
  r2_per_pc <- sapply(1:(ncol(pca_plot_df)-6), function(i) {
    # Fit a linear model to the phenotype
    model <- lm(pca_plot_df[[phenotype]] ~ pca_plot_df[, i])
    # Get the r-squared of that model
    summary(model)$r.squared
  })
  
  # Get the names of each PC, ranked by r-squared, and return the top 10
  names(r2_per_pc) <- colnames(pca_plot_df)[1:length(r2_per_pc)]
  r2_sorted <- sort(r2_per_pc, decreasing = TRUE)
  return (head(r2_sorted,10))
}

getTopPC("FT_norm")
getTopPC("PH_norm")
# Results: FT: PC3, PH: PC6

# Iterate through all 4 combinations of PCs
for (PCX in c(1, 3)) {
  for (PCY in c(2, 6)) {
    # Plot the PCA, color-coded by family
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
    
    # Plot the PCA, color-coded by suitability
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
    
    # Render each of the individuals, not on a surface, in space according to suitability
    pi <- render_individuals(pca_plot_df, plot_families, founders, PCX, PCY)
    pi
    htmlwidgets::saveWidget(as_widget(pi), file.path(output_dir, paste0("individuals_", PCX, "_", PCY, ".html")))
    
    # Create a smoothed version of the landscape with 14 iterations of a gaussian filter
    smoothed <- generate_landscape(pca_plot_df, PCX, PCY, thetas=c(1:14))
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
    
    # Interpolate a surface onto PCX and PCY
    pca_plot_df$Suitability <- fields::interp.surface(
      smoothed, 
      cbind(pca_plot_df[,PCX], pca_plot_df[,PCY])
    )
    # Add a slight offset to get the dots to render above the surface
    pca_plot_df$Suitability <- pca_plot_df$Suitability + 0.005
    
    # Recreate Wright's fitness landscape
    wright <- wright_landscape(smoothed, PCX, PCY, window=3)
    wright
    htmlwidgets::saveWidget(as_widget(wright), file.path(output_dir, "wright.html"))
    
    # Render the 3d fitness landscape with individuals overlaid
    p <- render_3d_landscape(smoothed, pca_plot_df, plot_families, founders, PCX, PCY)
    p
    htmlwidgets::saveWidget(as_widget(p), file.path(output_dir, paste0("landscape_", PCX, "_", PCY, ".html")))
    
  }
}

