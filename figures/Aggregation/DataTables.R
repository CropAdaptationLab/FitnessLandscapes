# Title: DATA TABLES
# Author: Ted Monyak
# Description: This script, which shold be called only from AggregateFigures.R,
# constructs a correlation matrix and summary tables of the metrics produced in
# QtlMonteCarlo.R

library(corrplot)
library(Hmisc)

theme <- theme_minimal(base_size = 12,
                       base_family="Helvetica") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = "black"),
    plot.margin= unit(c(1,1,1,1), unit="pt"),
    aspect.ratio = 1)

col_order <- c("Mean IE",
               "EV: Attained Traits",
               "EV: Breeding Fitness",
               "LOD Peaks: Attained Traits",
               "LOD Peaks: Breeding Fitness",
               "LOD Peaks: Desired Trait",
               "FST")

cor.df <- res.df %>%
  dplyr::filter(type == "Admixed") %>%
  dplyr::select(qtl, isoElite_Att, fst,
                ev_T1, ev_W,
                nLod_T1, nLod_T3, nLod_W) %>%
  dplyr::rename(
    "Mean IE" = isoElite_Att,
    "FST" = fst,
    "EV: Attained Traits"= ev_T1,
    "EV: Breeding Fitness"= ev_W,
    "LOD Peaks: Attained Traits" = nLod_T1,
    "LOD Peaks: Desired Trait" = nLod_T3,
    "LOD Peaks: Breeding Fitness"= nLod_W
  ) %>%
  dplyr::select(qtl, all_of(col_order))

qtl_vals <- c(10, 20, 50)

n_vars <- length(col_order)
n_qtls <- length(qtl_vals)

# Row order
row_names <- as.vector(outer(qtl_vals, col_order,
                           FUN = function(q, v) paste0(v, " L = ", q)))

# Store correlations and p-values
r_asym <- matrix(NA, nrow = n_vars * n_qtls, ncol = n_vars,
                 dimnames = list(row_names, col_order))
p_asym <- matrix(NA, nrow = n_vars * n_qtls, ncol = n_vars,
                 dimnames = list(row_names, col_order))

# Calculate a correlation matrix for each value of qtl_vals
for (qi in seq_along(qtl_vals)) {
  q <- qtl_vals[qi]
  # Subset the correlation matrix
  sub <- cor.df %>%
    dplyr::filter(qtl == q) %>%
    dplyr::select(all_of(col_order))
  
  # Calculate the correlation
  res <- rcorr(as.matrix(sub))
  diag(res$P) <- 1
  
  # For each column, populate the correlation and the p-value
  for (vi in seq_along(col_order)) {
    row_name <- paste0(col_order[vi], " L = ", q)
    r_asym[row_name, ] <- res$r[col_order[vi], col_order]
    p_asym[row_name, ] <- res$P[col_order[vi], col_order]
  }
}

# Separator lines between QTL blocks after row 3, 6, 9, 12, 15, 18
sep_rows <- seq(n_qtls, n_vars * n_qtls - n_qtls, by = n_qtls)  

# Plot an asymmetric correlation matrix
plot_cor_asym <- function() {
  corrplot(r_asym,
           method = "square",
           type = "full",
           is.corr = TRUE,
           p.mat = p_asym,
           sig.level = 0.05,
           insig = "blank",
           tl.cex = 0.55,
           tl.col = "black",
           tl.srt = 45,
           col = COL2("RdBu"),
           diag = TRUE,
           addCoef.col = "black",
           number.cex = 0.35,
           cl.pos = "n")

  # Draw horizontal lines separating QTL blocks
  # corrplot draws rows top-to-bottom; abline coords are in matrix row space
  abline(h = nrow(r_asym) - sep_rows + 0.5, col = "grey40", lwd = 1.2, lty = 2)
}

jpeg(file.path(output_dir, "cor.jpg"), width = 4, height = 5,
     units = "in", res = 600)
plot_cor_asym()
dev.off()

pdf(file.path(output_dir, "cor.pdf"), width = 4, height = 5)
plot_cor_asym()
dev.off()

# Determine the mean trait architecture values for the admixed RIL family
# Write a summary table
res.df %>%
  dplyr::filter(type=="Admixed") %>%
  dplyr::group_by(qtl) %>%
  dplyr::summarize(meanInitA=mean(initA, na.rm=TRUE),
                   meanInitVar=mean(initVar, na.rm=TRUE),
                   meanRilA=mean(rilA, na.rm=TRUE),
                   meanRilVar=mean(rilVar, na.rm=TRUE),
                   meanRank=mean(relRankMean, na.rm=TRUE)) %>%
  write.csv(file.path(output_dir, "aggregatedTraitArchData.csv"),quote=FALSE)

# Determine the mean linkage mapping and genetic divergence metrics
# Write a summary table
summary_table <- res.df %>%
  dplyr::group_by(qtl,type) %>%
  dplyr::summarize(meanIsoelitenessAtt=mean(isoElite_Att, na.rm=TRUE),
                   meanIsoelitenessT3=mean(isoElite_T3, na.rm=TRUE),
                   meanExcessVarAtt=mean(ev_Att, na.rm=TRUE),
                   meanExcessVarSuit=mean(ev_Suit, na.rm=TRUE),
                   meanExcessVarT3=mean(ev_T3, na.rm=TRUE),
                   meanExcessVarW=mean(ev_W, na.rm=TRUE),
                   meanLodPeaksAtt=mean(nLod_Att, na.rm=TRUE),
                   meanLodPeaksSuit=mean(nLod_Suit, na.rm=TRUE),
                   meanLodPeaksT3=mean(nLod_T3, na.rm=TRUE),
                   meanLodPeaksW=mean(nLod_W, na.rm=TRUE),
                   meanLodPeaksInt=mean(nLod_Int, na.rm=TRUE))
summary_table %>%
  write.csv(file.path(output_dir, "aggregatedQtlData.csv"),quote=FALSE)

# Calculate the percentage increase from unadmixed to admixed for each 'mean' metric
pct_diff_table <- summary_table %>%
  tidyr::pivot_longer(cols = starts_with("mean"), 
                      names_to = "metric", 
                      values_to = "value") %>%
  tidyr::pivot_wider(names_from = type, values_from = value) %>%
  dplyr::mutate(pct_diff = (Admixed - Unadmixed) / Unadmixed) %>%
  dplyr::select(-Admixed, -Unadmixed) %>%
  tidyr::pivot_wider(names_from = metric, values_from = pct_diff)
pct_diff_table %>%
  write.csv(file.path(output_dir, "aggregatedQtlPercentDiff.csv"),quote=FALSE)