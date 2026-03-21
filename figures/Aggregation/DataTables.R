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

non_paired_cols <- names(res.df)[!grepl("_T[12]$", names(res.df))]

cor.df <- res.df %>%
  dplyr::filter(type=="Admixed") %>%
  dplyr::select(isoElite_T1, isoElite_T2, isoElite_Att, hamm_T1, hamm_T2, hamm_Att, fst,
                ev_T1, ev_T2, ev_Suit, ev_W, ev_T3,
                nLod_T1, nLod_T2, nLod_Suit, nLod_W, nLod_Int, nLod_T3) %>%
                #emergentLod_W, gwpR_T1, gwpR_T2, gwpR_W) %>%
  dplyr::mutate(.row_id = row_number()) %>%
  pivot_longer(
    cols = matches("_T[12]$"),
    names_to = c(".value", "pair_num"),
    names_pattern = "^(.+)_T([12])$"
  ) %>%
  mutate(across(
    .cols = any_of(non_paired_cols),
    .fns = ~ if_else(pair_num == "2", NA, .x)
  )) %>%
  dplyr::select(-.row_id, -pair_num) %>%
  dplyr::rename("Isoeliteness: Attained Trait"=isoElite,
         "Isoeliteness: Mean"=isoElite_Att,
         "Hamming Distance: Attained Traits"=hamm,
         "Hamming Distance: Mean"=hamm_Att,
         "FST"=fst,
         "Excess Variance: Attained Traits"=ev,
         "Excess Variance: Suitability"=ev_Suit,
         "Excess Variance: Desired Trait"=ev_T3,
         "Excess Variance: Breeding Fitness"=ev_W,
         "LOD Peaks: Attained Traits"=nLod,
         "LOD Peaks: Suitability"=nLod_Suit,
         "LOD Peaks: Desired Trait"=nLod_T3,
         "LOD Peaks: Breeding Fitness"=nLod_W,
         "LOD Peaks: Fitness Interactions"=nLod_Int)
         #"Interaction LOD Peaks: Breeding Fitness"=emergentLod_W)
         #"GWP Accuracy: Attained Traits"=gwpR,
         #"GWP Accuracy: Breeding Fitness"=gwpR_W)

cor.df %>% 
  filter(!complete.cases(.)) %>%
  summarise(across(everything(), ~sum(is.na(.))))
# --- Compute correlations ---
cor_res <- rcorr(as.matrix(cor.df))
r_mat   <- cor_res$r
p_mat   <- cor_res$P

# --- Plot ---
jpeg(file.path(output_dir, "cor.jpg"), width = 7, height = 5, units = "in", res = 600)
corrplot(r_mat,
         method       = "square",
         type         = "lower",
         order        = "AOE",
         hclust.method = "complete",
         p.mat        = p_mat,
         sig.level    = 0.05,
         insig        = "blank",
         tl.cex       = 0.6,
         tl.col       = "black",
         tl.srt       = 45,
         col          = COL2("RdBu"),
         diag         = FALSE, 
         addCoef.col  = "black",
         number.cex   = 0.4,
         cl.pos="n"
)
dev.off()
pdf(file.path(output_dir, "cor.pdf"), width = 7, height = 5)
corrplot(r_mat,
         method       = "square",
         type         = "lower",
         order        = "AOE",
         hclust.method = "complete",
         p.mat        = p_mat,
         sig.level    = 0.05,
         insig        = "blank",
         tl.cex       = 0.6,
         tl.col       = "black",
         tl.srt       = 45,
         col          = COL2("RdBu"),
         diag         = FALSE, 
         addCoef.col  = "black",
         number.cex   = 0.4,
         cl.pos="n"
)
dev.off()

res.df %>%
  dplyr::filter(type=="Admixed") %>%
  dplyr::group_by(qtl) %>%
  dplyr::summarize(meanInitA=mean(initA, na.rm=TRUE),
                   meanInitVar=mean(initVar, na.rm=TRUE),
                   meanRilA=mean(rilA, na.rm=TRUE),
                   meanRilVar=mean(rilVar, na.rm=TRUE),
                   meanRank=mean(relRankMean, na.rm=TRUE)) %>%
  write.csv(file.path(output_dir, "aggregatedTraitArchData.csv"),quote=FALSE)

summary_table <- res.df %>%
  dplyr::group_by(qtl,type) %>%
  dplyr::mutate(gwpR_W=mean(c(gwpR_W, gwpR_W_Pop2), na.rm=TRUE)) %>%
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
                   #meanGWPR=mean(gwpR_W, na.rm=TRUE))
summary_table %>%
  write.csv(file.path(output_dir, "aggregatedQtlData.csv"),quote=FALSE)

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