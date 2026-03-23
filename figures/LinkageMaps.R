# Title: LINKAGE MAPS
# Author: Ted Monyak
# Description: This file contains functions for plotting linkage maps
library(ggnewscale)
library(ggplot2)

# Converts a scan1 output from R/qtl2 to a scanone output from R/qtl1
# s1.output: an output from scan1
# map: the result of a 'insert_pseudomarkers' call
# Returns: a dataframe with 2 + # phenotypes columns, with each rowname as the marker
# and 'chr', 'pos', and (phenotypes) columns following
qtl2toqtl1 <- function(s1.output, map) {
  # Flatten the map so it can be joined with the scan1 output
  flat_map <- as.data.frame(map_dfr(seq_along(map), function(i) {
    chr_name <- names(map)[i]
    chr_map <- map[[i]]
    
    tibble(
      marker = names(chr_map),
      chr = chr_name,
      pos = as.numeric(chr_map)
    )
  }))
  
  # Inner join the two dataframes
  s1.output <- merge(s1.output, flat_map, by.x="row.names", by.y="marker", all=FALSE) %>%
    column_to_rownames("Row.names") %>%
    relocate(chr, pos)
}

# Generates a ggplot with the effect sizes of the real QTl overlaid
# RIL: a recombinant inbred line family
# s1.output: the output of a scanone function
# s1.perm: the output of permutation test of scanone
# ril_dir: where to write the resulting plot
# Returns: nothing
plotLinkageMap <- function(RIL, s1.output, s1.perm, ril_dir) {
  
  # Get a dataframe of effect sizes per trait (including fitness)
  eff_sizes <- getQtlEffectSizes(RIL)
  
  # Create a tidy version of the scanone output
  plot_data <- as.data.frame(s1.output) %>%
    mutate(chr = as.numeric(chr)) %>%
    pivot_longer(
      cols = starts_with("pheno"),
      names_to = "phenotype",
      values_to = "lod",
    ) %>%
    mutate(phenotype = case_when(
      phenotype == "pheno.Trait1" ~ "Trait 1",
      phenotype == "pheno.Trait2" ~ "Trait 2",
      phenotype == "pheno.Trait3" ~ "Trait 3",
      phenotype == "pheno.Suitability" ~ "Suitability",
      phenotype == "pheno.W" ~ "W",
      TRUE ~ phenotype
    ))
  
  # Calculate cumulative positions for x-axis
  chr_lengths <- plot_data %>%
    dplyr::mutate(chr = as.numeric(chr)) %>%
    dplyr::arrange(chr) %>%
    dplyr::group_by(chr) %>%
    dplyr::summarize(max_pos = n.genMapLen, .groups = 'drop') %>%
    dplyr::mutate(
      cumul_start = cumsum(c(0, head(max_pos, -1))),
      cumul_mid = cumul_start + max_pos / 2
    )
  
  # Add cumulative positions to plot data
  plot_data <- plot_data %>%
    dplyr::left_join(chr_lengths %>% dplyr::select(chr, cumul_start), by = "chr") %>%
    dplyr::mutate(cumul_pos = pos + cumul_start)
  
  # Add cumulative positions to effect sizes
  eff_data <- eff_sizes %>%
    dplyr::left_join(chr_lengths %>% dplyr::select(chr, cumul_start), by = "chr") %>%
    dplyr::mutate(
      cumul_pos = pos + cumul_start,
      phenotype = trait  # Rename trait to phenotype for matching
    ) %>%
    dplyr::select(chr, pos, cumul_start, cumul_pos, phenotype, eff_size)
  
  # Calculate LOD range
  lod_min <- min(plot_data$lod, na.rm = TRUE)
  lod_max <- max(plot_data$lod, na.rm = TRUE)
  lod_range <- lod_max - lod_min
  
  # Calculate the effect size range
  if (nrow(eff_data) == 0) {
    eff_min <- 0
    eff_range <- 1
  } else {
    eff_min <- min(eff_data$eff_size, na.rm = TRUE)
    eff_max <- max(eff_data$eff_size, na.rm = TRUE)
    eff_range <- eff_max - eff_min
  }

  # Calculate a scaling multiple for the effect sizes
  scale_factor <- lod_max / eff_max
  
  # Apply global scaling to all effect sizes
  eff_data <- eff_data %>%
    dplyr::mutate(scaled_effect = eff_size * scale_factor)
  
  # Determine the significance thresholds for each trait
  perm_thresholds <- summary(s1.perm, alpha=0.05)
  thresholds <- data.frame(
    phenotype = c("Trait 1",
                  "Trait 2",
                  "Trait 3",
                  "Suitability",
                  "W"),
    threshold = c(perm_thresholds[1,1], perm_thresholds[1,2], perm_thresholds[1,3], perm_thresholds[1,4], perm_thresholds[1,5])
  )
  
  facet_labels <- c(
    "Trait 1"    = "Attained Trait 1 (Oligogenic; e.g. Maturity)",
    "Trait 2"    = "Attained Trait 2 (Oligogenic; e.g. Height)",
    "Trait 3"    = "Desired Trait (Polygenic; e.g. Yield)",
    "Suitability" = "Suitability",
    "W" = "Breeding Fitness"
  )
  
  # Builds a linkage map plot for a given phenotype
  # pheno: Name of the phenotype
  # eff_phenos: QTL for each of these phenotypes will be shown
  # show_eff: Whether to show the effect size label
  # show_qtl_guide: Whether to show the QTL legend
  # show_x: Whether to show the x-axis
  make_lod_plot <- function(pheno, eff_phenos=c(), show_eff=TRUE, show_qtl_guide=TRUE, show_x=FALSE) {
    # filter by phenotype
    ggplot(plot_data[plot_data$phenotype == pheno, ],
           aes(x = cumul_pos, y = lod, color = factor(chr), group=chr)) +
      geom_line(linewidth = 0.7) +
      scale_color_manual(values = rep(c("black", "grey"),
                                      length.out = nrow(chr_lengths)),
                                      guide="none") +
      new_scale_color() +
      
      # only print the QTL from eff_phenos
      geom_segment(data = eff_data[eff_data$phenotype %in% eff_phenos, ], 
                   aes(x = cumul_pos, xend = cumul_pos, y = 0, yend = scaled_effect), 
                   color = "gray50", alpha=0.8, linewidth = 0.5, inherit.aes = FALSE) +
      
      # only print the QTL from eff_phenos
      geom_point(data = eff_data[eff_data$phenotype %in% eff_phenos, ],
                 aes(x = cumul_pos, y = scaled_effect, color=phenotype),
                 shape=18, size = 2.5, alpha=0.8, inherit.aes = FALSE) +
      scale_color_manual(values=c("purple", "limegreen"),
                         name="QTL",
                         labels=c("QTL", "QTL"),
                         limits=c("Trait 1", "Trait 2"),
                         guide = if (show_qtl_guide) guide_legend() else "none") +
      
      # The significance threshold
      geom_hline(data = thresholds[thresholds$phenotype == pheno, ],
                 aes(yintercept = threshold),
                 linetype = "dashed", color = "gray30", linewidth = 0.5) +
      scale_x_continuous(
        breaks = chr_lengths$cumul_mid,
        labels = chr_lengths$chr,
        expand = c(0.01, 0)
      ) +
      scale_y_continuous(
        name = "LOD",
        expand=expansion(mult=c(0,0.2)),
        sec.axis = sec_axis(
          transform = ~ . / scale_factor,
          name = "Effect Size"
        )
      ) +
      coord_cartesian(ylim = c(lod_min, lod_max), clip = "off") +
      labs(x = "Genome Position", title = facet_labels[pheno]) +
      theme_minimal(base_family="Helvetica", base_size=10) +
      theme(
        panel.grid.major.y = element_line(color = "gray90", linewidth = 0.3),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.title.x = if (show_x) element_text() else element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", color = "black"),
        axis.title.y.left = element_text(margin = margin(r = 10)),
        axis.title.y.right = if (show_eff) element_text(margin = margin(l = 10), size=10) else element_blank(),
        axis.text.y.right = if (show_eff) element_text() else element_blank(),
        axis.ticks.y.right = if (!show_eff) element_blank(),
        strip.background = element_rect(fill = "gray95", color = "gray60"),
        legend.position = "none"
      )
  }
  
  # Build each linkage map plot
  p1 <- make_lod_plot("Trait 1", eff_phenos = "Trait 1", show_eff=TRUE, show_qtl_guide=FALSE, show_x=FALSE)
  p2 <- make_lod_plot("Trait 2", eff_phenos = "Trait 2", show_eff=TRUE, show_qtl_guide=FALSE, show_x=FALSE)
  pS <- make_lod_plot("Suitability", eff_phenos = c("Trait 1", "Trait 2"), show_eff=TRUE, show_qtl_guide=FALSE, show_x=FALSE)
  p3 <- make_lod_plot("Trait 3", show_eff=FALSE, show_qtl_guide=FALSE, show_x=FALSE)
  pW <- make_lod_plot("W", eff_phenos = c("Trait 1", "Trait 2"), show_eff=TRUE, show_qtl_guide=TRUE, show_x=TRUE)
  p <- (p1 / p2 / pS / p3 / pW) +
    plot_layout(guides="collect") &
    theme(legend.position="bottom",
          plot.margin = margin(t = 0, r = 5, b = 0, l = 5) ) & 
    plot_annotation(tag_levels='a')
  p
  ggplot2::ggsave(filename = file.path(ril_dir, "linkagemap.jpg"),
                  device = "jpg",
                  width=6,
                  height=8,
                  dpi=600)
  
  ggplot2::ggsave(filename = file.path(ril_dir, "linkagemap.pdf"),
                  device = "pdf",
                  width=6,
                  height=8)

}

# Generates a linkage map from a scantwo output, with effect sizes of actual QTL overlaid
# RIL: the population from which to get effect sizes
# s2.output: the result of a scantwo() call (r/qtl)
# ril_dir: the directory to write the result
# Returns: nothing
plot2DLinkageMap <- function(RIL, s2.output, ril_dir) {
  # Extract the interaction LOD scores matrix
  lod_int <- s2.output$lod
  
  # The interaction effects are calculated by subtracting the upper triangle (additive model)
  # from the lower triangle (full model)
  lod_int[lower.tri(lod_int)] <- lod_int[lower.tri(lod_int)] - t(lod_int)[lower.tri(lod_int)]
  
  # Get map information
  map <- as.data.frame(s2.output$map) %>%
    dplyr::select(chr, pos) %>%
    dplyr::arrange(chr, pos) %>%
    dplyr::group_by(chr) %>%
    dplyr::mutate(idx = row_number()) %>%
    dplyr::ungroup()
  
  # Find the marker closest to a QTL snp location
  find_closest_idx <- function(qtl_chr, qtl_pos) {
    chr_mars <- map[map$chr == qtl_chr, ]
    chr_mars$idx[which.min(abs(chr_mars$pos - qtl_pos))]
  }
  
  # Get a dataframe of effect sizes per trait
  eff_sizes <- getQtlEffectSizes(RIL) %>%
    dplyr::filter(trait %in% c("Trait 1", "Trait 2")) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(idx = find_closest_idx(chr, pos)) %>%
    dplyr::ungroup()
  
  # Get the number of markers per chromosome
  mar_per_chr <- map %>%
    dplyr::group_by(chr) %>%
    dplyr::summarize(n_markers = n())
  
  # Total number of markers
  nmar <- nrow(map)
  # The length of the long dataframe representing each pairwise interaction
  len <- tail(cumsum(seq(nmar)-1),1)
  lod_scores <- data.frame(chr1=numeric(len),
                           pos1=numeric(len),
                           idx1=numeric(len),
                           chr2=numeric(len),
                           pos2=numeric(len),
                           idx2=numeric(len),
                           lod=numeric(len))
  # Create a data frame with all pairwise positions in the lower triangle
  # Lower triangle stores the interaction LOD scores
  idx <- 1
  for (i in 1:nmar) {
    for (j in 1:nmar) {
      if (i > j) {  # Lower triangle only
        lod_scores$chr1[idx] <- map$chr[i]
        lod_scores$pos1[idx] <- map$pos[i]
        lod_scores$idx1[idx] <- map$idx[i]
        lod_scores$chr2[idx] <- map$chr[j]
        lod_scores$pos2[idx] <- map$pos[j]
        lod_scores$idx2[idx] <- map$idx[j]
        lod_scores$lod[idx] <- lod_int[i,j]
        idx <- idx + 1
      }
    }
  }

  # Get unique chromosomes and calculate cumulative positions
  chr_info <- data.frame(chr = seq(1:n.chr)) %>%
    dplyr::mutate(
      chr_length = mar_per_chr$n_markers,
      cumul_start = cumsum(c(0, head(chr_length, -1))),
      cumul_mid = cumul_start + chr_length/ 2
    )
  
  # Add cumulative positions to interaction data
  lod_scores <- lod_scores %>%
    dplyr::left_join(chr_info %>% dplyr::select(chr, cumul_start), by = c("chr1" = "chr")) %>%
    dplyr::rename(cumul_start1 = cumul_start) %>%
    dplyr::left_join(chr_info %>% dplyr::select(chr, cumul_start), by = c("chr2" = "chr")) %>%
    dplyr::rename(cumul_start2 = cumul_start) %>%
    dplyr::mutate(
      cumul_pos1 = idx1 + cumul_start1,
      cumul_pos2 = idx2 + cumul_start2
    )
  
  # Add cumulative positions to effect sizes
  eff_sizes <- eff_sizes %>%
    dplyr::left_join(chr_info %>% dplyr::select(chr, cumul_start), by = "chr") %>%
    dplyr::mutate(cumul_pos = idx + cumul_start)
  
  # Get the range of the plot
  plot_min <- 0
  plot_max <- nmar
  # Create the plot
  ggplot() +
    # Add heatmap of interaction LOD scores
    geom_tile(data = lod_scores, 
              aes(x = cumul_pos1, y=cumul_pos2, fill = lod),
              alpha = 1.0) +
    
    # Add white polygon to cover upper triangle
    geom_polygon(data = data.frame(
      x = c(plot_min, plot_min, plot_max),
      y = c(plot_min, plot_max, plot_max)
    ), aes(x = x, y = y), fill = "white", color = NA) +
    
    # Add chromosome boundary lines (vertical) - only in lower triangle
    geom_segment(data = chr_info[-1, ],
                 aes(x = cumul_start, xend = cumul_start,
                     y = plot_min, yend = cumul_start),
                 color = "white", linewidth = 0.3, alpha = 0.5) +
    
    # Add chromosome boundary lines (horizontal) - only in lower triangle
    geom_segment(data = chr_info[-1, ],
                 aes(x = cumul_start, xend = plot_max,
                     y = cumul_start, yend = cumul_start),
                 color = "white", linewidth = 0.3, alpha = 0.5) +
    
    # Add vertical lines for QTL (color by trait, width by effect size)
    geom_segment(data = eff_sizes,
                 aes(x = cumul_pos, xend = cumul_pos, 
                     y = plot_min, yend = cumul_pos,
                     color = trait,
                     linewidth = abs(eff_size)),
                 alpha = 0.3) +
    
    # Add horizontal lines for QTL (color by trait, width by effect size)
    geom_segment(data = eff_sizes,
                 aes(x = cumul_pos, xend = plot_max,
                     y = cumul_pos, yend = cumul_pos,
                     color = trait,
                     linewidth = abs(eff_size)),
                 alpha = 0.3) +
    
    # Color scale for lines (traits)
    scale_color_manual(values=c("purple", "limegreen"),
                       name="QTL",
                       labels=c("Attained Trait 1",
                                "Attained Trait 2")) +
    
    # Fill scale for heatmap (LOD scores)
    scale_fill_viridis_c(option="inferno",
                         name = "Interaction\nLOD") +
    
    # QTLs across x-axis
    geom_point(data = eff_sizes,
               aes(x = cumul_pos, y = plot_min-10, color=trait),
               size = 3, shape = 18,
               inherit.aes = FALSE) +
    
    # QTLs across y-axis
    geom_point(data = eff_sizes,
               aes(x = plot_max+10, y = cumul_pos, color=trait),
               size = 3, shape = 18,
               inherit.aes = FALSE) +
    
    # Size scale for line thickness
    scale_linewidth_continuous(range = c(0.1, 1), name = "Effect Size") +
    
    # Add chromosome labels
    scale_x_continuous(
      breaks = chr_info$cumul_mid,
      labels = chr_info$chr,
      expand = c(0.01, 0),
      position = "bottom",
    ) +
    scale_y_continuous(
      breaks = chr_info$cumul_mid,
      labels = chr_info$chr,
      expand = c(0.01, 0),
      position = "right",
    ) +
    labs(
      x = "Genome Position",
      y = "Genome Position"
    ) +
    theme_minimal(base_size = 14,
                  base_family="Helvetica") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
      legend.position = "right",
      legend.key = element_rect(fill = "white", color="white"),
      legend.text = element_text(size=16),
      legend.title = element_text(size=16),
      aspect.ratio = 1,
      axis.text.x = element_text(angle=0, hjust=1),
      axis.text.y = element_text(angle=45, hjust=1),
      axis.title.y = element_text(angle=270),
      axis.ticks.y = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid=element_blank(),
    ) +
    coord_fixed()  # Equal scaling on both axes

  ggplot2::ggsave(filename = file.path(ril_dir, "2D_linkagemap.jpg"),
                  device = "jpg",
                  width=9,
                  height=8.5,
                  dpi=600)
  
  ggplot2::ggsave(filename = file.path(ril_dir, "2D_linkagemap.pdf"),
                  device = "pdf",
                  width=9,
                  height=8.5)
}
