# Theresa Xiao
# DESeq2_deer_v1
# Aim 3: Differential abundance across prion pathology states (sex-stratified)

# 1. Load packages
library(phyloseq)
library(DESeq2)
library(tidyverse)

# 2. Set working directory & load phyloseq objects
# setwd("/Users/theresaxiao/Desktop/2025W1/475/475_project2")

load("deer_male_phyloseq.RData")    # object: deer_male
load("deer_female_phyloseq.RData")  # object: deer_female

ps_male   <- deer_male
ps_female <- deer_female

# 3. Clean pathology variable

clean_pathology_levels <- function(ps) {
  sd <- data.frame(sample_data(ps))
  
  if (!"pathology" %in% colnames(sd)) {
    stop("No 'pathology' column in sample_data.")
  }
  
  # Normalize spelling
  path <- as.character(sd$pathology)
  path <- trimws(path)
  path[path == "BR &LN POS"] <- "BR & LN POS"
  path[path == "B + LN POS"] <- "BR & LN POS"
  path[path == "LN+ POS"]    <- "LN POS"
  
  # Keep only the 3 main pathology states
  keep_levels <- c("Neg", "LN POS", "BR & LN POS")
  keep <- path %in% keep_levels
  
  ps   <- prune_samples(keep, ps)
  path <- path[keep]
  
  sample_data(ps)$pathology <- factor(path, levels = keep_levels)
  ps
}

ps_male   <- clean_pathology_levels(ps_male)
ps_female <- clean_pathology_levels(ps_female)

# 4. Run DESeq2
run_deseq_pathology <- function(ps,
                                path_var = "pathology",
                                alpha = 0.05,
                                sex_label = NA_character_) {
  
  sd <- data.frame(sample_data(ps))
  sd[[path_var]] <- factor(sd[[path_var]])
  sample_data(ps)[[path_var]] <- sd[[path_var]]
  
  # Convert phyloseq -> DESeqDataSet
  design_formula <- as.formula(paste("~", path_var))
  dds <- phyloseq_to_deseq2(ps, design = design_formula)
  
  # Size factors robust to many zeros
  geoMeans <- apply(counts(dds), 1, function(x) {
    if (all(x == 0)) 0 else exp(mean(log(x[x > 0])))
  })
  dds <- estimateSizeFactors(dds, geoMeans = geoMeans)
  
  # Run DESeq2
  dds <- DESeq(dds)
  
  # All pairwise contrasts
  levs  <- levels(colData(dds)[[path_var]])
  combs <- combn(levs, 2, simplify = FALSE)
  
  all_contrasts <- lapply(combs, function(cn) {
    ref <- cn[1]
    alt <- cn[2]
    
    res <- results(dds,
                   contrast = c(path_var, alt, ref),
                   alpha = alpha) %>%
      as.data.frame() %>%
      rownames_to_column(var = "ASV") %>%
      mutate(
        contrast = paste0(alt, "_vs_", ref),
        path_ref = ref,
        path_alt = alt,
        sex      = sex_label
      )
    
    res
  })
  
  res_all <- bind_rows(all_contrasts)
  
  list(dds = dds, results = res_all)
}

# Run DESeq2 for male and female deer
deseq_male <- run_deseq_pathology(ps_male,
                                  path_var  = "pathology",
                                  alpha     = 0.05,
                                  sex_label = "male")

deseq_female <- run_deseq_pathology(ps_female,
                                    path_var  = "pathology",
                                    alpha     = 0.05,
                                    sex_label = "female")

readLines("df_male_rel_EXTREME_ONLY.csv", n = 10)

# 5. Filter significantly different ASVs
filter_sig_asvs <- function(res_table,
                            contrast_label,
                            lfc_cutoff = 1,
                            padj_cutoff = 0.05) {
  res_table %>%
    filter(
      contrast == contrast_label,
      !is.na(padj),
      padj <= padj_cutoff,
      abs(log2FoldChange) >= lfc_cutoff
    ) %>%
    arrange(padj)
}

# NOTE: contrast labels use "Neg" (not "NEG")
sig_male_LN_vs_NEG <- filter_sig_asvs(deseq_male$results,
                                      contrast_label = "LN POS_vs_Neg")

sig_male_BR_LN_vs_NEG <- filter_sig_asvs(deseq_male$results,
                                         contrast_label = "BR & LN POS_vs_Neg")

sig_female_LN_vs_NEG <- filter_sig_asvs(deseq_female$results,
                                        contrast_label = "LN POS_vs_Neg")

sig_female_BR_LN_vs_NEG <- filter_sig_asvs(deseq_female$results,
                                           contrast_label = "BR & LN POS_vs_Neg")

# 6. Save results tables

write_csv(deseq_male$results,
          "DESeq2_male_pathology_allContrasts.csv")

write_csv(deseq_female$results,
          "DESeq2_female_pathology_allContrasts.csv")

write_csv(sig_male_LN_vs_NEG,
          "DESeq2_male_LN_POS_vs_Neg_sig.csv")

write_csv(sig_male_BR_LN_vs_NEG,
          "DESeq2_male_BR_LN_POS_vs_Neg_sig.csv")

write_csv(sig_female_LN_vs_NEG,
          "DESeq2_female_LN_POS_vs_Neg_sig.csv")

write_csv(sig_female_BR_LN_vs_NEG,
          "DESeq2_female_BR_LN_POS_vs_Neg_sig.csv")

# 6b. Helper: count significant ASVs by direction --------------------------

count_sig_by_direction <- function(res_table,
                                   contrast_label,
                                   lfc_cutoff = 1,
                                   padj_cutoff = 0.05) {
  df <- res_table %>%
    filter(contrast == contrast_label,
           !is.na(padj),
           padj <= padj_cutoff,
           abs(log2FoldChange) >= lfc_cutoff)
  
  tibble(
    contrast    = contrast_label,
    n_total_sig = nrow(df),
    n_up        = sum(df$log2FoldChange > 0),
    n_down      = sum(df$log2FoldChange < 0)
  )
}

# Example: get counts for all main contrasts
sig_counts <- bind_rows(
  count_sig_by_direction(deseq_male$results,   "LN POS_vs_Neg"),
  count_sig_by_direction(deseq_male$results,   "BR & LN POS_vs_Neg"),
  count_sig_by_direction(deseq_female$results, "LN POS_vs_Neg"),
  count_sig_by_direction(deseq_female$results, "BR & LN POS_vs_Neg")
)

# (Optional) save these counts
write_csv(sig_counts, "DESeq2_sig_counts_by_contrast.csv")

# 7. Volcano plots --------------------------------------------------------
# Helper function to make a volcano plot for one contrast
plot_volcano <- function(res_table, contrast_label, plot_title) {
  df <- res_table %>%
    filter(contrast == contrast_label) %>%
    drop_na(padj) %>%
    mutate(sig = padj <= 0.05 & abs(log2FoldChange) >= 1) %>%
    filter(sig)   # <- keep only significant points
  
  ggplot(df, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(color = "red", alpha = 0.8, size = 1.8) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    labs(
      title = plot_title,
      x = "log2 fold change",
      y = "-log10(adjusted p-value)"
    ) +
    theme_minimal()
}

# 7a. Male: LN POS vs Neg
plot_volcano(
  res_table      = deseq_male$results,
  contrast_label = "LN POS_vs_Neg",
  plot_title     = "Male deer: LN POS vs Neg"
)
ggsave("volcano_male_LN_POS_vs_Neg.png", width = 7, height = 5, dpi = 300)

# 7b. Male: BR & LN POS vs Neg
plot_volcano(
  res_table      = deseq_male$results,
  contrast_label = "BR & LN POS_vs_Neg",
  plot_title     = "Male deer: BR & LN POS vs Neg"
)
ggsave("volcano_male_BR_LN_POS_vs_Neg.png", width = 7, height = 5, dpi = 300)

# 7c. Male: BR & LN POS vs LN POS
plot_volcano(
  res_table      = deseq_male$results,
  contrast_label = "BR & LN POS_vs_LN POS",
  plot_title     = "Male deer: BR & LN POS vs LN POS"
)
ggsave("volcano_male_BR_LN_POS_vs_LN_POS.png",   width = 7, height = 5, dpi = 300)

# 7d. Female: LN POS vs Neg
plot_volcano(
  res_table      = deseq_female$results,
  contrast_label = "LN POS_vs_Neg",
  plot_title     = "Female deer: LN POS vs Neg"
)
ggsave("volcano_female_LN_POS_vs_Neg.png", width = 7, height = 5, dpi = 300)

# 7e. Female: BR & LN POS vs Neg
plot_volcano(
  res_table      = deseq_female$results,
  contrast_label = "BR & LN POS_vs_Neg",
  plot_title     = "Female deer: BR & LN POS vs Neg"
)
ggsave("volcano_female_BR_LN_POS_vs_Neg.png", width = 7, height = 5, dpi = 300)

# 7f. Female: BR & LN POS vs LN POS
plot_volcano(
  res_table      = deseq_female$results,
  contrast_label = "BR & LN POS_vs_LN POS",
  plot_title     = "Female deer: BR & LN POS vs LN POS"
)
ggsave("volcano_female_BR_LN_POS_vs_LN_POS.png", width = 7, height = 5, dpi = 300)


# 8. Bar plots of top differentially abundant taxa
# 8a. Attach taxonomy to DESeq results

tax_male <- as.data.frame(tax_table(ps_male))
tax_male$ASV <- rownames(tax_male)

tax_female <- as.data.frame(tax_table(ps_female))
tax_female$ASV <- rownames(tax_female)

# Join taxonomy into DESeq results
res_male_tax <- deseq_male$results %>%
  left_join(tax_male, by = "ASV")

res_female_tax <- deseq_female$results %>%
  left_join(tax_female, by = "ASV")

# Fix NAs, etc in Genus column
fix_genus = function(df){
  new = df %>% 
    mutate(Family = ifelse(is.na(Family) | Family == "" | Family == "f__Incertae_Sedis" | 
                            Family == "Unclassified" | str_detect(Family,'f__UCG') | str_detect(Family,'f__CAG'), 
                          Order, Family)) %>% 
    mutate(Genus = ifelse(is.na(Genus) | Genus == "" | Genus == "g__Incertae_Sedis" | 
                            Genus == "Unclassified" | str_detect(Genus,'g__UCG') | str_detect(Genus,'g__CAG'), 
                           paste(Family,'g__uncl',sep='.'), Genus))
  return(new)
}

res_male_tax_fixed = fix_genus(res_male_tax)
res_female_tax_fixed = fix_genus(res_female_tax)

# 8b. Function to generate bar plots of top taxa -------------------------

phyla_cols <- c(
  "p__Actinomycetota"    = "#E41A1C",
  "p__Bacillota"         = "#4DAF4A",
  "p__Bacteroidota"      = "#377EB8",
  "p__Spirochaetota"     = "#984EA3",
  "p__Verrucomicrobiota" = "#FF7F00"
)

plot_top_taxa_bar <- function(res_with_tax,
                              contrast_label,
                              n_top = 20,
                              lfc_cutoff = 1,
                              padj_cutoff = 0.05,
                              title_prefix = "") {
  
  df <- res_with_tax %>%
    filter(
      contrast == contrast_label,
      !is.na(padj),
      padj <= padj_cutoff,
      abs(log2FoldChange) >= lfc_cutoff
    ) %>%
    arrange(padj) %>%
    slice_head(n = n_top) %>%
    mutate(
      BarID = paste0(Genus, "_", ASV),
      GenusShort = ifelse(
        is.na(Genus) | Genus == "",
        "Unclassified",
        Genus
      )
    )
  
  df = df %>% group_by(GenusShort) %>% 
    mutate(genus_number = row_number()) %>%
    mutate(GenusShort = paste(GenusShort, genus_number, sep = ".")) %>%
    ungroup() %>%
    select(-genus_number)
  
  if (nrow(df) == 0) {
    message("No taxa pass filters for contrast: ", contrast_label)
    return(NULL)
  }
  
  # order bars by LFC
  df$BarID <- factor(df$BarID, levels = df$BarID[order(df$log2FoldChange)])
  axis_labels <- setNames(df$GenusShort, df$BarID)
  
  # 🔹 Zoom the log2FC axis to just cover the data (+ a bit of padding)
  lfc_min <- min(df$log2FoldChange, na.rm = TRUE)
  lfc_max <- max(df$log2FoldChange, na.rm = TRUE)
  pad     <- 0.1 * (lfc_max - lfc_min)
  lfc_lim <- c(lfc_min - pad, lfc_max + pad)
  
  ggplot(df, aes(x = BarID, y = log2FoldChange, fill = Phylum)) +
    geom_col() +
    coord_flip() +
    scale_y_continuous(limits = lfc_lim, expand = expansion(mult = 0.02)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_x_discrete(labels = axis_labels) +
    scale_fill_manual(values = phyla_cols) +   # <- no drop=FALSE
    labs(
      title = paste0(title_prefix, "Top ", n_top, " taxa in ", contrast_label),
      x = "Genus",
      y = "log2 fold change",
      fill = "Phylum"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.y = element_text(size = 11),
      axis.text.x = element_text(size = 11),
      plot.title  = element_text(size = 16, face = "bold"),
      legend.position = "right"
    )
}

# 8c. Generate/Save bar plots for the four key contrasts ----------------------

# Male: LN POS vs Neg
plot_top_taxa_bar(
  res_with_tax   = res_male_tax_fixed,
  contrast_label = "LN POS_vs_Neg",
  title_prefix   = "Male deer – "
)
ggsave("bar_male_LN_POS_vs_Neg_top20.png", width = 11, height = 9, dpi = 300)

# Male: BR & LN POS vs Neg
plot_top_taxa_bar(
  res_with_tax   = res_male_tax,
  contrast_label = "BR & LN POS_vs_Neg",
  title_prefix   = "Male deer – "
)
ggsave("bar_male_BR_LN_POS_vs_Neg_top20.png", width = 11, height = 9, dpi = 300)

# Male: BR & LN POS vs LN POS
plot_top_taxa_bar(
  res_with_tax   = res_male_tax,
  contrast_label = "BR & LN POS_vs_LN POS",
  title_prefix   = "Male deer – "
)
ggsave("bar_male_BR_LN_POS_vs_LN_POS_top20.png", width = 11, height = 9, dpi = 300)

# Female: LN POS vs Neg
plot_top_taxa_bar(
  res_with_tax   = res_female_tax,
  contrast_label = "LN POS_vs_Neg",
  title_prefix   = "Female deer – "
)
ggsave("bar_female_LN_POS_vs_Neg_top20.png", width = 11, height = 9, dpi = 300)

# Female: BR & LN POS vs Neg
plot_top_taxa_bar(
  res_with_tax   = res_female_tax,
  contrast_label = "BR & LN POS_vs_Neg",
  title_prefix   = "Female deer – "
)
ggsave("bar_female_BR_LN_POS_vs_Neg_top20.png", width = 11, height = 9, dpi = 300)

# Female: BR & LN POS vs LN POS
plot_top_taxa_bar(
  res_with_tax   = res_female_tax,
  contrast_label = "BR & LN POS_vs_LN POS",
  title_prefix   = "Female deer – "
)
ggsave("bar_female_BR_LN_POS_vs_LN_POS_top20.png", width = 11, height = 9, dpi = 300)


# 9. Validate DESeq2 results using relative-abundance boxplots --------------

# Helper: take a phyloseq object & return melted relative-abundance dataframe
prep_rel_abund_df <- function(ps) {
  ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
  psmelt(ps_rel)
}

df_male_rel   <- prep_rel_abund_df(ps_male)
df_female_rel <- prep_rel_abund_df(ps_female)

# Directory for saving plots
dir.create("validation_plots", showWarnings = FALSE)

# Helper: get ASVs with huge effect sizes (e.g., |log2FC| ≥ 25)
get_extreme_asvs <- function(deseq_res, threshold = 25) {
  deseq_res %>%
    filter(
      !is.na(log2FoldChange),
      !is.na(padj),
      padj <= 0.05,
      abs(log2FoldChange) >= threshold
    ) %>%
    arrange(desc(abs(log2FoldChange))) %>%
    pull(ASV) %>%
    unique()
}

extreme_female_asvs <- get_extreme_asvs(deseq_female$results, threshold = 25)
extreme_male_asvs   <- get_extreme_asvs(deseq_male$results,   threshold = 25)

# Boxplot function for one ASV
plot_rel_abundance_box <- function(df, asv_id, title_prefix = "") {
  df_asv <- df %>% filter(OTU == asv_id)
  
  if (nrow(df_asv) == 0) {
    message("No rows found for ASV ", asv_id, " – skipping.")
    return(NULL)
  }
  
  ggplot(df_asv, aes(x = pathology, y = Abundance, fill = pathology)) +
    geom_boxplot(alpha = 0.7, outlier.size = 1.5) +
    geom_jitter(width = 0.15, alpha = 0.7, size = 1) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 0.01)) +
    labs(
      title = paste0(title_prefix, "ASV ", asv_id,
                     " – relative abundance by pathology"),
      x = "Pathology",
      y = "Relative abundance (%)"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")
}

# 9a. Female extreme ASVs --------------------------------------------------

for (asv in extreme_female_asvs) {
  p <- plot_rel_abundance_box(df_female_rel, asv, "Female deer – ")
  if (!is.null(p)) {
    print(p)
    ggsave(
      filename = paste0("validation_plots/female_relAbund_ASV_", asv, ".png"),
      plot     = p,
      width    = 7,
      height   = 5,
      dpi      = 300
    )
  }
}

# 9b. Male extreme ASVs ----------------------------------------------------

for (asv in extreme_male_asvs) {
  p <- plot_rel_abundance_box(df_male_rel, asv, "Male deer – ")
  if (!is.null(p)) {
    print(p)
    ggsave(
      filename = paste0("validation_plots/male_relAbund_ASV_", asv, ".png"),
      plot     = p,
      width    = 7,
      height   = 5,
      dpi      = 300
    )
  }
}
