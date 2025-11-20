# Theresa Xiao
# DESeq2_deer_v1
# Aim 3: Differential abundance across prion pathology states (sex-stratified)

# 1. Load packages
install.packages("BiocManager")
BiocManager::install("DESeq2")
library(phyloseq)
library(DESeq2)
library(tidyverse)

# 2. Set working directory & load phyloseq objects
setwd("/Users/theresaxiao/Desktop/2025W1/475/475_project2")

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

# 7. Volcano plots --------------------------------------------------------
# Helper function to make a volcano plot for one contrast
plot_volcano <- function(res_table, contrast_label, plot_title) {
  df <- res_table %>%
    filter(contrast == contrast_label) %>%
    drop_na(padj) %>%
    mutate(sig = padj <= 0.05 & abs(log2FoldChange) >= 1)
  
  ggplot(df, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "red")) +
    labs(
      title = plot_title,
      x = "log2 fold change",
      y = "-log10(adjusted p-value)",
      color = "Significant\n(padj ≤ 0.05 & |LFC| ≥ 1)"
    ) +
    theme_minimal()
}

# 7a. Male: LN POS vs Neg
plot_volcano(
  res_table      = deseq_male$results,
  contrast_label = "LN POS_vs_Neg",
  plot_title     = "Male deer: LN POS vs Neg"
)

# 7b. Male: BR & LN POS vs Neg
plot_volcano(
  res_table      = deseq_male$results,
  contrast_label = "BR & LN POS_vs_Neg",
  plot_title     = "Male deer: BR & LN POS vs Neg"
)

# 7c. Male: BR & LN POS vs LN POS
plot_volcano(
  res_table      = deseq_male$results,
  contrast_label = "BR & LN POS_vs_LN POS",
  plot_title     = "Male deer: BR & LN POS vs LN POS"
)

# 7d. Female: LN POS vs Neg
plot_volcano(
  res_table      = deseq_female$results,
  contrast_label = "LN POS_vs_Neg",
  plot_title     = "Female deer: LN POS vs Neg"
)

# 7e. Female: BR & LN POS vs Neg
plot_volcano(
  res_table      = deseq_female$results,
  contrast_label = "BR & LN POS_vs_Neg",
  plot_title     = "Female deer: BR & LN POS vs Neg"
)

# 7f. Female: BR & LN POS vs LN POS
plot_volcano(
  res_table      = deseq_female$results,
  contrast_label = "BR & LN POS_vs_LN POS",
  plot_title     = "Female deer: BR & LN POS vs LN POS"
)

# Save plots
ggsave("volcano_male_LN_POS_vs_Neg.png", width = 7, height = 5, dpi = 300)
ggsave("volcano_male_BR_LN_POS_vs_Neg.png", width = 7, height = 5, dpi = 300)
ggsave("volcano_male_BR_LN_POS_vs_LN_POS.png",   width = 7, height = 5, dpi = 300)
ggsave("volcano_female_LN_POS_vs_Neg.png", width = 7, height = 5, dpi = 300)
ggsave("volcano_female_BR_LN_POS_vs_Neg.png", width = 7, height = 5, dpi = 300)
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


# 8b. Function to generate bar plots of top taxa -------------------------

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
      GenusLabel = ifelse(
        is.na(Genus) | Genus == "",
        ASV,
        paste0(Genus, " (", ASV, ")")  # ensures uniqueness
      ),
      GenusLabel = factor(GenusLabel, levels = GenusLabel[order(log2FoldChange)])
    )
  
  ggplot(df, aes(x = GenusLabel, y = log2FoldChange, fill = Phylum)) +
    geom_col() +
    coord_flip() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(
      title = paste0(title_prefix, "Top ", n_top, " taxa in ", contrast_label),
      x = "Genus (or ASV ID if unclassified)",
      y = "log2 fold change",
      fill = "Phylum"
    ) +
    theme_minimal()
}


# 8c. Generate bar plots for the four key contrasts ----------------------

# Male: LN POS vs Neg
plot_top_taxa_bar(
  res_with_tax   = res_male_tax,
  contrast_label = "LN POS_vs_Neg",
  title_prefix   = "Male deer – "
)

# Male: BR & LN POS vs Neg
plot_top_taxa_bar(
  res_with_tax   = res_male_tax,
  contrast_label = "BR & LN POS_vs_Neg",
  title_prefix   = "Male deer – "
)

# Male: BR & LN POS vs LN POS
plot_top_taxa_bar(
  res_with_tax   = res_male_tax,
  contrast_label = "BR & LN POS_vs_LN POS",
  title_prefix   = "Male deer – "
)

# Female: LN POS vs Neg
plot_top_taxa_bar(
  res_with_tax   = res_female_tax,
  contrast_label = "LN POS_vs_Neg",
  title_prefix   = "Female deer – "
)

# Female: BR & LN POS vs Neg
plot_top_taxa_bar(
  res_with_tax   = res_female_tax,
  contrast_label = "BR & LN POS_vs_Neg",
  title_prefix   = "Female deer – "
)

# Female: BR & LN POS vs LN POS
plot_top_taxa_bar(
  res_with_tax   = res_female_tax,
  contrast_label = "BR & LN POS_vs_LN POS",
  title_prefix   = "Female deer – "
)

# Save figures
ggsave("bar_male_LN_POS_vs_Neg_top20.png", width = 7, height = 6, dpi = 300)
ggsave("bar_male_BR_LN_POS_vs_Neg_top20.png", width = 7, height = 6, dpi = 300)
ggsave("bar_male_BR_LN_POS_vs_LN_POS_top20.png",   width = 7, height = 6, dpi = 300)
ggsave("bar_female_LN_POS_vs_Neg_top20.png", width = 7, height = 6, dpi = 300)
ggsave("bar_female_BR_LN_POS_vs_Neg_top20.png", width = 7, height = 6, dpi = 300)
ggsave("bar_female_BR_LN_POS_vs_LN_POS_top20.png", width = 7, height = 6, dpi = 300)


