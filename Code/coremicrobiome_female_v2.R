## Set working directory
setwd("C:/Users/ys_cl/OneDrive - UBC/UBC W25/MICB 475/Project 2/coremicrobiome")

## Load packages
library(phyloseq)
library(microbiome)
library(tidyverse)
library(ggVennDiagram)

## Load phyloseq object
load("../phyloseq/deer_female_phyloseq.RData")

## Convert phyloseq object to relative abundance 
phyloseq_female_RA <- transform_sample_counts(deer_female, fun=function(x) x/sum(x))

## Create phyloseq objects that contain only samples that belong to each pathology group

# Females: 
female_neg <- subset_samples(phyloseq_female_RA, pathology=="Neg")
female_ln_pos <- subset_samples(phyloseq_female_RA, pathology=="LN POS")
female_br_ln_pos <- subset_samples(phyloseq_female_RA, pathology=="BR & LN POS")

## Identify core members in each group
# Set choice of detection (abundance) = 0.001
# Set detection threshold (prevalance) = 0.2
female_neg_ASVs <- core_members(female_neg, detection = 0.001, prevalence = 0.2)
female_ln_pos_ASVs <- core_members(female_ln_pos, detection = 0.001, prevalence = 0.2)
female_br_ln_pos_ASVs <- core_members(female_br_ln_pos, detection = 0.001, prevalence = 0.2)

pathology_list_full <- list(Negative = female_neg_ASVs, LN_POS = female_ln_pos_ASVs, BR_LN_POS = female_br_ln_pos_ASVs)

## Create Venn diagram to visualize similarities between core members in each group
female_venn <- ggVennDiagram(
  x = pathology_list_full,
  label_alpha = 0,      # remove black background behind counts
  label_size = 6        # slightly increase font size of counts
) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(hjust = 0.5)  # centers the title
  ) +
  labs(title = "Core Microbiome Overlap in Female Deer by Pathology")

## Save the plot
ggsave(
  filename = "venn_female_pathology_v2.png",
  plot = female_venn,
  width = 7,
  height = 6,
  dpi = 300
)




