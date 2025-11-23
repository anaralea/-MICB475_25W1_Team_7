library(phyloseq)
library(ape) 
library(tidyverse)
library(vegan)
library(picante)

# Load in the phyloseq object
load("deer_female_phyloseq.RData")
load("deer_male_phyloseq.RData")

# Rarefy the female and males at 10,000
deer_female_r <- rarefy_even_depth(deer_female, rngseed = 1, sample.size = 10000)
deer_female_rare <- subset_samples(deer_female_r, !is.na(pathology))
deer_male_rare <- rarefy_even_depth(deer_male, rngseed = 1, sample.size = 10000)


#### Alpha diversity ######
# Alpha diversity for female
plot_richness(deer_female_rare) 

plot_richness(deer_female_rare, measures = c("Shannon","Observed")) 

gg_richness_female <- plot_richness(deer_female_rare, x = "pathology", measures = c("Shannon","Observed")) +
  xlab("pathology_female") +
  geom_boxplot()
gg_richness

ggsave(filename = "plot_richness_female.png"
       , gg_richness_female
       , height=4, width=6)

estimate_richness(deer_female_rare)

# Alpha diversity for males
gg_richness_male <- plot_richness(deer_male_rare, x = "pathology", measures = c("Shannon","Observed")) +
  xlab("pathology_male") +
  geom_boxplot()
gg_richness_male

ggsave(filename = "plot_richness_male.png"
       , gg_richness_male
       , height=4, width=6)

# phylogenetic diversity

# calculate Faith's phylogenetic diversity as PD

# Faith's female
phylo_dist_female <- pd(t(otu_table(deer_female_rare)), phy_tree(deer_female_rare),
                 include.root=F) 

# add PD to metadata table
sample_data(deer_female_rare)$PD <- phylo_dist_female$PD

# plot any metadata category against the PD
plot.pd <- ggplot(sample_data(deer_female_rare), aes(pathology, PD)) + 
  geom_boxplot() +
  xlab("Pathology_female") +
  ylab("Phylogenetic Diversity")

# view plot
plot.pd

# Faith's male
phylo_dist_male <- pd(t(otu_table(deer_male_rare)), phy_tree(deer_male_rare),
                        include.root=F) 

# add PD to metadata table
sample_data(deer_male_rare)$PD <- phylo_dist_male$PD

# plot any metadata category against the PD
plot.pd <- ggplot(sample_data(deer_male_rare), aes(pathology, PD)) + 
  geom_boxplot() +
  xlab("Pathology_male") +
  ylab("Phylogenetic Diversity")

# view plot
plot.pd

#### Beta diversity #####
# Beta diversity female
bc_dm <- distance(deer_female_rare, method="bray")
unweighted_dm <- distance(deer_female_rare, method = "unifrac")
weighted_dm <- distance(deer_female_rare, method = "wunifrac")

pcoa_bc <- ordinate(deer_female_rare, method="PCoA", distance=bc_dm)
pcoa_unweighted <-ordinate(deer_female_rare, method = "PCoA", distance = unweighted_dm)
pcoa_weighted <-ordinate(deer_female_rare, method = "PCoA", distance = weighted_dm)

plot_ordination(deer_female_rare, pcoa_bc, color = "pathology", shape="pathology")
plot_ordination(deer_female_rare, pcoa_unweighted, color = "pathology", shape="pathology")
plot_ordination(deer_female_rare, pcoa_weighted, color = "pathology", shape="pathology")

gg_pcoa_bc_female <- plot_ordination(deer_female_rare, pcoa_bc, color = "pathology", shape="pathology") +
  labs(col = "pathology") + stat_ellipse()
gg_pcoa_bc_female
  
ggsave("plot_pcoa_bc_female.png"
       , gg_pcoa_bc_female
       , height=4, width=5)

gg_pcoa_unweighted_female <- plot_ordination(deer_female_rare, pcoa_unweighted, color = "pathology", shape="pathology") +
  labs(col = "pathology") + stat_ellipse()
gg_pcoa_unweighted_female

ggsave("plot_pcoa_unweighted_female.png"
       , gg_pcoa_unweighted_female
       , height=4, width=5)

gg_pcoa_weighted_female <- plot_ordination(deer_female_rare, pcoa_weighted, color = "pathology", shape="pathology") +
  labs(col = "pathology") + stat_ellipse()
gg_pcoa_weighted_female

ggsave("plot_pcoa_weighted_female.png"
       , gg_pcoa_weighted_female
       , height=4, width=5)
# Beta diversity male
bc_dm_male <- distance(deer_male_rare, method="bray")
unweighted_dm_male <- distance(deer_male_rare, method = "unifrac")
weighted_dm_male <- distance(deer_male_rare, method = "wunifrac")

pcoa_bc_male <- ordinate(deer_male_rare, method="PCoA", distance=bc_dm_male)
pcoa_unweighted_male <-ordinate(deer_male_rare, method = "PCoA", distance = unweighted_dm_male)
pcoa_weighted_male <-ordinate(deer_male_rare, method = "PCoA", distance = weighted_dm_male)

plot_ordination(deer_male_rare, pcoa_bc_male, color = "pathology", shape="pathology")
plot_ordination(deer_male_rare, pcoa_unweighted_male, color = "pathology", shape="pathology")
plot_ordination(deer_male_rare, pcoa_weighted_male, color = "pathology", shape="pathology")

gg_pcoa_bc_male <- plot_ordination(deer_male_rare, pcoa_bc_male, color = "pathology", shape="pathology") +
  labs(col = "pathology") + stat_ellipse()
gg_pcoa_bc_male

ggsave("plot_pcoa_bc_male.png"
       , gg_pcoa_bc_male
       , height=4, width=5)

gg_pcoa_unweighted_male <- plot_ordination(deer_male_rare, pcoa_unweighted_male, color = "pathology", shape="pathology") +
  labs(col = "pathology") + stat_ellipse()
gg_pcoa_unweighted_male

ggsave("plot_pcoa_unweighted_male.png"
       , gg_pcoa_unweighted_male
       , height=4, width=5)

gg_pcoa_weighted_male <- plot_ordination(deer_male_rare, pcoa_weighted_male, color = "pathology", shape="pathology") +
  labs(col = "pathology") + stat_ellipse()
gg_pcoa_weighted_male

ggsave("plot_pcoa_weighted_male.png"
       , gg_pcoa_weighted_male
       , height=4, width=5)

# alpha diversity stats 
# female
dat_f <- as(sample_data(deer_female_rare), "data.frame")
kruskal.test(PD ~ pathology, data = dat_f)

#male
dat_m <- as(sample_data(deer_male_rare), "data.frame")
kruskal.test(PD ~ pathology, data = dat_m)

# beta diversity stats
# female
dm_unifrac_female <- UniFrac(deer_female_rare, weighted=TRUE) # Weighted UniFrac
dm_braycurtis_female <- vegdist(t(otu_table(deer_female_rare)), method="bray") # Bray-curtis
dm_jaccard_female <- vegdist(t(otu_table(deer_female_rare)), method="jaccard") # Jaccard

dat_female <- data.frame (sample_data (deer_female_rare))
str(dat_female)
adonis2(dm_unifrac_female ~ pathology, data = dat_female, by = "terms")
adonis2(dm_braycurtis_female ~ pathology, data=dat_female, by = "terms")
adonis2(dm_jaccard_female ~ pathology, data=dat_female, by = "terms")

# male
dm_unifrac_male <- UniFrac(deer_male_rare, weighted=TRUE) # Weighted UniFrac
dm_braycurtis_male <- vegdist(t(otu_table(deer_male_rare)), method="bray") # Bray-curtis
dm_jaccard_male <- vegdist(t(otu_table(deer_male_rare)), method="jaccard") # Jaccard

dat_male <- data.frame (sample_data (deer_male_rare))
adonis2(dm_unifrac_male ~ pathology, data = dat_male, by = "terms")
adonis2(dm_braycurtis_male ~ pathology, data=dat_male, by = "terms")
adonis2(dm_jaccard_male ~ pathology, data=dat_male, by = "terms")

