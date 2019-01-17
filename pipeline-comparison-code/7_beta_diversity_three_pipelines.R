rm(list=ls())
getwd()
setwd("C:/Users/sierr/Documents/ThesisProject_Pipelines/")

#packages
library(ggforce)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(phyloseq)
library(plyr)
library(PMCMR)
library(tidyverse)
library(vegan)

tax_group <- read.csv("tax_group.csv")
tax_group <- tax_group %>% group_by(Kingdom, Phylum, Class, Order, Family, Genus) %>% summarise_all(funs(sum))
tax_group <- tax_group %>%  ungroup()

regroup_physeq_object <-function(table) {
  tax <- table %>% select(Kingdom, Phylum, Class, Order, Family, Genus)
  tax <- as.matrix(tax)
  otu <- table %>% select(contains("WCME"))
  OTU=otu_table(otu, taxa_are_rows=TRUE)
  TAX=tax_table(tax)
  taxa_names(TAX)=row.names(OTU)
  physeq_all=phyloseq(OTU,TAX)
  return(physeq_all)
}

physeq_all <- regroup_physeq_object(tax_group)

#import metadata and combine
metadata=(read.csv("Metadata/HPMMMeta_r_merge.csv",header=TRUE))
sampdat=sample_data(metadata)
sample_names(sampdat)=metadata$SampleID
physeq=merge_phyloseq(physeq_all, sampdat)
physeq

#of samples
physeq_rec <- subset_samples(physeq, Sample_Area == "Rectum")
physeq_mou <- subset_samples(physeq, Sample_Area == "Mouth")

physeq_Q_r <- subset_samples(physeq_rec, Pipeline == "QIIME2")
physeq_M_r <- subset_samples(physeq_rec, Pipeline == "mothur")
physeq_G_r <- subset_samples(physeq_rec, Pipeline == "MG-RAST")

physeq_QM_r <- merge_phyloseq(physeq_Q_r, physeq_M_r)
physeq_QG_r <- merge_phyloseq(physeq_Q_r, physeq_G_r)
physeq_GM_r <- merge_phyloseq(physeq_M_r, physeq_G_r)

physeq_Q_m <- subset_samples(physeq_mou, Pipeline == "QIIME2")
physeq_M_m <- subset_samples(physeq_mou, Pipeline == "mothur")
physeq_G_m <- subset_samples(physeq_mou, Pipeline == "MG-RAST")

physeq_QM_m <- merge_phyloseq(physeq_Q_m, physeq_M_m)
physeq_QG_m <- merge_phyloseq(physeq_Q_m, physeq_G_m)
physeq_GM_m <- merge_phyloseq(physeq_M_m, physeq_G_m)

#beta diversity
#distance= jaccard 

theme_set(theme_classic(base_size = 14))

#total
ord = ordinate(physeq, method="PCoA", distance="(A+B-2*J)/(A+B-J)")
ordplot=plot_ordination(physeq, ord, color="Pipeline")
ordplot
tiff("total_beta_jaccard_three_pipelines.TIF", width = 800, height = 600)
ordplot+ geom_point(size = 2, aes(color = sample_data(physeq)$Pipeline)) + scale_color_manual(values = c("#787878", "#ffb31a", "#800000", "#5c5c8a", "#918151")) +
  stat_ellipse(alpha = 1, size=2, aes(color= Sample_Area))
dev.off()

ord = ordinate(physeq_rec, method="PCoA", distance="(A+B-2*J)/(A+B-J)")
ordplot=plot_ordination(physeq_rec, ord, color="Pipeline")
ordplot
tiff("rec_beta_jaccard_three_pipelines.TIF", width = 1400, height = 1200)
ordplot+ geom_point(size = 5, aes(color = sample_data(physeq_rec)$Pipeline)) + 
  scale_color_manual(values = c("#787878", "#ffb31a", "#5c5c8a")) + 
  stat_ellipse(alpha = 1, size=3, aes(fill = Pipeline)) +
  scale_fill_manual(values = c("#787878", "#ffb31a", "#5c5c8a")) + 
  theme_bw(base_size = 40) + theme(panel.grid.major = element_blank(), 
                                   panel.grid.minor = element_blank(), 
                                   axis.line = element_line(colour = "black"))
dev.off()

ord = ordinate(physeq_mou, method="PCoA", distance="(A+B-2*J)/(A+B-J)")
ordplot=plot_ordination(physeq_mou, ord, color="Pipeline")
ordplot
tiff("mou_beta_jaccard_three_pipelines.TIF",  width = 1400, height = 1200)
ordplot+ geom_point(size = 5, aes(color = sample_data(physeq_mou)$Pipeline)) + 
  scale_color_manual(values = c("#787878", "#ffb31a", "#5c5c8a")) + 
  stat_ellipse(alpha = 1, size=3, aes(fill = Pipeline)) +
  scale_fill_manual(values = c("#787878", "#ffb31a", "#5c5c8a")) + 
  theme_bw(base_size = 40) + theme(panel.grid.major = element_blank(), 
                                   panel.grid.minor = element_blank(), 
                                   axis.line = element_line(colour = "black"))
dev.off()

beta_diversity_calc <- function(physeq) {
  GPdist=phyloseq::distance(physeq, "(A+B-2*J)/(A+B-J)")
  sampledf <- data.frame(sample_data(physeq))
  print(return(adonis(GPdist ~ Pipeline, data = sampledf)))
}

beta_dispersion_calc <- function(physeq) {
  GPdist=phyloseq::distance(physeq, "(A+B-2*J)/(A+B-J)")
  sampledf <- data.frame(sample_data(physeq))
  beta <- betadisper(GPdist, sampledf$Pipeline)
  print(return(permutest(beta)))
}

beta_diversity_calc(physeq)
beta_dispersion_calc(physeq)

beta_diversity_calc(physeq_rec)
beta_dispersion_calc(physeq_rec)

beta_diversity_calc(physeq_mou)
beta_dispersion_calc(physeq_mou)

beta_diversity_calc(physeq_QG_r)
beta_dispersion_calc(physeq_QG_r)

beta_diversity_calc(physeq_QM_r)
beta_dispersion_calc(physeq_QM_r)

beta_diversity_calc(physeq_GM_r)
beta_dispersion_calc(physeq_GM_r)

beta_diversity_calc(physeq_QG_m)
beta_dispersion_calc(physeq_QG_m)

beta_diversity_calc(physeq_QM_m)
beta_dispersion_calc(physeq_QM_m)

beta_diversity_calc(physeq_GM_m)
beta_dispersion_calc(physeq_GM_m)
