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

#make a table in excel for rare data
tax_group_nr <- read.csv("tax_group_norare.csv")
tax_group_nr <- tax_group_nr %>% group_by(Kingdom, Phylum, Class, Order, Family, Genus) %>% summarise_all(funs(sum))
tax_group_nr <- tax_group_nr %>%  ungroup()

tax_group_7 <- read.csv("tax_group_7000.csv")
tax_group_7 <- tax_group_7 %>% group_by(Kingdom, Phylum, Class, Order, Family, Genus) %>% summarise_all(funs(sum))
tax_group_7 <- tax_group_7 %>%  ungroup()

tax_group_1 <- read.csv("tax_group_1000.csv")
tax_group_1 <- tax_group_1 %>% group_by(Kingdom, Phylum, Class, Order, Family, Genus) %>% summarise_all(funs(sum))
tax_group_1 <- tax_group_1 %>%  ungroup()

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

physeq_nr <- regroup_physeq_object(tax_group_nr)
physeq_7 <- regroup_physeq_object(tax_group_7)
physeq_1 <- regroup_physeq_object(tax_group_1)

#import metadata and combine
metadata=(read.csv("Metadata/HPMMMeta_r_merge_rare.csv",header=TRUE))
sampdat=sample_data(metadata)
sample_names(sampdat)=metadata$SampleID
sampdat$Subsample <- factor(sampdat$Subsample, levels = c('60', '120', '188'))

physeq_nr=merge_phyloseq(physeq_nr, sampdat)
physeq_7=merge_phyloseq(physeq_7, sampdat)
physeq_1=merge_phyloseq(physeq_1, sampdat)

#testing subsample
physeq_rec_nr <- subset_samples(physeq_nr, Sample_Area == "Rectum")
physeq_mou_nr <- subset_samples(physeq_nr, Sample_Area == "Mouth")

physeq_rec_7 <- subset_samples(physeq_7, Sample_Area == "Rectum")
physeq_mou_7 <- subset_samples(physeq_7, Sample_Area == "Mouth")

physeq_rec_1 <- subset_samples(physeq_1, Sample_Area == "Rectum")
physeq_mou_1 <- subset_samples(physeq_1, Sample_Area == "Mouth")


#pairwise subsample
pairwise_subsample <- function(physeq, physeq1, physeq7) {
  physeq_rec_nr <- subset_samples(physeq, Sample_Area == "Rectum")
  physeq_mou_nr <- subset_samples(physeq, Sample_Area == "Mouth")
  physeq_rec_7 <- subset_samples(physeq7, Sample_Area == "Rectum")
  physeq_mou_7 <- subset_samples(physeq7, Sample_Area == "Mouth")
  physeq_rec_1 <- subset_samples(physeq1, Sample_Area == "Rectum")
  physeq_mou_1 <- subset_samples(physeq1, Sample_Area == "Mouth")
  physeq_rec_nr_188 <- subset_samples(physeq_rec_nr, Subsample=='188')
  physeq_rec_nr_120 <- subset_samples(physeq_rec_nr, Subsample=='120')
  physeq_rec_nr_60 <- subset_samples(physeq_rec_nr, Subsample=='60')
  physeq_rec_7_188 <- subset_samples(physeq_rec_7, Subsample=='188')
  physeq_rec_7_120 <- subset_samples(physeq_rec_7, Subsample=='120')
  physeq_rec_7_60 <- subset_samples(physeq_rec_7, Subsample=='60')
  physeq_rec_1_188 <- subset_samples(physeq_rec_1, Subsample=='188')
  physeq_rec_1_120 <- subset_samples(physeq_rec_1, Subsample=='120')
  physeq_rec_1_60 <- subset_samples(physeq_rec_1, Subsample=='60')
  physeq_mou_nr_188 <- subset_samples(physeq_mou_nr, Subsample=='188')
  physeq_mou_nr_120 <- subset_samples(physeq_mou_nr, Subsample=='120')
  physeq_mou_nr_60 <- subset_samples(physeq_mou_nr, Subsample=='60')
  physeq_mou_7_188 <- subset_samples(physeq_mou_7, Subsample=='188')
  physeq_mou_7_120 <- subset_samples(physeq_mou_7, Subsample=='120')
  physeq_mou_7_60 <- subset_samples(physeq_mou_7, Subsample=='60')
  physeq_mou_1_188 <- subset_samples(physeq_mou_1, Subsample=='188')
  physeq_mou_1_120 <- subset_samples(physeq_mou_1, Subsample=='120')
  physeq_mou_1_60 <- subset_samples(physeq_mou_1, Subsample=='60')
  return(list(physeq_rec_nr_60120 <- merge_phyloseq(physeq_rec_nr_60, physeq_rec_nr_120),
              physeq_rec_nr_60188 <- merge_phyloseq(physeq_rec_nr_60, physeq_rec_nr_188),
              physeq_rec_nr_120188 <- merge_phyloseq(physeq_rec_nr_188, physeq_rec_nr_120),
              physeq_rec_7_60120 <- merge_phyloseq(physeq_rec_7_60, physeq_rec_7_120),
              physeq_rec_7_60188 <- merge_phyloseq(physeq_rec_7_60, physeq_rec_7_188),
              physeq_rec_7_120188 <- merge_phyloseq(physeq_rec_7_188, physeq_rec_7_120),
              physeq_rec_1_60120 <- merge_phyloseq(physeq_rec_1_60, physeq_rec_1_120),
              physeq_rec_1_60188 <- merge_phyloseq(physeq_rec_1_60, physeq_rec_1_188),
              physeq_rec_1_120188 <- merge_phyloseq(physeq_rec_1_188, physeq_rec_1_120),
              physeq_mou_nr_60120 <- merge_phyloseq(physeq_mou_nr_60, physeq_mou_nr_120),
              physeq_mou_nr_60188 <- merge_phyloseq(physeq_mou_nr_60, physeq_mou_nr_188),
              physeq_mou_nr_120188 <- merge_phyloseq(physeq_mou_nr_188, physeq_mou_nr_120),
              physeq_mou_7_60120 <- merge_phyloseq(physeq_mou_7_60, physeq_mou_7_120),
              physeq_mou_7_60188 <- merge_phyloseq(physeq_mou_7_60, physeq_mou_7_188),
              physeq_mou_7_120188 <- merge_phyloseq(physeq_mou_7_188, physeq_mou_7_120),
              physeq_mou_1_60120 <- merge_phyloseq(physeq_mou_1_60, physeq_mou_1_120),
              physeq_mou_1_60188 <- merge_phyloseq(physeq_mou_1_60, physeq_mou_1_188),
              physeq_mou_1_120188 <- merge_phyloseq(physeq_mou_1_188, physeq_mou_1_120)))
}

#testing rarefaction 
sample_names(physeq_nr) <- paste("Norare_", sample_names(physeq_nr), sep="")
sample_names(physeq_7) <- paste("7000_", sample_names(physeq_7), sep="")
sample_names(physeq_1) <- paste("1000_", sample_names(physeq_1), sep="")

merger= merge_phyloseq(physeq_nr, physeq_7, physeq_1)
metadata_rare=(read.csv("Metadata/HPMMMeta_r_merge_rare_levels.csv",header=TRUE))
sampdat=sample_data(metadata_rare)
sample_names(sampdat)=metadata_rare$SampleID
sampdat$Subsample <- factor(sampdat$Subsample, levels = c('60', '120', '188'))
merger = merge_phyloseq(merger, sampdat)

physeq_rec_188 <- subset_samples(merger, Sample_Area == "Rectum")
physeq_rec_188 <- subset_samples(physeq_rec_188, Subsample=='188')
physeq_mou_188 <- subset_samples(merger, Sample_Area == "Mouth")
physeq_mou_188 <- subset_samples(physeq_mou_188, Subsample=='188')
physeq_rec_120 <- subset_samples(merger, Sample_Area == "Rectum")
physeq_rec_120 <- subset_samples(physeq_rec_120, Subsample=='120')
physeq_mou_120 <- subset_samples(merger, Sample_Area == "Mouth")
physeq_mou_120 <- subset_samples(physeq_mou_120, Subsample=='120')
physeq_rec_60 <- subset_samples(merger, Sample_Area == "Rectum")
physeq_rec_60 <- subset_samples(physeq_rec_60, Subsample=='60')
physeq_mou_60 <- subset_samples(merger, Sample_Area == "Mouth")
physeq_mou_60 <- subset_samples(physeq_mou_60, Subsample=='60')

#rarefy pairwise
rarefy_pairwise <- function(physeq) {
  physeq_rec <- subset_samples(merger, Sample_Area == "Rectum")
  physeq_mou <- subset_samples(merger, Sample_Area == "Mouth")
  physeq_rec_188 <- subset_samples(physeq_rec, Subsample=='188')
  physeq_mou_188 <- subset_samples(physeq_mou, Subsample=='188')
  physeq_rec_120 <- subset_samples(physeq_rec, Subsample=='120')
  physeq_mou_120 <- subset_samples(physeq_mou, Subsample=='120')
  physeq_rec_60 <- subset_samples(physeq_rec, Subsample=='60')
  physeq_mou_60 <- subset_samples(physeq_mou, Subsample=='60')
  physeq_rec_188_nr <- subset_samples(physeq_rec_188, Rarefy == 'No Rarefaction')
  physeq_rec_188_7 <- subset_samples(physeq_rec_188, Rarefy == '7000')
  physeq_rec_188_1 <- subset_samples(physeq_rec_188, Rarefy == '1000')
  physeq_rec_120_nr <- subset_samples(physeq_rec_120, Rarefy == 'No Rarefaction')
  physeq_rec_120_7 <- subset_samples(physeq_rec_120, Rarefy == '7000')
  physeq_rec_120_1 <- subset_samples(physeq_rec_120, Rarefy == '1000')
  physeq_rec_60_nr <- subset_samples(physeq_rec_60, Rarefy == 'No Rarefaction')
  physeq_rec_60_7 <- subset_samples(physeq_rec_60, Rarefy == '7000')
  physeq_rec_60_1 <- subset_samples(physeq_rec_60, Rarefy == '1000')
  physeq_mou_188_nr <- subset_samples(physeq_mou_188, Rarefy == 'No Rarefaction')
  physeq_mou_188_7 <- subset_samples(physeq_mou_188, Rarefy == '7000')
  physeq_mou_188_1 <- subset_samples(physeq_mou_188, Rarefy == '1000')
  physeq_mou_120_nr <- subset_samples(physeq_mou_120, Rarefy == 'No Rarefaction')
  physeq_mou_120_7 <- subset_samples(physeq_mou_120, Rarefy == '7000')
  physeq_mou_120_1 <- subset_samples(physeq_mou_120, Rarefy == '1000')
  physeq_mou_60_nr <- subset_samples(physeq_mou_60, Rarefy == 'No Rarefaction')
  physeq_mou_60_7 <- subset_samples(physeq_mou_60, Rarefy == '7000')
  physeq_mou_60_1 <- subset_samples(physeq_mou_60, Rarefy == '1000')
  return(list(physeq_rec_188_nr7 <- merge_phyloseq(physeq_rec_188_nr, physeq_rec_188_7),
  physeq_rec_188_nr1 <- merge_phyloseq(physeq_rec_188_nr, physeq_rec_188_1),
  physeq_rec_188_71 <- merge_phyloseq(physeq_rec_188_1, physeq_rec_188_7),
  physeq_rec_120_nr7 <- merge_phyloseq(physeq_rec_120_nr, physeq_rec_120_7),
  physeq_rec_120_nr1 <- merge_phyloseq(physeq_rec_120_nr, physeq_rec_120_1),
  physeq_rec_120_71 <- merge_phyloseq(physeq_rec_120_1, physeq_rec_120_7),
  physeq_rec_60_nr7 <- merge_phyloseq(physeq_rec_60_nr, physeq_rec_60_7),
  physeq_rec_60_nr1 <- merge_phyloseq(physeq_rec_60_nr, physeq_rec_60_1),
  physeq_rec_60_71 <- merge_phyloseq(physeq_rec_60_1, physeq_rec_60_7),
  physeq_mou_188_nr7 <- merge_phyloseq(physeq_mou_188_nr, physeq_mou_188_7),
  physeq_mou_188_nr1 <- merge_phyloseq(physeq_mou_188_nr, physeq_mou_188_1),
  physeq_mou_188_71 <- merge_phyloseq(physeq_mou_188_1, physeq_mou_188_7),
  physeq_mou_120_nr7 <- merge_phyloseq(physeq_mou_120_nr, physeq_mou_120_7),
  physeq_mou_120_nr1 <- merge_phyloseq(physeq_mou_120_nr, physeq_mou_120_1),
  physeq_mou_120_71 <- merge_phyloseq(physeq_mou_120_1, physeq_mou_120_7),
  physeq_mou_60_nr7 <- merge_phyloseq(physeq_mou_60_nr, physeq_mou_60_7),
  physeq_mou_60_nr1 <- merge_phyloseq(physeq_mou_60_nr, physeq_mou_60_1),
  physeq_mou_60_71 <- merge_phyloseq(physeq_mou_60_1, physeq_mou_60_7)))
}

ord = ordinate(physeq_rec_188, method="PCoA", distance="(A+B-2*J)/(A+B-J)")
ordplot=plot_ordination(physeq_rec_188, ord, color="Rarefy")
ordplot
theme_set(theme_bw(base_size = 14))
tiff("rec_188_beta_jaccard_.TIF", width = 800, height = 600)
ordplot+ geom_point(size = 2, aes(color = factor(sample_data(physeq_rec_188)$Rarefy))) + scale_color_manual(values = c('purple', 'lightgreen', 'brown'))
dev.off()

ord = ordinate(physeq_rec_120, method="PCoA", distance="(A+B-2*J)/(A+B-J)")
ordplot=plot_ordination(physeq_rec_120, ord, color="Rarefy")
ordplot
theme_set(theme_bw(base_size = 14))
tiff("rec_120_beta_jaccard_.TIF", width = 800, height = 600)
ordplot+ geom_point(size = 2, aes(color = factor(sample_data(physeq_rec_120)$Rarefy))) + scale_color_manual(values = c('purple', 'lightgreen', 'brown'))
dev.off()

ord = ordinate(physeq_rec_60, method="PCoA", distance="(A+B-2*J)/(A+B-J)")
ordplot=plot_ordination(physeq_rec_60, ord, color="Rarefy")
ordplot
theme_set(theme_bw(base_size = 14))
tiff("rec_60_beta_jaccard_.TIF", width = 800, height = 600)
ordplot+ geom_point(size = 2, aes(color = factor(sample_data(physeq_rec_60)$Rarefy))) + scale_color_manual(values = c('purple', 'lightgreen', 'brown'))
dev.off()

ord = ordinate(physeq_mou_188, method="PCoA", distance="(A+B-2*J)/(A+B-J)")
ordplot=plot_ordination(physeq_mou_188, ord, color="Rarefy")
ordplot
theme_set(theme_bw(base_size = 14))
tiff("mou_188_beta_jaccard_.TIF", width = 800, height = 600)
ordplot+ geom_point(size = 2, aes(color = factor(sample_data(physeq_mou_188)$Rarefy))) + scale_color_manual(values = c('purple', 'lightgreen', 'brown'))
dev.off()

ord = ordinate(physeq_mou_120, method="PCoA", distance="(A+B-2*J)/(A+B-J)")
ordplot=plot_ordination(physeq_mou_120, ord, color="Rarefy")
ordplot
theme_set(theme_bw(base_size = 14))
tiff("mou_120_beta_jaccard_.TIF", width = 800, height = 600)
ordplot+ geom_point(size = 2, aes(color = factor(sample_data(physeq_mou_120)$Rarefy))) + scale_color_manual(values = c('purple', 'lightgreen', 'brown'))
dev.off()

ord = ordinate(physeq_mou_60, method="PCoA", distance="(A+B-2*J)/(A+B-J)")
ordplot=plot_ordination(physeq_mou_60, ord, color="Rarefy")
ordplot
theme_set(theme_bw(base_size = 14))
tiff("mou_60_beta_jaccard_.TIF", width = 800, height = 600)
ordplot+ geom_point(size = 2, aes(color = factor(sample_data(physeq_mou_60)$Rarefy))) + scale_color_manual(values = c('purple', 'lightgreen', 'brown'))
dev.off()

beta_diversity_calc_rarefy <- function(physeq) {
  GPdist=phyloseq::distance(physeq, "(A+B-2*J)/(A+B-J)")
  sampledf <- data.frame(sample_data(physeq))
  print(return(adonis(GPdist ~ Rarefy, data = sampledf)))
}

beta_dispersion_calc_rarefy <- function(physeq) {
  GPdist=phyloseq::distance(physeq, "(A+B-2*J)/(A+B-J)")
  sampledf <- data.frame(sample_data(physeq))
  beta <- betadisper(GPdist, sampledf$Rarefy)
  print(return(permutest(beta)))
}

rarefy_list <- list(physeq_rec_188, physeq_rec_120, physeq_rec_60, 
                    physeq_mou_188, physeq_mou_120, physeq_mou_60 )

for(i in 1:length(rarefy_list)) {
  print(beta_diversity_calc_rarefy(rarefy_list[[i]]))
  print(beta_dispersion_calc_rarefy(rarefy_list[[i]]))
}

pairwise_rare_list <- rarefy_pairwise(merger)

for(i in 1:length(pairwise_rare_list)) {
  print(beta_diversity_calc_rarefy(pairwise_rare_list[[i]]))
  print(beta_dispersion_calc_rarefy(pairwise_rare_list[[i]]))
}

beta_diversity_calc_subsample <- function(physeq) {
  GPdist=phyloseq::distance(physeq, "(A+B-2*J)/(A+B-J)")
  sampledf <- data.frame(sample_data(physeq))
  print(return(adonis(GPdist ~ Subsample, data = sampledf)))
}

beta_dispersion_calc_subsample <- function(physeq) {
  GPdist=phyloseq::distance(physeq, "(A+B-2*J)/(A+B-J)")
  sampledf <- data.frame(sample_data(physeq))
  beta <- betadisper(GPdist, sampledf$Subsample)
  print(return(permutest(beta)))
}

subsample_list <- list(physeq_rec_nr, physeq_rec_7, physeq_rec_1,
                       physeq_mou_nr, physeq_mou_7, physeq_mou_1)

for(i in 1:length(subsample_list)) {
  print(beta_diversity_calc_subsample(subsample_list[[i]]))
  print(beta_dispersion_calc_subsample(subsample_list[[i]]))
}
