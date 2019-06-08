# Import Data (Rarefying/Subsample) -------------------------------------------------------------

#packages
library(ggforce)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(phyloseq)
library(plyr)
library(PMCMR)
library(randomForest)
library(tidyverse)


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

tax_group_rare <- read.csv("tax_group_rare.csv")
tax_group_rare <- tax_group_rare %>% group_by(Kingdom, Phylum, Class, Order, Family, Genus) %>% summarise_all(funs(sum))
tax_group_rare <- tax_group_rare %>%  ungroup()

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
physeq_rare <- regroup_physeq_object(tax_group_rare)

#import metadata and combine
metadata=(read.csv("Metadata/HPMMMeta_r_merge_rare.csv",header=TRUE))
sampdat=sample_data(metadata)
sample_names(sampdat)=metadata$SampleID

metadata_rare=(read.csv("Metadata/HPMMMeta_r_merge_rare_levels.csv",header=TRUE))
sampdat_rare=sample_data(metadata_rare)
sampdat_rare$Rarefy <- factor(sampdat_rare$Rarefy, levels = c("No Rarefaction", "7000", "1000"))
sample_names(sampdat_rare)=metadata_rare$SampleID

merger = merge_phyloseq(physeq_rare, sampdat_rare)
sample_data(merger)$Rarefy <- factor(sample_data(merger)$Rarefy, levels = c("No Rarefaction", "7000", "1000"))

physeq_nr=merge_phyloseq(physeq_nr, sampdat)
physeq_7=merge_phyloseq(physeq_7, sampdat)
physeq_1=merge_phyloseq(physeq_1, sampdat)

# Table S14 --------------------------------------------------------------------
unclassified_sequences <- function(physeq) {
  physeq_rec <- subset_samples(physeq, Sample_Area == "Rectum")
  physeq_mou <- subset_samples(physeq, Sample_Area == "Mouth")
  physeq_rec_60 <- subset_samples(physeq_rec, Subsample == '60')
  physeq_rec_120 <- subset_samples(physeq_rec, Subsample == '120')
  physeq_rec_188 <- subset_samples(physeq_rec, Subsample == '188')
  physeq_mou_60 <- subset_samples(physeq_mou, Subsample == '60')
  physeq_mou_120 <- subset_samples(physeq_mou, Subsample == '120')
  physeq_mou_188 <- subset_samples(physeq_mou, Subsample == '188')
  physeq_Q_r_60_unc = subset_taxa(physeq_rec_60, Phylum=="Unclassified")
  physeq_Q_r_120_unc = subset_taxa(physeq_rec_120, Phylum=="Unclassified")
  physeq_Q_r_188_unc = subset_taxa(physeq_rec_188, Phylum=="Unclassified")
  physeq_Q_m_60_unc = subset_taxa(physeq_mou_60, Phylum=="Unclassified")
  physeq_Q_m_120_unc = subset_taxa(physeq_mou_120, Phylum=="Unclassified")
  physeq_Q_m_188_unc = subset_taxa(physeq_mou_188, Phylum=="Unclassified")
  physeq_Q_r_60_uncf = subset_taxa(physeq_rec_60, Family=="Unclassified")
  physeq_Q_r_120_uncf = subset_taxa(physeq_rec_120, Family=="Unclassified")
  physeq_Q_r_188_uncf = subset_taxa(physeq_rec_188, Family=="Unclassified")
  physeq_Q_m_60_uncf = subset_taxa(physeq_mou_60, Family=="Unclassified")
  physeq_Q_m_120_uncf = subset_taxa(physeq_mou_120, Family=="Unclassified")
  physeq_Q_m_188_uncf = subset_taxa(physeq_mou_188, Family=="Unclassified")
  a <- sum(sample_sums(physeq_rec_60))
  b <- sum(sample_sums(physeq_rec_120))
  c <- sum(sample_sums(physeq_rec_188))
  d <- sum(sample_sums(physeq_mou_60))
  e <- sum(sample_sums(physeq_mou_120))
  f <- sum(sample_sums(physeq_mou_188))
  aa <- sum(sample_sums(physeq_Q_r_60_unc))
  bb <- sum(sample_sums(physeq_Q_r_120_unc))
  cc <- sum(sample_sums(physeq_Q_r_188_unc))
  dd <- sum(sample_sums(physeq_Q_m_60_unc))
  ee <- sum(sample_sums(physeq_Q_m_120_unc))
  ff <- sum(sample_sums(physeq_Q_m_188_unc))
  aaa <- sum(sample_sums(physeq_Q_r_60_uncf))
  bbb <- sum(sample_sums(physeq_Q_r_120_uncf))
  ccc <- sum(sample_sums(physeq_Q_r_188_uncf))
  ddd <- sum(sample_sums(physeq_Q_m_60_uncf))
  eee <- sum(sample_sums(physeq_Q_m_120_uncf))
  fff <- sum(sample_sums(physeq_Q_m_188_uncf))
  df <- data.frame('Total' = c(a,b,c,d,e,f), 'Phylum' = c(aa,bb,cc,dd,ee,ff),
                   'Family' = c(aaa,bbb,ccc,ddd,eee,fff))
}

df <- unclassified_sequences(physeq_nr)
dfa <- unclassified_sequences(physeq_7)
dfc <- unclassified_sequences(physeq_1)
