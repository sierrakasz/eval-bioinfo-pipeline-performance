## Table 4
# Table 4 is a summary of the differential effect of minimum library size and sample size
# on filtered and unclassified sequences. 

#load required packages
#phyloseq requires BiocManager 
# code for phyloseq install: if(!requireNamespace("BiocManager")){
#install.packages("BiocManager")
#}
#BiocManager::install("phyloseq")

#packages
library(httr)
library(HTSSIP)
library(gridExtra)
library(ggforce)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(PMCMR)
library(phyloseq)
library(plyr)
library(tidyverse)
library(vegan)
library(HTSSIP)
library(gridExtra)
library(ggforce)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(PMCMR)
library(phyloseq)
library(plyr)
library(tidyverse)
library(vegan)


# import data ------------------------------------------------------------

###tax group is a file that combines the OTU file and taxonomy file. 
#no minimum library size
tax_group_nr <- read.csv(text=GET('https://github.com/sierrakasz/eval-bioinfo-pipeline-performance/blob/master/min_lib_sample_size/tax_group_norare.csv'))
tax_group_nr <- tax_group_nr %>% group_by(Kingdom, Phylum, Class, Order, Family, Genus) %>% summarise_all(funs(sum))
tax_group_nr <- tax_group_nr %>%  ungroup()
#Minimum library size: 7000
tax_group_7 <- read.csv(text=GET("https://github.com/sierrakasz/eval-bioinfo-pipeline-performance/blob/master/min_lib_sample_size/tax_group_7000.csv"))
tax_group_7 <- tax_group_7 %>% group_by(Kingdom, Phylum, Class, Order, Family, Genus) %>% summarise_all(funs(sum))
tax_group_7 <- tax_group_7 %>%  ungroup()
#Minimum library: 1000
tax_group_1 <- read.csv(text=GET("https://github.com/sierrakasz/eval-bioinfo-pipeline-performance/blob/master/min_lib_sample_size/tax_group_1000.csv"))
tax_group_1 <- tax_group_1 %>% group_by(Kingdom, Phylum, Class, Order, Family, Genus) %>% summarise_all(funs(sum))
tax_group_1 <- tax_group_1 %>%  ungroup()
#All of the combined together
#So, same samples among minimum library size, but with different minimum libraries
tax_group_rare <- read.csv(text=GET("https://github.com/sierrakasz/eval-bioinfo-pipeline-performance/blob/master/min_lib_sample_size/tax_group_rare.csv"))
tax_group_rare <- tax_group_rare %>% group_by(Kingdom, Phylum, Class, Order, Family, Genus) %>% summarise_all(funs(sum))
tax_group_rare <- tax_group_rare %>%  ungroup()

#this function specifically collapses the OTU table and taxonomy file 
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

#phyloseq objects moving forward
physeq_nr <- regroup_physeq_object(tax_group_nr)
physeq_7 <- regroup_physeq_object(tax_group_7)
physeq_1 <- regroup_physeq_object(tax_group_1)
physeq_rare <- regroup_physeq_object(tax_group_rare)

#import metadata
metadata=(text=GET("https://github.com/sierrakasz/eval-bioinfo-pipeline-performance/blob/master/min_lib_sample_size/HPMMMeta_r_merge_rare.csv"))
sampdat=sample_data(metadata)
sample_names(sampdat)=metadata$SampleID
#Change subsample level into a factor instead of numeric
sampdat$Subsample <- factor(sampdat$Subsample, levels = c('60', '120', '188'))
#combine into one phyloseq object
physeq_nr=merge_phyloseq(physeq_nr, sampdat)
physeq_7=merge_phyloseq(physeq_7, sampdat)
physeq_1=merge_phyloseq(physeq_1, sampdat)

#separate metadata and phyloseq object for having all minimum library sizes at largest subsample size (188)
metadata_rare=(text=GET("https://github.com/sierrakasz/eval-bioinfo-pipeline-performance/blob/master/min_lib_sample_size/HPMMMeta_r_merge_rare_levels.csv"))
sampdat_rare=sample_data(metadata_rare)
#change to factors with level
sampdat_rare$Rarefy <- factor(sampdat_rare$Rarefy, levels = c("No Rarefaction", "7000", "1000"))
sample_names(sampdat_rare)=metadata_rare$SampleID
merger = merge_phyloseq(physeq_rare, sampdat_rare)
sample_data(merger)$Rarefy <- factor(sample_data(merger)$Rarefy, levels = c("No Rarefaction", "7000", "1000"))
#final phyloseq object
physeq_188 <- subset_samples(merger, Subsample == '188')

##unclassified sequences function
#subset by sample area, sample size, and minimum library size
# subset out Unclassified taxa at two taxonomic levels: phylum and family
# putting together all information into a dataframe
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

#all information was added to a table in excel for formatting into the manuscript. 
