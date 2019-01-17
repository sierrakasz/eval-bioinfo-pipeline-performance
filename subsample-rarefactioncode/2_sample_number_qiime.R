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

physeq_nr=merge_phyloseq(physeq_nr, sampdat)
physeq_7=merge_phyloseq(physeq_7, sampdat)
physeq_1=merge_phyloseq(physeq_1, sampdat)

sample_num_function <- function(physeq) {
  physeq_rec <- subset_samples(physeq, Sample_Area == "Rectum")
  physeq_mou <- subset_samples(physeq, Sample_Area == "Mouth")
  physeq_rec_60 <- subset_samples(physeq_rec, Subsample == '60')
  physeq_rec_120 <- subset_samples(physeq_rec, Subsample == '120')
  physeq_rec_188 <- subset_samples(physeq_rec, Subsample == '188')
  physeq_mou_60 <- subset_samples(physeq_mou, Subsample == '60')
  physeq_mou_120 <- subset_samples(physeq_mou, Subsample == '120')
  physeq_mou_188 <- subset_samples(physeq_mou, Subsample == '188')
  QR60 <- length(sample_names(physeq_rec_60))
  QR120 <- length(sample_names(physeq_rec_120))
  QR188 <- length(sample_names(physeq_rec_188))
  QM60 <- length(sample_names(physeq_mou_60))
  QM120 <- length(sample_names(physeq_mou_120))
  QM188 <- length(sample_names(physeq_mou_188))
  QT60 <- QR60 + QM60
  QT120 <- QR120 + QM120
  QT188 <- QR188 + QM188
  sample_num_vec <- c(QT60,QT120,QT188,QM60,QM120,QM188,QR60,QR120,QR188)
  return(sample_num_vec)
}

nr <- sample_num_function(physeq_nr)
write.csv(nr, 'samplenumrare.csv', append=T)

sev <- sample_num_function(physeq_7)
write.csv(sev, 'sample7.csv', append=T)

one <- sample_num_function(physeq_1)
write.csv(one, 'sample1.csv', append=T)

