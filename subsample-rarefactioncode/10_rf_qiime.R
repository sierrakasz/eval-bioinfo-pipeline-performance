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

physeq_nr_phylum <- tax_glom(physeq_nr, taxrank = 'Phylum')
physeq_nr_family <- tax_glom(physeq_nr, taxrank = 'Family')
physeq_7_phylum <- tax_glom(physeq_7, taxrank = 'Phylum')
physeq_7_family <- tax_glom(physeq_7, taxrank = 'Family')
physeq_1_phylum <- tax_glom(physeq_1, taxrank = 'Phylum')
physeq_1_family <- tax_glom(physeq_1, taxrank = 'Family')

samples_rel_abundance <- function(physeq) {
  return(list(
    physeq_60 <- subset_samples(physeq, Subsample == '60'),
    physeq_120 <- subset_samples(physeq, Subsample == '120'),
    physeq_188 <- subset_samples(physeq, Subsample == '188')))
}

physeq_nr_phy_list <- samples_rel_abundance(physeq_nr_phylum)
physeq_7_phy_list <- samples_rel_abundance(physeq_7_phylum)
physeq_1_phy_list <- samples_rel_abundance(physeq_1_phylum)
physeq_nr_fam_list <- samples_rel_abundance(physeq_nr_family)
physeq_7_fam_list <- samples_rel_abundance(physeq_7_family)
physeq_1_fam_list <- samples_rel_abundance(physeq_1_family)

#testing within subsample and rarefy 
#which one is better, for predicting sample area and MoD?

random_foresting_sample_area <- function(data_phy) {
  fd = data_phy
  predictors=t(otu_table(fd))
  dim(predictors)
  resp <- as.factor(sample_data(fd)$Sample_Area)
  rf.data <- data.frame(resp, predictors)
  Forest <- randomForest(resp~., data=rf.data, ntree=1000)
  return(Forest)
}

for1 <- random_foresting_sample_area(physeq_nr_phy_list[[1]])
for2 <- random_foresting_sample_area(physeq_nr_phy_list[[2]])
for3 <- random_foresting_sample_area(physeq_nr_phy_list[[3]])
for4 <- random_foresting_sample_area(physeq_nr_fam_list[[1]])
for5 <- random_foresting_sample_area(physeq_nr_fam_list[[2]])
for6 <- random_foresting_sample_area(physeq_nr_fam_list[[3]])
for7 <- random_foresting_sample_area(physeq_7_phy_list[[1]])
for8 <- random_foresting_sample_area(physeq_7_phy_list[[2]])
for9 <- random_foresting_sample_area(physeq_7_phy_list[[3]])
for10 <- random_foresting_sample_area(physeq_7_fam_list[[1]])
for11 <- random_foresting_sample_area(physeq_7_fam_list[[2]])
for12 <- random_foresting_sample_area(physeq_7_fam_list[[3]])
for13 <- random_foresting_sample_area(physeq_1_phy_list[[1]])
for14 <- random_foresting_sample_area(physeq_1_phy_list[[2]])
for15 <- random_foresting_sample_area(physeq_1_phy_list[[3]])
for16 <- random_foresting_sample_area(physeq_1_fam_list[[1]])
for17 <- random_foresting_sample_area(physeq_1_fam_list[[2]])
for18 <- random_foresting_sample_area(physeq_1_fam_list[[3]])


random_foresting_MoD <- function(data_phy) {
  fd = data_phy
  predictors=t(otu_table(fd))
  dim(predictors)
  resp <- as.factor(sample_data(fd)$MoD)
  rf.data <- data.frame(resp, predictors)
  Forest <- randomForest(resp~., data=rf.data, ntree=1000)
  return(Forest)
}


for19 <- random_foresting_MoD(physeq_nr_phy_list[[1]])
for20 <- random_foresting_MoD(physeq_nr_phy_list[[2]])
for21 <- random_foresting_MoD(physeq_nr_phy_list[[3]])
for22 <- random_foresting_MoD(physeq_nr_fam_list[[1]])
for23 <- random_foresting_MoD(physeq_nr_fam_list[[2]])
for24 <- random_foresting_MoD(physeq_nr_fam_list[[3]])
for25 <- random_foresting_MoD(physeq_7_phy_list[[1]])
for26 <- random_foresting_MoD(physeq_7_phy_list[[2]])
for27 <- random_foresting_MoD(physeq_7_phy_list[[3]])
for28 <- random_foresting_MoD(physeq_7_fam_list[[1]])
for29 <- random_foresting_MoD(physeq_7_fam_list[[2]])
for30 <- random_foresting_MoD(physeq_7_fam_list[[3]])
for31 <- random_foresting_MoD(physeq_1_phy_list[[1]])
for32 <- random_foresting_MoD(physeq_1_phy_list[[2]])
for33 <- random_foresting_MoD(physeq_1_phy_list[[3]])
for34 <- random_foresting_MoD(physeq_1_fam_list[[1]])
for35 <- random_foresting_MoD(physeq_1_fam_list[[2]])
for36 <- random_foresting_MoD(physeq_1_fam_list[[3]])

forest_predictors_nr <- function(forest) {
  imp <- importance(forest)
  imp <- data.frame(predictors = rownames(imp), imp)
  imp.sort <- arrange(imp, desc(MeanDecreaseGini))
  imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)
  imp.20 <- imp.sort[1:20, ]
  tax <- data.frame(tax_table(physeq_nr))
  head(tax)
  tax <- tax %>% select(Kingdom, Phylum, Class, Order, Family, Genus)
  tax$predictors <- rownames(tax)
  imp.20 <- merge(imp.20, tax)
  write.table(imp.20, sep=",", "random_forest_predictors_qiime_nr.csv", append = TRUE)
  return(imp.20)
}

for_list_nr <- list(for1, for2, for3, for4, for5, for6, 
                    for19, for20, for21, for22, for23, for24)

print(for_list_nr)

for(i in 1:length(for_list_nr)) {
  forest_predictors_nr(for_list_nr[[i]])
}

forest_predictors_7 <- function(forest) {
  imp <- importance(forest)
  imp <- data.frame(predictors = rownames(imp), imp)
  imp.sort <- arrange(imp, desc(MeanDecreaseGini))
  imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)
  imp.20 <- imp.sort[1:20, ]
  tax <- data.frame(tax_table(physeq_7))
  head(tax)
  tax <- tax %>% select(Kingdom, Phylum, Class, Order, Family, Genus)
  tax$predictors <- rownames(tax)
  imp.20 <- merge(imp.20, tax)
  write.table(imp.20, sep=",", "random_forest_predictors_qiime_7.csv", append = TRUE)
  return(imp.20)
}

for_list_7 <- list(for7, for8, for9, for10, for11, for12, 
                    for25, for26, for27, for28, for29, for30)

print(for_list_7)

for(i in 1:length(for_list_7)) {
  forest_predictors_7(for_list_7[[i]])
}

forest_predictors_1 <- function(forest) {
  imp <- importance(forest)
  imp <- data.frame(predictors = rownames(imp), imp)
  imp.sort <- arrange(imp, desc(MeanDecreaseGini))
  imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)
  imp.20 <- imp.sort[1:20, ]
  tax <- data.frame(tax_table(physeq_1))
  head(tax)
  tax <- tax %>% select(Kingdom, Phylum, Class, Order, Family, Genus)
  tax$predictors <- rownames(tax)
  imp.20 <- merge(imp.20, tax)
  write.table(imp.20, sep=",", "random_forest_predictors_qiime_1.csv", append = TRUE)
  return(imp.20)
}

for_list_1 <- list(for13, for14, for15, for16, for17, for18, 
                   for31, for32, for33, for34, for35, for36)

print(for_list_1)

for(i in 1:length(for_list_1)) {
  forest_predictors_1(for_list_1[[i]])
}
