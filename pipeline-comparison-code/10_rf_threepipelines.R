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

random_forest_data <- function(physeq) {
  return(list(pipeline_mouth <- subset_samples(physeq, Sample_Area == 'Mouth'),
  pipeline_rectum <- subset_samples(physeq, Sample_Area == 'Rectum'),
  sample_area_mot <- subset_samples(physeq, Pipeline == 'mothur'),
  sample_area_mg <- subset_samples(physeq, Pipeline == 'MG-RAST'),
  sample_area_qim <- subset_samples(physeq, Pipeline == 'QIIME2'),
  rec_mot_mod <- subset_samples(sample_area_mot, Sample_Area=='Rectum'),
  rec_mg_mod <- subset_samples(sample_area_mg, Sample_Area=='Rectum'),
  rec_qim_mod <- subset_samples(sample_area_qim, Sample_Area=='Rectum'),
  mou_mot_mod <- subset_samples(sample_area_mot, Sample_Area=='Mouth'),
  mou_mg_mod <- subset_samples(sample_area_mg, Sample_Area=='Mouth'),
  mou_qim_mod <- subset_samples(sample_area_qim, Sample_Area=='Mouth')))
}

rf_data_list <- random_forest_data(physeq)

#random forest
#functions
random_foresting_pipeline <- function(data_phy) {
  fd = data_phy
  predictors=t(otu_table(fd))
  dim(predictors)
  resp <- as.factor(sample_data(fd)$Pipeline)
  rf.data <- data.frame(resp, predictors)
  Forest <- randomForest(resp~., data=rf.data, ntree=1000)
  return(Forest)
}

for1 <- random_foresting_pipeline(rf_data_list[[1]])
for2 <- random_foresting_pipeline(rf_data_list[[2]])

random_foresting_sample_area <- function(data_phy) {
  fd = data_phy
  predictors=t(otu_table(fd))
  dim(predictors)
  resp <- as.factor(sample_data(fd)$Sample_Area)
  rf.data <- data.frame(resp, predictors)
  Forest <- randomForest(resp~., data=rf.data, ntree=1000)
  return(Forest)
}

for3 <- random_foresting_sample_area(rf_data_list[[3]])
for4 <- random_foresting_sample_area(rf_data_list[[4]])
for5 <- random_foresting_sample_area(rf_data_list[[5]])

random_foresting_MoD <- function(data_phy) {
  fd = data_phy
  predictors=t(otu_table(fd))
  dim(predictors)
  resp <- as.factor(sample_data(fd)$MoD)
  rf.data <- data.frame(resp, predictors)
  Forest <- randomForest(resp~., data=rf.data, ntree=1000)
  return(Forest)
}

for6 <- random_foresting_MoD(rf_data_list[[6]])
for7 <- random_foresting_MoD(rf_data_list[[7]])
for8 <- random_foresting_MoD(rf_data_list[[8]])
for9 <- random_foresting_MoD(rf_data_list[[9]])
for10 <- random_foresting_MoD(rf_data_list[[10]])
for11 <- random_foresting_MoD(rf_data_list[[11]])

forest_predictors <- function(forest) {
  imp <- importance(forest)
  imp <- data.frame(predictors = rownames(imp), imp)
  imp.sort <- arrange(imp, desc(MeanDecreaseGini))
  imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)
  imp.20 <- imp.sort[1:20, ]
  tax <- data.frame(tax_table(physeq))
  head(tax)
  tax <- tax %>% select(Kingdom, Phylum, Class, Order, Family, Genus)
  tax$predictors <- rownames(tax)
  imp.20 <- merge(imp.20, tax)
  write.table(imp.20, sep=",", "random_forest_predictors.csv", append = TRUE)
  return(imp.20)
}

for_list <- list(for1, for2, for3, for4, for5, for6, for7, for8, for9, for10, for11)

for(i in 1:length(for_list)) {
  forest_predictors(for_list[[i]])
}


