##Table 3
# Table 3 is a summary random forest classification error among pipelines, taxonomic levels. 
# Both sample area (mouth and rectum) and manner of death (homicide, suicide, accident, natural)
# were classified. 

#load required packages
#phyloseq requires BiocManager 
# code for phyloseq install: if(!requireNamespace("BiocManager")){
#install.packages("BiocManager")
#}
#BiocManager::install("phyloseq")
library(ggforce)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(HMP)
library(phyloseq)
library(plyr)
library(PMCMR)
library(randomForest)
library(rsample)
library(tidyverse)
library(phyloseq)
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

#subsample out minimum library sizes
physeq_nr_phylum <- tax_glom(physeq_nr, taxrank = 'Phylum')
physeq_nr_family <- tax_glom(physeq_nr, taxrank = 'Family')
physeq_7_phylum <- tax_glom(physeq_7, taxrank = 'Phylum')
physeq_7_family <- tax_glom(physeq_7, taxrank = 'Family')
physeq_1_phylum <- tax_glom(physeq_1, taxrank = 'Phylum')
physeq_1_family <- tax_glom(physeq_1, taxrank = 'Family')

#within minimum library size, subsample out by sample size
samples_rel_abundance <- function(physeq) {
  return(list(
    physeq_60 <- subset_samples(physeq, Subsample == '60'),
    physeq_120 <- subset_samples(physeq, Subsample == '120'),
    physeq_188 <- subset_samples(physeq, Subsample == '188')))
}

#make lists with each sample size within each minimum library size and taxonomic level
physeq_nr_phy_list <- samples_rel_abundance(physeq_nr_phylum)
physeq_7_phy_list <- samples_rel_abundance(physeq_7_phylum)
physeq_1_phy_list <- samples_rel_abundance(physeq_1_phylum)
physeq_nr_fam_list <- samples_rel_abundance(physeq_nr_family)
physeq_7_fam_list <- samples_rel_abundance(physeq_7_family)
physeq_1_fam_list <- samples_rel_abundance(physeq_1_family)

#function for running random forest classification among sample areas
#output is the classification error matrix and overall classification error
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

#same as above, except classifying manner of death
# random forest classification error and overall error
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

#outputs are summarized in supplemental tables, and reformatted into excel for manuscript table