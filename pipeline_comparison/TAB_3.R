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

# import data -------------------------------------------------------------
###tax group is a file that combines the OTU file and taxonomy file. 
##combining the files made it easier for collapsing differences among the pipeline outputs
#rarefied to 1,000 sequences
tax_group <- read.csv(text=GET('https://github.com/sierrakasz/eval-bioinfo-pipeline-performance/blob/master/pipeline_comparison/tax_group.csv'), header=T)

#collapsing down same taxa names to get true abundance sum
tax_group <- tax_group %>% group_by(Kingdom, Phylum, Class, Order, Family, Genus) %>% summarise_all(funs(sum))
tax_group <- tax_group %>%  ungroup()

#this function specifically collapses the OTU table and taxonomy file 
#when the taxonomy is the same name among pipelines
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

#import metadata and combine into one phyloseq object
metadata=(read.csv("Metadata/HPMMMeta_r_merge.csv",header=TRUE))
sampdat=sample_data(metadata)
sample_names(sampdat)=metadata$SampleID
physeq=merge_phyloseq(physeq_all, sampdat)
physeq

# subsetting samples among pipelines, and taxonomic levels
# q - QIIME2, m - mothur, g - MG-RAST, p - Phylum, f - family
physeq_q <- subset_samples(physeq, Pipeline == 'QIIME2')
physeq_m <- subset_samples(physeq, Pipeline == 'mothur')
physeq_mg <- subset_samples(physeq, Pipeline == 'MG-RAST')

physeq_q_p <- tax_glom(physeq_q, taxrank = 'Phylum')
physeq_m_p <- tax_glom(physeq_m, taxrank = 'Phylum')
physeq_mg_p <- tax_glom(physeq_mg, taxrank = 'Phylum')
physeq_q_f <- tax_glom(physeq_q, taxrank = 'Family')
physeq_m_f <- tax_glom(physeq_m, taxrank = 'Family')
physeq_mg_f <- tax_glom(physeq_mg, taxrank = 'Family')

#function for running random forest classification among sample areas
#output is the classification error matrix and overall classification error
random_foresting_sample_area <- function(data_phy) {
  fd = data_phy
  predictors=t(otu_table(fd))
  dim(predictors)
  resp <- as.factor(sample_data(fd)$Sample_Area)
  rf.data <- data.frame(resp, predictors)
  Forest <- randomForest(resp~., data=rf.data)
  return(Forest)
}

forest_q_p <- random_foresting_sample_area(physeq_q_p)
forest_m_p <- random_foresting_sample_area(physeq_m_p)
forest_mg_p <- random_foresting_sample_area(physeq_mg_p)
forest_q_f <- random_foresting_sample_area(physeq_q_f)
forest_m_f <- random_foresting_sample_area(physeq_m_f)
forest_mg_f <- random_foresting_sample_area(physeq_mg_f)

#function for determing indicator taxa for random forest classification among sample areas
#output is top 20 predictors based on mean decrease in Gini
#function writes predictors to a csv file
forest_predictors_sample_area <- function(forest) {
  imp <- importance(forest)
  imp <- data.frame(predictors = rownames(imp), imp)
  imp.sort <- arrange(imp, desc(MeanDecreaseGini))
  imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)
  imp.20 <- imp.sort[1:20, ]
  tax <- data.frame(tax_table(physeq))
  head(tax)
  tax <- tax %>% select(Kingdom, Phylum, Class, Order, Family)
  tax$predictors <- rownames(tax)
  imp.20 <- merge(imp.20, tax)
  write.table(imp.20, sep=",", "random_forest_predictors_pipe_sa.csv", append = TRUE)
  return(imp.20)
}


pred_q_f <- forest_predictors_sample_area(forest_q_f)
pred_m_f <- forest_predictors_sample_area(forest_m_f)
pred_mg_f <- forest_predictors_sample_area(forest_mg_f)

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

forest_q_p <- random_foresting_MoD(physeq_q_p)
forest_m_p <- random_foresting_MoD(physeq_m_p)
forest_mg_p <- random_foresting_MoD(physeq_mg_p)
forest_q_f <- random_foresting_MoD(physeq_q_f)
forest_m_f <- random_foresting_MoD(physeq_m_f)
forest_mg_f <- random_foresting_MoD(physeq_mg_f)

#indicator taxa for manner of death classification 
forest_predictors_MoD <- function(forest) {
  imp <- importance(forest)
  imp <- data.frame(predictors = rownames(imp), imp)
  imp.sort <- arrange(imp, desc(MeanDecreaseGini))
  imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)
  imp.20 <- imp.sort[1:20, ]
  tax <- data.frame(tax_table(physeq))
  head(tax)
  tax <- tax %>% select(Kingdom, Phylum, Class, Order, Family)
  tax$predictors <- rownames(tax)
  imp.20 <- merge(imp.20, tax)
  write.table(imp.20, sep=",", "random_forest_predictors_pipe_md.csv", append = TRUE)
  return(imp.20)
}


pred_q_f <- forest_predictors_MoD(forest_q_f)
pred_m_f <- forest_predictors_MoD(forest_m_f)
pred_mg_f <- forest_predictors_MoD(forest_mg_f)

