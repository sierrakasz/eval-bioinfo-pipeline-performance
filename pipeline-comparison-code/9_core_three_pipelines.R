rm(list=ls())
getwd()
setwd("C:/Users/sierr/Documents/ThesisProject_Pipelines/")

#packages
library(ggforce)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(PMCMR)
library(phyloseq)
library(plyr)
library(tidyverse)
library(exactRankTests)
library(nlme)
library(UpSetR)

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

#core
preparing_data_subsample_upset <- function(physeq) {
  otus <- data.frame(otu_table(physeq))
  totus <- data.frame(t(otus))
  totus$SampleID <- rownames(totus)
  met <- metadata[,c('SampleID', 'Pipeline')]
  mtotus <- merge(totus, met)
  mtotus <- mtotus[,-1]
  total <- as.vector(colSums(Filter(is.numeric, mtotus)))
  new_df <- mtotus %>% group_by(Pipeline) %>% summarise_all(funs(sum))
  new_df <- data.frame(t(new_df))
  colnames(new_df) <- as.character(unlist(new_df[1,]))
  new_df = new_df[-1, ]
  new_df$OTU <- rownames(new_df)
  rownames(new_df) <- NULL
  Upset <- cbind(new_df, total)
  return(upset <- Upset[,c(4,1,2,3,5)])
}

physeq_upset <- preparing_data_subsample_upset(physeq)

#export to excel to change format
#Change to binary per OTU
#OTU number,(binary), Total
write.csv(physeq_upset,"pipeline_upset.csv", row.names=FALSE)
pipeCore=read.csv("pipeline_upset.csv",header=TRUE)

tiff("upset_pipeline.TIF", width = 1200, height = 1200)
upset(pipeCore, 
      order.by = c("degree", 'freq'), 
      point.size = 8, 
      line.size = 4,  
      empty.intersections = "on", 
      mainbar.y.label = "Number of Taxa", 
      sets.x.label = "Core Taxa",
      queries = list(list(query = intersects, params = list("MG.RAST","mothur"), active = T, color= "#F2AD00"),
                     list(query = intersects, params = list("MG.RAST"), active = T, color = "#FF0000"),
                     list(query = intersects, params = list("mothur"), active = T, color = '#00A08A'),
                     list(query = intersects, params = list("QIIME2"), active = T, color = '#663399'),
                     list(query = intersects, params = list("QIIME2","mothur"), active = T, color = "#5BBCD6"),
                     list(query = intersects, params = list("QIIME2","MG.RAST"), active = T, color = "#F98400"),
                     list(query = intersects, params = list("QIIME2","MG.RAST", "mothur"), active = T, color = "#ECCBAE")),
      sets.bar.color = "#56B4E9",
      text.scale = 6)
dev.off()
