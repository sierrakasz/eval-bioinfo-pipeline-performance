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
QR <- length(sample_names(physeq_Q_r))

physeq_M_r <- subset_samples(physeq_rec, Pipeline == "mothur")
MR <- length(sample_names(physeq_M_r))

physeq_G_r <- subset_samples(physeq_rec, Pipeline == "MG-RAST")
GR <- length(sample_names(physeq_G_r))

physeq_Q_m <- subset_samples(physeq_mou, Pipeline == "QIIME2")
QM <- length(sample_names(physeq_Q_m))

physeq_M_m <- subset_samples(physeq_mou, Pipeline == "mothur")
MM <- length(sample_names(physeq_M_m))

physeq_G_m <- subset_samples(physeq_mou, Pipeline == "MG-RAST")
GM <- length(sample_names(physeq_G_m))

QT <- QR + QM
MT <- MR + MM
GT <- GR + GM

sample_num_vec <- c(QT, QM, QR, MT, MM, MR, GT, GM, GR)
sample_num_vec
tiff('sample_number_three_pipelines.TIF', width = 800, height = 600)
sam_num_plot <- barplot(sample_num_vec, names.arg = c("Total", "Mouth", "Rectum", 
                                      "Total", "Mouth", "Rectum","Total", "Mouth", "Rectum"),
        col = c("#5c5c8a", "#5c5c8a", "#5c5c8a","#ffb31a", "#ffb31a", "#ffb31a", 
                "#787878", "#787878", "#787878"), xlab = "Sample Area", ylab = "Number of Samples",
        ylim = c(0,375), cex.axis = 1, cex.names = 1)
text(x = sam_num_plot, y = sample_num_vec, label = sample_num_vec, pos = 3, cex = 0.8, col = "black")
legend(x='topright', legend = c("QIIME2", "mothur", "MG-RAST"), pch = 15, pt.cex = 3, cex=1, col = c("#5c5c8a", "#ffb31a", "#787878"))
dev.off()
