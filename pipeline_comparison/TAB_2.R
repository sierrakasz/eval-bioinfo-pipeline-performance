##Table 2. 
# Table 2 is a summary of the differential effect of filtering among pipelines on sample number
# and unclassified sequences. 

#load required packages
#phyloseq requires BiocManager 
# code for phyloseq install: if(!requireNamespace("BiocManager")){
#install.packages("BiocManager")
#}
#BiocManager::install("phyloseq")

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


# sample number -----------------------------------------------------------

#number of samples that remain. Subset among sample areas and pipelines.
#Length function for sample names of subsetted phyloseq object gives number of samples
#Q - QIIME2, M - mothur, G - MG-RAST, R - rectum, Mouth - m
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

#adding rectum and mouth subset samples among pipelines to get total sample number
QT <- QR + QM
MT <- MR + MM
GT <- GR + GM

#add them into all one vector
sample_num_vec <- c(QT, QM, QR, MT, MM, MR, GT, GM, GR)
sample_num_vec


# unclassified sequences --------------------------------------------------


##unclassified sequences function
#subset by sample area and pipeline
# subset out Unclassified taxa at two taxonomic levels: phylum and family
# putting together all information into a dataframe
unclassified_sequences <- function(physeq) {
  physeq_rec <- subset_samples(physeq, Sample_Area == "Rectum")
  physeq_mou <- subset_samples(physeq, Sample_Area == "Mouth")
  physeq_Q_r <- subset_samples(physeq_rec, Pipeline == "QIIME2")
  physeq_M_r <- subset_samples(physeq_rec, Pipeline == "mothur")
  physeq_G_r <- subset_samples(physeq_rec, Pipeline == "MG-RAST")
  physeq_Q_m <- subset_samples(physeq_mou, Pipeline == "QIIME2")
  physeq_M_m <- subset_samples(physeq_mou, Pipeline == "mothur")
  physeq_G_m <- subset_samples(physeq_mou, Pipeline == "MG-RAST")
  physeq_Q_r_unc = subset_taxa(physeq_Q_r, Phylum=="Unclassified")
  physeq_M_r_unc = subset_taxa(physeq_M_r, Phylum=="Unclassified")
  physeq_G_r_unc = subset_taxa(physeq_G_r, Phylum=="Unclassified")
  physeq_Q_m_unc = subset_taxa(physeq_Q_m, Phylum=="Unclassified")
  physeq_M_m_unc = subset_taxa(physeq_M_m, Phylum=="Unclassified")
  physeq_G_m_unc = subset_taxa(physeq_G_m, Phylum=="Unclassified")
  physeq_Q_r_uncf = subset_taxa(physeq_Q_r, Family=="Unclassified")
  physeq_M_r_uncf = subset_taxa(physeq_M_r, Family=="Unclassified")
  physeq_G_r_uncf = subset_taxa(physeq_G_r, Family=="Unclassified")
  physeq_Q_m_uncf = subset_taxa(physeq_Q_m, Family=="Unclassified")
  physeq_M_m_uncf = subset_taxa(physeq_M_m, Family=="Unclassified")
  physeq_G_m_uncf = subset_taxa(physeq_G_m, Family=="Unclassified")
  a <- sum(sample_sums(physeq_Q_r))
  b <- sum(sample_sums(physeq_M_r))
  c <- sum(sample_sums(physeq_G_r))
  d <- sum(sample_sums(physeq_Q_m))
  e <- sum(sample_sums(physeq_M_m))
  f <- sum(sample_sums(physeq_G_m))
  aa <- sum(sample_sums(physeq_Q_r_unc))
  bb <- sum(sample_sums(physeq_M_r_unc))
  cc <- sum(sample_sums(physeq_G_r_unc))
  dd <- sum(sample_sums(physeq_Q_m_unc))
  ee <- sum(sample_sums(physeq_M_m_unc))
  ff <- sum(sample_sums(physeq_G_m_unc))
  aaa <- sum(sample_sums(physeq_Q_r_uncf))
  bbb <- sum(sample_sums(physeq_M_r_uncf))
  ccc <- sum(sample_sums(physeq_G_r_uncf))
  ddd <- sum(sample_sums(physeq_Q_m_uncf))
  eee <- sum(sample_sums(physeq_M_m_uncf))
  fff <- sum(sample_sums(physeq_G_m_uncf))
  df <- data.frame('Total' = c(a,b,c,d,e,f), 'Phylum' = c(aa,bb,cc,dd,ee,ff), 
                   'Family' = c(aaa,bbb,ccc,ddd,eee,fff))
}

df <- unclassified_sequences(physeq)

#all information was added to a table in excel for formatting into the manuscript. 
