rm(list=ls())
getwd()
setwd("C:/Users/sierr/Documents/ThesisProject_Pipelines/")

#packages
library(phyloseq)
library(plyr)
library(tidyverse)

#import mothur
import_mothur_files <- function(shared, tax) {
  physeq_mothur <- import_mothur(mothur_shared_file = shared,
                                 mothur_constaxonomy_file = tax)
  colnames(tax_table(physeq_mothur)) <- c("Kingdom", "Phylum", "Class", 
                                          "Order", "Family",  "Genus")
  return(physeq_mothur)
}

physeq_mothur <- import_mothur_files("mothur biom files/mothur.all.shared", 
                                     "mothur biom files/mothur.all.taxonomy")
sample_names(physeq_mothur) <- paste("M_188_", sample_names(physeq_mothur), sep="")

#import qiime
import_qiime_files <- function(biom) {
  biom_qiime <- import_biom(biom)
  a <- data.frame(tax_table(biom_qiime))
  b <- a %>% 
    mutate(Rank1 = str_replace(Rank1, "D_0__", "")) %>% 
    mutate(Rank2 = str_replace(Rank2, "D_1__", "")) %>% 
    mutate(Rank3 = str_replace(Rank3, "D_2__", "")) %>% 
    mutate(Rank4 = str_replace(Rank4, "D_3__", "")) %>% 
    mutate(Rank5 = str_replace(Rank5, "D_4__", "")) %>% 
    mutate(Rank6 = str_replace(Rank6, "D_5__", ""))
  c <- b[,c("Rank1", "Rank2", "Rank3", "Rank4", "Rank5", "Rank6")]
  c <- as.matrix(c)
  d <- data.frame(otu_table(biom_qiime))
  OTUq=otu_table(d, taxa_are_rows=TRUE)
  TAXq=tax_table(c)
  taxa_names(TAXq)=row.names(OTUq)
  physeq_qiime = merge_phyloseq(OTUq, TAXq)
  colnames(tax_table(physeq_qiime)) <- c("Kingdom", "Phylum", "Class",
                                         "Order","Family", "Genus")
  return(physeq_qiime)
}

physeq_qiime <- import_qiime_files("qiime biom files/table-with-taxonomy.biom")
sample_names(physeq_qiime) <- paste("Q_188_", sample_names(physeq_qiime), sep="")

#import MG-RAST
import_MGRAST <- function(file) {
  file_name <- read.csv(file)
  otu <- select(file_name, contains('WCME'))
  tax <-select(file_name, domain, phylum, className, order, family, genus)
  tax <- as.matrix(tax)
  OTU=otu_table(otu, taxa_are_rows=TRUE)
  TAX=tax_table(tax)
  taxa_names(TAX)=row.names(OTU)
  physeq_MG=phyloseq(OTU,TAX)
  colnames(tax_table(physeq_MG)) <- c("Kingdom", "Phylum", "Class", 
                                      "Order", "Family",  "Genus")
  return(physeq_MG)
}
physeq_MG <- import_MGRAST("MG-Rast biom files\\HPMMSAll.csv")
sample_names(physeq_MG) <- paste("G_188_", sample_names(physeq_MG), sep="")

#merging phyloseq objects together
merge = merge_phyloseq(physeq_qiime,physeq_mothur,physeq_MG)

#rarefaction
physeq_rare <- rarefy_even_depth(merge, 
                                 rngseed = 711, sample.size = 1000)

#tidying data to insure that taxa are consistent 
tidying_data <- function(physeq) {
  x <- data.frame(otu_table(physeq))
  y <- data.frame(tax_table(physeq))
  x$names <- rownames(x)
  y$names <- rownames(y)
  z <- merge(x,y)
  z <- z %>% select(-contains("names"))
  tax_group <- z %>% group_by(Kingdom, Phylum, Class, Order, Family, Genus) %>% summarise_all(funs(sum))
  tax_group <- tax_group %>%  ungroup()
  return(tax_group)
}
tax_group <- tidying_data(physeq_rare)
write.csv(tax_group, "tax_group.csv")