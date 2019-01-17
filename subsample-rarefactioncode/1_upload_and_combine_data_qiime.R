rm(list=ls())
getwd()
setwd("C:/Users/sierr/Documents/ThesisProject_Pipelines/")

#packages
library(phyloseq)
library(plyr)
library(tidyverse)

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

physeq_qiime_120 <- import_qiime_files("qiime biom files/table-with-taxonomy_120.biom")
sample_names(physeq_qiime_120) <- paste("Q_120_", sample_names(physeq_qiime_120), sep="")

physeq_qiime_60 <- import_qiime_files("qiime biom files/table-with-taxonomy_60.biom")
sample_names(physeq_qiime_60) <- paste("Q_60_", sample_names(physeq_qiime_60), sep="")

#merging phyloseq objects together
merge = merge_phyloseq(physeq_qiime,physeq_qiime_120, physeq_qiime_60)

#rarefaction
physeq_norare <- merge

physeq_rare_7000 <- rarefy_even_depth(merge, 
                                 rngseed = 711, sample.size = 7000)

physeq_rare_1000 <- rarefy_even_depth(merge, 
                                 rngseed = 711, sample.size = 1000)

sample_names(physeq_nr) <- paste("Norare_", sample_names(physeq_nr), sep="")
sample_names(physeq_7) <- paste("S7000_", sample_names(physeq_7), sep="")
sample_names(physeq_1) <- paste("O1000_", sample_names(physeq_1), sep="")

merger= merge_phyloseq(physeq_nr, physeq_7, physeq_1)
metadata_rare=(read.csv("Metadata/HPMMMeta_r_merge_rare_levels.csv",header=TRUE))
sampdat_rare=sample_data(metadata_rare)
sampdat_rare$Rarefy <- factor(sampdat_rare$Rarefy, levels = c("No Rarefaction", "7000", "1000"))
sample_names(sampdat_rare)=metadata_rare$SampleID
merger = merge_phyloseq(merger, sampdat_rare)
sample_data(merger)$Rarefy <- factor(sample_data(merger)$Rarefy, levels = c("No Rarefaction", "7000", "1000"))

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

tax_group <- tidying_data(physeq_norare)
write.csv(tax_group, "tax_group_norare.csv")

tax_group_a <- tidying_data(physeq_rare_7000)
write.csv(tax_group_a, "tax_group_7000.csv")

tax_group_c <- tidying_data(physeq_rare_1000)
write.csv(tax_group_c, "tax_group_1000.csv")

tax_group_merge <- tidying_data(merger)
write.csv(tax_group_merge, "tax_group_rare.csv")
