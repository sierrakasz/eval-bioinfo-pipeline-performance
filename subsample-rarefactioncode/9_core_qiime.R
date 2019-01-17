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
library(wesanderson)

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

tax_group_rare <- read.csv("tax_group_rare.csv")
tax_group_rare <- tax_group_1 %>% group_by(Kingdom, Phylum, Class, Order, Family, Genus) %>% summarise_all(funs(sum))
tax_group_rare <- tax_group_1 %>%  ungroup()

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
physeq_rare <- regroup_physeq_object(tax_group_rare)

#import metadata and combine
metadata=(read.csv("Metadata/HPMMMeta_r_merge_rare.csv",header=TRUE))
sampdat=sample_data(metadata)
sample_names(sampdat)=metadata$SampleID
sampdat$Subsample <- factor(sampdat$Subsample, levels = c('60', '120', '188'))

physeq_nr=merge_phyloseq(physeq_nr, sampdat)
physeq_7=merge_phyloseq(physeq_7, sampdat)
physeq_1=merge_phyloseq(physeq_1, sampdat)

#core
preparing_data_subsample_upset <- function(physeq) {
  otus <- data.frame(otu_table(physeq))
  totus <- data.frame(t(otus))
  totus$SampleID <- rownames(totus)
  met <- metadata[,c('SampleID', 'Subsample')]
  mtotus <- merge(totus, met)
  mtotus <- mtotus[,-1]
  mtotus$Subsample <- factor(mtotus$Subsample, levels = c('60', '120', '188'))
  total <- as.vector(colSums(Filter(is.numeric, mtotus)))
  mtotus$Subsample <- factor(mtotus$Subsample, levels = c('60', '120', '188'))
  new_df <- mtotus %>% group_by(Subsample) %>% summarise_all(funs(sum))
  new_df <- data.frame(t(new_df))
  colnames(new_df) <- as.character(unlist(new_df[1,]))
  new_df = new_df[-1, ]
  new_df$OTU <- rownames(new_df)
  rownames(new_df) <- NULL
  Upset <- cbind(new_df, total)
  return(upset <- Upset[,c(4,1,2,3,5)])
}

physeq_nr_upset <- preparing_data_subsample_upset(physeq_nr)
physeq_7_upset <- preparing_data_subsample_upset(physeq_7)
physeq_1_upset <- preparing_data_subsample_upset(physeq_1)

#export to excel to change format
#Change to binary per OTU
#OTU number,(binary), Total
write.csv(physeq_nr_upset,"upset_nr.csv", row.names=FALSE)
write.csv(physeq_7_upset,"upset_7.csv", row.names=FALSE)
write.csv(physeq_1_upset,"upset_1.csv", row.names=FALSE)

#Add two columns sequences and samples 
samples_andseqs_upset <- function(physeq) {
  otu_df <- data.frame(otu_table(physeq))
  sequences_df <- rowSums(otu_df)
  colnames(otu_df)
  for(i in 1:length(rownames(otu_df))){
    for(j in 1:length(colnames(otu_df))) {
      if((otu_df[i,j] > 0) == TRUE ) {
        otu_df[i,j] <- 1
      }
    }
  }
  sample_df <- rowSums(otu_df)
  return(list(sample_df, sequences_df))
}

physeq_nr_sampandseq <- samples_andseqs_upset(physeq_nr)
physeq_7_sampandseq <- samples_andseqs_upset(physeq_7)
physeq_1_sampandseq <- samples_andseqs_upset(physeq_1)

pipeCore_nr=read.csv("upset_nr.csv",header=TRUE)
pipeCore_7=read.csv("upset_7.csv",header=TRUE)
pipeCore_1=read.csv("upset_1.csv",header=TRUE)

pipeCore_nr <- cbind(pipeCore_nr, physeq_nr_sampandseq[[1]], physeq_nr_sampandseq[[2]])
names(pipeCore_nr)[names(pipeCore_nr) == 'physeq_nr_sampandseq[[1]]'] <- 'Samples'
names(pipeCore_nr)[names(pipeCore_nr) == 'physeq_nr_sampandseq[[2]]'] <- 'Sequences'

pipeCore_7 <- cbind(pipeCore_7, physeq_7_sampandseq[[1]], physeq_7_sampandseq[[2]])
names(pipeCore_7)[names(pipeCore_7) == 'physeq_7_sampandseq[[1]]'] <- 'Samples'
names(pipeCore_7)[names(pipeCore_7) == 'physeq_7_sampandseq[[2]]'] <- 'Sequences'

pipeCore_1 <- cbind(pipeCore_1, physeq_1_sampandseq[[1]], physeq_1_sampandseq[[2]])
names(pipeCore_1)[names(pipeCore_1) == 'physeq_1_sampandseq[[1]]'] <- 'Samples'
names(pipeCore_1)[names(pipeCore_1) == 'physeq_1_sampandseq[[2]]'] <- 'Sequences'

tiff("upset_nr.TIF", width = 550, height = 750)
theme_set(theme_bw(base_size = 14))
myplot <- function(mydata, x, y) {
  plot <- (ggplot(data = mydata, aes_string(x = x, y = y, colour = "color")) + 
             geom_point() + scale_color_identity() + 
             theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))) +
    ylim(0,1500) +
    xlim(0,20) +
    theme_bw()
}

upset(pipeCore_nr, 
      order.by = c("degree", 'freq'), 
      point.size = 4, 
      line.size = 2,  
      empty.intersections = "on", 
      mainbar.y.label = "Number of Taxa", 
      sets.x.label = "Core Taxa",
      queries = list(list(query = intersects, params = list("Subsample_188","Subsample_120"), active = T, color= "#5BBCD6"),
                     list(query = intersects, params = list("Subsample_120"), active = T, color = "#FF0000"),
                     list(query = intersects, params = list("Subsample_188"), active = T, color = '#00A08A'),
                     list(query = intersects, params = list("Subsample_60","Subsample_120"), active = T, color = "#F2AD00"),
                     list(query = intersects, params = list("Subsample_60","Subsample_188"), active = T, color = "#F98400"),
                     list(query = intersects, params = list("Subsample_60","Subsample_120", "Subsample_188"), active = T, color = "#ECCBAE")),
      sets.bar.color = "#56B4E9",
      text.scale = 1.5,
      attribute.plots=list(gridrows = 60, ncols = 1, 
                           plots = list(list(plot = myplot, x = "Samples", y = "Sequences", queries = TRUE))))

dev.off()


#rarefaction
metadata_rare=(read.csv("Metadata/HPMMMeta_r_merge_rare_levels.csv",header=TRUE))
sampdat_rare=sample_data(metadata_rare)
sampdat_rare$Rarefy <- factor(sampdat_rare$Rarefy, levels = c("No Rarefaction", "7000", "1000"))
sample_names(sampdat_rare)=metadata_rare$SampleID
merger = merge_phyloseq(tax_group_rare, sampdat_rare)
sample_data(merger)$Rarefy <- factor(sample_data(merger)$Rarefy, levels = c("No Rarefaction", "7000", "1000"))

#core
preparing_data_rarefy_upset <- function(physeq) {
  otus <- data.frame(otu_table(physeq))
  totus <- data.frame(t(otus))
  totus$SampleID <- rownames(totus)
  met <- metadata_rare[,c('SampleID', 'Rarefy')]
  mtotus <- merge(totus, met)
  mtotus <- mtotus[,-1]
  mtotus$Rarefy <- factor(mtotus$Rarefy, levels = c('No Rarefaction', '7000', '1000'))
  total <- as.vector(colSums(Filter(is.numeric, mtotus)))
  mtotus$Rarefy <- factor(mtotus$Rarefy, levels = c('No Rarefaction', '7000', '1000'))
  new_df <- mtotus %>% group_by(Rarefy) %>% summarise_all(funs(sum))
  new_df <- data.frame(t(new_df))
  colnames(new_df) <- as.character(unlist(new_df[1,]))
  new_df = new_df[-1, ]
  new_df$OTU <- rownames(new_df)
  rownames(new_df) <- NULL
  Upset <- cbind(new_df, total)
  return(upset <- Upset[,c(4,1,2,3,5)])
}

merger_upset <- preparing_data_rarefy_upset(merger)

#export to excel to change format
#Change to binary per OTU
#OTU number,(binary), Total
write.csv(merge_upset,"merge_upset.csv", row.names=FALSE)

#Add two columns sequences and samples 
samples_andseqs_upset <- function(physeq) {
  otu_df <- data.frame(otu_table(physeq))
  sequences_df <- rowSums(otu_df)
  colnames(otu_df)
  for(i in 1:length(rownames(otu_df))){
    for(j in 1:length(colnames(otu_df))) {
      if((otu_df[i,j] > 0) == TRUE ) {
        otu_df[i,j] <- 1
      }
    }
  }
  sample_df <- rowSums(otu_df)
  return(list(sample_df, sequences_df))
}

merge_sampandseq <- samples_andseqs_upset(merge)

pipeCore_merge=read.csv("merge_upset.csv",header=TRUE)

pipeCore_merge <- cbind(pipeCore_merge, merge_sampandseq[[1]], merge_sampandseq[[2]])
names(pipeCore_merge)[names(pipeCore_merge) == 'merge_sampandseq[[1]]'] <- 'Samples'
names(pipeCore_merge)[names(pipeCore_merge) == 'merge_sampandseq[[2]]'] <- 'Sequences'

tiff("merge_upset.TIF", width = 550, height = 750)
theme_set(theme_bw(base_size = 14))
myplot <- function(mydata, x, y) {
  plot <- (ggplot(data = mydata, aes_string(x = x, y = y, colour = "color")) + 
             geom_point() + scale_color_identity() + 
             theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))) +
    ylim(0,100000) +
    xlim(0,400) +
    theme_bw()
}

upset(pipeCore_merge, 
      order.by = c("degree", 'freq'), 
      point.size = 4, 
      line.size = 2,  
      empty.intersections = "on", 
      mainbar.y.label = "Number of Taxa", 
      sets.x.label = "Core Taxa",
      queries = list(list(query = intersects, params = list("No.Rarefaction","Seqs_7000"), active = T),
                     list(query = intersects, params = list("No.Rarefaction","Seqs_1000"), active = T),
                     list(query = intersects, params = list("Seqs_1000","Seqs_7000"), active = T),
                     list(query = intersects, params = list("No.Rarefaction"), active = T),
                     list(query = intersects, params = list("Seqs_1000"), active = T),
                     list(query = intersects, params = list("No.Rarefaction","Seqs_7000", "Seqs_1000"), active = T),
                     list(query = intersects, params = list("Seqs_7000"), active = T)),
      sets.bar.color = "#56B4E9",
      text.scale = 1.5,
      attribute.plots=list(gridrows = 60, ncols = 1, 
                           plots = list(list(plot = myplot, x = "Samples", y = "Sequences", queries = TRUE))))

dev.off()

