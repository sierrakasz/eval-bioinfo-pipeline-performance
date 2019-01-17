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

samples_rel_abundance <- function(physeq) {
  return(list(
  physeq_60 <- subset_samples(physeq, Subsample == '60'),
  physeq_120 <- subset_samples(physeq, Subsample == '120'),
  physeq_188 <- subset_samples(physeq, Subsample == '188')))
}

physeq_nr_list <- samples_rel_abundance(physeq_nr)
physeq_7_list <- samples_rel_abundance(physeq_7)
physeq_1_list <- samples_rel_abundance(physeq_1)

relative_abundance_setup_rec <- function(physeq) {
  physeq_rec <- subset_samples(physeq, Sample_Area == "Rectum")
  GPrPhylum=tax_glom(physeq_rec, "Phylum")
  PhylumLevel_rec = transform_sample_counts(GPrPhylum, function(x) x / sum(x))
  PhylumLevel_rec = filter_taxa(PhylumLevel_rec, function(x) mean(x) > 0.0001, TRUE) 
  df <- psmelt(PhylumLevel_rec) 
  df$Abundance=df$Abundance*100
  Trtdata <- ddply(df, c("Phylum", "Pipeline", 'Subsample'), summarise,
                   N    = length(Abundance),
                   mean = mean(Abundance),
                   sd   = sd(Abundance),
                   se   = sd / sqrt(N)
  )
}

df_rec_nr <- list()
for(i in 1:length(physeq_nr_list)) {
  df_rec_nr[[i]] <- relative_abundance_setup_rec(physeq_nr_list[[i]])
}

sink("rel_abundance_phy_rec_norare.csv")
print(df_rec_nr)
sink()

df_rec_graph_nr <- relative_abundance_setup_rec(physeq_nr)
df_rec_graph_nr$Subsample <- factor(df_rec_graph_nr$Subsample, levels= c('60', '120', '188'))
theme_set(theme_classic(base_size = 14))
a=ggplot(df_rec_graph_nr, aes(x=Subsample, y=mean, fill=Subsample))+
  geom_bar(stat="identity")+
  ylab("Relative Abundance Phylum Level (> 1%)")+facet_wrap(~Phylum, scales="free")+
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="black",width=1.0)
a
tiff("rel_abund_rec__phy_nr_188.TIF", width = 800, height = 600)
a + scale_fill_manual(values = c("#006699", "#FF6600", "#FFCC00")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  xlab("")
dev.off()


df_rec_7 <- list()
for(i in 1:length(physeq_7_list)) {
  df_rec_7[[i]] <- relative_abundance_setup_rec(physeq_7_list[[i]])
}

sink("rel_abundance_phy_rec_7000.csv")
print(df_rec_7)
sink()

df_rec_graph_7 <- relative_abundance_setup_rec(physeq_7)
df_rec_graph_7$Subsample <- factor(df_rec_graph_7$Subsample, levels= c('60', '120', '188'))
theme_set(theme_classic(base_size = 14))
a=ggplot(df_rec_graph_7, aes(x=Subsample, y=mean, fill=Subsample))+
  geom_bar(stat="identity")+
  ylab("Relative Abundance Phylum Level (> 1%)")+facet_wrap(~Phylum, scales="free")+
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="black",width=1.0)
a
tiff("rel_abund_rec__phy_7_188.TIF", width = 800, height = 600)
a + scale_fill_manual(values = c("#006699", "#FF6600", "#FFCC00")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  xlab("")
dev.off()

df_rec_1 <- list()
for(i in 1:length(physeq_1_list)) {
  df_rec_1[[i]] <- relative_abundance_setup_rec(physeq_1_list[[i]])
}

sink("rel_abundance_phy_rec_1000.csv")
print(df_rec_1)
sink()

df_rec_graph_1 <- relative_abundance_setup_rec(physeq_1)
df_rec_graph_1$Subsample <- factor(df_rec_graph_1$Subsample, levels= c('60', '120', '188'))
theme_set(theme_classic(base_size = 14))
a=ggplot(df_rec_graph_1, aes(x=Subsample, y=mean, fill=Subsample))+
  geom_bar(stat="identity")+
  ylab("Relative Abundance Phylum Level (> 1%)")+facet_wrap(~Phylum, scales="free")+
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="black",width=1.0)
a
tiff("rel_abund_rec__phy_1_188.TIF", width = 800, height = 600)
a + scale_fill_manual(values = c("#006699", "#FF6600", "#FFCC00")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  xlab("")
dev.off()

relative_abundance_setup_mou <- function(physeq) {
  physeq_mou <- subset_samples(physeq, Sample_Area == "Mouth")
  GPrPhylum_a=tax_glom(physeq_mou, "Phylum")
  PhylumLevel_mou = transform_sample_counts(GPrPhylum_a, function(x) x / sum(x))
  PhylumLevel_mou = filter_taxa(PhylumLevel_mou, function(x) mean(x) > 0.0001, TRUE) 
  df <- psmelt(PhylumLevel_mou) 
  df$Abundance=df$Abundance*100
  Trtdata <- ddply(df, c("Phylum", "Pipeline", 'Subsample'), summarise,
                   N    = length(Abundance),
                   mean = mean(Abundance),
                   sd   = sd(Abundance),
                   se   = sd / sqrt(N)
  )
}

df_mou_nr <- list()
for(i in 1:length(physeq_nr_list)) {
  df_mou_nr[[i]] <- relative_abundance_setup_mou(physeq_nr_list[[i]])
}

sink("rel_abundance_phy_mou_norare.csv")
print(df_mou_nr)
sink()

df_mou_graph_nr <- relative_abundance_setup_mou(physeq_nr)
df_mou_graph_nr$Subsample <- factor(df_mou_graph_nr$Subsample, levels= c('60', '120', '188'))
theme_set(theme_classic(base_size = 14))
a=ggplot(df_mou_graph_nr, aes(x=Subsample, y=mean, fill=Subsample))+
  geom_bar(stat="identity")+
  ylab("Relative Abundance Phylum Level (> 1%)")+facet_wrap(~Phylum, scales="free")+
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="black",width=1.0)
a
tiff("rel_abund_mou__phy_nr_188.TIF", width = 800, height = 600)
a + scale_fill_manual(values = c("#006699", "#FF6600", "#FFCC00")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  xlab("")
dev.off()

df_mou_7 <- list()
for(i in 1:length(physeq_7_list)) {
  df_mou_7[[i]] <- relative_abundance_setup_mou(physeq_7_list[[i]])
}

sink("rel_abundance_phy_mou_7000.csv")
print(df_mou_7)
sink()

df_mou_graph_7 <- relative_abundance_setup_mou(physeq_7)
df_mou_graph_7$Subsample <- factor(df_mou_graph_7$Subsample, levels= c('60', '120', '188'))
theme_set(theme_classic(base_size = 14))
a=ggplot(df_mou_graph_7, aes(x=Subsample, y=mean, fill=Subsample))+
  geom_bar(stat="identity")+
  ylab("Relative Abundance Phylum Level (> 1%)")+facet_wrap(~Phylum, scales="free")+
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="black",width=1.0)
a
tiff("rel_abund_mou__phy_7_188.TIF", width = 800, height = 600)
a + scale_fill_manual(values = c("#006699", "#FF6600", "#FFCC00")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  xlab("")
dev.off()

df_mou_1 <- list()
for(i in 1:length(physeq_1_list)) {
  df_mou_1[[i]] <- relative_abundance_setup_mou(physeq_1_list[[i]])
}

sink("rel_abundance_phy_mou_1000.csv")
print(df_mou_1)
sink()

df_mou_graph_1 <- relative_abundance_setup_mou(physeq_1)
df_mou_graph_1$Subsample <- factor(df_mou_graph_1$Subsample, levels= c('60', '120', '188'))
theme_set(theme_classic(base_size = 14))
a=ggplot(df_mou_graph_1, aes(x=Subsample, y=mean, fill=Subsample))+
  geom_bar(stat="identity")+
  ylab("Relative Abundance Phylum Level (> 1%)")+facet_wrap(~Phylum, scales="free")+
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="black",width=1.0)
a
tiff("rel_abund_mou__phy_1_188.TIF", width = 800, height = 600)
a + scale_fill_manual(values = c("#006699", "#FF6600", "#FFCC00")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  xlab("")
dev.off()

relative_abundance_setup_rec_fam <- function(physeq) {
  physeq_rec <- subset_samples(physeq, Sample_Area == "Rectum")
  GPrFamily=tax_glom(physeq_rec, "Family")
  FamilyLevel_rec = transform_sample_counts(GPrFamily, function(x) x / sum(x))
  FamilyLevel_rec = filter_taxa(FamilyLevel_rec, function(x) mean(x) > 0.0001, TRUE) 
  df <- psmelt(FamilyLevel_rec) 
  df$Abundance=df$Abundance*100
  Trtdata <- ddply(df, c("Family", "Pipeline", 'Subsample'), summarise,
                   N    = length(Abundance),
                   mean = mean(Abundance),
                   sd   = sd(Abundance),
                   se   = sd / sqrt(N)
  )
}

df_rec_nr_fam <- list()
for(i in 1:length(physeq_nr_list)) {
  df_rec_nr_fam[[i]] <- relative_abundance_setup_rec_fam(physeq_nr_list[[i]])
}

sink("rel_abundance_fam_rec_nr.csv")
print(df_rec_nr_fam)
sink()

df_rec_7_fam <- list()
for(i in 1:length(physeq_7_list)) {
  df_rec_7_fam[[i]] <- relative_abundance_setup_rec_fam(physeq_7_list[[i]])
}

sink("rel_abundance_fam_rec_7000.csv")
print(df_rec_7_fam)
sink()

df_rec_1_fam <- list()
for(i in 1:length(physeq_1_list)) {
  df_rec_1_fam[[i]] <- relative_abundance_setup_rec_fam(physeq_1_list[[i]])
}

sink("rel_abundance_fam_rec_1000.csv")
print(df_rec_1_fam)
sink()

relative_abundance_setup_mou_fam <- function(physeq) {
  physeq_mou <- subset_samples(physeq, Sample_Area == "Mouth")
  GPrFamily=tax_glom(physeq_mou, "Family")
  FamilyLevel_mou = transform_sample_counts(GPrFamily, function(x) x / sum(x))
  FamilyLevel_mou = filter_taxa(FamilyLevel_mou, function(x) mean(x) > 0.0001, TRUE)  
  df <- psmelt(FamilyLevel_mou) 
  df$Abundance=df$Abundance*100
  Trtdata <- ddply(df, c("Family", "Pipeline", "Subsample"), summarise,
                   N    = length(Abundance),
                   mean = mean(Abundance),
                   sd   = sd(Abundance),
                   se   = sd / sqrt(N)
  )
}

df_mou_nr_fam <- list()
for(i in 1:length(physeq_nr_list)) {
  df_mou_nr_fam[[i]] <- relative_abundance_setup_mou_fam(physeq_nr_list[[i]])
}

sink("rel_abundance_fam_mou_nr.csv")
print(df_mou_nr_fam)
sink()

df_mou_7_fam <- list()
for(i in 1:length(physeq_7_list)) {
  df_mou_7_fam[[i]] <- relative_abundance_setup_mou_fam(physeq_7_list[[i]])
}

sink("rel_abundance_fam_mou_7000.csv")
print(df_mou_7_fam)
sink()

df_mou_1_fam <- list()
for(i in 1:length(physeq_1_list)) {
  df_mou_1_fam[[i]] <- relative_abundance_setup_mou_fam(physeq_1_list[[i]])
}

sink("rel_abundance_fam_mou_1000.csv")
print(df_mou_1_fam)
sink()