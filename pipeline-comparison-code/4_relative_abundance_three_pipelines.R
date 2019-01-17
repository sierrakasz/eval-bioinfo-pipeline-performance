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


#Tax glom
GPrPhylum=tax_glom(physeq_rec, "Phylum")
PhylumLevel_rec = transform_sample_counts(GPrPhylum, function(x) x / sum(x))
PhylumLevel_rec = filter_taxa(PhylumLevel_rec, function(x) mean(x) > 0.0001, TRUE) 

GPrPhylum=tax_glom(physeq_mou, "Phylum")
PhylumLevel_mou = transform_sample_counts(GPrPhylum, function(x) x / sum(x))
PhylumLevel_mou = filter_taxa(PhylumLevel_mou, function(x) mean(x) > 0.0001, TRUE) 

#Relative abundance
#relative abundance by variable
df <- psmelt(PhylumLevel_mou) 
df$Abundance=df$Abundance*100
Trtdata <- ddply(df, c("Phylum", "Pipeline"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)
write.csv(Trtdata, "rel_abundance_phy_mou.csv")

theme_set(theme_classic(base_size = 40))
a=ggplot(Trtdata, aes(x=Pipeline, y=mean, fill=Pipeline))+
  geom_bar(stat="identity")+
  ylab("Relative Abundance Phylum Level (> 1%)")+facet_wrap(~Phylum, scales="free")+
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="black",width=1.0)
a
tiff("rel_abund_mou_phy.TIF", width = 1800, height = 1200)
a + scale_fill_manual(values = c("#787878", "#ffb31a", "#5c5c8a")) + 
  xlab("Mouth") +
  scale_x_discrete(labels=c('MG-RAST' = 'MG', 'mothur' = 'M', "QIIME2" = 'Q'))
dev.off()

df <- psmelt(PhylumLevel_rec) 
df$Abundance=df$Abundance*100
Trtdata <- ddply(df, c("Phylum", "Pipeline"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)
write.csv(Trtdata, "rel_abundance_phy_rec.csv")

theme_set(theme_classic(base_size = 14))
a=ggplot(Trtdata, aes(x=Pipeline, y=mean, fill=Pipeline))+
  geom_bar(stat="identity")+
  ylab("Relative Abundance Phylum Level (> 1%)")+facet_wrap(~Phylum, scales="free")+
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="black",width=1.0)
a
tiff("rel_abund_rec_phy.TIF", width = 1850, height = 1200)
a + scale_fill_manual(values = c("#787878", "#ffb31a", "#5c5c8a")) + 
  xlab("Rectum") +
  scale_x_discrete(labels=c('MG-RAST' = 'MG', 'mothur' = 'M', "QIIME2" = 'Q'))
dev.off()

#relative abundance family level
#Tax glom
GPrFamily=tax_glom(physeq_rec, "Family")
FamilyLevel_rec = transform_sample_counts(GPrFamily, function(x) x / sum(x))
FamilyLevel_rec = filter_taxa(FamilyLevel_rec, function(x) mean(x) > 0.0001, TRUE) 

GPrFamily=tax_glom(physeq_mou, "Family")
FamilyLevel_mou = transform_sample_counts(GPrFamily, function(x) x / sum(x))
FamilyLevel_mou = filter_taxa(FamilyLevel_mou, function(x) mean(x) > 0.0001, TRUE) 

#Relative abundance
#realtive abundance by variable
df <- psmelt(FamilyLevel_rec) 
df$Abundance=df$Abundance*100
Trtdata <- ddply(df, c("Family", "Pipeline"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)
write.csv(Trtdata, "rel_abundance_rec_fam.csv")

df <- psmelt(FamilyLevel_mou) 
df$Abundance=df$Abundance*100
Trtdata <- ddply(df, c("Family", "Pipeline"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)
write.csv(Trtdata, "rel_abundance_mou_fam.csv")

