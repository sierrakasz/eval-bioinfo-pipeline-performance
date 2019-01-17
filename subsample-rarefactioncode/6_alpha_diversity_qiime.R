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
library(vegan)

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

erich <- estimate_richness(physeq_nr, measures = c("Observed", "Shannon", "InvSimpson"))
write.table(erich,"alpdiv_qiime_nr.txt",sep="\t",row.names=TRUE)
erich <- add_rownames(erich, "SampleID")
erich$Rarefy <- 'No Rarefaction'

erich7 <- estimate_richness(physeq_7, measures = c("Observed", "Shannon", "InvSimpson"))
write.table(erich7,"alpdiv_qiime_7000.txt",sep="\t",row.names=TRUE)
erich7 <- add_rownames(erich7, "SampleID")
erich7$Rarefy <- '7000'

erich1 <- estimate_richness(physeq_1, measures = c("Observed", "Shannon", "InvSimpson"))
write.table(erich1,"alpdiv_qiime_1000.txt",sep="\t",row.names=TRUE)
erich1 <- add_rownames(erich1, "SampleID")
erich1$Rarefy <- '1000'

#make data tidy 
erich <- erich %>%
  gather(Index, Observation, c("Observed", "Shannon", "InvSimpson"), na.rm = TRUE)

erich7 <- erich7 %>%
  gather(Index, Observation, c("Observed", "Shannon", "InvSimpson"), na.rm = TRUE)

erich1 <- erich1 %>%
  gather(Index, Observation, c("Observed", "Shannon", "InvSimpson"), na.rm = TRUE)

#add metadata
rich = merge(erich, metadata)
rich <- rich %>% select(SampleID, Index, Observation, SampleName, Pipeline,
                        Sample_Area, Subsample, Rarefy)
rich$Index <- factor(rich$Index, levels = c("Observed", "Shannon", "InvSimpson"))
rich$Subsample <- factor(rich$Subsample, levels = c('60', '120', '188'))

theme_set(theme_bw(base_size = 14))
tiff("alphadiv_qiime_pipelines_nr.TIF", width = 800, height = 600)
p <- ggplot(rich, aes(x=Subsample, y=Observation, fill=Subsample)) +
  geom_boxplot()
p + facet_grid(Index~Sample_Area, scales="free") + scale_fill_manual(values = c("#006699", "#FF6600", "#FFCC00"))
dev.off()

rich7 = merge(erich7, metadata)
rich7 <- rich7 %>% select(SampleID, Index, Observation, SampleName, Pipeline,
                        Sample_Area, Subsample, Rarefy)
rich7$Index <- factor(rich7$Index, levels = c("Observed", "Shannon", "InvSimpson"))
rich7$Subsample <- factor(rich7$Subsample, levels = c('60', '120', '188'))

theme_set(theme_bw(base_size = 14))
tiff("alphadiv_qiime_7000_pipelines.TIF", width = 800, height = 600)
p <- ggplot(rich7, aes(x=Subsample, y=Observation, fill=Subsample)) +
  geom_boxplot()
p + facet_grid(Index~Sample_Area, scales="free") + scale_fill_manual(values = c("#006699", "#FF6600", "#FFCC00"))
dev.off()

rich1 = merge(erich1, metadata)
rich1 <- rich1 %>% select(SampleID, Index, Observation, SampleName, Pipeline,
                        Sample_Area, Subsample, Rarefy)
rich1$Index <- factor(rich1$Index, levels = c("Observed", "Shannon", "InvSimpson"))
rich1$Subsample <- factor(rich1$Subsample, levels = c('60', '120', '188'))

theme_set(theme_bw(base_size = 14))
tiff("alphadiv_qiime_1000_pipelines.TIF", width = 800, height = 600)
p <- ggplot(rich1, aes(x=Subsample, y=Observation, fill=Subsample)) +
  geom_boxplot()
p + facet_grid(Index~Sample_Area, scales="free") + scale_fill_manual(values = c("#006699", "#FF6600", "#FFCC00"))
dev.off()

rare_rich <- rbind(rich, rich7, rich1)
rare_rich$Rarefy <- factor(rare_rich$Rarefy, levels = c("No Rarefaction", '7000', '1000'))
rare_rich$Rarefy <- as.factor(rare_rich$Rarefy)

#Subsamples
prepare_samples_kw <- function(rich) {
  return(list(rich_qim_mou_obs <- rich %>% filter(Index == "Observed", Sample_Area == 'Mouth', Pipeline == 'QIIME2'),
  rich_qim_mou_obs <- rich %>% filter(Index == "Observed", Sample_Area == 'Mouth', Pipeline == 'QIIME2'),
  rich_qim_mou_sha <- rich %>% filter(Index == "Shannon", Sample_Area == 'Mouth', Pipeline == 'QIIME2'),
  rich_qim_mou_sha <- rich %>% filter(Index == "Shannon", Sample_Area == 'Mouth', Pipeline == 'QIIME2'),
  rich_qim_mou_inv <- rich %>% filter(Index == "InvSimpson", Sample_Area == 'Mouth', Pipeline == 'QIIME2'),
  rich_qim_mou_inv <- rich %>% filter(Index == "InvSimpson", Sample_Area == 'Mouth', Pipeline == 'QIIME2'),
  rich_qim_rec_obs <- rich %>% filter(Index == "Observed", Sample_Area == 'Rectum', Pipeline == 'QIIME2'),
  rich_qim_rec_obs <- rich %>% filter(Index == "Observed", Sample_Area == 'Rectum', Pipeline == 'QIIME2'),
  rich_qim_rec_sha <- rich %>% filter(Index == "Shannon", Sample_Area == 'Rectum', Pipeline == 'QIIME2'),
  rich_qim_rec_sha <- rich %>% filter(Index == "Shannon", Sample_Area == 'Rectum', Pipeline == 'QIIME2'),
  rich_qim_rec_inv <- rich %>% filter(Index == "InvSimpson", Sample_Area == 'Rectum', Pipeline == 'QIIME2'),
  rich_qim_rec_inv <- rich %>% filter(Index == "InvSimpson", Sample_Area == 'Rectum', Pipeline == 'QIIME2')))
}

kw_values <- prepare_samples_kw(rich)

for(i in 1:length(kw_values)) {
  print(kruskal.test(Observation ~ Subsample, data = kw_values[[i]]))
  out <- posthoc.kruskal.nemenyi.test(x=kw_values[[i]]$Observation, g=kw_values[[i]]$Subsample, dist='Tukey', p.adjust.method = 'bonf' )
  print(out$p.value)
}

kw_values <- prepare_samples_kw(rich7)

for(i in 1:length(kw_values)) {
  print(kruskal.test(Observation ~ Subsample, data = kw_values[[i]]))
  out <- posthoc.kruskal.nemenyi.test(x=kw_values[[i]]$Observation, g=kw_values[[i]]$Subsample, dist='Tukey', p.adjust.method = 'bonf' )
  print(out$p.value)
}

kw_values <- prepare_samples_kw(rich1)

for(i in 1:length(kw_values)) {
  print(kruskal.test(Observation ~ Subsample, data = kw_values[[i]]))
  out <- posthoc.kruskal.nemenyi.test(x=kw_values[[i]]$Observation, g=kw_values[[i]]$Subsample, dist='Tukey', p.adjust.method = 'bonf' )
  print(out$p.value)
}


kw_values_rare <- list(rich_obs_rec_60_qim <- rare_rich%>% filter(Index == "Observed", Sample_Area == 'Rectum', Subsample=='60', Pipeline== 'QIIME2'),
  rich_obs_mou_60_qim <- rare_rich%>% filter(Index == "Observed", Sample_Area == 'Mouth', Subsample=='60', Pipeline== 'QIIME2'),
  rich_sha_rec_60_qim <- rare_rich%>% filter(Index == "Shannon", Sample_Area == 'Rectum', Subsample=='60', Pipeline== 'QIIME2'),
  rich_sha_mou_60_qim <- rare_rich%>% filter(Index == "Shannon", Sample_Area == 'Mouth', Subsample=='60', Pipeline== 'QIIME2'),
  rich_inv_rec_60_qim <- rare_rich%>% filter(Index == "InvSimpson", Sample_Area == 'Rectum', Subsample=='60', Pipeline== 'QIIME2'),
  rich_inv_mou_60_qim <- rare_rich%>% filter(Index == "InvSimpson", Sample_Area == 'Mouth', Subsample=='60', Pipeline== 'QIIME2'),
  rich_obs_rec_120_qim <- rare_rich%>% filter(Index == "Observed", Sample_Area == 'Rectum', Subsample=='120', Pipeline== 'QIIME2'),
  rich_obs_mou_120_qim <- rare_rich%>% filter(Index == "Observed", Sample_Area == 'Mouth', Subsample=='120', Pipeline== 'QIIME2'),
  rich_sha_rec_120_qim <- rare_rich%>% filter(Index == "Shannon", Sample_Area == 'Rectum', Subsample=='120', Pipeline== 'QIIME2'),
  rich_sha_mou_120_qim <- rare_rich%>% filter(Index == "Shannon", Sample_Area == 'Mouth', Subsample=='120', Pipeline== 'QIIME2'),
  rich_inv_rec_120_qim <- rare_rich%>% filter(Index == "InvSimpson", Sample_Area == 'Rectum', Subsample=='120', Pipeline== 'QIIME2'),
  rich_inv_mou_120_qim <- rare_rich%>% filter(Index == "InvSimpson", Sample_Area == 'Mouth', Subsample=='120', Pipeline== 'QIIME2'),
  rich_obs_rec_188_qim <- rare_rich%>% filter(Index == "Observed", Sample_Area == 'Rectum', Subsample=='188', Pipeline== 'QIIME2'),
  rich_obs_mou_188_qim <- rare_rich%>% filter(Index == "Observed", Sample_Area == 'Mouth', Subsample=='188', Pipeline== 'QIIME2'),
  rich_sha_rec_188_qim <- rare_rich%>% filter(Index == "Shannon", Sample_Area == 'Rectum', Subsample=='188', Pipeline== 'QIIME2'),
  rich_sha_mou_188_qim <- rare_rich%>% filter(Index == "Shannon", Sample_Area == 'Mouth', Subsample=='188', Pipeline== 'QIIME2'),
  rich_inv_rec_188_qim <- rare_rich%>% filter(Index == "InvSimpson", Sample_Area == 'Rectum', Subsample=='188', Pipeline== 'QIIME2'),
  rich_inv_mou_188_qim <- rare_rich%>% filter(Index == "InvSimpson", Sample_Area == 'Mouth', Subsample=='188', Pipeline== 'QIIME2'))

for(i in 1:length(kw_values_rare)) {
  print(kruskal.test(Observation ~ Rarefy, data = kw_values_rare[[i]]))
  out <- posthoc.kruskal.nemenyi.test(x=kw_values_rare[[i]]$Observation, g=kw_values_rare[[i]]$Rarefy, dist='Tukey', p.adjust.method = 'bonf' )
  print(out$p.value)
}

