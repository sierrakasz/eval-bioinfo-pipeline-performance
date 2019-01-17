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

#Alpha Diversity
erich <- estimate_richness(physeq, measures = c("Observed", "Shannon", "InvSimpson"))
write.table(erich,"alpdiv_three_pipelines.txt",sep="\t",row.names=TRUE)
erich <- add_rownames(erich, "SampleID")

#make data tidy 
erich2 <- erich %>%
  gather(Index, Observation, c("Observed", "Shannon", "InvSimpson"), na.rm = TRUE)

#add metadata
rich = merge(erich2, metadata)
rich <- rich %>% select(SampleID, Index, Observation, SampleName, Pipeline,
                        Sample_Area)

theme_set(theme_bw(base_size = 40))
tiff("alphadiv_three_pipelines.TIF", width = 1400, height = 1200)
rich$Index <- factor(rich$Index, levels = c("Observed", "Shannon", "InvSimpson"))
p <- ggplot(rich, aes(x=Pipeline, y=Observation, fill=Pipeline)) +
  geom_boxplot(outlier.size = 3)
p + facet_grid(Index~Sample_Area, scales="free") + scale_fill_manual(values = c("#787878", "#ffb31a", "#5c5c8a")) +
  scale_x_discrete(labels=c('MG-RAST' = 'MG', 'mothur' = 'M', "QIIME2" = 'Q')) +
  xlab("") + theme_bw(base_size = 40) + theme(panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

#kruskal wallis test
#nemenyi  test

#kruskal.test(alpha_div~pipeline, data = alpha)
#out <- posthoc.kruskal.nemenyi.test(x=alpha_div, g=pipeline, dist="Tukey", p.adjust.method = 'bonf')
#print(otu$statistic)
#https://cran.r-project.org/web/packages/PMCMR/vignettes/PMCMR.pdf

rich_obs_rec <- rich %>% filter(Index == "Observed", Sample_Area == 'Rectum')
rich_obs_mou <- rich %>% filter(Index == "Observed", Sample_Area == 'Mouth')
rich_sha_rec <- rich %>% filter(Index == "Shannon", Sample_Area == 'Rectum')
rich_sha_mou <- rich %>% filter(Index == "Shannon", Sample_Area == 'Mouth')
rich_inv_rec <- rich %>% filter(Index == "InvSimpson", Sample_Area == 'Rectum')
rich_inv_mou <- rich %>% filter(Index == "InvSimpson", Sample_Area == 'Mouth')

kruskal.test(Observation ~ Pipeline, data = rich_obs_rec)
out <- posthoc.kruskal.nemenyi.test(x=rich_obs_rec$Observation, g=rich_obs_rec$Pipeline, dist='Tukey', p.adjust.method = 'bonf' )
print(out$p.value)

kruskal.test(Observation ~ Pipeline, data = rich_obs_mou)
out <- posthoc.kruskal.nemenyi.test(x=rich_obs_mou$Observation, g=rich_obs_mou$Pipeline, dist='Tukey', p.adjust.method = 'bonf' )
print(out$p.value)

kruskal.test(Observation ~ Pipeline, data = rich_sha_rec)
out <- posthoc.kruskal.nemenyi.test(x=rich_sha_rec$Observation, g=rich_sha_rec$Pipeline, dist='Tukey', p.adjust.method = 'bonf' )
print(out$p.value)

kruskal.test(Observation ~ Pipeline, data = rich_sha_mou)
out <- posthoc.kruskal.nemenyi.test(x=rich_sha_mou$Observation, g=rich_sha_mou$Pipeline, dist='Tukey', p.adjust.method = 'bonf' )
print(out$p.value)

kruskal.test(Observation ~ Pipeline, data = rich_inv_rec)
out <- posthoc.kruskal.nemenyi.test(x=rich_inv_rec$Observation, g=rich_inv_rec$Pipeline, dist='Tukey', p.adjust.method = 'bonf' )
print(out$p.value)

kruskal.test(Observation ~ Pipeline, data = rich_inv_mou)
out <- posthoc.kruskal.nemenyi.test(x=rich_inv_mou$Observation, g=rich_inv_mou$Pipeline, dist='Tukey', p.adjust.method = 'bonf' )
print(out$p.value)


