##Supplementary Material
# Each supplemental table is separated according to the order of tables in the manuscript

# Import Data  -------------------------------------------------------------

#load required packages
#phyloseq requires BiocManager 
# code for phyloseq install: if(!requireNamespace("BiocManager")){
#install.packages("BiocManager")
#}
#BiocManager::install("phyloseq")

#packages
library(httr)
library(HTSSIP)
library(gridExtra)
library(ggforce)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(PMCMR)
library(phyloseq)
library(plyr)
library(tidyverse)
library(vegan)
library(HTSSIP)
library(gridExtra)
library(ggforce)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(PMCMR)
library(phyloseq)
library(plyr)
library(tidyverse)
library(vegan)


###tax group is a file that combines the OTU file and taxonomy file. 
#no minimum library size
tax_group_nr <- read.csv(text=GET('https://github.com/sierrakasz/eval-bioinfo-pipeline-performance/blob/master/min_lib_sample_size/tax_group_norare.csv'))
tax_group_nr <- tax_group_nr %>% group_by(Kingdom, Phylum, Class, Order, Family, Genus) %>% summarise_all(funs(sum))
tax_group_nr <- tax_group_nr %>%  ungroup()
#Minimum library size: 7000
tax_group_7 <- read.csv(text=GET("https://github.com/sierrakasz/eval-bioinfo-pipeline-performance/blob/master/min_lib_sample_size/tax_group_7000.csv"))
tax_group_7 <- tax_group_7 %>% group_by(Kingdom, Phylum, Class, Order, Family, Genus) %>% summarise_all(funs(sum))
tax_group_7 <- tax_group_7 %>%  ungroup()
#Minimum library: 1000
tax_group_1 <- read.csv(text=GET("https://github.com/sierrakasz/eval-bioinfo-pipeline-performance/blob/master/min_lib_sample_size/tax_group_1000.csv"))
tax_group_1 <- tax_group_1 %>% group_by(Kingdom, Phylum, Class, Order, Family, Genus) %>% summarise_all(funs(sum))
tax_group_1 <- tax_group_1 %>%  ungroup()
#All of the combined together
#So, same samples among minimum library size, but with different minimum libraries
tax_group_rare <- read.csv(text=GET("https://github.com/sierrakasz/eval-bioinfo-pipeline-performance/blob/master/min_lib_sample_size/tax_group_rare.csv"))
tax_group_rare <- tax_group_rare %>% group_by(Kingdom, Phylum, Class, Order, Family, Genus) %>% summarise_all(funs(sum))
tax_group_rare <- tax_group_rare %>%  ungroup()

#this function specifically collapses the OTU table and taxonomy file 
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

#phyloseq objects moving forward
physeq_nr <- regroup_physeq_object(tax_group_nr)
physeq_7 <- regroup_physeq_object(tax_group_7)
physeq_1 <- regroup_physeq_object(tax_group_1)
physeq_rare <- regroup_physeq_object(tax_group_rare)

#import metadata
metadata=(text=GET("https://github.com/sierrakasz/eval-bioinfo-pipeline-performance/blob/master/min_lib_sample_size/HPMMMeta_r_merge_rare.csv"))
sampdat=sample_data(metadata)
sample_names(sampdat)=metadata$SampleID
#Change subsample level into a factor instead of numeric
sampdat$Subsample <- factor(sampdat$Subsample, levels = c('60', '120', '188'))
#combine into one phyloseq object
physeq_nr=merge_phyloseq(physeq_nr, sampdat)
physeq_7=merge_phyloseq(physeq_7, sampdat)
physeq_1=merge_phyloseq(physeq_1, sampdat)

#separate metadata and phyloseq object for having all minimum library sizes at largest subsample size (188)
metadata_rare=(text=GET("https://github.com/sierrakasz/eval-bioinfo-pipeline-performance/blob/master/min_lib_sample_size/HPMMMeta_r_merge_rare_levels.csv"))
sampdat_rare=sample_data(metadata_rare)
#change to factors with level
sampdat_rare$Rarefy <- factor(sampdat_rare$Rarefy, levels = c("No Rarefaction", "7000", "1000"))
sample_names(sampdat_rare)=metadata_rare$SampleID
merger = merge_phyloseq(physeq_rare, sampdat_rare)
sample_data(merger)$Rarefy <- factor(sample_data(merger)$Rarefy, levels = c("No Rarefaction", "7000", "1000"))
#final phyloseq object
physeq_188 <- subset_samples(merger, Subsample == '188')


# Table S17 -----------------------------------------------------------------
# S17 Table. Mean relative abundances of phyla greater than 1% relative abundance among sample 
# sizes, minimum library sizes, and sample areas

GPrPhylum=tax_glom(merger, "Phylum")
PhylumLevel = transform_sample_counts(GPrPhylum, function(x) x / sum(x))
PhylumLevel = filter_taxa(PhylumLevel, function(x) mean(x) > 0.0001, TRUE) 

df <- psmelt(PhylumLevel) 
df$Abundance=df$Abundance*100
Trtdata <- ddply(df, c('Pipeline', 'Sample Area', 'Rarefaction Level', 'Subsample'
), summarise,
N    = length(Abundance),
mean = mean(Abundance),
sd   = sd(Abundance),
se   = sd / sqrt(N)
)


# Table S18 ---------------------------------------------------------------
#S18 Table. Results from permutational multivariate analysis of variance (PERMANOVA) tests 
#on weighted Jaccard distance matrix including 

#pairwise subsample
pairwise_subsample <- function(physeq, physeq1, physeq7) {
  physeq_rec_nr <- subset_samples(physeq, Sample_Area == "Rectum")
  physeq_mou_nr <- subset_samples(physeq, Sample_Area == "Mouth")
  physeq_rec_7 <- subset_samples(physeq7, Sample_Area == "Rectum")
  physeq_mou_7 <- subset_samples(physeq7, Sample_Area == "Mouth")
  physeq_rec_1 <- subset_samples(physeq1, Sample_Area == "Rectum")
  physeq_mou_1 <- subset_samples(physeq1, Sample_Area == "Mouth")
  physeq_rec_nr_188 <- subset_samples(physeq_rec_nr, Subsample=='188')
  physeq_rec_nr_120 <- subset_samples(physeq_rec_nr, Subsample=='120')
  physeq_rec_nr_60 <- subset_samples(physeq_rec_nr, Subsample=='60')
  physeq_rec_7_188 <- subset_samples(physeq_rec_7, Subsample=='188')
  physeq_rec_7_120 <- subset_samples(physeq_rec_7, Subsample=='120')
  physeq_rec_7_60 <- subset_samples(physeq_rec_7, Subsample=='60')
  physeq_rec_1_188 <- subset_samples(physeq_rec_1, Subsample=='188')
  physeq_rec_1_120 <- subset_samples(physeq_rec_1, Subsample=='120')
  physeq_rec_1_60 <- subset_samples(physeq_rec_1, Subsample=='60')
  physeq_mou_nr_188 <- subset_samples(physeq_mou_nr, Subsample=='188')
  physeq_mou_nr_120 <- subset_samples(physeq_mou_nr, Subsample=='120')
  physeq_mou_nr_60 <- subset_samples(physeq_mou_nr, Subsample=='60')
  physeq_mou_7_188 <- subset_samples(physeq_mou_7, Subsample=='188')
  physeq_mou_7_120 <- subset_samples(physeq_mou_7, Subsample=='120')
  physeq_mou_7_60 <- subset_samples(physeq_mou_7, Subsample=='60')
  physeq_mou_1_188 <- subset_samples(physeq_mou_1, Subsample=='188')
  physeq_mou_1_120 <- subset_samples(physeq_mou_1, Subsample=='120')
  physeq_mou_1_60 <- subset_samples(physeq_mou_1, Subsample=='60')
  return(list(physeq_rec_nr_60120 <- merge_phyloseq(physeq_rec_nr_60, physeq_rec_nr_120),
              physeq_rec_nr_60188 <- merge_phyloseq(physeq_rec_nr_60, physeq_rec_nr_188),
              physeq_rec_nr_120188 <- merge_phyloseq(physeq_rec_nr_188, physeq_rec_nr_120),
              physeq_rec_7_60120 <- merge_phyloseq(physeq_rec_7_60, physeq_rec_7_120),
              physeq_rec_7_60188 <- merge_phyloseq(physeq_rec_7_60, physeq_rec_7_188),
              physeq_rec_7_120188 <- merge_phyloseq(physeq_rec_7_188, physeq_rec_7_120),
              physeq_rec_1_60120 <- merge_phyloseq(physeq_rec_1_60, physeq_rec_1_120),
              physeq_rec_1_60188 <- merge_phyloseq(physeq_rec_1_60, physeq_rec_1_188),
              physeq_rec_1_120188 <- merge_phyloseq(physeq_rec_1_188, physeq_rec_1_120),
              physeq_mou_nr_60120 <- merge_phyloseq(physeq_mou_nr_60, physeq_mou_nr_120),
              physeq_mou_nr_60188 <- merge_phyloseq(physeq_mou_nr_60, physeq_mou_nr_188),
              physeq_mou_nr_120188 <- merge_phyloseq(physeq_mou_nr_188, physeq_mou_nr_120),
              physeq_mou_7_60120 <- merge_phyloseq(physeq_mou_7_60, physeq_mou_7_120),
              physeq_mou_7_60188 <- merge_phyloseq(physeq_mou_7_60, physeq_mou_7_188),
              physeq_mou_7_120188 <- merge_phyloseq(physeq_mou_7_188, physeq_mou_7_120),
              physeq_mou_1_60120 <- merge_phyloseq(physeq_mou_1_60, physeq_mou_1_120),
              physeq_mou_1_60188 <- merge_phyloseq(physeq_mou_1_60, physeq_mou_1_188),
              physeq_mou_1_120188 <- merge_phyloseq(physeq_mou_1_188, physeq_mou_1_120)))
}

#testing rarefaction 
physeq_rec_188 <- subset_samples(merger, Sample_Area == "Rectum")
physeq_188 <- subset_samples(merger, Subsample=='188')
physeq_mou_188 <- subset_samples(merger, Sample_Area == "Mouth")
physeq_mou_188 <- subset_samples(physeq_mou_188, Subsample=='188')
physeq_rec_120 <- subset_samples(merger, Sample_Area == "Rectum")
physeq_rec_120 <- subset_samples(physeq_rec_120, Subsample=='120')
physeq_mou_120 <- subset_samples(merger, Sample_Area == "Mouth")
physeq_mou_120 <- subset_samples(physeq_mou_120, Subsample=='120')
physeq_rec_60 <- subset_samples(merger, Sample_Area == "Rectum")
physeq_rec_60 <- subset_samples(physeq_rec_60, Subsample=='60')
physeq_mou_60 <- subset_samples(merger, Sample_Area == "Mouth")
physeq_mou_60 <- subset_samples(physeq_mou_60, Subsample=='60')

#rarefy pairwise
rarefy_pairwise <- function(physeq) {
  physeq_rec <- subset_samples(merger, Sample_Area == "Rectum")
  physeq_mou <- subset_samples(merger, Sample_Area == "Mouth")
  physeq_rec_188 <- subset_samples(physeq_rec, Subsample=='188')
  physeq_mou_188 <- subset_samples(physeq_mou, Subsample=='188')
  physeq_rec_120 <- subset_samples(physeq_rec, Subsample=='120')
  physeq_mou_120 <- subset_samples(physeq_mou, Subsample=='120')
  physeq_rec_60 <- subset_samples(physeq_rec, Subsample=='60')
  physeq_mou_60 <- subset_samples(physeq_mou, Subsample=='60')
  physeq_rec_188_nr <- subset_samples(physeq_rec_188, Rarefy == 'No Rarefaction')
  physeq_rec_188_7 <- subset_samples(physeq_rec_188, Rarefy == '7000')
  physeq_rec_188_1 <- subset_samples(physeq_rec_188, Rarefy == '1000')
  physeq_rec_120_nr <- subset_samples(physeq_rec_120, Rarefy == 'No Rarefaction')
  physeq_rec_120_7 <- subset_samples(physeq_rec_120, Rarefy == '7000')
  physeq_rec_120_1 <- subset_samples(physeq_rec_120, Rarefy == '1000')
  physeq_rec_60_nr <- subset_samples(physeq_rec_60, Rarefy == 'No Rarefaction')
  physeq_rec_60_7 <- subset_samples(physeq_rec_60, Rarefy == '7000')
  physeq_rec_60_1 <- subset_samples(physeq_rec_60, Rarefy == '1000')
  physeq_mou_188_nr <- subset_samples(physeq_mou_188, Rarefy == 'No Rarefaction')
  physeq_mou_188_7 <- subset_samples(physeq_mou_188, Rarefy == '7000')
  physeq_mou_188_1 <- subset_samples(physeq_mou_188, Rarefy == '1000')
  physeq_mou_120_nr <- subset_samples(physeq_mou_120, Rarefy == 'No Rarefaction')
  physeq_mou_120_7 <- subset_samples(physeq_mou_120, Rarefy == '7000')
  physeq_mou_120_1 <- subset_samples(physeq_mou_120, Rarefy == '1000')
  physeq_mou_60_nr <- subset_samples(physeq_mou_60, Rarefy == 'No Rarefaction')
  physeq_mou_60_7 <- subset_samples(physeq_mou_60, Rarefy == '7000')
  physeq_mou_60_1 <- subset_samples(physeq_mou_60, Rarefy == '1000')
  return(list(physeq_rec_188_nr7 <- merge_phyloseq(physeq_rec_188_nr, physeq_rec_188_7),
              physeq_rec_188_nr1 <- merge_phyloseq(physeq_rec_188_nr, physeq_rec_188_1),
              physeq_rec_188_71 <- merge_phyloseq(physeq_rec_188_1, physeq_rec_188_7),
              physeq_rec_120_nr7 <- merge_phyloseq(physeq_rec_120_nr, physeq_rec_120_7),
              physeq_rec_120_nr1 <- merge_phyloseq(physeq_rec_120_nr, physeq_rec_120_1),
              physeq_rec_120_71 <- merge_phyloseq(physeq_rec_120_1, physeq_rec_120_7),
              physeq_rec_60_nr7 <- merge_phyloseq(physeq_rec_60_nr, physeq_rec_60_7),
              physeq_rec_60_nr1 <- merge_phyloseq(physeq_rec_60_nr, physeq_rec_60_1),
              physeq_rec_60_71 <- merge_phyloseq(physeq_rec_60_1, physeq_rec_60_7),
              physeq_mou_188_nr7 <- merge_phyloseq(physeq_mou_188_nr, physeq_mou_188_7),
              physeq_mou_188_nr1 <- merge_phyloseq(physeq_mou_188_nr, physeq_mou_188_1),
              physeq_mou_188_71 <- merge_phyloseq(physeq_mou_188_1, physeq_mou_188_7),
              physeq_mou_120_nr7 <- merge_phyloseq(physeq_mou_120_nr, physeq_mou_120_7),
              physeq_mou_120_nr1 <- merge_phyloseq(physeq_mou_120_nr, physeq_mou_120_1),
              physeq_mou_120_71 <- merge_phyloseq(physeq_mou_120_1, physeq_mou_120_7),
              physeq_mou_60_nr7 <- merge_phyloseq(physeq_mou_60_nr, physeq_mou_60_7),
              physeq_mou_60_nr1 <- merge_phyloseq(physeq_mou_60_nr, physeq_mou_60_1),
              physeq_mou_60_71 <- merge_phyloseq(physeq_mou_60_1, physeq_mou_60_7)))
}

beta_diversity_calc_rarefy <- function(physeq) {
  GPdist=phyloseq::distance(physeq, "(A+B-2*J)/(A+B-J)")
  sampledf <- data.frame(sample_data(physeq))
  print(return(adonis(GPdist ~ Rarefy, data = sampledf)))
}

beta_dispersion_calc_rarefy <- function(physeq) {
  GPdist=phyloseq::distance(physeq, "(A+B-2*J)/(A+B-J)")
  sampledf <- data.frame(sample_data(physeq))
  beta <- betadisper(GPdist, sampledf$Rarefy)
  print(return(permutest(beta)))
}

rarefy_list <- list(physeq_rec_188, physeq_rec_120, physeq_rec_60, 
                    physeq_mou_188, physeq_mou_120, physeq_mou_60 )

for(i in 1:length(rarefy_list)) {
  print(beta_diversity_calc_rarefy(rarefy_list[[i]]))
  print(beta_dispersion_calc_rarefy(rarefy_list[[i]]))
}

pairwise_rare_list <- rarefy_pairwise(merger)

for(i in 1:length(pairwise_rare_list)) {
  print(beta_diversity_calc_rarefy(pairwise_rare_list[[i]]))
  print(beta_dispersion_calc_rarefy(pairwise_rare_list[[i]]))
}

beta_diversity_calc_subsample <- function(physeq) {
  GPdist=phyloseq::distance(physeq, "(A+B-2*J)/(A+B-J)")
  sampledf <- data.frame(sample_data(physeq))
  print(return(adonis(GPdist ~ Subsample, data = sampledf)))
}

beta_dispersion_calc_subsample <- function(physeq) {
  GPdist=phyloseq::distance(physeq, "(A+B-2*J)/(A+B-J)")
  sampledf <- data.frame(sample_data(physeq))
  beta <- betadisper(GPdist, sampledf$Subsample)
  print(return(permutest(beta)))
}

subsample_list <- list(physeq_rec_nr, physeq_rec_7, physeq_rec_1,
                       physeq_mou_nr, physeq_mou_7, physeq_mou_1)

for(i in 1:length(subsample_list)) {
  print(beta_diversity_calc_subsample(subsample_list[[i]]))
  print(beta_dispersion_calc_subsample(subsample_list[[i]]))
}


# Table S19 ---------------------------------------------------------------
# S19 Table. Alpha diversity metrics for each sample among sample size and minimum library size 
# including: observed richness, Chao1, Shannon diversity (1-D), and Inverse Simpson diversity (1/D) 

erichnr <- estimate_richness(merger, measures = c("Observed", 'Chao1', "Shannon", "InvSimpson"))
erichnr <- add_rownames(erichnr, "SampleID")

#make data tidy

erichnr <- erichnr %>%
  gather(Index, Observation, c("Observed", 'Chao1', "Shannon", "InvSimpson"), na.rm = TRUE)

#add metadata
rich = merge(erichnr, metadata_rare)
rich <- rich %>% select(SampleID, Index, Observation, SampleName, Rarefy,
                        Subsample, Sample_Area, Pipeline)
rich$Index <- factor(rich$Index, levels = c("Observed", 'Chao1', "Shannon", "InvSimpson"))
rich$Subsample <- factor(rich$Subsample, levels = c('60', '120', '188'))
rich$Rarefy <- factor(rich$Rarefy, levels = c("No Rarefaction", "7000", "1000"))
write.csv(rich, file = 'alpha_chao1.csv')


# Table S20 ---------------------------------------------------------------
# S20 Table. Results of A) Kruskal-Wallis and B) pairwise post-hoc Nemenyi tests of alpha-diversity 
# metrics among sample size and minimum library size. 

#Subsamples
prepare_samples_kw <- function(rich) {
  return(list(rich_qim_mou_obs_nr <- rich %>% filter(Index == "Observed", Sample_Area == 'Mouth', Rarefy == 'No Rarefaction', Pipeline == 'QIIME2'),
              rich_qim_mou_sha_nr <- rich %>% filter(Index == "Shannon", Sample_Area == 'Mouth', Rarefy == 'No Rarefaction', Pipeline == 'QIIME2'),
              rich_qim_mou_inv_nr <- rich %>% filter(Index == "InvSimpson", Sample_Area == 'Mouth', Rarefy == 'No Rarefaction', Pipeline == 'QIIME2'),
              rich_qim_mou_cha_nr <- rich %>% filter(Index == 'Chao1', Sample_Area == 'Mouth', Rarefy == 'No Rarefaction', Pipeline == 'QIIME2'),
              rich_qim_rec_obs_nr <- rich %>% filter(Index == "Observed", Sample_Area == 'Rectum', Rarefy == 'No Rarefaction', Pipeline == 'QIIME2'),
              rich_qim_rec_sha_nr <- rich %>% filter(Index == "Shannon", Sample_Area == 'Rectum', Rarefy == 'No Rarefaction', Pipeline == 'QIIME2'),
              rich_qim_rec_inv_nr <- rich %>% filter(Index == "InvSimpson", Sample_Area == 'Rectum', Rarefy == 'No Rarefaction', Pipeline == 'QIIME2'),
              rich_qim_rec_cha_nr <- rich %>% filter(Index == 'Chao1', Sample_Area == 'Rectum', Rarefy == 'No Rarefaction', Pipeline == 'QIIME2'),
              rich_qim_mou_obs_7 <- rich %>% filter(Index == "Observed", Sample_Area == 'Mouth', Rarefy == '7000', Pipeline == 'QIIME2'),
              rich_qim_mou_sha_7 <- rich %>% filter(Index == "Shannon", Sample_Area == 'Mouth', Rarefy == '7000', Pipeline == 'QIIME2'),
              rich_qim_mou_inv_7 <- rich %>% filter(Index == "InvSimpson", Sample_Area == 'Mouth', Rarefy == '7000', Pipeline == 'QIIME2'),
              rich_qim_mou_cha_7 <- rich %>% filter(Index == 'Chao1', Sample_Area == 'Mouth', Rarefy == '7000', Pipeline == 'QIIME2'),
              rich_qim_rec_obs_7 <- rich %>% filter(Index == "Observed", Sample_Area == 'Rectum', Rarefy == '7000', Pipeline == 'QIIME2'),
              rich_qim_rec_sha_7 <- rich %>% filter(Index == "Shannon", Sample_Area == 'Rectum', Rarefy == '7000', Pipeline == 'QIIME2'),
              rich_qim_rec_inv_7 <- rich %>% filter(Index == "InvSimpson", Sample_Area == 'Rectum', Rarefy == '7000', Pipeline == 'QIIME2'),
              rich_qim_rec_cha_7 <- rich %>% filter(Index == 'Chao1', Sample_Area == 'Rectum', Rarefy == '7000', Pipeline == 'QIIME2'),
              rich_qim_mou_obs_1 <- rich %>% filter(Index == "Observed", Sample_Area == 'Mouth', Rarefy == '1000', Pipeline == 'QIIME2'),
              rich_qim_mou_sha_1 <- rich %>% filter(Index == "Shannon", Sample_Area == 'Mouth', Rarefy == '1000', Pipeline == 'QIIME2'),
              rich_qim_mou_inv_1 <- rich %>% filter(Index == "InvSimpson", Sample_Area == 'Mouth', Rarefy == '1000', Pipeline == 'QIIME2'),
              rich_qim_mou_cha_1 <- rich %>% filter(Index == 'Chao1', Sample_Area == 'Mouth', Rarefy == '1000', Pipeline == 'QIIME2'),
              rich_qim_rec_obs_1 <- rich %>% filter(Index == "Observed", Sample_Area == 'Rectum', Rarefy == '1000', Pipeline == 'QIIME2'),
              rich_qim_rec_sha_1 <- rich %>% filter(Index == "Shannon", Sample_Area == 'Rectum', Rarefy == '1000', Pipeline == 'QIIME2'),
              rich_qim_rec_inv_1 <- rich %>% filter(Index == "InvSimpson", Sample_Area == 'Rectum', Rarefy == '1000', Pipeline == 'QIIME2'),
              rich_qim_rec_cha_1 <- rich %>% filter(Index == 'Chao1', Sample_Area == 'Rectum', Rarefy == '1000', Pipeline == 'QIIME2')
  ))
  
}

kw_values <- prepare_samples_kw(rich)

for(i in 1:length(kw_values)) {
  print(kruskal.test(Observation ~ Subsample, data = kw_values[[i]]))
  out <- posthoc.kruskal.nemenyi.test(x=kw_values[[i]]$Observation, g=kw_values[[i]]$Subsample, dist='Tukey', p.adjust.method = 'bonf' )
  print(out$p.value)
}

#Rarefy
kw_values_rare <- list(rich_obs_rec_60_qim <- rich%>% filter(Index == "Observed", Sample_Area == 'Rectum', Subsample=='60', Pipeline== 'QIIME2'),
                       rich_obs_mou_60_qim <- rich%>% filter(Index == "Observed", Sample_Area == 'Mouth', Subsample=='60', Pipeline== 'QIIME2'),
                       rich_sha_rec_60_qim <- rich%>% filter(Index == "Shannon", Sample_Area == 'Rectum', Subsample=='60', Pipeline== 'QIIME2'),
                       rich_sha_mou_60_qim <- rich%>% filter(Index == "Shannon", Sample_Area == 'Mouth', Subsample=='60', Pipeline== 'QIIME2'),
                       rich_inv_rec_60_qim <- rich%>% filter(Index == "InvSimpson", Sample_Area == 'Rectum', Subsample=='60', Pipeline== 'QIIME2'),
                       rich_inv_mou_60_qim <- rich%>% filter(Index == "InvSimpson", Sample_Area == 'Mouth', Subsample=='60', Pipeline== 'QIIME2'),
                       rich_cha_rec_60_qim <- rich%>% filter(Index == "Chao1", Sample_Area == 'Rectum', Subsample=='60', Pipeline== 'QIIME2'),
                       rich_cha_mou_60_qim <- rich%>% filter(Index == "Chao1", Sample_Area == 'Mouth', Subsample=='60', Pipeline== 'QIIME2'),
                       rich_obs_rec_120_qim <- rich%>% filter(Index == "Observed", Sample_Area == 'Rectum', Subsample=='120', Pipeline== 'QIIME2'),
                       rich_obs_mou_120_qim <- rich%>% filter(Index == "Observed", Sample_Area == 'Mouth', Subsample=='120', Pipeline== 'QIIME2'),
                       rich_sha_rec_120_qim <- rich%>% filter(Index == "Shannon", Sample_Area == 'Rectum', Subsample=='120', Pipeline== 'QIIME2'),
                       rich_sha_mou_120_qim <- rich%>% filter(Index == "Shannon", Sample_Area == 'Mouth', Subsample=='120', Pipeline== 'QIIME2'),
                       rich_inv_rec_120_qim <- rich%>% filter(Index == "InvSimpson", Sample_Area == 'Rectum', Subsample=='120', Pipeline== 'QIIME2'),
                       rich_inv_mou_120_qim <- rich%>% filter(Index == "InvSimpson", Sample_Area == 'Mouth', Subsample=='120', Pipeline== 'QIIME2'),
                       rich_cha_rec_120_qim <- rich%>% filter(Index == "Chao1", Sample_Area == 'Rectum', Subsample=='120', Pipeline== 'QIIME2'),
                       rich_cha_mou_120_qim <- rich%>% filter(Index == "Chao1", Sample_Area == 'Mouth', Subsample=='120', Pipeline== 'QIIME2'),
                       rich_obs_rec_188_qim <- rich%>% filter(Index == "Observed", Sample_Area == 'Rectum', Subsample=='188', Pipeline== 'QIIME2'),
                       rich_obs_mou_188_qim <- rich%>% filter(Index == "Observed", Sample_Area == 'Mouth', Subsample=='188', Pipeline== 'QIIME2'),
                       rich_sha_rec_188_qim <- rich%>% filter(Index == "Shannon", Sample_Area == 'Rectum', Subsample=='188', Pipeline== 'QIIME2'),
                       rich_sha_mou_188_qim <- rich%>% filter(Index == "Shannon", Sample_Area == 'Mouth', Subsample=='188', Pipeline== 'QIIME2'),
                       rich_inv_rec_188_qim <- rich%>% filter(Index == "InvSimpson", Sample_Area == 'Rectum', Subsample=='188', Pipeline== 'QIIME2'),
                       rich_inv_mou_188_qim <- rich%>% filter(Index == "InvSimpson", Sample_Area == 'Mouth', Subsample=='188', Pipeline== 'QIIME2'),
                       rich_cha_rec_188_qim <- rich%>% filter(Index == "Chao1", Sample_Area == 'Rectum', Subsample=='188', Pipeline== 'QIIME2'),
                       rich_cha_mou_188_qim <- rich%>% filter(Index == "Chao1", Sample_Area == 'Mouth', Subsample=='188', Pipeline== 'QIIME2'))

for(i in 1:length(kw_values_rare)) {
  print(kruskal.test(Observation ~ Rarefy, data = kw_values_rare[[i]]))
  out <- posthoc.kruskal.nemenyi.test(x=kw_values_rare[[i]]$Observation, g=kw_values_rare[[i]]$Rarefy, dist='Tukey', p.adjust.method = 'bonf' )
  print(out$p.value)
}




# Table S21 ---------------------------------------------------------------

#S21 Table. Core microbiome among minimum library size. 

#rarefaction
metadata_rare=(read.csv("Metadata/HPMMMeta_r_merge_rare_levels.csv",header=TRUE))
sampdat_rare=sample_data(metadata_rare)
sampdat_rare$Rarefy <- factor(sampdat_rare$Rarefy, levels = c("No Rarefaction", "7000", "1000"))
sample_names(sampdat_rare)=metadata_rare$SampleID
merger = merge_phyloseq(physeq_rare, sampdat_rare)
sample_data(merger)$Rarefy <- factor(sample_data(merger)$Rarefy, levels = c("No Rarefaction", "7000", "1000"))
merger <- subset_samples(merger, Subsample = '188')

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
dfr <- data.frame(sapply(merger_upset, function(x) as.numeric(as.character(x))))
dfr <- dfr[,-1]

colnames(dfr) <- c("No", "S7000", "S1000", "total")
dfr$logAbundance <- log(dfr$total)


dfr$Rarefy[dfr$S1000 > 0 & dfr$S7000 > 0 & dfr$No > 0] <- 'All'
dfr$Rarefy[dfr$S1000 == 0 & dfr$S7000 > 0 & dfr$No > 0] <- '7000_No'
dfr$Rarefy[dfr$S1000 == 0 & dfr$S7000 == 0 & dfr$No > 0] <- 'No'
dfr <- dfr[!is.na(dfr$Rarefy),]
dfr$OTU <- c(1:length(dfr$total))
dfr$Rarefy <- factor(dfr$Rarefy, levels = c("All", "7000_No", "No"))

length(grep('All', dfr$Rarefy))
length(grep('7000', dfr$Rarefy))
length(grep('No', dfr$Rarefy))






# Table S22 ---------------------------------------------------------------
#S22 Table. Confusion matrices of random forest classification of 1,000 decision trees. 
# Classification out-of-box error rate of sample area and manner of death were determined among
# minimum library size, sample size, and taxonomic levels.  

physeq_nr_phylum <- tax_glom(physeq_nr, taxrank = 'Phylum')
physeq_nr_family <- tax_glom(physeq_nr, taxrank = 'Family')
physeq_7_phylum <- tax_glom(physeq_7, taxrank = 'Phylum')
physeq_7_family <- tax_glom(physeq_7, taxrank = 'Family')
physeq_1_phylum <- tax_glom(physeq_1, taxrank = 'Phylum')
physeq_1_family <- tax_glom(physeq_1, taxrank = 'Family')

samples_rel_abundance <- function(physeq) {
  return(list(
    physeq_60 <- subset_samples(physeq, Subsample == '60'),
    physeq_120 <- subset_samples(physeq, Subsample == '120'),
    physeq_188 <- subset_samples(physeq, Subsample == '188')))
}

physeq_nr_phy_list <- samples_rel_abundance(physeq_nr_phylum)
physeq_7_phy_list <- samples_rel_abundance(physeq_7_phylum)
physeq_1_phy_list <- samples_rel_abundance(physeq_1_phylum)
physeq_nr_fam_list <- samples_rel_abundance(physeq_nr_family)
physeq_7_fam_list <- samples_rel_abundance(physeq_7_family)
physeq_1_fam_list <- samples_rel_abundance(physeq_1_family)

#testing within subsample and rarefy 
#which one is better, for predicting sample area and MoD?

random_foresting_sample_area <- function(data_phy) {
  fd = data_phy
  predictors=t(otu_table(fd))
  dim(predictors)
  resp <- as.factor(sample_data(fd)$Sample_Area)
  rf.data <- data.frame(resp, predictors)
  Forest <- randomForest(resp~., data=rf.data, ntree=1000)
  return(Forest)
}

for1 <- random_foresting_sample_area(physeq_nr_phy_list[[1]])
for2 <- random_foresting_sample_area(physeq_nr_phy_list[[2]])
for3 <- random_foresting_sample_area(physeq_nr_phy_list[[3]])
for4 <- random_foresting_sample_area(physeq_nr_fam_list[[1]])
for5 <- random_foresting_sample_area(physeq_nr_fam_list[[2]])
for6 <- random_foresting_sample_area(physeq_nr_fam_list[[3]])
for7 <- random_foresting_sample_area(physeq_7_phy_list[[1]])
for8 <- random_foresting_sample_area(physeq_7_phy_list[[2]])
for9 <- random_foresting_sample_area(physeq_7_phy_list[[3]])
for10 <- random_foresting_sample_area(physeq_7_fam_list[[1]])
for11 <- random_foresting_sample_area(physeq_7_fam_list[[2]])
for12 <- random_foresting_sample_area(physeq_7_fam_list[[3]])
for13 <- random_foresting_sample_area(physeq_1_phy_list[[1]])
for14 <- random_foresting_sample_area(physeq_1_phy_list[[2]])
for15 <- random_foresting_sample_area(physeq_1_phy_list[[3]])
for16 <- random_foresting_sample_area(physeq_1_fam_list[[1]])
for17 <- random_foresting_sample_area(physeq_1_fam_list[[2]])
for18 <- random_foresting_sample_area(physeq_1_fam_list[[3]])


random_foresting_MoD <- function(data_phy) {
  fd = data_phy
  predictors=t(otu_table(fd))
  dim(predictors)
  resp <- as.factor(sample_data(fd)$MoD)
  rf.data <- data.frame(resp, predictors)
  Forest <- randomForest(resp~., data=rf.data, ntree=1000)
  return(Forest)
}


for19 <- random_foresting_MoD(physeq_nr_phy_list[[1]])
for20 <- random_foresting_MoD(physeq_nr_phy_list[[2]])
for21 <- random_foresting_MoD(physeq_nr_phy_list[[3]])
for22 <- random_foresting_MoD(physeq_nr_fam_list[[1]])
for23 <- random_foresting_MoD(physeq_nr_fam_list[[2]])
for24 <- random_foresting_MoD(physeq_nr_fam_list[[3]])
for25 <- random_foresting_MoD(physeq_7_phy_list[[1]])
for26 <- random_foresting_MoD(physeq_7_phy_list[[2]])
for27 <- random_foresting_MoD(physeq_7_phy_list[[3]])
for28 <- random_foresting_MoD(physeq_7_fam_list[[1]])
for29 <- random_foresting_MoD(physeq_7_fam_list[[2]])
for30 <- random_foresting_MoD(physeq_7_fam_list[[3]])
for31 <- random_foresting_MoD(physeq_1_phy_list[[1]])
for32 <- random_foresting_MoD(physeq_1_phy_list[[2]])
for33 <- random_foresting_MoD(physeq_1_phy_list[[3]])
for34 <- random_foresting_MoD(physeq_1_fam_list[[1]])
for35 <- random_foresting_MoD(physeq_1_fam_list[[2]])
for36 <- random_foresting_MoD(physeq_1_fam_list[[3]])

# Table S23 ----------------------------------------------------------------
# S23 Table. Top twenty predictor taxa based on mean decrease in Gini for random forest classifying manner of death and sample area.  

forest_predictors_nr <- function(forest) {
  imp <- importance(forest)
  imp <- data.frame(predictors = rownames(imp), imp)
  imp.sort <- arrange(imp, desc(MeanDecreaseGini))
  imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)
  imp.20 <- imp.sort[1:20, ]
  tax <- data.frame(tax_table(physeq_nr))
  head(tax)
  tax <- tax %>% select(Kingdom, Phylum, Class, Order, Family, Genus)
  tax$predictors <- rownames(tax)
  imp.20 <- merge(imp.20, tax)
  write.table(imp.20, sep=",", "random_forest_predictors_qiime_nr.csv", append = TRUE)
  return(imp.20)
}

for_list_nr <- list(for1, for2, for3, for4, for5, for6, 
                    for19, for20, for21, for22, for23, for24)

print(for_list_nr)

for(i in 1:length(for_list_nr)) {
  forest_predictors_nr(for_list_nr[[i]])
}

forest_predictors_7 <- function(forest) {
  imp <- importance(forest)
  imp <- data.frame(predictors = rownames(imp), imp)
  imp.sort <- arrange(imp, desc(MeanDecreaseGini))
  imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)
  imp.20 <- imp.sort[1:20, ]
  tax <- data.frame(tax_table(physeq_7))
  head(tax)
  tax <- tax %>% select(Kingdom, Phylum, Class, Order, Family, Genus)
  tax$predictors <- rownames(tax)
  imp.20 <- merge(imp.20, tax)
  write.table(imp.20, sep=",", "random_forest_predictors_qiime_7.csv", append = TRUE)
  return(imp.20)
}

for_list_7 <- list(for7, for8, for9, for10, for11, for12, 
                   for25, for26, for27, for28, for29, for30)

print(for_list_7)

for(i in 1:length(for_list_7)) {
  forest_predictors_7(for_list_7[[i]])
}

forest_predictors_1 <- function(forest) {
  imp <- importance(forest)
  imp <- data.frame(predictors = rownames(imp), imp)
  imp.sort <- arrange(imp, desc(MeanDecreaseGini))
  imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)
  imp.20 <- imp.sort[1:20, ]
  tax <- data.frame(tax_table(physeq_1))
  head(tax)
  tax <- tax %>% select(Kingdom, Phylum, Class, Order, Family, Genus)
  tax$predictors <- rownames(tax)
  imp.20 <- merge(imp.20, tax)
  write.table(imp.20, sep=",", "random_forest_predictors_qiime_1.csv", append = TRUE)
  return(imp.20)
}

for_list_1 <- list(for13, for14, for15, for16, for17, for18, 
                   for31, for32, for33, for34, for35, for36)

print(for_list_1)

for(i in 1:length(for_list_1)) {
  forest_predictors_1(for_list_1[[i]])
}




# Table S24 --------------------------------------------------------------
# S24 Table. Number of samples that remained after filtering (including chimera detection and 
# rarefaction) among sample area, sample size and minimum library size.

sample_num_function <- function(physeq) {
  physeq_rec <- subset_samples(physeq, Sample_Area == "Rectum")
  physeq_mou <- subset_samples(physeq, Sample_Area == "Mouth")
  physeq_rec_60 <- subset_samples(physeq_rec, Subsample == '60')
  physeq_rec_120 <- subset_samples(physeq_rec, Subsample == '120')
  physeq_rec_188 <- subset_samples(physeq_rec, Subsample == '188')
  physeq_mou_60 <- subset_samples(physeq_mou, Subsample == '60')
  physeq_mou_120 <- subset_samples(physeq_mou, Subsample == '120')
  physeq_mou_188 <- subset_samples(physeq_mou, Subsample == '188')
  QR60 <- length(sample_names(physeq_rec_60))
  QR120 <- length(sample_names(physeq_rec_120))
  QR188 <- length(sample_names(physeq_rec_188))
  QM60 <- length(sample_names(physeq_mou_60))
  QM120 <- length(sample_names(physeq_mou_120))
  QM188 <- length(sample_names(physeq_mou_188))
  QT60 <- QR60 + QM60
  QT120 <- QR120 + QM120
  QT188 <- QR188 + QM188
  sample_num_vec <- c(QT60,QT120,QT188,QM60,QM120,QM188,QR60,QR120,QR188)
  return(sample_num_vec)
}


# Table S25----------------------------------------------------------------
# S25 Table. Core microbiome among sample size. 

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
df <- data.frame(sapply(physeq_nr_upset, function(x) as.numeric(as.character(x))))
df <- df[,-1]

colnames(df) <- c("S60", "S120", "S188", "total")
df$logAbundance <- log(df$total)


df$Subsample[df$S60 > 0 & df$S120 > 0 & df$S188 > 0] <- 'All'
df$Subsample[df$S60 > 0 & df$S120 > 0 & df$S188 == 0] <- '60_120'
df$Subsample[df$S60 > 0 & df$S120 == 0 & df$S188 > 0] <- '60_188'
df$Subsample[df$S60 == 0 & df$S120 > 0 & df$S188 > 0] <- '120_188'
df$Subsample[df$S60 == 0 & df$S120 == 0 & df$S188 > 0] <- '188'
df <- df[!is.na(df$Subsample),]
df$OTU <- c(1:length(df$total))
df$Subsample <- factor(df$Subsample, levels = c("All", "60_120", "60_188", "120_188", "188"))

length(grep('All', df$Subsample))
length(grep('60_188', df$Subsample))
length(grep('60_120', df$Subsample))
length(grep('120_188', df$Subsample))
length(grep('188', df$Subsample))

