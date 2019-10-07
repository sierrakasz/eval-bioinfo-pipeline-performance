## Figure 5
# Multipanel figure showing core taxa among sample sizes and minimum library size

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
library(UpsetR)
library(vegan)


# import data ------------------------------------------------------------

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

#function for calculating core taxa among sample sizes
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

#Figure 5, panel a - core taxa among sample size
a <- ggplot(data = df, aes(x = OTU, y = logAbundance, color = factor(Subsample, 
                                                                labels = c("All", "60 & 120", "60 & 188", "120 & 188", "188")))) +
  geom_point(size = 4, alpha= 1/2) +
  geom_hline(yintercept = mean(df$logAbundance), color="black") +
  labs(color = 'Subsample') +
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#999999", "#D638E0", "#38E0A2"))


#calculating core among minimum library sizes
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

#Figure 5, panel B - core taxa among minimum library sizes
b <- ggplot(data = dfr, aes(x = OTU, y = logAbundance, color = factor(Rarefy, 
                                                                labels = c("All", "7000 & No Rarefaction", "No Rarefaction")))) +
  geom_point(size = 4, alpha= 1/2) +
  geom_hline(yintercept = mean(dfr$logAbundance), color="black") +
  labs(color = 'Rarefy') +
  scale_color_manual(values=c("#D8513C", "#5878D7", "#793C0E"))

#Final Figure
theme_set(theme_classic(base_size = 16))
tiff("Fig5.TIF", width = 4000, height = 2000, res=300)
ggarrange(a,b, 
          labels = c("A", "B"),
          ncol = 2)
dev.off()





