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

#import metadata and combine
metadata=(read.csv("Metadata/HPMMMeta_r_merge_rare.csv",header=TRUE))
sampdat=sample_data(metadata)
sample_names(sampdat)=metadata$SampleID
sampdat$Subsample <- factor(sampdat$Subsample, levels = c('60', '120', '188'))

physeq_nr=merge_phyloseq(physeq_nr, sampdat)

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


a <- ggplot(data = df, aes(x = OTU, y = logAbundance, color = factor(Subsample, 
                                                                labels = c("All", "60 & 120", "60 & 188", "120 & 188", "188")))) +
  geom_point(size = 4) +
  geom_hline(yintercept = mean(df$logAbundance), color="black") +
  labs(color = 'Subsample') +
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#999999", "#D638E0", "#38E0A2"))

tax_group_rare <- read.csv("tax_group_rare.csv")
tax_group_rare <- tax_group_rare %>% group_by(Kingdom, Phylum, Class, Order, Family, Genus) %>% summarise_all(funs(sum))
tax_group_rare <- tax_group_rare %>%  ungroup()

physeq_rare <- regroup_physeq_object(tax_group_rare)

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


b <- ggplot(data = dfr, aes(x = OTU, y = logAbundance, color = factor(Rarefy, 
                                                                labels = c("All", "7000 & No Rarefaction", "No Rarefaction")))) +
  geom_point(size = 4) +
  geom_hline(yintercept = mean(dfr$logAbundance), color="black") +
  labs(color = 'Rarefy') +
  scale_color_manual(values=c("#D8513C", "#5878D7", "#793C0E"))

#Final Figure
theme_set(theme_classic(base_size = 16))
tiff("Fig5.TIF", width = 1200, height = 600)
ggarrange(a,b, 
          labels = c("A", "B"),
          ncol = 2)
dev.off()





