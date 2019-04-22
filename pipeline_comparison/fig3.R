#packages
library(gridExtra)
library(ggforce)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(phyloseq)
library(plyr)
library(tidyverse)

tax_group <- read.csv("tax_group_insilico.csv")
tax_group <- tax_group %>% group_by(Kingdom, Phylum, Class, Order, Family, Genus) %>% summarise_all(funs(sum))
tax_group <- tax_group %>%  ungroup()

regroup_physeq_object <-function(table) {
  tax <- table %>% select(Kingdom, Phylum, Class, Order, Family, Genus)
  tax <- as.matrix(tax)
  otu <- table %>% select(contains("HMP"))
  OTU=otu_table(otu, taxa_are_rows=TRUE)
  TAX=tax_table(tax)
  taxa_names(TAX)=row.names(OTU)
  physeq_all=phyloseq(OTU,TAX)
  return(physeq_all)
}

physeq_all <- regroup_physeq_object(tax_group)

#import metadata and combine
metadata=read.csv("insilico/mock3/mock3_meta_merge.csv",header=TRUE)
sampdat=sample_data(metadata)
sample_names(sampdat)=metadata$SampleID
physeq=merge_phyloseq(physeq_all, sampdat)
physeq=rarefy_even_depth(physeq)

#rarefied level
otu <- as.data.frame(otu_table(physeq))
seq_num <- sum(otu[,3])

df <- read.csv("rel_abundance_insilico_asis.csv", header=T)
new_df <- df[,c("Sample", "Abundance", "Pipeline", "Kingdom", "Phylum", "Class", 
                "Order", "Family")]
new_df <- new_df %>% group_by(Sample, Pipeline, Kingdom, Phylum, Class, Order, Family) %>% summarise_all(funs(sum))

#excel sheet from mockrobiota website
#https://github.com/caporaso-lab/mockrobiota
taxa_answers <- read.csv('insilico/taxa_answers_nogenus.csv', header=T)
taxa_answers <- taxa_answers %>%  mutate(Expected = Expected * seq_num)

mg <- new_df %>% filter(Pipeline == 'MG-RAST')
mot <- new_df %>% filter(Pipeline == 'mothur')
qim <- new_df %>% filter(Pipeline == 'QIIME2')

compare_mg <- merge(mg, taxa_answers, all= T)
compare_mot <- merge(mot, taxa_answers, all= T)
compare_qim <- merge(qim, taxa_answers, all= T)
compare_mg[is.na(compare_mg)] <- 0
compare_mot[is.na(compare_mot)] <- 0
compare_qim[is.na(compare_qim)] <- 0

compare_mg_even1 <- compare_mg %>% filter(Sample == 'HMPMockV1.1.Even1')
compare_mg_even2 <- compare_mg %>% filter(Sample == 'HMPMockV1.1.Even2')
compare_mg_stag1 <- compare_mg %>% filter(Sample == 'HMPMockV1.2.Staggered1')
compare_mg_stag2 <- compare_mg %>% filter(Sample == 'HMPMockV1.2.Staggered2')

compare_mot_even1 <- compare_mot %>% filter(Sample == 'HMPMockV1.1.Even1')
compare_mot_even2 <- compare_mot %>% filter(Sample == 'HMPMockV1.1.Even2')
compare_mot_stag1 <- compare_mot %>% filter(Sample == 'HMPMockV1.2.Staggered1')
compare_mot_stag2 <- compare_mot %>% filter(Sample == 'HMPMockV1.2.Staggered2')

compare_qim_even1 <- compare_qim %>% filter(Sample == 'HMPMockV1.1.Even1')
compare_qim_even2 <- compare_qim %>% filter(Sample == 'HMPMockV1.1.Even2')
compare_qim_stag1 <- compare_qim %>% filter(Sample == 'HMPMockV1.2.Staggered1')
compare_qim_stag2 <- compare_qim %>% filter(Sample == 'HMPMockV1.2.Staggered2')

find_taxa_errors <- function(dataframe) {
  correct_taxa <- c()
  false_positive <- c()
  false_negative <- c()
  for(i in 1:nrow(dataframe)) {
    if (dataframe$Abundance[[i]] > 0 & dataframe$Expected[[i]] > 0) {
      correct_taxa<- append(correct_taxa, 1)
    } else if (dataframe$Abundance[[i]] == 0 & dataframe$Expected[[i]] > 0) {
      false_negative <-append(false_negative, 1)
    } else if (dataframe$Abundance[[i]] > 0 & dataframe$Expected[[i]]== 0) {
      false_positive <- append(false_positive, 1)
    }
  }
  a <- c(sum(correct_taxa), sum(false_negative), sum(false_positive))
  return(a)
}

mg_even1 <- find_taxa_errors(compare_mg_even1)
mg_even2 <- find_taxa_errors(compare_mg_even2)
mg_stag1 <- find_taxa_errors(compare_mg_stag1)
mg_stag2 <- find_taxa_errors(compare_mg_stag2)

mot_even1 <- find_taxa_errors(compare_mot_even1)
mot_even2 <- find_taxa_errors(compare_mot_even2)
mot_stag1 <- find_taxa_errors(compare_mot_stag1)
mot_stag2 <- find_taxa_errors(compare_mot_stag2)

qim_even1 <- find_taxa_errors(compare_qim_even1)
qim_even2 <- find_taxa_errors(compare_qim_even2)
qim_stag1 <- find_taxa_errors(compare_qim_stag1)
qim_stag2 <- find_taxa_errors(compare_qim_stag2)

counts <- c(mg_even1, mg_even2, mg_stag1, mg_stag2, mot_even1, mot_even2, mot_stag1,
                   mot_stag2, qim_even1, qim_even2, qim_stag1, qim_stag2)
names <- c(rep("mg_even1", 3), rep("mg_even2",3), rep("mg_stag1",3), rep("mg_stag2",3), 
           rep("mot_even1", 3), rep("mot_even2",3), rep("mot_stag1",3), rep("mot_stag2",3), 
           rep("qim_even1", 3), rep("qim_even2",3), rep("qim_stag1",3), rep("qim_stag2",3))
taxa_Data <- data.frame(cbind(names, counts))
taxa_Data$Pipeline <- 'QIIME2'
taxa_Data$Pipeline[grep('mot', taxa_Data$names)] = 'mothur'
taxa_Data$Pipeline[grep('mg', taxa_Data$names)] = 'MG-RAST'
taxa_Data$Community <- 'Even'
taxa_Data$Community[grep('stag', taxa_Data$names)] = 'Staggered'
taxa_Data$Sample_number <- factor(1:36)
taxa_Data$Id <- rep(c("Correct Taxa", "False Negatives", "False Positives"), 12)
taxa_Data$counts <- as.numeric(as.character(taxa_Data$counts))

a <- ggplot(data = taxa_Data, aes(x = Sample_number, y = counts, fill = Pipeline)) + 
  geom_bar(stat = 'identity') + facet_wrap(~Id, scales = 'free') + xlab("") +
  scale_fill_manual(values=c("#787878", "#ffb31a", "#5c5c8a")) + ylab("Counts")
a <- a  + theme(axis.text.x = element_blank()) +
  theme(legend.position="top")
a

changing_false_negatives <- function(dataframe) {
  for(i in 1:nrow(dataframe)) {
    if(dataframe$Abundance[[i]] > 0 & dataframe$Expected[[i]] == 0) {
      dataframe$Abundance[[i]] <- -dataframe$Abundance[[i]]
    }
  }
  return(dataframe)
}

compare_mg <- changing_false_negatives(compare_mg)

ggplotRegression <- function(fit){
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point(color = 'grey', size = 4) + 
    geom_abline(slope =1, intercept = 0, linetype = 3) +
    geom_abline(slope =0, intercept = 0, color='black') +
    geom_vline(xintercept = 0, color='black') +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste('MG-RAST')) +
    annotate('text', label = "Adjusted R-squared = 0.238",
                                              x = 1500, y = 200)
}

b <- ggplotRegression(lm(Expected ~ Abundance, data = compare_mg))
b

compare_mot <- changing_false_negatives(compare_mot)

ggplotRegression <- function(fit){
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point(color = "#ffb31a", size = 4) + 
    geom_abline(slope =1, intercept = 0, linetype = 3) +
    geom_abline(slope =0, intercept = 0, color='black') +
    geom_vline(xintercept = 0, color='black') +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste('mothur')) +
    annotate('text', label = "Adjusted R-squared = 0.083",
             x = -1000, y = 800)
    
}

c <- ggplotRegression(lm(Expected ~ Abundance, data = compare_mot))
c


compare_qim <- changing_false_negatives(compare_qim)

ggplotRegression <- function(fit){
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point(color = "#5c5c8a", size = 4) + 
    geom_abline(slope =1, intercept = 0, linetype=3) +
    geom_abline(slope =0, intercept = 0, color='black') +
    geom_vline(xintercept = 0, color='black') +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste('QIIME2')) +
    annotate('text', label = "Adjusted R-squared = 0.267",
             x = 1200, y = 200)
}

d <- ggplotRegression(lm(Expected ~ Abundance, data = compare_qim))
d

theme_set(theme_classic(base_size = 18))
tiff("Fig3.TIF", width = 1000, height = 800)
ggarrange(a,b,c,d , 
          labels = c("A", "B", "C", "D"),
          nrow = 2, ncol = 2)
dev.off()
