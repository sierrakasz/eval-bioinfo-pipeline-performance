## FIGURE 2.
# Figure two is a multi-panel comparison of error rate, false positives, and false negatives among pipelines. Using 'mockrobiota' sequences, the true taxa and 
# abundanes are known. Raw sequence files were run through each pipeline following parameters used for the postmortem microbiome dataset and assessed


#load required packages
#phyloseq requires BiocManager 
# code for phyloseq install: if(!requireNamespace("BiocManager")){
#install.packages("BiocManager")
#}
#BiocManager::install("phyloseq")

library(httr)
library(gridExtra)
library(ggforce)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(phyloseq)
library(plyr)
library(tidyverse)

#pulling in a csv file from GitHub
#tax_group is a tax/otu table of 'mockrobiota' mock-3 communities run through QIIME2, mothur, and MG-RAST
#https://github.com/caporaso-lab/mockrobiota/tree/master/data/mock-3
tax_group <- read.csv(text=GET('https://github.com/sierrakasz/eval-bioinfo-pipeline-performance/blob/master/pipeline_comparison/tax_group_insilico.csv'), header=T)

#collapsing down same taxa names to get true abundance sum
tax_group <- tax_group %>% group_by(Kingdom, Phylum, Class, Order, Family, Genus) %>% summarise_all(funs(sum))
tax_group <- tax_group %>%  ungroup()

#formatting tax_group table into a phyloseq object
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

#import metadata and combine into final phyloseq object
metadata=read.csv(text=GET('https://github.com/sierrakasz/eval-bioinfo-pipeline-performance/blob/master/pipeline_comparison/mock3_meta_merge.csv'),header=TRUE)
sampdat=sample_data(metadata)
sample_names(sampdat)=metadata$SampleID
physeq=merge_phyloseq(physeq_all, sampdat)

#rarefy among pipelines to get comparable results
physeq=rarefy_even_depth(physeq)

#check what minimum library size it rarefied to
otu <- as.data.frame(otu_table(physeq))
seq_num <- sum(otu[,3])

#file was exported to excel
#formatting was fixed
#add metadata, abundance, and taxonomy in one table
#load in file for downstream use
#file for comparing to the 'true' taxonomy and abundances from mockrobiota website
df <- read.csv(text=GET('https://github.com/sierrakasz/eval-bioinfo-pipeline-performance/blob/master/pipeline_comparison/rel_abundance_insilico_asis.csv'), header=T)
new_df <- df[,c("Sample", "Abundance", "Pipeline", "Kingdom", "Phylum", "Class", 
                "Order", "Family")]
new_df <- new_df %>% group_by(Sample, Pipeline, Kingdom, Phylum, Class, Order, Family) %>% summarise_all(funs(sum))

#excel sheet from mockrobiota website
#true abundances to compare to pipeline outputs
taxa_answers <- read.csv(text=GET('https://github.com/sierrakasz/eval-bioinfo-pipeline-performance/blob/master/pipeline_comparison/taxa_answers_nogenus.csv'), header=T)
# mulitpling by minimum library size to properly compare to pipelines
taxa_answers <- taxa_answers %>%  mutate(Expected = Expected * seq_num)

#separate out analyses for each pipeline based on metadata
mg <- new_df %>% filter(Pipeline == 'MG-RAST')
mot <- new_df %>% filter(Pipeline == 'mothur')
qim <- new_df %>% filter(Pipeline == 'QIIME2')

#merge together pipeline outputs and true abundances
compare_mg <- merge(mg, taxa_answers, all= T)
compare_mot <- merge(mot, taxa_answers, all= T)
compare_qim <- merge(qim, taxa_answers, all= T)
#remove NAs and make them zeros. False postives and negatives will not properly merge and show up as NAs
compare_mg[is.na(compare_mg)] <- 0
compare_mot[is.na(compare_mot)] <- 0
compare_qim[is.na(compare_qim)] <- 0

#separate out each sample within pipelines
#there are even and staggered samples 1 and 2. Even and staggered refers to the abundances of taxa
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

#this function will determine the number of false postives, false negatives, and correct taxa for each samples and pipeline
#false negatives: samples in mockrobiota but not in pipeline output
#false postive: samples in pipeline output but not in mockrobiota
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

#resulting vector from function for each sample and pipeline
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

#take separate vectors and but them into one dataframe
#add counts and sample names and important metadata (pipeline, staggered vs. even)
#in the order of correct, false negatives, false positives 
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

#Figure 2, part A
#simple bar chart showing the number of correct taxa, false postives, and false negatives among samples and pipelines
a <- ggplot(data = taxa_Data, aes(x = Sample_number, y = counts, fill = Pipeline)) + 
  geom_bar(stat = 'identity') + facet_wrap(~Id, scales = 'free') + xlab("") +
  scale_fill_manual(values=c("#787878", "#ffb31a", "#5c5c8a")) + ylab("Counts")
a <- a  + theme(axis.text.x = element_blank()) +
  theme(legend.position="top")
a

#making false negatives negative abundance values for an abundance comparison 
changing_false_negatives <- function(dataframe) {
  for(i in 1:nrow(dataframe)) {
    if(dataframe$Abundance[[i]] > 0 & dataframe$Expected[[i]] == 0) {
      dataframe$Abundance[[i]] <- -dataframe$Abundance[[i]]
    }
  }
  return(dataframe)
}

#starting with MG-RAST, changing false negatives to negative abundances
compare_mg <- changing_false_negatives(compare_mg)

#plot regression for MG-RAST expected and actual abundances based on mockrobiota abundances and pipeline outputs 
#R squared values ~ error rate
#Adjusted R Squared was found in regression list and printed on figure
ggplotRegression <- function(fit){
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point(color = 'grey', size = 4, alpha = 1/2) + 
    geom_abline(slope =0, intercept = 0, color='black') +
    geom_vline(xintercept = 0, color='black') +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste('MG-RAST')) +
    annotate('text', label = "Adjusted R-squared = 0.238",
                                              x = 1500, y = 200)
}

#Figure 2, part B
#Error of MG-RAST
b <- ggplotRegression(lm(Expected ~ Abundance, data = compare_mg))
b

#same as above, except with mothur
compare_mot <- changing_false_negatives(compare_mot)

ggplotRegression <- function(fit){
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point(color = "#ffb31a", size = 4, alpha = 1/2) + 
    geom_abline(slope =0, intercept = 0, color='black') +
    geom_vline(xintercept = 0, color='black') +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste('mothur')) +
    annotate('text', label = "Adjusted R-squared = 0.083",
             x = -1000, y = 800)
    
}

#Figure 2, part C
#Error of mothur
c <- ggplotRegression(lm(Expected ~ Abundance, data = compare_mot))
c

#same as above except for QIIME2
compare_qim <- changing_false_negatives(compare_qim)

ggplotRegression <- function(fit){
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point(color = "#5c5c8a", size = 4, alpha = 1/2) + 
    geom_abline(slope =0, intercept = 0, color='black') +
    geom_vline(xintercept = 0, color='black') +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste('QIIME2')) +
    annotate('text', label = "Adjusted R-squared = 0.267",
             x = 1200, y = 200)
}

#Figure 2, part c
#Error of QIIME2
d <- ggplotRegression(lm(Expected ~ Abundance, data = compare_qim))
d

#make the final figure
theme_set(theme_classic(base_size = 18))
tiff("FIG2.TIF", width = 4000, height = 3500, res=300)
ggarrange(a,b,c,d , 
          labels = c("A", "B", "C", "D"),
          nrow = 2, ncol = 2)
dev.off()
