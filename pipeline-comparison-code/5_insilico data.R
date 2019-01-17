rm(list=ls())
getwd()
setwd("C:/Users/sierr/Documents/ThesisProject_Pipelines/")

#samples ran through the pipelines came from here
#https://github.com/caporaso-lab/mockrobiota

#packages
library(phyloseq)
library(plyr)
library(tidyverse)

#import mothur
import_mothur_files <- function(shared, tax) {
  physeq_mothur <- import_mothur(mothur_shared_file = shared,
                                 mothur_constaxonomy_file = tax)
  colnames(tax_table(physeq_mothur)) <- c("Kingdom", "Phylum", "Class", 
                                          "Order", "Family",  "Genus")
  return(physeq_mothur)
}

physeq_mothur <- import_mothur_files("insilico/mock3/mothur_mock3.shared", 
                                     "insilico/mock3/mothur_mock3.taxonomy")
sample_names(physeq_mothur) <- c('M_HMPMockV1.1.Even1','M_HMPMockV1.1.Even2',
                                 'M_HMPMockV1.2.Staggered1','M_HMPMockV1.2.Staggered2'
)

#import qiime
import_qiime_files <- function(biom) {
  biom_qiime <- import_biom(biom)
  a <- data.frame(tax_table(biom_qiime))
  b <- a %>% 
    mutate(Rank1 = str_replace(Rank1, "D_0__", "")) %>% 
    mutate(Rank2 = str_replace(Rank2, "D_1__", "")) %>% 
    mutate(Rank3 = str_replace(Rank3, "D_2__", "")) %>% 
    mutate(Rank4 = str_replace(Rank4, "D_3__", "")) %>% 
    mutate(Rank5 = str_replace(Rank5, "D_4__", "")) %>% 
    mutate(Rank6 = str_replace(Rank6, "D_5__", ""))
  c <- b[,c("Rank1", "Rank2", "Rank3", "Rank4", "Rank5", "Rank6")]
  c <- as.matrix(c)
  d <- data.frame(otu_table(biom_qiime))
  OTUq=otu_table(d, taxa_are_rows=TRUE)
  TAXq=tax_table(c)
  taxa_names(TAXq)=row.names(OTUq)
  physeq_qiime = merge_phyloseq(OTUq, TAXq)
  colnames(tax_table(physeq_qiime)) <- c("Kingdom", "Phylum", "Class",
                                         "Order","Family", "Genus")
  return(physeq_qiime)
}

physeq_qiime <- import_qiime_files("insilico/mock3/table-with-taxonomy.biom")
sample_names(physeq_qiime) <- paste("Q_", sample_names(physeq_qiime), sep="")

#import MG-RAST
import_MGRAST <- function(file) {
  file_name <- read.csv(file)
  otu <- select(file_name, contains('HMP'))
  tax <-select(file_name, domain, phylum, className, order, family, genus)
  tax <- as.matrix(tax)
  OTU=otu_table(otu, taxa_are_rows=TRUE)
  TAX=tax_table(tax)
  taxa_names(TAX)=row.names(OTU)
  physeq_MG=phyloseq(OTU,TAX)
  colnames(tax_table(physeq_MG)) <- c("Kingdom", "Phylum", "Class", 
                                      "Order", "Family",  "Genus")
  return(physeq_MG)
}
physeq_MG <- import_MGRAST("insilico\\mock3\\Mock3_MG-RAST.csv")
sample_names(physeq_MG) <- paste("G_", sample_names(physeq_MG), sep="")

#merging phyloseq objects together
merge = merge_phyloseq(physeq_qiime,physeq_mothur,physeq_MG)

#tidying data to insure that taxa are consistent 
tidying_data <- function(physeq) {
  x <- data.frame(otu_table(physeq))
  y <- data.frame(tax_table(physeq))
  x$names <- rownames(x)
  y$names <- rownames(y)
  z <- merge(x,y)
  z <- z %>% select(-contains("names"))
  tax_group <- z %>% group_by(Kingdom, Phylum, Class, Order, Family, Genus) %>% summarise_all(funs(sum))
  tax_group <- tax_group %>%  ungroup()
  return(tax_group)
}
tax_group <- tidying_data(merge)
write.csv(tax_group, "tax_group_insilico.csv")
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
physeq = rarefy_even_depth(physeq)

#rarefyed level
otu <- as.data.frame(otu_table(physeq))
seq_num <- sum(otu[,3])

df <- psmelt(physeq)
write.csv(df, "rel_abundance_insilico_asis.csv")
#remove zeros in files in excel
#fix names in excel 'sample' column - take out prefix
#use this excel formula: =RIGHT(C2,(LEN(C2)-2))
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

taxa_Data <- rbind(mg_even1, mg_even2, mg_stag1, mg_stag2, mot_even1, mot_even2, mot_stag1,
                   mot_stag2, qim_even1, qim_even2, qim_stag1, qim_stag2)
colnames(taxa_Data) <- c("Correct Taxa", "False Negatives", 'False Positives')

tiff('taxa_error_three_pipelines.TIF', width = 500, height = 400)
taxa_error_plot <- barplot(taxa_Data, beside = T, col = c("#787878", "#787878", "#787878", "#787878",
                                       "#ffb31a", "#ffb31a", "#ffb31a", "#ffb31a",
                                       "#5c5c8a", "#5c5c8a", "#5c5c8a", "#5c5c8a"), ylim = c(0,25),
                           ylab = ('Number of Taxa at Family Level'), density = c(500,500,10,10))
text(x = taxa_error_plot, y = taxa_Data, label = taxa_Data, pos = 3, cex = 0.75, col = "black")
legend(x='topright', legend = c("MG-RAST", "mothur", "QIIME2", 'Even' , "Staggered"), pch = c(15,15,15,15,12), 
       pt.cex = 3, cex=1, col = c("#787878",  "#ffb31a", "#5c5c8a", 'black', 'black'))
dev.off()

changing_false_negatives <- function(dataframe) {
  for(i in 1:nrow(dataframe)) {
  if(dataframe$Abundance[[i]] > 0 & dataframe$Expected[[i]] == 0) {
    dataframe$Abundance[[i]] <- -dataframe$Abundance[[i]]
   }
  }
  return(dataframe)
}

theme_set(theme_classic(base_size = 14))
compare_mg_even1_seq <- changing_false_negatives(compare_mg_even1)
compare_mg_even2_seq <- changing_false_negatives(compare_mg_even2)
compare_mg_stag1_seq <- changing_false_negatives(compare_mg_stag1)
compare_mg_stag2_seq <- changing_false_negatives(compare_mg_stag2)

tiff('sequence_error_MG_even1.TIF', width = 500, height = 400)
ggplot(data=compare_mg_even1_seq,aes(x=Family))+
  geom_bar(aes(y=Expected),stat="identity",position ="identity",alpha=.3,fill='white', color = 'black') +
  geom_bar(aes(y=Abundance),stat="identity",position ="identity",alpha=.8,fill="#787878",color="#787878") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab('Sequence Count')
dev.off()

tiff('sequence_error_MG_even2.TIF', width = 500, height = 400)
ggplot(data=compare_mg_even2_seq,aes(x=Family))+
  geom_bar(aes(y=Expected),stat="identity",position ="identity",alpha=.3,fill='white', color = 'black') +
  geom_bar(aes(y=Abundance),stat="identity",position ="identity",alpha=.8,fill="#787878",color="#787878") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab('Sequence Count')
dev.off()

tiff('sequence_error_MG_stag1.TIF', width = 500, height = 400)
ggplot(data=compare_mg_stag1_seq,aes(x=Family))+
  geom_bar(aes(y=Expected),stat="identity",position ="identity",alpha=.3,fill='white', color = 'black') +
  geom_bar(aes(y=Abundance),stat="identity",position ="identity",alpha=.8,fill="#787878",color="#787878") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab('Sequence Count')
dev.off()

tiff('sequence_error_MG_stag2.TIF', width = 500, height = 400)
ggplot(data=compare_mg_stag2_seq,aes(x=Family))+
  geom_bar(aes(y=Expected),stat="identity",position ="identity",alpha=.3,fill='white', color = 'black') +
  geom_bar(aes(y=Abundance),stat="identity",position ="identity",alpha=.8,fill="#787878",color="#787878") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab('Sequence Count')
dev.off()

compare_mot_even1_seq <- changing_false_negatives(compare_mot_even1)
compare_mot_even2_seq <- changing_false_negatives(compare_mot_even2)
compare_mot_stag1_seq <- changing_false_negatives(compare_mot_stag1)
compare_mot_stag2_seq <- changing_false_negatives(compare_mot_stag2)

tiff('sequence_error_mot_even1.TIF', width = 500, height = 400)
ggplot(data=compare_mot_even1_seq,aes(x=Family))+
  geom_bar(aes(y=Expected),stat="identity",position ="identity",alpha=.3,fill='white', color = 'black') +
  geom_bar(aes(y=Abundance),stat="identity",position ="identity",alpha=.8,fill="#ffb31a",color="#ffb31a") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab('Sequence Count')
dev.off()

tiff('sequence_error_mot_even2.TIF', width = 500, height = 400)
ggplot(data=compare_mot_even2_seq,aes(x=Family))+
  geom_bar(aes(y=Expected),stat="identity",position ="identity",alpha=.3,fill='white', color = 'black') +
  geom_bar(aes(y=Abundance),stat="identity",position ="identity",alpha=.8,fill="#ffb31a",color="#ffb31a") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab('Sequence Count')
dev.off()

tiff('sequence_error_mot_stag1.TIF', width = 500, height = 400)
ggplot(data=compare_mot_stag1_seq,aes(x=Family))+
  geom_bar(aes(y=Expected),stat="identity",position ="identity",alpha=.3,fill='white', color = 'black') +
  geom_bar(aes(y=Abundance),stat="identity",position ="identity",alpha=.8,fill="#ffb31a",color="#ffb31a") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab('Sequence Count')
dev.off()

tiff('sequence_error_mot_stag2.TIF', width = 500, height = 400)
ggplot(data=compare_mot_stag2_seq,aes(x=Family))+
  geom_bar(aes(y=Expected),stat="identity",position ="identity",alpha=.3,fill='white', color = 'black') +
  geom_bar(aes(y=Abundance),stat="identity",position ="identity",alpha=.8,fill="#ffb31a",color="#ffb31a") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab('Sequence Count')
dev.off()

compare_qim_even1_seq <- changing_false_negatives(compare_qim_even1)
compare_qim_even2_seq <- changing_false_negatives(compare_qim_even2)
compare_qim_stag1_seq <- changing_false_negatives(compare_qim_stag1)
compare_qim_stag2_seq <- changing_false_negatives(compare_qim_stag2)

tiff('sequence_error_qim_even1.TIF', width = 500, height = 400)
ggplot(data=compare_qim_even1_seq,aes(x=Family))+
  geom_bar(aes(y=Expected),stat="identity",position ="identity",alpha=.3,fill='white', color = 'black') +
  geom_bar(aes(y=Abundance),stat="identity",position ="identity",alpha=.8,fill="#5c5c8a",color="#5c5c8a") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab('Sequence Count')
dev.off()

tiff('sequence_error_qim_even2.TIF', width = 500, height = 400)
ggplot(data=compare_qim_even2_seq,aes(x=Family))+
  geom_bar(aes(y=Expected),stat="identity",position ="identity",alpha=.3,fill='white', color = 'black') +
  geom_bar(aes(y=Abundance),stat="identity",position ="identity",alpha=.8,fill="#5c5c8a",color="#5c5c8a") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab('Sequence Count')
dev.off()

tiff('sequence_error_qim_stag1.TIF', width = 500, height = 400)
ggplot(data=compare_qim_stag1_seq,aes(x=Family))+
  geom_bar(aes(y=Expected),stat="identity",position ="identity",alpha=.3,fill='white', color = 'black') +
  geom_bar(aes(y=Abundance),stat="identity",position ="identity",alpha=.8,fill="#5c5c8a",color="#5c5c8a") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab('Sequence Count')
dev.off()

tiff('sequence_error_qim_stag2.TIF', width = 500, height = 400)
ggplot(data=compare_qim_stag2_seq,aes(x=Family))+
  geom_bar(aes(y=Expected),stat="identity",position ="identity",alpha=.3,fill='white', color = 'black') +
  geom_bar(aes(y=Abundance),stat="identity",position ="identity",alpha=.8,fill="#5c5c8a",color="#5c5c8a") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab('Sequence Count')
dev.off()





