## Figure 4
# Multipanel figure with microbial community metrics among minimum library size and sample size
# including: relative abundances, alpha-diversity, and beta-diversity

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


# relative abundance ------------------------------------------------------


##relative abundance
#tax glom is a function that collapses taxonomy into phylum
#also removing very rare taxa (less than 1%)
GPrPhylum=tax_glom(physeq_nr, "Phylum")
PhylumLevel = transform_sample_counts(GPrPhylum, function(x) x / sum(x))
PhylumLevel = filter_taxa(PhylumLevel, function(x) mean(x) > 0.0001, TRUE) 

#melting into dataframe to calculate abundances 
df <- psmelt(PhylumLevel) 
df$Abundance=df$Abundance*100
Trtdata <- ddply(df, c("Phylum", "Subsample"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)

#pull out top five phyla
Topfive_trtdata <- filter(Trtdata, Phylum %in% c("Firmicutes", "Bacteroidetes",
                                                 "Actinobacteria", "Proteobacteria",
                                                 "Fusobacteria")) 
Topfive_trtdata$Subsample <- factor(Topfive_trtdata$Subsample, levels = c('60', '120', '188') )

#figure 4, panel A - among sample size
a <- ggplot(Topfive_trtdata, aes(x=Subsample, y=mean, fill=Phylum)) + 
  geom_bar(stat = 'identity', position = position_stack()) + ylab("Relative Abundance (> 1%)") +
  scale_x_discrete(labels=c('MG-RAST' = 'MG', 'mothur' = 'M', "QIIME2" = 'Q')) +
  scale_fill_manual(values=  c("#FF0000", "#006699", "#FF6600", "#FFCC00", "#00A08A"))
a

#repeat steps but for minimum library sizes
GPrPhylum=tax_glom(physeq_188, "Phylum")
PhylumLevel = transform_sample_counts(GPrPhylum, function(x) x / sum(x))
PhylumLevel = filter_taxa(PhylumLevel, function(x) mean(x) > 0.0001, TRUE) 

df <- psmelt(PhylumLevel) 
df$Abundance=df$Abundance*100
Trtdata <- ddply(df, c("Phylum", "Rarefy"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)

Topfive_trtdata <- filter(Trtdata, Phylum %in% c("Firmicutes", "Bacteroidetes",
                                                 "Actinobacteria", "Proteobacteria",
                                                 "Fusobacteria")) 

#figure 4, panel b - among minimum library sizes 
b <- ggplot(Topfive_trtdata, aes(x=Rarefy, y=mean, fill=Phylum)) + 
  geom_bar(stat = 'identity', position = position_stack()) + ylab("Relative Abundance (> 1%)") +
  scale_fill_manual(values=  c("#FF0000", "#006699", "#FF6600", "#FFCC00", "#00A08A")) +
  scale_x_discrete(labels=c("No Rarefaction" = 'NR'))
b


# alp div -----------------------------------------------------------------


## calculate alpha diversity
erichnr <- estimate_richness(merger, measures = c("Observed", 'Chao1', "Shannon", "InvSimpson"))
erichnr <- add_rownames(erichnr, "SampleID")

#make data tidy 
erichnr <- erichnr %>%
  gather(Index, Observation, c("Observed", 'Chao1', "Shannon", "InvSimpson"), na.rm = TRUE)

#add metadata
rich = merge(erichnr, metadata_rare)
rich <- rich %>% select(SampleID, Index, Observation, SampleName, Rarefy,
                        Subsample, Rarefy)
#make everything a factor and order it for the subsequent plot
rich$Index <- factor(rich$Index, levels = c("Observed", 'Chao1', "Shannon", "InvSimpson"))
rich$Subsample <- factor(rich$Subsample, levels = c('60', '120', '188'))
rich$Rarefy <- factor(rich$Rarefy, levels = c("No Rarefaction", "7000", "1000"))

#figure 4, panel C - among sample size
p <- ggplot(rich, aes(x=Subsample, y=Observation, fill=Subsample)) +
  geom_boxplot() + ylab("Counts") + xlab("")
c <- p + facet_wrap(~Index, scales="free") + scale_fill_manual(values = c('#273746', '#F5B041', 'brown')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
c 

#figure 4, panel d - among minimum library size
p <- ggplot(rich, aes(x=Rarefy, y=Observation, fill=Rarefy)) +
  geom_boxplot() + ylab("Counts") + xlab("Rarefaction Level")
d <- p + facet_wrap(~Index, scales="free") + scale_fill_manual(values = c("#512465", "#197C75", "#DAF7A6")) +
  scale_x_discrete(labels=c("No Rarefaction" = 'NR')) + xlab("")  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
d


# beta div ----------------------------------------------------------------


##beta diversity
#Beta-div values using Jaccard index and PCoA plots
ord = ordinate(physeq_nr, method="PCoA", distance="(A+B-2*J)/(A+B-J)")
ordplot=plot_ordination(physeq_nr, ord, color="Subsample")
ordplot
#Figure 4, panel e - PCoA among sample sizes
e <- ordplot+ geom_point(size = 3, alpha = 1/5, aes(color = factor(sample_data(physeq_nr)$Subsample))) + scale_color_manual(values = c('#273746', '#F5B041', 'brown'))
e

ord = ordinate(physeq_188, method="PCoA", distance="(A+B-2*J)/(A+B-J)")
ordplot=plot_ordination(physeq_188, ord, color="Rarefy")
ordplot 
#Figure 4, panel f - PCoA among minimum library size 
f <- ordplot+ geom_point(size = 3, alpha=1/5, aes(color = factor(sample_data(physeq_188)$Rarefy))) + scale_color_manual(values = c("#512465", "#197C75", "#DAF7A6"))
f


# Figure ------------------------------------------------------------------


#Final Figure
theme_set(theme_classic(base_size = 20))
tiff("Fig4.TIF", width = 4700, height = 5500, res=300)
ggarrange(a,b,c,d,e,f + rremove("x.text"), 
          labels = c("A", "B", "C", "D", "E", "F"),
          nrow = 3, ncol = 2)
dev.off()

