### Making Figure 4

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


# import data -------------------------------------------------------------


###tax group files that combine the OTU file and taxonomy file. 
#tax.group is available on github
#tax groups nr, 7000, and 1,000 include all samples at different rarefaction levels
#tax group rare includes everything
tax_group_nr <- read.csv("tax_group_norare.csv")
tax_group_nr <- tax_group_nr %>% group_by(Kingdom, Phylum, Class, Order, Family, Genus) %>% summarise_all(funs(sum))
tax_group_nr <- tax_group_nr %>%  ungroup()

tax_group_7 <- read.csv("tax_group_7000.csv")
tax_group_7 <- tax_group_7 %>% group_by(Kingdom, Phylum, Class, Order, Family, Genus) %>% summarise_all(funs(sum))
tax_group_7 <- tax_group_7 %>%  ungroup()

tax_group_1 <- read.csv("tax_group_1000.csv")
tax_group_1 <- tax_group_1 %>% group_by(Kingdom, Phylum, Class, Order, Family, Genus) %>% summarise_all(funs(sum))
tax_group_1 <- tax_group_1 %>%  ungroup()

tax_group_rare <- read.csv("tax_group_rare.csv")
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

physeq_nr <- regroup_physeq_object(tax_group_nr)
physeq_7 <- regroup_physeq_object(tax_group_7)
physeq_1 <- regroup_physeq_object(tax_group_1)
physeq_rare <- regroup_physeq_object(tax_group_rare)

#import metadata and combine into one phyloseq object
metadata=(read.csv("Metadata/HPMMMeta_r_merge_rare.csv",header=TRUE))
sampdat=sample_data(metadata)
sample_names(sampdat)=metadata$SampleID
sampdat$Subsample <- factor(sampdat$Subsample, levels = c('60', '120', '188'))

physeq_nr=merge_phyloseq(physeq_nr, sampdat)
physeq_7=merge_phyloseq(physeq_7, sampdat)
physeq_1=merge_phyloseq(physeq_1, sampdat)

metadata_rare=(read.csv("Metadata/HPMMMeta_r_merge_rare_levels.csv",header=TRUE))
sampdat_rare=sample_data(metadata_rare)
sampdat_rare$Rarefy <- factor(sampdat_rare$Rarefy, levels = c("No Rarefaction", "7000", "1000"))
sample_names(sampdat_rare)=metadata_rare$SampleID
merger = merge_phyloseq(physeq_rare, sampdat_rare)
sample_data(merger)$Rarefy <- factor(sample_data(merger)$Rarefy, levels = c("No Rarefaction", "7000", "1000"))

physeq_188 <- subset_samples(merger, Subsample == '188')


# relative abundance ------------------------------------------------------


##relative abundance
#tax glom is a function that collapses taxonomy into phylum
GPrPhylum=tax_glom(physeq_nr, "Phylum")
PhylumLevel = transform_sample_counts(GPrPhylum, function(x) x / sum(x))
PhylumLevel = filter_taxa(PhylumLevel, function(x) mean(x) > 0.0001, TRUE) 

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

#figure 4, panel A
a <- ggplot(Topfive_trtdata, aes(x=Subsample, y=mean, fill=Phylum)) + 
  geom_bar(stat = 'identity', position = position_stack()) + ylab("Relative Abundance (> 1%)") +
  scale_x_discrete(labels=c('MG-RAST' = 'MG', 'mothur' = 'M', "QIIME2" = 'Q')) +
  scale_fill_manual(values=  c("#FF0000", "#006699", "#FF6600", "#FFCC00", "#00A08A"))
a

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


b <- ggplot(Topfive_trtdata, aes(x=Rarefy, y=mean, fill=Phylum)) + 
  geom_bar(stat = 'identity', position = position_stack()) + ylab("Relative Abundance (> 1%)") +
  scale_fill_manual(values=  c("#FF0000", "#006699", "#FF6600", "#FFCC00", "#00A08A")) +
  scale_x_discrete(labels=c("No Rarefaction" = 'NR'))
b


# alp div -----------------------------------------------------------------


##alpha diversity
erichnr <- estimate_richness(merger, measures = c("Observed", "Shannon", "InvSimpson"))
erichnr <- add_rownames(erichnr, "SampleID")

#make data tidy 
erichnr <- erichnr %>%
  gather(Index, Observation, c("Observed", "Shannon", "InvSimpson"), na.rm = TRUE)

#add metadata
rich = merge(erichnr, metadata_rare)
rich <- rich %>% select(SampleID, Index, Observation, SampleName, Rarefy,
                        Subsample, Rarefy)
rich$Index <- factor(rich$Index, levels = c("Observed", "Shannon", "InvSimpson"))
rich$Subsample <- factor(rich$Subsample, levels = c('60', '120', '188'))
rich$Rarefy <- factor(rich$Rarefy, levels = c("No Rarefaction", "7000", "1000"))

#figure 4, panel C
p <- ggplot(rich, aes(x=Subsample, y=Observation, fill=Subsample)) +
  geom_boxplot() + ylab("Counts") + xlab("")
c <- p + facet_wrap(~Index, scales="free") + scale_fill_manual(values = c('#273746', '#F5B041', 'brown')) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
c 

#figure 4, panel d
p <- ggplot(rich, aes(x=Rarefy, y=Observation, fill=Rarefy)) +
  geom_boxplot() + ylab("Counts") + xlab("Rarefaction Level")
d <- p + facet_wrap(~Index, scales="free") + scale_fill_manual(values = c("#512465", "#197C75", "#DAF7A6")) +
  scale_x_discrete(labels=c("No Rarefaction" = 'NR')) + xlab("")  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
d


# beta div ----------------------------------------------------------------


##beta diversity
ord = ordinate(physeq_nr, method="PCoA", distance="(A+B-2*J)/(A+B-J)")
ordplot=plot_ordination(physeq_nr, ord, color="Subsample")
ordplot
e <- ordplot+ geom_point(size = 3, aes(color = factor(sample_data(physeq_nr)$Subsample))) + scale_color_manual(values = c('#273746', '#F5B041', 'brown'))
e

ord = ordinate(physeq_188, method="PCoA", distance="(A+B-2*J)/(A+B-J)")
ordplot=plot_ordination(physeq_188, ord, color="Rarefy")
ordplot 
f <- ordplot+ geom_point(size = 3, aes(color = factor(sample_data(physeq_188)$Rarefy))) + scale_color_manual(values = c("#512465", "#197C75", "#DAF7A6"))
f


# Figure ------------------------------------------------------------------


#Final Figure
theme_set(theme_classic(base_size = 20))
tiff("Fig4.TIF", width = 4700, height = 5000, res=300)
ggarrange(a,b,c,d,e,f + rremove("x.text"), 
          labels = c("A", "B", "C", "D", "E", "F"),
          nrow = 3, ncol = 2)
dev.off()

