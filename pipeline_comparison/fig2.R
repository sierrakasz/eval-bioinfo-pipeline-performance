### Making Figure 2

#packages
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
##combining the files made it easier for collapsing differences among the pipeline outputs
#tax.group is available on github

tax_group <- read.csv("tax_group.csv")
tax_group <- tax_group %>% group_by(Kingdom, Phylum, Class, Order, Family, Genus) %>% summarise_all(funs(sum))
tax_group <- tax_group %>%  ungroup()

#this function specifically collapses the OTU table and taxonomy file 
#when the taxonomy is the same name among pipelines
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

#import metadata and combine into one phyloseq object
metadata=(read.csv("HPMMMeta_r_merge.csv",header=TRUE))
sampdat=sample_data(metadata)
sample_names(sampdat)=metadata$SampleID
physeq=merge_phyloseq(physeq_all, sampdat)
physeq

#subsetting among both body sites compared
physeq_rec <- subset_samples(physeq, Sample_Area == "Rectum")
physeq_mou <- subset_samples(physeq, Sample_Area == "Mouth")


#Tax glom is a function that collapses the taxonomic levels into phylum
#also are removing very rare taxa (less than 1%)
GPrPhylum=tax_glom(physeq_rec, "Phylum")
PhylumLevel_rec = transform_sample_counts(GPrPhylum, function(x) x / sum(x))
PhylumLevel_rec = filter_taxa(PhylumLevel_rec, function(x) mean(x) > 0.0001, TRUE) 

GPrPhylum=tax_glom(physeq_mou, "Phylum")
PhylumLevel_mou = transform_sample_counts(GPrPhylum, function(x) x / sum(x))
PhylumLevel_mou = filter_taxa(PhylumLevel_mou, function(x) mean(x) > 0.0001, TRUE) 

#Relative abundance
#split by sample area
#figures include top five phyla, indciated by mean determined in Trtdata
df <- psmelt(PhylumLevel_mou) 
df$Abundance=df$Abundance*100
Trtdata <- ddply(df, c("Phylum", "Pipeline"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)

Topfive_trtdata <- filter(Trtdata, Phylum %in% c("Firmicutes", "Bacteroidetes",
                                                 "Actinobacteria", "Proteobacteria",
                                                 "Fusobacteria")) 
#figure 2, panel A
a <- ggplot(Topfive_trtdata, aes(x=Pipeline, y=mean, fill=Phylum)) + 
  geom_bar(stat = 'identity', position = position_stack()) +
  xlab("Mouth") + ylab("Relative Abundance Phylum Level (> 1%)") +
  scale_x_discrete(labels=c('MG-RAST' = 'MG', 'mothur' = 'M', "QIIME2" = 'Q')) +
  scale_fill_manual(values=  c("#FF0000", "#006699", "#FF6600", "#FFCC00", "#00A08A")) +
  theme(legend.position="left")

df <- psmelt(PhylumLevel_rec) 
df$Abundance=df$Abundance*100
Trtdata <- ddply(df, c("Phylum", "Pipeline"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)

Topfive_trtdata <- filter(Trtdata, Phylum %in% c("Firmicutes", "Bacteroidetes",
                                                 "Actinobacteria", "Proteobacteria",
                                                 "Fusobacteria")) 
                  
#figure 2, panel B                                                                                "Fusobacteria")) 
b <- ggplot(Topfive_trtdata, aes(x=Pipeline, y=mean, fill=Phylum)) + 
  geom_bar(stat = 'identity', position = position_stack()) +
  xlab("Rectum") + ylab("Relative Abundance Phylum Level (> 1%)") +
  scale_x_discrete(labels=c('MG-RAST' = 'MG', 'mothur' = 'M', "QIIME2" = 'Q')) +
  scale_fill_manual(values=  c("#FF0000", "#006699", "#FF6600", "#FFCC00", "#00A08A")) +
  theme(legend.position="left")

#Alpha Diversity
erich <- estimate_richness(physeq, measures = c("Observed", "Shannon", "InvSimpson"))
erich <- add_rownames(erich, "SampleID")

#make data tidy 
erich2 <- erich %>%
  gather(Index, Observation, c("Observed", "Shannon", "InvSimpson"), na.rm = TRUE)

#add metadata
rich = merge(erich2, metadata)
rich <- rich %>% select(SampleID, Index, Observation, SampleName, Pipeline,
                        Sample_Area)

rich$Index <- factor(rich$Index, levels = c("Observed", "Shannon", "InvSimpson"))
#figure 2, panel C
p <- ggplot(rich, aes(x=Pipeline, y=Observation, fill=Pipeline)) +
  geom_boxplot(outlier.size = 3)
c <- p + facet_grid(Index~Sample_Area, scales="free") + scale_fill_manual(values = c("#787878", "#ffb31a", "#5c5c8a")) +
  scale_x_discrete(labels=c('MG-RAST' = 'MG', 'mothur' = 'M', "QIIME2" = 'Q')) +
  xlab("") + ylab('Counts') + theme_bw(base_size = 20) + theme(panel.grid.major = element_blank(),
                                                         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.position = "top")

#beta diversity
#ellipses indicate sample area
ord = ordinate(physeq, method="PCoA", distance="(A+B-2*J)/(A+B-J)")
ordplot=plot_ordination(physeq, ord, color="Pipeline")
#figure 2, panel D
d <- ordplot+ geom_point(size = 3,aes(color = sample_data(physeq)$Pipeline)) + scale_color_manual(values = c("#787878", "#ffb31a", "black", "#5c5c8a", "#918151")) +
  stat_ellipse(alpha = 1, size=2, aes(color= Sample_Area)) +
  theme(legend.position = 'top')
d


#Final Figure
theme_set(theme_classic(base_size = 16))
tiff("Fig2.TIF", width = 1000, height = 1000)
ggarrange(a,b,c,d + rremove("x.text"), 
          labels = c("A", "B", "C", "D"),
          nrow = 2, ncol = 2)
dev.off()

