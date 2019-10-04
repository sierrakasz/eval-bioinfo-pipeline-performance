##Supplementary Material
# Each supplemental table is separated according to the order of tables in the manuscript


# Import Data (Pipeline comparison) -------------------------------------------------------------
#load required packages
#phyloseq requires BiocManager 
# code for phyloseq install: if(!requireNamespace("BiocManager")){
#install.packages("BiocManager")
#}
#BiocManager::install("phyloseq")
library(HMP)
library(phyloseq)
library(plyr)
library(PMCMR)
library(randomForest)
library(rsample)
library(tidyverse)
library(phyloseq)
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
metadata=(read.csv("Metadata/HPMMMeta_r_merge.csv",header=TRUE))
sampdat=sample_data(metadata)
sample_names(sampdat)=metadata$SampleID
physeq=merge_phyloseq(physeq_all, sampdat)
physeq



# Table S1 ----------------------------------------------------------------

#table was constructed in excel
# Table S2 ----------------------------------------------------------------

#table was constructed in excel
# Table S3 ----------------------------------------------------------------
#Table S3: Kullback-Leibler Divergence matrices among pipelines

#subset phyloseq objects by pipelines
physeq_q <- subset_samples(physeq, Pipeline == 'QIIME2')
physeq_m <- subset_samples(physeq, Pipeline == 'mothur')
physeq_mg <- subset_samples(physeq, Pipeline == 'MG-RAST')

#separate out otu tables from each phyloseq object as dataframes so we can work with them
#in the HMP package
data_mg <- as.data.frame(t(otu_table(physeq_mg)))
data_m <- as.data.frame(t(otu_table(physeq_m)))
data_q <- as.data.frame(t(otu_table(physeq_q)))

#group together - format for HMP package
group.data <- list(data_mg, data_m, data_q)

#run Kullback Leibler distance
kl <- Kullback.Leibler(group.data, plot = FALSE)
kl

# Table S4 --------------------------------------------------------
# Table S4: List of sample names that remained after filtering for both sample areas. 

#subset samples by sample area
physeq_rec <- subset_samples(physeq, Sample_Area == "Rectum")
physeq_mou <- subset_samples(physeq, Sample_Area == "Mouth")

#subset samples by pipelines
#print out list of sample names within each pipeline and sample area
#added to final table
physeq_Q_r <- subset_samples(physeq_rec, Pipeline == "QIIME2")
QR <- sample_names(physeq_Q_r)
View(QR)

physeq_M_r <- subset_samples(physeq_rec, Pipeline == "mothur")
MR <- sample_names(physeq_M_r)
View(MR)

physeq_G_r <- subset_samples(physeq_rec, Pipeline == "MG-RAST")
GR <- sample_names(physeq_G_r)
View(GR)

physeq_Q_m <- subset_samples(physeq_mou, Pipeline == "QIIME2")
QM <- sample_names(physeq_Q_m)
View(QM)

physeq_M_m <- subset_samples(physeq_mou, Pipeline == "mothur")
MM <- sample_names(physeq_M_m)
View(MM)

physeq_G_m <- subset_samples(physeq_mou, Pipeline == "MG-RAST")
GM <- sample_names(physeq_G_m)
View(GM)


# Table S5 ----------------------------------------------------------------

#data was analyzed using basic excel functions. 
# Table S6 --------------------------------------------------------
#Table S6: Confusion matrices of random forest classification of 1,000 decision trees

# subsetting samples among pipelines, and taxonomic levels
# q - QIIME2, m - mothur, g - MG-RAST, p - Phylum, f - family
physeq_q <- subset_samples(physeq, Pipeline == 'QIIME2')
physeq_m <- subset_samples(physeq, Pipeline == 'mothur')
physeq_mg <- subset_samples(physeq, Pipeline == 'MG-RAST')

physeq_q_p <- tax_glom(physeq_q, taxrank = 'Phylum')
physeq_m_p <- tax_glom(physeq_m, taxrank = 'Phylum')
physeq_mg_p <- tax_glom(physeq_mg, taxrank = 'Phylum')
physeq_q_f <- tax_glom(physeq_q, taxrank = 'Family')
physeq_m_f <- tax_glom(physeq_m, taxrank = 'Family')
physeq_mg_f <- tax_glom(physeq_mg, taxrank = 'Family')

#function for running random forest classification among sample areas
#output is the classification error matrix and overall classification error
random_foresting_sample_area <- function(data_phy) {
  fd = data_phy
  predictors=t(otu_table(fd))
  dim(predictors)
  resp <- as.factor(sample_data(fd)$Sample_Area)
  rf.data <- data.frame(resp, predictors)
  Forest <- randomForest(resp~., data=rf.data)
  return(Forest)
}

random_foresting_sample_area(physeq_q_p)
random_foresting_sample_area(physeq_m_p)
random_foresting_sample_area(physeq_mg_p)
random_foresting_sample_area(physeq_q_f)
random_foresting_sample_area(physeq_m_f) 
random_foresting_sample_area(physeq_mg_f)

#same as above, except classifying manner of death
# random forest classification error and overall error
random_foresting_MoD <- function(data_phy) {
  fd = data_phy
  predictors=t(otu_table(fd))
  dim(predictors)
  resp <- as.factor(sample_data(fd)$MoD)
  rf.data <- data.frame(resp, predictors)
  Forest <- randomForest(resp~., data=rf.data, ntree=1000)
  return(Forest)
}

random_foresting_MoD(physeq_q_p)
random_foresting_MoD(physeq_m_p)
random_foresting_MoD(physeq_mg_p)
random_foresting_MoD(physeq_q_f)
random_foresting_MoD(physeq_m_f)
random_foresting_MoD(physeq_mg_f)

# Table S7 ----------------------------------------------------------------
#Table S7:Mean relative abundances of phyla greater than 1% relative abundance among pipelines and sample areas
#Tax glom function cutoffs to specified taxonomic level
#Separated out sample area
GPrPhylum=tax_glom(physeq_rec, "Phylum")
#transforms data from sequence counts to relative abundances
PhylumLevel_rec = transform_sample_counts(GPrPhylum, function(x) x / sum(x))
#removes phyla with relative abundances level than .01% of average relative abundance
PhylumLevel_rec = filter_taxa(PhylumLevel_rec, function(x) mean(x) > 0.0001, TRUE) 

GPrPhylum=tax_glom(physeq_mou, "Phylum")
PhylumLevel_mou = transform_sample_counts(GPrPhylum, function(x) x / sum(x))
PhylumLevel_mou = filter_taxa(PhylumLevel_mou, function(x) mean(x) > 0.0001, TRUE) 

#Phyloseq object melted into dataframe for easier use
#caluclating relative abundance as overall percentage
#summary statistics at phylum level among pipeline including mean and standard deviation 
#for rectum
df <- psmelt(PhylumLevel_rec) 
df$Abundance=df$Abundance*100
Trtdata <- ddply(df, c("Phylum", "Pipeline"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)

#and mouth
df <- psmelt(PhylumLevel_mou) 
df$Abundance=df$Abundance*100
Trtdata <- ddply(df, c("Phylum", "Pipeline"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)

# Table S8 --------------------------------------------------------
# Table S8: Mean relative abundances of families greater than 1% relative abundance among pipelines and sample areas.

#See Table S7 for more description
#Repeat same steps as above, however specifying family level instead of phylum
GPrFamily=tax_glom(physeq_rec, "Family")
FamilyLevel_rec = transform_sample_counts(GPrFamily, function(x) x / sum(x))
FamilyLevel_rec = filter_taxa(FamilyLevel_rec, function(x) mean(x) > 0.0001, TRUE) 

GPrFamily=tax_glom(physeq_mou, "Family")
FamilyLevel_mou = transform_sample_counts(GPrFamily, function(x) x / sum(x))
FamilyLevel_mou = filter_taxa(FamilyLevel_mou, function(x) mean(x) > 0.0001, TRUE) 

df <- psmelt(FamilyLevel_rec) 
df$Abundance=df$Abundance*100
Trtdata <- ddply(df, c("Family", "Pipeline"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)

df <- psmelt(FamilyLevel_mou) 
df$Abundance=df$Abundance*100
Trtdata <- ddply(df, c("Family", "Pipeline"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)
# Table S9 ---------------------------------------------------------------
# Table S9: Analysis of the composition of microbiomes (ANCOM) result for differentially abundant phyla among pipelines and sample areas. 

#code from https://sites.google.com/site/siddharthamandal1985/research
ancom.W = function(otu_data,var_data,
                   adjusted,repeated,
                   main.var,adj.formula,
                   repeat.var,long,rand.formula,
                   multcorr,sig){
  
  n_otu=dim(otu_data)[2]-1
  
  otu_ids=colnames(otu_data)[-1]
  
  if(repeated==F){
    data_comp=data.frame(merge(otu_data,var_data,by="Sample.ID",all.y=T),row.names=NULL)
    #data_comp=data.frame(merge(otu_data,var_data[,c("Sample.ID",main.var)],by="Sample.ID",all.y=T),row.names=NULL)
  }else if(repeated==T){
    data_comp=data.frame(merge(otu_data,var_data,by="Sample.ID"),row.names=NULL)
    # data_comp=data.frame(merge(otu_data,var_data[,c("Sample.ID",main.var,repeat.var)],by="Sample.ID"),row.names=NULL)
  }
  
  base.formula = paste0("lr ~ ",main.var)
  if(repeated==T){
    repeat.formula = paste0(base.formula," | ", repeat.var)
  }
  if(adjusted==T){
    adjusted.formula = paste0(base.formula," + ", adj.formula)
  }
  
  if( adjusted == F & repeated == F ){
    fformula  <- formula(base.formula)
  } else if( adjusted == F & repeated == T & long == T ){
    fformula  <- formula(base.formula)   
  }else if( adjusted == F & repeated == T & long == F ){
    fformula  <- formula(repeat.formula)   
  }else if( adjusted == T & repeated == F  ){
    fformula  <- formula(adjusted.formula)   
  }else if( adjusted == T & repeated == T  ){
    fformula  <- formula(adjusted.formula)   
  }else{
    stop("Problem with data. Dataset should contain OTU abundances, groups, 
         and optionally an ID for repeated measures.")
  }
  
  
  
  if( repeated==FALSE & adjusted == FALSE){
    if( length(unique(data_comp[,which(colnames(data_comp)==main.var)]))==2 ){
      tfun <- exactRankTests::wilcox.exact
    } else{
      tfun <- stats::kruskal.test
    }
  }else if( repeated==FALSE & adjusted == TRUE){
    tfun <- stats::aov
  }else if( repeated== TRUE & adjusted == FALSE & long == FALSE){
    tfun <- stats::friedman.test
  }else if( repeated== TRUE & adjusted == FALSE & long == TRUE){
    tfun <- nlme::lme
  }else if( repeated== TRUE & adjusted == TRUE){
    tfun <- nlme::lme
  }
  
  logratio.mat <- matrix(NA, nrow=n_otu, ncol=n_otu)
  for(ii in 1:(n_otu-1)){
    for(jj in (ii+1):n_otu){
      data.pair <- data_comp[,which(colnames(data_comp)%in%otu_ids[c(ii,jj)])]
      lr <- log((1+as.numeric(data.pair[,1]))/(1+as.numeric(data.pair[,2])))
      
      lr_dat <- data.frame( lr=lr, data_comp,row.names=NULL )
      
      if(adjusted==FALSE&repeated==FALSE){  ## Wilcox, Kruskal Wallis
        logratio.mat[ii,jj] <- tfun( formula=fformula, data = lr_dat)$p.value
      }else if(adjusted==FALSE&repeated==TRUE&long==FALSE){ ## Friedman's 
        logratio.mat[ii,jj] <- tfun( formula=fformula, data = lr_dat)$p.value
      }else if(adjusted==TRUE&repeated==FALSE){ ## ANOVA
        model=tfun(formula=fformula, data = lr_dat,na.action=na.omit)   
        picker=which(gsub(" ","",row.names(summary(model)[[1]]))==main.var)  
        logratio.mat[ii,jj] <- summary(model)[[1]][["Pr(>F)"]][picker]
      }else if(repeated==TRUE&long==TRUE){ ## GEE
        model=tfun(fixed=fformula,data = lr_dat,
                   random = formula(rand.formula),
                   correlation=corAR1(),
                   na.action=na.omit)   
        picker=which(gsub(" ","",row.names(anova(model)))==main.var)
        logratio.mat[ii,jj] <- anova(model)[["p-value"]][picker]
      }
      
    }
  } 
  
  ind <- lower.tri(logratio.mat)
  logratio.mat[ind] <- t(logratio.mat)[ind]
  
  
  logratio.mat[which(is.finite(logratio.mat)==FALSE)] <- 1
  
  mc.pval <- t(apply(logratio.mat,1,function(x){
    s <- p.adjust(x, method = "BH")
    return(s)
  }))
  
  a <- logratio.mat[upper.tri(logratio.mat,diag=FALSE)==TRUE]
  
  b <- matrix(0,ncol=n_otu,nrow=n_otu)
  b[upper.tri(b)==T] <- p.adjust(a, method = "BH")
  diag(b)  <- NA
  ind.1    <- lower.tri(b)
  b[ind.1] <- t(b)[ind.1]
  
  #########################################
  ### Code to extract surrogate p-value
  surr.pval <- apply(mc.pval,1,function(x){
    s0=quantile(x[which(as.numeric(as.character(x))<sig)],0.95)
    # s0=max(x[which(as.numeric(as.character(x))<alpha)])
    return(s0)
  })
  #########################################
  ### Conservative
  if(multcorr==1){
    W <- apply(b,1,function(x){
      subp <- length(which(x<sig))
    })
    ### Moderate
  } else if(multcorr==2){
    W <- apply(mc.pval,1,function(x){
      subp <- length(which(x<sig))
    })
    ### No correction
  } else if(multcorr==3){
    W <- apply(logratio.mat,1,function(x){
      subp <- length(which(x<sig))
    })
  }
  
  return(W)
  }



ANCOM.main = function(OTUdat,Vardat,
                      adjusted,repeated,
                      main.var,adj.formula,
                      repeat.var,longitudinal,
                      random.formula,
                      multcorr,sig,
                      prev.cut){
  
  p.zeroes=apply(OTUdat[,-1],2,function(x){
    s=length(which(x==0))/length(x)
  })
  
  zeroes.dist=data.frame(colnames(OTUdat)[-1],p.zeroes,row.names=NULL)
  colnames(zeroes.dist)=c("Taxon","Proportion_zero")
  
  zero.plot = ggplot(zeroes.dist, aes(x=Proportion_zero)) + 
    geom_histogram(binwidth=0.1,colour="black",fill="white") + 
    xlab("Proportion of zeroes") + ylab("Number of taxa") +
    theme_bw()
  
  #print(zero.plot)
  
  OTUdat.thinned=OTUdat
  OTUdat.thinned=OTUdat.thinned[,c(1,1+which(p.zeroes<prev.cut))]
  
  otu.names=colnames(OTUdat.thinned)[-1]
  
  W.detected   <- ancom.W(OTUdat.thinned,Vardat,
                          adjusted,repeated,
                          main.var,adj.formula,
                          repeat.var,longitudinal,random.formula,
                          multcorr,sig)
  
  W_stat       <- W.detected
  
  
  ### Bubble plot
  
  W_frame = data.frame(otu.names,W_stat,row.names=NULL)
  W_frame = W_frame[order(-W_frame$W_stat),]
  
  W_frame$detected_0.9=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.8=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.7=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.6=rep(FALSE,dim(W_frame)[1])
  
  W_frame$detected_0.9[which(W_frame$W_stat>0.9*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.8[which(W_frame$W_stat>0.8*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.7[which(W_frame$W_stat>0.7*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.6[which(W_frame$W_stat>0.6*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  
  final_results=list(W_frame,zero.plot)
  names(final_results)=c("W.taxa","PLot.zeroes")
  return(final_results)
}

#Specify tax cut off as phylum
physeq_phy <- tax_glom(physeq, taxrank = 'Phylum')

#subset data for comparisons
comparisons_for_ancom <- function(physeq) {
  rec_phy <- subset_samples(physeq, Sample_Area == 'Rectum')
  mou_phy <- subset_samples(physeq, Sample_Area == 'Mouth')
  rec_phy_G <- subset_samples(rec_phy, Pipeline == 'MG-RAST')
  rec_phy_Q <- subset_samples(rec_phy, Pipeline == 'QIIME2')
  rec_phy_M <- subset_samples(rec_phy, Pipeline == 'mothur')
  mou_phy_G <- subset_samples(mou_phy, Pipeline == 'MG-RAST')
  mou_phy_Q <- subset_samples(mou_phy, Pipeline == 'QIIME2')
  mou_phy_M <- subset_samples(mou_phy, Pipeline == 'mothur')
  return(list(rec_phy_GM <-merge_phyloseq(rec_phy_G, rec_phy_M),
              rec_phy_QM <-merge_phyloseq(rec_phy_Q, rec_phy_M),
              rec_phy_GQ <-merge_phyloseq(rec_phy_G, rec_phy_Q),
              mou_phy_GM <-merge_phyloseq(mou_phy_G, mou_phy_M),
              mou_phy_QM <-merge_phyloseq(mou_phy_Q, mou_phy_M),
              mou_phy_GQ <-merge_phyloseq(mou_phy_G, mou_phy_Q)))
}

#add each subsetted phyloseq object into list
list_for_ancom_phy <- comparisons_for_ancom(physeq_phy)

#puts otu table into correct format required for ANCOM function
otu_ancom_make <- function(physeq) {
  otu_ancom <- data.frame(otu_table(physeq))
  otu_ancom <- data.frame(t(otu_ancom))
  Sample.ID <- rownames(otu_ancom)
  rownames(otu_ancom) <- NULL
  otu_ancom <- cbind(Sample.ID, otu_ancom)
  return(otu_ancom)
}

#metadata in correct format
metadata_ancom <- metadata
colnames(metadata_ancom)[1] <- 'Sample.ID'

#make empty list to add results to
ancom_results_phy <- list()

#for each subsetted phyloseq object, find differentially abundant phyla among pipelines
for(i in 1:length(list_for_ancom_phy)) {
  otu_ancom <- otu_ancom_make(list_for_ancom_phy[[i]])
  comparison_test <- ANCOM.main(OTUdat = otu_ancom,
                                Vardat = metadata_ancom,
                                adjusted = FALSE,
                                repeated = F,
                                main.var = "Pipeline",
                                adj.formula = NULL,
                                repeat.var=NULL,
                                longitudinal=FALSE,
                                random.formula=NULL,
                                multcorr=2,
                                sig=0.05,
                                prev.cut=0.90)
  w_values <- data.frame(comparison_test$W.taxa)
  tax <- data.frame(tax_table(physeq))
  tax <- tax %>% select(Kingdom, Phylum)
  tax$otu.names <- rownames(tax)
  ancom_sign_taxa <- merge(w_values, tax, by='otu.names')
  ancom_sign_taxa <- ancom_sign_taxa[,-1]
  ancom_results_phy[[i]] <- ancom_sign_taxa
}


# Table S10 ---------------------------------------------------------------
#Table S10: . Analysis of the composition of microbiomes (ANCOM) results for differentially abundant familes among pipelines and sample areas

#See Table S9 for more description
#Same as above, only with family
list_for_ancom_fam <- comparisons_for_ancom(physeq_fam)

otu_ancom_make <- function(physeq) {
  otu_ancom <- data.frame(otu_table(physeq))
  otu_ancom <- data.frame(t(otu_ancom))
  Sample.ID <- rownames(otu_ancom)
  rownames(otu_ancom) <- NULL
  otu_ancom <- cbind(Sample.ID, otu_ancom)
  return(otu_ancom)
}

metadata_ancom <- metadata
colnames(metadata_ancom)[1] <- 'Sample.ID'

ancom_results_fam <- list()

for(i in 1:length(list_for_ancom_fam)) {
  otu_ancom <- otu_ancom_make(list_for_ancom_fam[[i]])
  comparison_test <- ANCOM.main(OTUdat = otu_ancom,
                                Vardat = metadata_ancom,
                                adjusted = FALSE,
                                repeated = F,
                                main.var = "Pipeline",
                                adj.formula = NULL,
                                repeat.var=NULL,
                                longitudinal=FALSE,
                                random.formula=NULL,
                                multcorr=2,
                                sig=0.05,
                                prev.cut=0.90)
  w_values <- data.frame(comparison_test$W.taxa)
  tax <- data.frame(tax_table(physeq))
  tax <- tax %>% select(Kingdom, Phylum, Class, Order, Family)
  tax$otu.names <- rownames(tax)
  ancom_sign_taxa <- merge(w_values, tax, by='otu.names')
  ancom_sign_taxa <- ancom_sign_taxa[,-1]
  ancom_results_fam[[i]] <- ancom_sign_taxa
}

# Table S11 ---------------------------------------------------------------

#data was analyzed in excel
# Table S12 --------------------------------------------------------------
#Table S12: Alpha-diversity metrics for each sample among pipeline and sample area including: observed richness, 
#Chao1, Shannon diversity (1-D), and Inverse Simpson diversity (1/D)

#calculate alpha-div metrics
erich <- estimate_richness(physeq, measures = c("Observed", 'Chao1', "Shannon", "InvSimpson"))
erich <- add_rownames(erich, "SampleID")

#make data tidy 
erich2 <- erich %>%
  gather(Index, Observation, c("Observed", 'Chao1', "Shannon", "InvSimpson"), na.rm = TRUE)

#add metadata
rich = merge(erich2, metadata)
rich <- rich %>% select(SampleID, Index, Observation, SampleName, Pipeline,
                        Sample_Area)
rich %>% group_by(Index, Observation, Pipeline) %>% summarise_all(funs(mean, sd))

# Table S13 -----------------------------------------------------------
#Table S13: Results of A) Kruskal-Wallis and 
#B) pairwise post-hoc Nemenyi tests of alpha diversity metrics among pipeline and sample area

#subset data to make comparions for each metric among pipeline within each sample area
rich_obs_rec <- rich %>% filter(Index == "Observed", Sample_Area == 'Rectum')
rich_obs_mou <- rich %>% filter(Index == "Observed", Sample_Area == 'Mouth')
rich_cha_rec <- rich %>% filter(Index == "Chao1", Sample_Area == 'Rectum')
rich_cha_mou <- rich %>% filter(Index == "Chao1", Sample_Area == 'Mouth')
rich_sha_rec <- rich %>% filter(Index == "Shannon", Sample_Area == 'Rectum')
rich_sha_mou <- rich %>% filter(Index == "Shannon", Sample_Area == 'Mouth')
rich_inv_rec <- rich %>% filter(Index == "InvSimpson", Sample_Area == 'Rectum')
rich_inv_mou <- rich %>% filter(Index == "InvSimpson", Sample_Area == 'Mouth')

#kruskal-wallis test and post hoc nemenyi test with bonferroni correction 
#could have done this in a loop, but hey! Can't do everything perfect
kruskal.test(Observation ~ Pipeline, data = rich_obs_rec)
out <- posthoc.kruskal.nemenyi.test(x=rich_obs_rec$Observation, g=rich_obs_rec$Pipeline, dist='Tukey', p.adjust.method = 'bonf' )
print(out$p.value)

kruskal.test(Observation ~ Pipeline, data = rich_obs_mou)
out <- posthoc.kruskal.nemenyi.test(x=rich_obs_mou$Observation, g=rich_obs_mou$Pipeline, dist='Tukey', p.adjust.method = 'bonf' )
print(out$p.value)

kruskal.test(Observation ~ Pipeline, data = rich_cha_rec)
out <- posthoc.kruskal.nemenyi.test(x=rich_cha_rec$Observation, g=rich_cha_rec$Pipeline, dist='Tukey', p.adjust.method = 'bonf' )
print(out$p.value)

kruskal.test(Observation ~ Pipeline, data = rich_cha_mou)
out <- posthoc.kruskal.nemenyi.test(x=rich_cha_mou$Observation, g=rich_cha_mou$Pipeline, dist='Tukey', p.adjust.method = 'bonf' )
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


# Table S14 ------------------------------------------------------------
#Table S14: Results from permutational multivariate analysis of variance (PERMANOVA) tests on weighted Jaccard 
#distance matrix including A) beta diversity and B) beta dispersion for 999 permutations.  

#subset out by pipeline to get beta-div calculations
physeq_q <- subset_samples(physeq, Pipeline == 'QIIME2')
physeq_m <- subset_samples(physeq, Pipeline == 'mothur')
physeq_mg <- subset_samples(physeq, Pipeline == 'MG-RAST')

#subset out within sample area and pipeline for beta-div and beta-disp pairwise comparisons
physeq_rec <- subset_samples(physeq, Sample_Area == "Rectum")
physeq_mou <- subset_samples(physeq, Sample_Area == "Mouth")

physeq_Q_r <- subset_samples(physeq_rec, Pipeline == "QIIME2")
physeq_M_r <- subset_samples(physeq_rec, Pipeline == "mothur")
physeq_G_r <- subset_samples(physeq_rec, Pipeline == "MG-RAST")

physeq_QM_r <- merge_phyloseq(physeq_Q_r, physeq_M_r)
physeq_QG_r <- merge_phyloseq(physeq_Q_r, physeq_G_r)
physeq_GM_r <- merge_phyloseq(physeq_M_r, physeq_G_r)

physeq_Q_m <- subset_samples(physeq_mou, Pipeline == "QIIME2")
physeq_M_m <- subset_samples(physeq_mou, Pipeline == "mothur")
physeq_G_m <- subset_samples(physeq_mou, Pipeline == "MG-RAST")

physeq_QM_m <- merge_phyloseq(physeq_Q_m, physeq_M_m)
physeq_QG_m <- merge_phyloseq(physeq_Q_m, physeq_G_m)
physeq_GM_m <- merge_phyloseq(physeq_M_m, physeq_G_m)

#calculating average jaccard distance for each pipeline
#included in text
otu <- data.frame(otu_table(physeq_m))
otu= otu[rowSums(otu)!=0,] 
beta <- vegdist(otu, method = "jaccard")
mean(beta)
sd(beta)

otu <- data.frame(otu_table(physeq_q))
otu= otu[rowSums(otu)!=0,] 
beta <- vegdist(otu, method = "jaccard")
mean(beta)
sd(beta)

otu <- data.frame(otu_table(physeq_mg))
otu= otu[rowSums(otu)!=0,] 
beta <- vegdist(otu, method = "jaccard")
mean(beta)
sd(beta)

#functions for testing beta-diversity and beta-dispersion PERMANOVA comparisons among pipeline
beta_diversity_calc <- function(physeq) {
  GPdist=phyloseq::distance(physeq, "(A+B-2*J)/(A+B-J)")
  sampledf <- data.frame(sample_data(physeq))
  print(return(adonis(GPdist ~ Pipeline, data = sampledf)))
}

beta_dispersion_calc <- function(physeq) {
  GPdist=phyloseq::distance(physeq, "(A+B-2*J)/(A+B-J)")
  sampledf <- data.frame(sample_data(physeq))
  beta <- betadisper(GPdist, sampledf$Pipeline)
  print(return(permutest(beta)))
}

#results for PERMANOVA
beta_diversity_calc(physeq)
beta_dispersion_calc(physeq)

beta_diversity_calc(physeq_rec)
beta_dispersion_calc(physeq_rec)

beta_diversity_calc(physeq_mou)
beta_dispersion_calc(physeq_mou)

beta_diversity_calc(physeq_QG_r)
beta_dispersion_calc(physeq_QG_r)

beta_diversity_calc(physeq_QM_r)
beta_dispersion_calc(physeq_QM_r)

beta_diversity_calc(physeq_GM_r)
beta_dispersion_calc(physeq_GM_r)

beta_diversity_calc(physeq_QG_m)
beta_dispersion_calc(physeq_QG_m)

beta_diversity_calc(physeq_QM_m)
beta_dispersion_calc(physeq_QM_m)

beta_diversity_calc(physeq_GM_m)
beta_dispersion_calc(physeq_GM_m)


# Table S15 ---------------------------------------------------------------
#Table S15: Top twenty predictor taxa based on mean decrease in Gini for random forest 
#classifying manner of death and sample area.  

# subsetting samples among pipelines, and taxonomic levels
# q - QIIME2, m - mothur, g - MG-RAST, p - Phylum, f - family
physeq_q <- subset_samples(physeq, Pipeline == 'QIIME2')
physeq_m <- subset_samples(physeq, Pipeline == 'mothur')
physeq_mg <- subset_samples(physeq, Pipeline == 'MG-RAST')

physeq_q_p <- tax_glom(physeq_q, taxrank = 'Phylum')
physeq_m_p <- tax_glom(physeq_m, taxrank = 'Phylum')
physeq_mg_p <- tax_glom(physeq_mg, taxrank = 'Phylum')
physeq_q_f <- tax_glom(physeq_q, taxrank = 'Family')
physeq_m_f <- tax_glom(physeq_m, taxrank = 'Family')
physeq_mg_f <- tax_glom(physeq_mg, taxrank = 'Family')

#function for running random forest classification among sample areas
#output is the classification error matrix and overall classification error
random_foresting_sample_area <- function(data_phy) {
  fd = data_phy
  predictors=t(otu_table(fd))
  dim(predictors)
  resp <- as.factor(sample_data(fd)$Sample_Area)
  rf.data <- data.frame(resp, predictors)
  Forest <- randomForest(resp~., data=rf.data)
  return(Forest)
}

forest_q_p <- random_foresting_sample_area(physeq_q_p)
forest_m_p <- random_foresting_sample_area(physeq_m_p)
forest_mg_p <- random_foresting_sample_area(physeq_mg_p)
forest_q_f <- random_foresting_sample_area(physeq_q_f)
forest_m_f <- random_foresting_sample_area(physeq_m_f)
forest_mg_f <- random_foresting_sample_area(physeq_mg_f)

#function for determing indicator taxa for random forest classification among sample areas
#output is top 20 predictors based on mean decrease in Gini
#function writes predictors to a csv file
forest_predictors_sample_area <- function(forest) {
  imp <- importance(forest)
  imp <- data.frame(predictors = rownames(imp), imp)
  imp.sort <- arrange(imp, desc(MeanDecreaseGini))
  imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)
  imp.20 <- imp.sort[1:20, ]
  tax <- data.frame(tax_table(physeq))
  head(tax)
  tax <- tax %>% select(Kingdom, Phylum, Class, Order, Family)
  tax$predictors <- rownames(tax)
  imp.20 <- merge(imp.20, tax)
  write.table(imp.20, sep=",", "random_forest_predictors_pipe_sa.csv", append = TRUE)
  return(imp.20)
}

#resulting output went into table
pred_q_f <- forest_predictors_sample_area(forest_q_f)
pred_m_f <- forest_predictors_sample_area(forest_m_f)
pred_mg_f <- forest_predictors_sample_area(forest_mg_f)

#same as above, except classifying manner of death
# random forest classification error and overall error
random_foresting_MoD <- function(data_phy) {
  fd = data_phy
  predictors=t(otu_table(fd))
  dim(predictors)
  resp <- as.factor(sample_data(fd)$MoD)
  rf.data <- data.frame(resp, predictors)
  Forest <- randomForest(resp~., data=rf.data, ntree=1000)
  return(Forest)
}

forest_q_p <- random_foresting_MoD(physeq_q_p)
forest_m_p <- random_foresting_MoD(physeq_m_p)
forest_mg_p <- random_foresting_MoD(physeq_mg_p)
forest_q_f <- random_foresting_MoD(physeq_q_f)
forest_m_f <- random_foresting_MoD(physeq_m_f)
forest_mg_f <- random_foresting_MoD(physeq_mg_f)

#indicator taxa for manner of death classification 
forest_predictors_MoD <- function(forest) {
  imp <- importance(forest)
  imp <- data.frame(predictors = rownames(imp), imp)
  imp.sort <- arrange(imp, desc(MeanDecreaseGini))
  imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)
  imp.20 <- imp.sort[1:20, ]
  tax <- data.frame(tax_table(physeq))
  head(tax)
  tax <- tax %>% select(Kingdom, Phylum, Class, Order, Family)
  tax$predictors <- rownames(tax)
  imp.20 <- merge(imp.20, tax)
  write.table(imp.20, sep=",", "random_forest_predictors_pipe_md.csv", append = TRUE)
  return(imp.20)
}

#resulting output went into table
pred_q_f <- forest_predictors_MoD(forest_q_f)
pred_m_f <- forest_predictors_MoD(forest_m_f)
pred_mg_f <- forest_predictors_MoD(forest_mg_f)

# Table S16 -----------------------------------------------------------
#Table S16: . The number of correct taxa, false negatives, and 
#false positives for in silico data among pipelines and community types (even and staggered).

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
#values were put into a table
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




