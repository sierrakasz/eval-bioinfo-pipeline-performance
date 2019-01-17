rm(list=ls())
getwd()
setwd("C:/Users/sierr/Documents/ThesisProject_Pipelines/")

#packages
library(ggforce)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(PMCMR)
library(phyloseq)
library(plyr)
library(tidyverse)
library(exactRankTests)
library(nlme)

#make a table in excel for rare data
tax_group_nr <- read.csv("tax_group_norare.csv")
tax_group_nr <- tax_group_nr %>% group_by(Kingdom, Phylum, Class, Order, Family, Genus) %>% summarise_all(funs(sum))
tax_group_nr <- tax_group_nr %>%  ungroup()

tax_group_7 <- read.csv("tax_group_7000.csv")
tax_group_7 <- tax_group_7 %>% group_by(Kingdom, Phylum, Class, Order, Family, Genus) %>% summarise_all(funs(sum))
tax_group_7 <- tax_group_7 %>%  ungroup()

tax_group_1 <- read.csv("tax_group_1000.csv")
tax_group_1 <- tax_group_1 %>% group_by(Kingdom, Phylum, Class, Order, Family, Genus) %>% summarise_all(funs(sum))
tax_group_1 <- tax_group_1 %>%  ungroup()

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

#import metadata and combine
metadata=(read.csv("Metadata/HPMMMeta_r_merge_rare.csv",header=TRUE))
sampdat=sample_data(metadata)
sample_names(sampdat)=metadata$SampleID
sampdat$Subsample <- factor(sampdat$Subsample, levels = c('60', '120', '188'))

physeq_nr=merge_phyloseq(physeq_nr, sampdat)
physeq_7=merge_phyloseq(physeq_7, sampdat)
physeq_1=merge_phyloseq(physeq_1, sampdat)

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

physeq_nr_phy <- tax_glom(physeq_nr, taxrank = 'Phylum')
physeq_nr_fam <- tax_glom(physeq_nr, taxrank = 'Family')
physeq_7_phy <- tax_glom(physeq_7, taxrank = 'Phylum')
physeq_7_fam <- tax_glom(physeq_7, taxrank = 'Family')
physeq_1_phy <- tax_glom(physeq_1, taxrank = 'Phylum')
physeq_1_fam <- tax_glom(physeq_1, taxrank = 'Family')

comparisons_for_ancom_subsample <- function(physeq) {
  rec_phy <- subset_samples(physeq, Sample_Area == 'Rectum')
  mou_phy <- subset_samples(physeq, Sample_Area == 'Mouth')
  physeq_rec_188 <- subset_samples(rec_phy, Subsample=='188')
  physeq_rec_120 <- subset_samples(rec_phy, Subsample=='120')
  physeq_rec_60 <- subset_samples(rec_phy, Subsample=='60')
  physeq_mou_188 <- subset_samples(mou_phy, Subsample=='188')
  physeq_mou_120 <- subset_samples(mou_phy, Subsample=='120')
  physeq_mou_60 <- subset_samples(mou_phy, Subsample=='60')
  return(list(physeq_rec_60120 <- merge_phyloseq(physeq_rec_60, physeq_rec_120),
              physeq_rec_60188 <- merge_phyloseq(physeq_rec_60, physeq_rec_188),
              physeq_rec_120188 <- merge_phyloseq(physeq_rec_188, physeq_rec_120),
              physeq_mou_60120 <- merge_phyloseq(physeq_mou_60, physeq_mou_120),
              physeq_mou_60188 <- merge_phyloseq(physeq_mou_60, physeq_mou_188),
              physeq_mou_120188 <- merge_phyloseq(physeq_mou_188, physeq_mou_120)))
}

list_for_ancom_phy_nr <- comparisons_for_ancom_subsample(physeq_nr_phy) 
list_for_ancom_fam_nr <- comparisons_for_ancom_subsample(physeq_nr_fam) 
list_for_ancom_phy_7 <- comparisons_for_ancom_subsample(physeq_7_phy) 
list_for_ancom_fam_7 <- comparisons_for_ancom_subsample(physeq_7_fam) 
list_for_ancom_phy_1 <- comparisons_for_ancom_subsample(physeq_1_phy) 
list_for_ancom_fam_1 <- comparisons_for_ancom_subsample(physeq_1_fam) 

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

ancom_results_phy_nr <- list()
ancom_results_fam_nr <- list()
ancom_results_phy_7 <- list()
ancom_results_fam_7 <- list()
ancom_results_phy_1 <- list()
ancom_results_fam_1 <- list()

for(i in 1:length(list_for_ancom_phy_nr)) {
  otu_ancom <- otu_ancom_make(list_for_ancom_phy_nr[[i]])
  comparison_test <- ANCOM.main(OTUdat = otu_ancom,
                                Vardat = metadata_ancom,
                                adjusted = FALSE,
                                repeated = F,
                                main.var = "Subsample",
                                adj.formula = NULL,
                                repeat.var=NULL,
                                longitudinal=FALSE,
                                random.formula=NULL,
                                multcorr=2,
                                sig=0.05,
                                prev.cut=0.90)
  w_values <- data.frame(comparison_test$W.taxa)
  tax <- data.frame(tax_table(physeq_nr))
  tax <- tax %>% select(Kingdom, Phylum)
  tax$otu.names <- rownames(tax)
  ancom_sign_taxa <- merge(w_values, tax, by='otu.names')
  ancom_sign_taxa <- ancom_sign_taxa[,-1]
  ancom_results_phy_nr[[i]] <- ancom_sign_taxa
}

for(i in 1:length(list_for_ancom_phy_7)) {
  otu_ancom <- otu_ancom_make(list_for_ancom_phy_7[[i]])
  comparison_test <- ANCOM.main(OTUdat = otu_ancom,
                                Vardat = metadata_ancom,
                                adjusted = FALSE,
                                repeated = F,
                                main.var = "Subsample",
                                adj.formula = NULL,
                                repeat.var=NULL,
                                longitudinal=FALSE,
                                random.formula=NULL,
                                multcorr=2,
                                sig=0.05,
                                prev.cut=0.90)
  w_values <- data.frame(comparison_test$W.taxa)
  tax <- data.frame(tax_table(physeq_7))
  tax <- tax %>% select(Kingdom, Phylum)
  tax$otu.names <- rownames(tax)
  ancom_sign_taxa <- merge(w_values, tax, by='otu.names')
  ancom_sign_taxa <- ancom_sign_taxa[,-1]
  ancom_results_phy_7[[i]] <- ancom_sign_taxa
}

for(i in 1:length(list_for_ancom_phy_1)) {
  otu_ancom <- otu_ancom_make(list_for_ancom_phy_1[[i]])
  comparison_test <- ANCOM.main(OTUdat = otu_ancom,
                                Vardat = metadata_ancom,
                                adjusted = FALSE,
                                repeated = F,
                                main.var = "Subsample",
                                adj.formula = NULL,
                                repeat.var=NULL,
                                longitudinal=FALSE,
                                random.formula=NULL,
                                multcorr=2,
                                sig=0.05,
                                prev.cut=0.90)
  w_values <- data.frame(comparison_test$W.taxa)
  tax <- data.frame(tax_table(physeq_1))
  tax <- tax %>% select(Kingdom, Phylum)
  tax$otu.names <- rownames(tax)
  ancom_sign_taxa <- merge(w_values, tax, by='otu.names')
  ancom_sign_taxa <- ancom_sign_taxa[,-1]
  ancom_results_phy_1[[i]] <- ancom_sign_taxa
}

#all have null results 
sink('phylum_ANCOM_qiime_nr.csv')
rec_phy_60120_nr <- ancom_results_phy_nr[[1]]
rec_phy_60120_nr <-rec_phy_60120_nr %>% filter(detected_0.9 == 'TRUE')
write.csv(rec_phy_60120_nr)
rec_phy_60188_nr <- ancom_results_phy_nr[[2]]
rec_phy_60188_nr <-rec_phy_60188_nr %>% filter(detected_0.9 == 'TRUE')
write.csv(rec_phy_60188_nr)
rec_phy_120188_nr <- ancom_results_phy_nr[[3]]
rec_phy_120188_nr <-rec_phy_120188_nr %>% filter(detected_0.9 == 'TRUE')
write.csv(rec_phy_120188_nr)
mou_phy_60120_nr <- ancom_results_phy_nr[[4]]
mou_phy_60120_nr <-mou_phy_60120_nr %>% filter(detected_0.9 == 'TRUE')
write.csv(mou_phy_60120_nr)
mou_phy_60188_nr <- ancom_results_phy_nr[[5]]
mou_phy_60188_nr <-mou_phy_60188_nr %>% filter(detected_0.9 == 'TRUE')
write.csv(mou_phy_60188_nr)
mou_phy_120188_nr <- ancom_results_phy_nr[[6]]
mou_phy_120188_nr <-mou_phy_120188_nr %>% filter(detected_0.9 == 'TRUE')
write.csv(mou_phy_120188_nr)
sink()

sink('phylum_ANCOM_qiime_7.csv')
rec_phy_60120_7 <- ancom_results_phy_7[[1]]
rec_phy_60120_7 <-rec_phy_60120_7 %>% filter(detected_0.9 == 'TRUE')
write.csv(rec_phy_60120_7)
rec_phy_60188_7 <- ancom_results_phy_7[[2]]
rec_phy_60188_7 <-rec_phy_60188_7 %>% filter(detected_0.9 == 'TRUE')
write.csv(rec_phy_60188_7)
rec_phy_120188_7 <- ancom_results_phy_7[[3]]
rec_phy_120188_7 <-rec_phy_120188_7 %>% filter(detected_0.9 == 'TRUE')
write.csv(rec_phy_120188_7)
mou_phy_60120_7 <- ancom_results_phy_7[[4]]
mou_phy_60120_7 <-mou_phy_60120_7 %>% filter(detected_0.9 == 'TRUE')
write.csv(mou_phy_60120_7)
mou_phy_60188_7 <- ancom_results_phy_7[[5]]
mou_phy_60188_7 <-mou_phy_60188_7 %>% filter(detected_0.9 == 'TRUE')
write.csv(mou_phy_60188_7)
mou_phy_120188_7 <- ancom_results_phy_7[[6]]
mou_phy_120188_7 <-mou_phy_120188_7 %>% filter(detected_0.9 == 'TRUE')
write.csv(mou_phy_120188_7)
sink()

sink('phylum_ANCOM_qiime_1.csv')
rec_phy_60120_1 <- ancom_results_phy_1[[1]]
rec_phy_60120_1 <-rec_phy_60120_1 %>% filter(detected_0.9 == 'TRUE')
write.csv(rec_phy_60120_1)
rec_phy_60188_1 <- ancom_results_phy_1[[2]]
rec_phy_60188_1 <-rec_phy_60188_1 %>% filter(detected_0.9 == 'TRUE')
write.csv(rec_phy_60188_1)
rec_phy_120188_1 <- ancom_results_phy_1[[3]]
rec_phy_120188_1 <-rec_phy_120188_1 %>% filter(detected_0.9 == 'TRUE')
write.csv(rec_phy_120188_1)
mou_phy_60120_1 <- ancom_results_phy_1[[4]]
mou_phy_60120_1 <-mou_phy_60120_1 %>% filter(detected_0.9 == 'TRUE')
write.csv(mou_phy_60120_1)
mou_phy_60188_1 <- ancom_results_phy_1[[5]]
mou_phy_60188_1 <-mou_phy_60188_1 %>% filter(detected_0.9 == 'TRUE')
write.csv(mou_phy_60188_1)
mou_phy_120188_1 <- ancom_results_phy_1[[6]]
mou_phy_120188_1 <-mou_phy_120188_1 %>% filter(detected_0.9 == 'TRUE')
write.csv(mou_phy_120188_1)
sink()

for(i in 1:length(list_for_ancom_fam_nr)) {
  otu_ancom <- otu_ancom_make(list_for_ancom_fam_nr[[i]])
  comparison_test <- ANCOM.main(OTUdat = otu_ancom,
                                Vardat = metadata_ancom,
                                adjusted = FALSE,
                                repeated = F,
                                main.var = "Subsample",
                                adj.formula = NULL,
                                repeat.var=NULL,
                                longitudinal=FALSE,
                                random.formula=NULL,
                                multcorr=2,
                                sig=0.05,
                                prev.cut=0.90)
  w_values <- data.frame(comparison_test$W.taxa)
  tax <- data.frame(tax_table(physeq_nr))
  tax <- tax %>% select(Kingdom, Phylum, Class, Order, Family)
  tax$otu.names <- rownames(tax)
  ancom_sign_taxa <- merge(w_values, tax, by='otu.names')
  ancom_sign_taxa <- ancom_sign_taxa[,-1]
  ancom_results_fam_nr[[i]] <- ancom_sign_taxa
}

for(i in 1:length(list_for_ancom_fam_7)) {
  otu_ancom <- otu_ancom_make(list_for_ancom_fam_7[[i]])
  comparison_test <- ANCOM.main(OTUdat = otu_ancom,
                                Vardat = metadata_ancom,
                                adjusted = FALSE,
                                repeated = F,
                                main.var = "Subsample",
                                adj.formula = NULL,
                                repeat.var=NULL,
                                longitudinal=FALSE,
                                random.formula=NULL,
                                multcorr=2,
                                sig=0.05,
                                prev.cut=0.90)
  w_values <- data.frame(comparison_test$W.taxa)
  tax <- data.frame(tax_table(physeq_7))
  tax <- tax %>% select(Kingdom, Phylum, Class, Order, Family)
  tax$otu.names <- rownames(tax)
  ancom_sign_taxa <- merge(w_values, tax, by='otu.names')
  ancom_sign_taxa <- ancom_sign_taxa[,-1]
  ancom_results_fam_7[[i]] <- ancom_sign_taxa
}

for(i in 1:length(list_for_ancom_fam_1)) {
  otu_ancom <- otu_ancom_make(list_for_ancom_fam_1[[i]])
  comparison_test <- ANCOM.main(OTUdat = otu_ancom,
                                Vardat = metadata_ancom,
                                adjusted = FALSE,
                                repeated = F,
                                main.var = "Subsample",
                                adj.formula = NULL,
                                repeat.var=NULL,
                                longitudinal=FALSE,
                                random.formula=NULL,
                                multcorr=2,
                                sig=0.05,
                                prev.cut=0.90)
  w_values <- data.frame(comparison_test$W.taxa)
  tax <- data.frame(tax_table(physeq_1))
  tax <- tax %>% select(Kingdom, Phylum, Class, Order, Family)
  tax$otu.names <- rownames(tax)
  ancom_sign_taxa <- merge(w_values, tax, by='otu.names')
  ancom_sign_taxa <- ancom_sign_taxa[,-1]
  ancom_results_fam_1[[i]] <- ancom_sign_taxa
}

#also null results 

sink('family_ANCOM_qiime_nr.csv')
rec_fam_60120_nr <- ancom_results_fam_nr[[1]]
rec_fam_60120_nr <-rec_fam_60120_nr %>% filter(detected_0.9 == 'TRUE')
write.csv(rec_fam_60120_nr)
rec_fam_60188_nr <- ancom_results_fam_nr[[2]]
rec_fam_60188_nr <-rec_fam_60188_nr %>% filter(detected_0.9 == 'TRUE')
write.csv(rec_fam_60188_nr)
rec_fam_120188_nr <- ancom_results_fam_nr[[3]]
rec_fam_120188_nr <-rec_fam_120188_nr %>% filter(detected_0.9 == 'TRUE')
write.csv(rec_fam_120188_nr)
mou_fam_60120_nr <- ancom_results_fam_nr[[4]]
mou_fam_60120_nr <-mou_fam_60120_nr %>% filter(detected_0.9 == 'TRUE')
write.csv(mou_fam_60120_nr)
mou_fam_60188_nr <- ancom_results_fam_nr[[5]]
mou_fam_60188_nr <-mou_fam_60188_nr %>% filter(detected_0.9 == 'TRUE')
write.csv(mou_fam_60188_nr)
mou_fam_120188_nr <- ancom_results_fam_nr[[6]]
mou_fam_120188_nr <-mou_fam_120188_nr %>% filter(detected_0.9 == 'TRUE')
write.csv(mou_fam_120188_nr)
sink()

sink('family_ANCOM_qiime_7.csv')
rec_fam_60120_7 <- ancom_results_fam_7[[1]]
rec_fam_60120_7 <-rec_fam_60120_7 %>% filter(detected_0.9 == 'TRUE')
write.csv(rec_fam_60120_7)
rec_fam_60188_7 <- ancom_results_fam_7[[2]]
rec_fam_60188_7 <-rec_fam_60188_7 %>% filter(detected_0.9 == 'TRUE')
write.csv(rec_fam_60188_7)
rec_fam_120188_7 <- ancom_results_fam_7[[3]]
rec_fam_120188_7 <-rec_fam_120188_7 %>% filter(detected_0.9 == 'TRUE')
write.csv(rec_fam_120188_7)
mou_fam_60120_7 <- ancom_results_fam_7[[4]]
mou_fam_60120_7 <-mou_fam_60120_7 %>% filter(detected_0.9 == 'TRUE')
write.csv(mou_fam_60120_7)
mou_fam_60188_7 <- ancom_results_fam_7[[5]]
mou_fam_60188_7 <-mou_fam_60188_7 %>% filter(detected_0.9 == 'TRUE')
write.csv(mou_fam_60188_7)
mou_fam_120188_7 <- ancom_results_fam_7[[6]]
mou_fam_120188_7 <-mou_fam_120188_7 %>% filter(detected_0.9 == 'TRUE')
write.csv(mou_fam_120188_7)
sink()

sink('family_ANCOM_qiime_1.csv')
rec_fam_60120_1 <- ancom_results_fam_1[[1]]
rec_fam_60120_1 <-rec_fam_60120_1 %>% filter(detected_0.9 == 'TRUE')
write.csv(rec_fam_60120_1)
rec_fam_60188_1 <- ancom_results_fam_1[[2]]
rec_fam_60188_1 <-rec_fam_60188_1 %>% filter(detected_0.9 == 'TRUE')
write.csv(rec_fam_60188_1)
rec_fam_120188_1 <- ancom_results_fam_1[[3]]
rec_fam_120188_1 <-rec_fam_120188_1 %>% filter(detected_0.9 == 'TRUE')
write.csv(rec_fam_120188_1)
mou_fam_60120_1 <- ancom_results_fam_1[[4]]
mou_fam_60120_1 <-mou_fam_60120_1 %>% filter(detected_0.9 == 'TRUE')
write.csv(mou_fam_60120_1)
mou_fam_60188_1 <- ancom_results_fam_1[[5]]
mou_fam_60188_1 <-mou_fam_60188_1 %>% filter(detected_0.9 == 'TRUE')
write.csv(mou_fam_60188_1)
mou_fam_120188_1 <- ancom_results_fam_1[[6]]
mou_fam_120188_1 <-mou_fam_120188_1 %>% filter(detected_0.9 == 'TRUE')
write.csv(mou_fam_120188_1)
sink()
