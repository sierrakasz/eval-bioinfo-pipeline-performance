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

tax_group <- read.csv("tax_group.csv")
tax_group <- tax_group %>% group_by(Kingdom, Phylum, Class, Order, Family, Genus) %>% summarise_all(funs(sum))
tax_group <- tax_group %>%  ungroup()

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

#import metadata and combine
metadata=(read.csv("Metadata/HPMMMeta_r_merge.csv",header=TRUE))
sampdat=sample_data(metadata)
sample_names(sampdat)=metadata$SampleID
physeq=merge_phyloseq(physeq_all, sampdat)
physeq

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

physeq_phy <- tax_glom(physeq, taxrank = 'Phylum')
physeq_fam <- tax_glom(physeq, taxrank = 'Family')

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

list_for_ancom_phy <- comparisons_for_ancom(physeq_phy)
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

ancom_results_phy <- list()
ancom_results_fam <- list()

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

sink('phylum_ANCOM_threepipelines.csv')
rec_phy_GM <- ancom_results_phy[[1]]
rec_phy_GM <-rec_phy_GM %>% filter(detected_0.9 == 'TRUE')
write.csv(rec_phy_GM)
rec_phy_QM <- ancom_results_phy[[2]]
rec_phy_QM <-rec_phy_QM %>% filter(detected_0.9 == 'TRUE')
write.csv(rec_phy_QM)
rec_phy_GQ <- ancom_results_phy[[3]]
rec_phy_GQ <-rec_phy_GQ %>% filter(detected_0.9 == 'TRUE')
write.csv(rec_phy_GQ)
mou_phy_GM <- ancom_results_phy[[4]]
mou_phy_GM <-mou_phy_GM %>% filter(detected_0.9 == 'TRUE')
write.csv(mou_phy_GM)
mou_phy_QM <- ancom_results_phy[[5]]
mou_phy_QM <-mou_phy_QM %>% filter(detected_0.9 == 'TRUE')
write.csv(mou_phy_QM)
mou_phy_GQ <- ancom_results_phy[[6]]
mou_phy_GQ <-mou_phy_GQ %>% filter(detected_0.9 == 'TRUE')
write.csv(mou_phy_GQ)
sink()

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

sink('family_ANCOM_threepipelines.csv')
rec_fam_GM <- ancom_results_fam[[1]]
rec_fam_GM <-rec_fam_GM %>% filter(detected_0.9 == 'TRUE')
write.csv(rec_fam_GM)
rec_fam_QM <- ancom_results_fam[[2]]
rec_fam_QM <-rec_fam_QM %>% filter(detected_0.9 == 'TRUE')
write.csv(rec_fam_QM)
rec_fam_GQ <- ancom_results_fam[[3]]
rec_fam_GQ <-rec_fam_GQ %>% filter(detected_0.9 == 'TRUE')
write.csv(rec_fam_GQ)
mou_fam_GM <- ancom_results_fam[[4]]
mou_fam_GM <-mou_fam_GM %>% filter(detected_0.9 == 'TRUE')
write.csv(mou_fam_GM)
mou_fam_QM <- ancom_results_fam[[5]]
mou_fam_QM <-mou_fam_QM %>% filter(detected_0.9 == 'TRUE')
write.csv(mou_fam_QM)
mou_fam_GQ <- ancom_results_fam[[6]]
mou_fam_GQ <-mou_fam_GQ %>% filter(detected_0.9 == 'TRUE')
write.csv(mou_fam_GQ)
sink()

MG_samp <- grep("G_", colnames(test2))
test3 <- test2[, c("Kingdom", "Phylum", "Class", "Order", "Family")]
test4 <- test2[,MG_samp]
test5 <- cbind(test3, test4)
sum(test5[1, 6:length(colnames(test5))])
blank <- c()

for(i in 6:length(colnames(test5))) {
  if(test5[1,i] > 0) {
    blank <- c(blank, 1)
  }
}

for(i in 1:length(row(test5))) {
  for(j in 6:length(colnames(test5))) {
    if(test5[i,j] > 0) {
      print(1)
    }
  }
}

rowSums(test5 > 0)
test5$col1 <- apply(test5, 1, function(i) sum(i > 0))
length(colnames(test5))

test5[1,]


#Figures
ANCOM_figures_setup <- function(ancom, physeq, results) {
  otu <- data.frame(otu_table(ancom))
  tax_table <- data.frame(tax_table(physeq))
  otu$names <- rownames(otu)
  tax_table$names <- rownames(tax_table)
  tax_table <- tax_table %>% select(-contains('Genus'))
  otu_tax <- merge(otu,tax_table)
  otu_tax <- otu_tax %>% select(-contains('Genus'))
  merging <- merge(otu_tax, results, by = c("Kingdom", "Phylum", "Class", "Order", "Family"))
  return(merging)
}

testing <- ANCOM_figures(list_for_ancom_fam[[1]], physeq, rec_fam_GM)
write.csv(testing, 'rec_fam_GM.csv')
testing <- ANCOM_figures(list_for_ancom_fam[[2]], physeq, rec_fam_QM)
write.csv(testing, 'rec_fam_QM.csv')
testing <- ANCOM_figures(list_for_ancom_fam[[3]], physeq, rec_fam_GQ)
write.csv(testing, 'rec_fam_GQ.csv')
testing <- ANCOM_figures(list_for_ancom_fam[[4]], physeq, mou_fam_GM)
write.csv(testing, 'mou_fam_GM.csv')
testing <- ANCOM_figures(list_for_ancom_fam[[5]], physeq, mou_fam_QM)
write.csv(testing, 'mou_fam_QM.csv')
testing <- ANCOM_figures(list_for_ancom_fam[[6]], physeq, mou_fam_GQ)
write.csv(testing, 'mou_fam_GQ.csv')

#Figures
#import form excel, just got too frusterated 
#just going to do table instead of figures


ANCOM_figures_numbers <- function(test_df) {
  MG_samp <- grep("G_", colnames(test_df))
  mot_samp <- grep("M_", colnames(test_df))
  qim_samp <- grep("Q_", colnames(test_df))
  MG_tax <- test_df[, c("Kingdom", "Phylum", "Class", "Order", "Family")]
  MGsamples <- test_df[,MG_samp]
  MG_df <- cbind(MG_tax, MGsamples)
  mg_vec <- c()
  for(i in 1:length(row(MG_df))) {
    mg_vec <- c(mg_vec, sum(MG_df[i, 6:length(colnames(MG_df))]))
  }
  print(mg_vec[1:7])
  samp_df <- data.frame()
  for(i in 1:length(row(MG_df))) {
    samp_df$i <- apply(MG_df, 1, function(i) sum(i > 0))
  }
}

first <- ANCOM_figures_numbers(testing)

#write a function that does the following:
  #merges results from ancom with tax and otu information
  #for each pipeline, counts the number of samples and sequences taxa is in
  #outputs taxa, W value, sequences and sample numbers into dataframe 


