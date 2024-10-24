## ---------------------------
##
##
## Purpose of script: polyA+/non-polyA classification
##
## Author: Lucia Lorenzi
##
## Date Created: 2018-Nov
##
## Email: lucialorenzi90@gmail.com
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

options(scipen = 999) 
## ---------------------------
## load up the packages we will need:  (uncomment as required)
require(tidyverse)
require(data.table)
## ---------------------------
#Load data
#datasets polyA and TotalRNA
# select matching samples
#require(data.table)
setwd("/user/RNA_atlas/analyses/new_polyadenylation_after_masking_mono_exonic/new_polyadenylation_November_2018_htseq_counts")

##load subsampled tpm data
polyA=read.table("/user/RNA_atlas/analyses/HTSeq_count/correction_libsizes/polyA_htseq_subsampled_counts.txt",sep = '\t',header = T,check.names = F)
TotalRNA=read.table("/user/RNA_atlas/analyses/HTSeq_count/correction_libsizes/TotalRNA_htseq_subsampled_counts.txt",sep = '\t',header = T,check.names = F)

#load known polyA(+) and polyA(-) genes
polyadenylated_genes=as.character(read.table("/user/RNA_atlas/analyses/new_polyadenylation_after_masking_mono_exonic/polyadenylated_genes.txt",header = F)[,1])
non_polyadenylated_genes=as.character(read.table("/user/RNA_atlas/analyses/new_polyadenylation_after_masking_mono_exonic/non_polyadenylated_genes.txt",header = F,sep = '\t')[,1])


###Note I will filter in all cases based on counts ####

samples=colnames(polyA)

count_cutoff=10

##Classification with htseq.counts

system("mkdir htseq.counts")
setwd("htseq.counts")
system("mkdir densities_with_cutoff")
system("mkdir accuracy_plots")
system("mkdir spec_vs_sens_plots")
system("mkdir densities")

htseq.counts_all_samples_call_mean_cutoff_accu=c()
htseq.counts_all_samples_ratios=c()
htseq.counts_all_samples_mean_cutoff_accu=c()
htseq.counts_all_samples_mean_cutoff_opt=c()
htseq.counts_all_samples_var_cutoff_accu=c()
htseq.counts_all_samples_var_cutoff_opt=c()
htseq.counts_all_samples_cutoffs_accu=c()
htseq.counts_all_samples_cutoffs_opt=c()
htseq.counts_all_samples_cutoffs_pROCR=c()
htseq.counts_all_samples_accuracies=list()
htseq.counts_all_samples_specificities=list()
htseq.counts_all_samples_sensitivities=list()
htseq.counts_all_samples_AUC=c()
htseq.counts_all_samples_indices=list()


set.seed(4)
j=0
library(pROC)
for(i in samples){
  tmp_data=data.frame(polyA=polyA[,i],TotalRNA=TotalRNA[,i])
  tmp_data$gene=rownames(polyA)
  tmp_data$ratio=tmp_data$polyA/tmp_data$TotalRNA
  
  #create tables
  tmp_data_call_mean_cutoff=matrix(data=NA,nrow = nrow(polyA),ncol=1)
  rownames(tmp_data_call_mean_cutoff)=rownames(polyA)
  colnames(tmp_data_call_mean_cutoff)=i
  
  tmp_data_ratios=tmp_data_call_mean_cutoff
  
  #fill in the call tables with values from -2 to 2
  #samples that have less than cutoff counts in TotalRNA are undetermined
  tmp_data_call_mean_cutoff[tmp_data$TotalRNA<count_cutoff]=0
  #samples that have 0 counts in polyA and more than cutoff in TotalRNA are non-polyadenylated
  tmp_data_call_mean_cutoff[tmp_data$polyA==0&tmp_data$TotalRNA>=count_cutoff]=-2
  
  
  #Separate data to classify
  tmp_data_to_classify=tmp_data[tmp_data$polyA!=0&tmp_data$TotalRNA>=count_cutoff,]
  rownames(tmp_data_to_classify)=tmp_data_to_classify$gene
  
  
  #Select the known genes to use in the classification
  #we will also require the polyadenylated genes to be expressed in 
  #at least cutoff counts in polyA libraries to be used as reference
  polyadenylated_genes_sample=tmp_data$gene[tmp_data$gene%in%polyadenylated_genes&tmp_data$polyA>=count_cutoff]
  
  #convert to log
  tmp_data_to_classify$polyA=log2(tmp_data_to_classify$polyA)
  tmp_data_to_classify$TotalRNA=log2(tmp_data_to_classify$TotalRNA)
  tmp_data_to_classify$ratio=log2(tmp_data_to_classify$ratio)
  
  #Include a plot for each sample here 
  
  pdf(paste0("densities/",i,"_density.pdf"))
  
  par(las=1)
  plot(density(tmp_data_to_classify$ratio[tmp_data_to_classify$gene%in%polyadenylated_genes_sample]),col="red",xlim=c(min(tmp_data_to_classify$ratio[tmp_data_to_classify$gene%in%non_polyadenylated_genes]),max(tmp_data_to_classify$ratio[tmp_data_to_classify$gene%in%polyadenylated_genes])),xlab="htseq.counts_ratio",main = i)
  points(density(tmp_data_to_classify$ratio[tmp_data_to_classify$gene%in%non_polyadenylated_genes]),col="blue",type = 's')
  
  dev.off()
  

  
  tmp_data_to_classify$class[tmp_data_to_classify$gene%in%polyadenylated_genes_sample]="polyadenylated"
  
  tmp_data_to_classify$class[tmp_data_to_classify$gene%in%non_polyadenylated_genes]="non_polyadenylated"
  #divide in known and unknown
  known=tmp_data_to_classify[tmp_data_to_classify$gene%in%c(polyadenylated_genes_sample,as.character(non_polyadenylated_genes)),]
  log.data=known[,c("ratio","class")]
  
  
  # classification 
  
  cutoffs_accu=c()
  cutoffs_opt=c()
  cutoffs_pROC=c()
  all_accuracy=c()
  all_TPRs=c()
  all_TNRs=c()
  all_AUCs=c()
  all_pROCs=list()
  all_indices=c()
  
  for(j in 1:100){
    training_data <- log.data
    training_data_polyA <- training_data[training_data$class=="polyadenylated", ]
    training_data_non_polyA <- training_data[training_data$class == "non_polyadenylated", ] 
    
    #get subsample to balance dataset
    indices <- sample(nrow(training_data_polyA), nrow(training_data_non_polyA))
    training_data_polyA <- training_data_polyA[indices, ]
    training_data <- rbind(training_data_polyA, training_data_non_polyA)
    accus=c()
    FPRs=c()
    TPRs=c()
    FNRs=c()
    TNRs=c()
    #cnf_matrices=c()
    for(k in seq(-5,2,0.01)){
      pred=rep(NA,nrow(training_data))
      pred[training_data$ratio<k]="non_polyadenylated"
      pred[training_data$ratio>=k]="polyadenylated"
      pred=factor(pred,levels=c("polyadenylated","non_polyadenylated"))
      
      res=table(pred,training_data$class)
      
      accu=(res[1,2]+res[2,1])/sum(res)
      FPR=res[1,1]/sum(res[,1])
      TPR=res[1,2]/sum(res[,2])
      FNR=res[2,2]/sum(res[,2])
      TNR=res[2,1]/sum(res[,1])
      
      
      accus=c(accus,accu)
      FPRs=c(FPRs,FPR)
      TPRs=c(TPRs,TPR)
      FNRs=c(FNRs,FNR)
      TNRs=c(TNRs,TNR)
      #cnf_matrices=rbind(cnf_matrices,as.vector(res))
      
    }
    all_indices=rbind(all_indices,indices)
    tmp_cutoff_accu=mean(seq(-5,2,0.01)[which(accus==max(accus))])
    cutoffs_accu=c(cutoffs_accu,tmp_cutoff_accu)
    d=(TPRs-1)^2 + (TNRs-1)^2
    tmp_cutoff_opt=mean(seq(-5,2,0.01)[which(d==min(d))])
    cutoffs_opt=c(cutoffs_opt,tmp_cutoff_opt)
    all_accuracy=rbind(all_accuracy,accus)
    all_TPRs=rbind(all_TPRs,TPRs)
    all_TNRs=rbind(all_TNRs,TNRs)
    
    pROC=roc(training_data$class,training_data$ratio)
    d2=(pROC$sensitivities-1)^2 +(pROC$specificities-1)^2
    pROC_cutoff=mean(pROC$thresholds[which(d2==min(d2))])
    cutoffs_pROC=c(cutoffs_pROC,pROC_cutoff)
    all_pROCs[[j]]=pROC
    all_AUCs=c(all_AUCs,as.numeric(pROC$auc))
    
  }
  
  ######Plot ROC#######
  pdf(paste0("spec_vs_sens_plots/",i,"_ROC_curves.pdf"))
  par(las=1)
  plot(all_pROCs[[1]])
  for(k in 2:100){
    plot(all_pROCs[[k]],add=T)
  }
  dev.off()
  
  ####Plot accuracies######
  pdf(paste0("accuracy_plots/",i,"_accu_curves.pdf"))
  par(las=1)
  plot(seq(-5,2,0.01),all_accuracy[1,],type = 's',ylab = "accuracy",xlab = "cutoff",main = i)
  for(k in 2:100){
    points(seq(-5,2,0.01),all_accuracy[k,],type = 's')
  }
  dev.off()
  
  var_cutoff_accu=var(cutoffs_accu)
  mean_cutoff_accu=mean(cutoffs_accu)
  var_cutoff_opt=var(cutoffs_opt)
  mean_cutoff_opt=mean(cutoffs_opt)
  mean_cutoff_pROC=mean(cutoffs_pROC)
  
  #classify with calculated cutoff
  tmp_data_call_mean_cutoff[rownames(tmp_data_call_mean_cutoff)%in%tmp_data_to_classify$gene[tmp_data_to_classify$ratio>=mean_cutoff_accu]]=1
  tmp_data_call_mean_cutoff[rownames(tmp_data_call_mean_cutoff)%in%tmp_data_to_classify$gene[tmp_data_to_classify$ratio<mean_cutoff_accu]]=-1
  
  tmp_data_ratios[,1]=tmp_data_to_classify$ratio[match(rownames(tmp_data_ratios),rownames(tmp_data_to_classify))]
  
  
  htseq.counts_all_samples_call_mean_cutoff_accu=cbind(htseq.counts_all_samples_call_mean_cutoff_accu,tmp_data_call_mean_cutoff)
  htseq.counts_all_samples_ratios=cbind(htseq.counts_all_samples_ratios,tmp_data_ratios)
  htseq.counts_all_samples_mean_cutoff_accu=c(htseq.counts_all_samples_mean_cutoff_accu,mean_cutoff_accu)
  htseq.counts_all_samples_var_cutoff_accu=c(htseq.counts_all_samples_var_cutoff_accu,var_cutoff_accu)
  htseq.counts_all_samples_mean_cutoff_opt=c(htseq.counts_all_samples_mean_cutoff_opt,mean_cutoff_opt)
  htseq.counts_all_samples_var_cutoff_opt=c(htseq.counts_all_samples_var_cutoff_opt,var_cutoff_opt)
  htseq.counts_all_samples_accuracies[[which(samples==i)]]=all_accuracy
  htseq.counts_all_samples_specificities[[which(samples==i)]]=all_TPRs
  htseq.counts_all_samples_sensitivities[[which(samples==i)]]=all_TNRs
  
  htseq.counts_all_samples_cutoffs_accu=rbind(htseq.counts_all_samples_cutoffs_accu,cutoffs_accu)
  htseq.counts_all_samples_cutoffs_opt=rbind(htseq.counts_all_samples_cutoffs_opt,cutoffs_opt)
  htseq.counts_all_samples_AUC=rbind(htseq.counts_all_samples_AUC,all_AUCs)
  htseq.counts_all_samples_cutoffs_pROCR=rbind(htseq.counts_all_samples_cutoffs_pROCR,cutoffs_pROC)
  htseq.counts_all_samples_indices[[which(samples==i)]]=all_indices
  
  pdf(paste0("densities_with_cutoffs/",i,"_density_with_cutoffs.pdf"))
  par(las=1)
  plot(density(tmp_data_to_classify$ratio[tmp_data_to_classify$gene%in%polyadenylated_genes]),xlim=c(min(tmp_data_to_classify$ratio[tmp_data_to_classify$gene%in%non_polyadenylated_genes])-1,max(tmp_data_to_classify$ratio[tmp_data_to_classify$gene%in%polyadenylated_genes])+1),main = i,col="red",xlab="htseq.counts_ratio")
  points(density(tmp_data_to_classify$ratio[tmp_data_to_classify$gene%in%non_polyadenylated_genes]),type = 's',col="blue")
  abline(v=mean_cutoff_accu)
  abline(v=mean_cutoff_opt,col="darkgreen")
  abline(v=mean_cutoff_pROC,col="blue")
  dev.off()
}

save.image("/user/RNA_atlas/analyses/new_polyadenylation_after_masking_mono_exonic/new_polyadenylation_November_2018_htseq_counts/polyA_classification_htseq.counts_November_2018.RData")

#################################################
##################################################

