# set up environment
# update 02032020: if apply bs from other datasets, fixed_b=a 30-digit numerical vector in *.RData. 
# the vector should be named as b.
args=commandArgs(trailingOnly = T)
input_rdata=args[1];xi=as.numeric(args[2]);g=as.numeric(args[3]);initialize_cluster_factor=as.numeric(args[4])
save=args[5];is_sampleCluster=args[6];fixed_b=args[7]
#final params:
#xi=1e+25;g=0.001;initialize_cluster_factor=6
#test:
#input_rdata='~/projects/scTCR/data/matched_exp_TCR/CD8_allinone/Natmed_bcc.RData'
#save='~/temp/testbcc';is_sampleCluster=TRUE;fixed_b='~/temp/testbcc_b.RData'
set.seed(123)
setwd("~/projects/scTCR/code/Tessa/")
#load('~/projects/scTCR/data/matched_exp_TCR/CD8_allinone/Natmed_bcc.RData')
load(input_rdata)
source('update.R')
source('initialization.R')
source('MCMC_control.R')
source('utility.R')
source('post_analysis.R')
library(MASS)
library(LaplacesDemon)
# the users need to provide these data: 
# The columns/rows/length of them should be matched up wherever applicable
# PCs of expression data, cells on columns, and PCs on rows (e). 
# e can be constructed by PCA or t-SNE, the first row is the first PC, the second row is the second PC, etc
# CDR3 sequences, a character vector (cdr3)
# encoded CDR3 values, cells on columns, and embeddings on rows (t)
# (optional) a vector of sample categories (sample_id) 
print(input_rdata)
print(xi)
print(g)
print(initialize_cluster_factor)
n=ncol(exp_data)
library(Rtsne)
tmp=apply(exp_data,1,sd)
e=t(Rtsne(t(exp_data[tmp>quantile(tmp,0.9),1:n]),dims=3)$Y) # Run TSNE
colnames(e)=colnames(exp_data)[1:n]
cdr3=contigs_encoded[1:n,1]
t=t(contigs_encoded[1:n,-1])
#08132019 test: sample_id
if(is.na(is_sampleCluster)){is_sampleCluster=FALSE}
if(is_sampleCluster){
  sample_id=contigs_refined$sample[1:n]
}else{
  sample_id=NULL
}
#02032020 test fixed b
if(!is.na(fixed_b)){
  load(fixed_b)
}
# the users need to provide these parameters, here are the suggested values
hyper_priors=list(lambda=mean(apply(t,1,var)),
  xi=xi,g=g,tau=100,u=0.1,v=0.1,initialize_cluster_factor=initialize_cluster_factor)
max_iter=1000
#save="~/projects/scTCR/data/Tessa_save"
# Tessa
tessa_results=Tessa(e,cdr3,t,hyper_priors,max_iter,sample_id,save,b)
#folder="~/projects/scTCR/data/Tessa_results"
#plot_tessa(tessa_results,folder)
t_new=t
cluster_new=predict_tessa(tessa_results,t_new)
library(mclust)
adjustedRandIndex(cluster_new,tessa_results$meta$cluster_number)
