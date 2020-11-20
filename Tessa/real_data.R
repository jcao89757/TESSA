args=commandArgs(trailingOnly = TRUE)
exp_file=args[2]
contigs_file=args[3]
cdr3_file=args[4]
save=args[5]
is_sampleCluster=as.logical(args[6])
fixed_b=args[7]
xi=1e+25
g=0.001
initialize_cluster_factor=6
source(paste(args[1],'update.R',sep='/'))
source(paste(args[1],'initialization.R',sep='/'))
source(paste(args[1],'MCMC_control.R',sep='/'))
source(paste(args[1],'utility.R',sep='/'))
source(paste(args[1],'post_analysis.R',sep='/'))
library(MASS)
library(LaplacesDemon)
library(Rtsne)
# the users need to provide these data: 
# The columns/rows/length of them should be matched up wherever applicable
# exp_file: expression data, cells on columns, and genes on rows.
#   e can be constructed by PCA or t-SNE, the first row is the first PC, the second row is the second PC, etc
# contigs_file: encoded CDR3 values, cells on columns, and embeddings on rows。
# cdr3: character vectors of CDR3 sequences
# save: a file dir to store tessa results
# (optional) sample_id: a column vector of sample categories. If is_sampleCluster=TRUE, users must provide an additional
#   column next to the cdr3 column.
# (optional) fixed_b: a vector of pre-defined b. The vector must be numerical and has the length of TCR embeddings.
exp_data=read.csv(exp_file,row.names=1,stringsAsFactors=F)
n=ncol(exp_data)
tmp=apply(exp_data,1,sd)
e=t(Rtsne(t(exp_data[tmp>quantile(tmp,0.9),1:n]),dims=3)$Y) # Run TSNE
colnames(e)=colnames(exp_data)[1:n]
contigs_encoded=read.csv(contigs_file,stringsAsFactors=F)
t=t(contigs_encoded[1:n,-1])
meta=read.csv(cdr3_file,header=TRUE,stringsAsFactors=F)
cdr3=meta$cdr3
if(is_sampleCluster){
  sample_id=meta$sample
}else{
  sample_id=NULL
}
if(fixed_b!='NA'){
  b=read.csv(fixed_b,header=TRUE,stringsAsFactors=F)$b
}else{
  b=NULL
}
# the users need to provide these parameters, here are the suggested values
hyper_priors=list(lambda=mean(apply(t,1,var)),
  xi=xi,g=g,tau=100,u=0.1,v=0.1,initialize_cluster_factor=initialize_cluster_factor)
max_iter=1000
#save="~/projects/scTCR/data/Tessa_save"
# Tessa
tessa_results=Tessa(e,cdr3,t,hyper_priors,max_iter,sample_id,save,b)
plot_tessaClsuters(tessa_results,save)
