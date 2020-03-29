######  set up simulation  ############
setwd('~/projects/scTCR/code/Tessa/')
source('update.R')
source('initialization.R')
source('MCMC_control.R')
source('utility.R')
source('post_analysis.R')
library(MASS)
library(LaplacesDemon)
library(mclust)

simuData=function(K,dist_between_cluster_factor,
                  cluster_size_max,e_dispersion_factor,shape=5,scale=2,seed_num){
  set.seed(seed_num)
  dim_t=30
  b=rinvgamma(dim_t, shape=shape, scale=scale)
  #Change the shape and the scale resulted in poor ARIs.
  t0_new=list()
  t=c()
  phi=c()
  cdr3=c()
  e=c()
  for (k in 1:K) {
    #real data size: ~100 groups (10 clusters) to ~10k groups (1k clusters)
    phi=c(phi,sample(1:cluster_size_max,1))
    #real data mean cluster size: ~ 10 groups
    names=paste("cluster",k,"cell",1:phi[k],sep="")
    cdr3=c(cdr3,names)
    t0_new[[k]]=rnorm(n=dim_t,rep(0,dim_t),mean(sqrt(b))*dist_between_cluster_factor)#random var
    #increase dist_between_cluster_factor increased the performance
    tmp=sapply(1:phi[k],function(i) rnorm(n=dim_t,t0_new[[k]],sqrt(b)))#t0=t/sqrt(b), assume t0 is known
    colnames(tmp)=names
    t=cbind(t,tmp)
    e_center=rnorm(n=3,rep(0,3),10)
    dist_t=colSums((tmp-t0_new[[k]])^2/b)
    tmp[]=sapply(1:phi[k],function(i) rnorm(n=dim_t,e_center,e_dispersion_factor*dist_t[i]^0.5))#linear correlation
    #increase e_dispersion_factor had minimal effect. 
    e=cbind(e,tmp)#e after PCA transformation, 3 dimentions
  }
  return(list(t=t,
              t0_new=t0_new,
              e=e,
              cdr3=cdr3,
              phi=phi))
}
# count_eudist_cor=function(exp_data,tessa_results){
#   dist_exp=as.matrix(dist(t(exp_data)))
#   dist_exp[lower.tri(dist_exp,diag=T)]=NA
#   dist_exp=as.vector(dist_exp)
#   dist_exp=dist_exp[!is.na(dist_exp)]
#   dist_exp=log(dist_exp+1)
#   contigs_new=t(tessa_results$t/sqrt(tessa_results$b))
#   dist_cdr=as.matrix(dist(contigs_new))
#   dist_cdr[lower.tri(dist_cdr,diag=T)]=NA
#   dist_cdr=as.vector(dist_cdr)
#   dist_cdr=dist_cdr[!is.na(dist_cdr)]
#   group_size=floor(length(dist_cdr)/170)
#   if(group_size<5){
#     group_size=floor(length(dist_cdr)/50)
#   }
#   dist_exp=dist_exp[order(dist_cdr,decreasing = F)]
#   dist_cdr=dist_cdr[order(dist_cdr,decreasing = F)]
#   aggrlist=rep(c(1:floor(length(dist_cdr)/group_size)),each=group_size)
#   aggrlist=c(aggrlist,rep(aggrlist[length(aggrlist)],(length(dist_cdr)-length(aggrlist))))
#   dist_exp=aggregate(dist_exp,by=list(aggrlist),mean)
#   dist_exp=dist_exp[,2]
#   dist_cdr=aggregate(dist_cdr,by=list(aggrlist),mean)
#   dist_cdr=dist_cdr[,2]
#   cor=cor.test(dist_exp,dist_cdr)$estimate
#   return(cor)
# }
###########  Tessa  ###############
args = commandArgs(trailingOnly=TRUE)
seed_num=as.numeric(args[1])
#K_list=c(1:10)*100
#dist_between_cluster_list=c(0.5,1,1.5,2,2.5,3)
#e_dispersion_list=c(0.01,0.1,0.5,1,1.5,2)
K=as.numeric(args[2])
cluster_size_max=as.numeric(args[3])
tdist=as.numeric(args[4])
edist=as.numeric(args[5])
save=paste("/work/SCCC/s421955/tessa_tmp/simu_results",seed_num,'_',K,'_',cluster_size_max,'_',tdist,'_',edist,sep='')
max_iter=1000
if(!file.exists(save)){dir.create(save)}
testdata=simuData(K,dist_between_cluster_factor = tdist,cluster_size_max=cluster_size_max,e_dispersion_factor=edist,shape=5,scale=2,seed_num=123)
# Tessa
hyper_priors=list(lambda=mean(apply(testdata$t,1,var)),
                  xi=1e25,g=0.001,tau=100,u=0.1,v=0.1,initialize_cluster_factor=6)#same as final conditions
tessa_results=Tessa(testdata$e,testdata$cdr3,testdata$t,hyper_priors,max_iter,save = save,seed_num=seed_num)
#############  validate Tessa results  ################
# check cluster assignment
ari=adjustedRandIndex(sub("cell.*","",testdata$cdr3),tessa_results$meta$cluster_number)
save(ari,testdata,file = paste(save,'/simu_results.RData',sep=''))




