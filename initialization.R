# find center of each cluster
find_center<-function(t_cluster,b)
{
  names(which.min(colSums((t_cluster-rowMeans(t_cluster))^2/b)))[1]
}  

# to initialize the assignment of clusters. A good initialization will save a lot of MCMC cycles
initialize_cluster<-function(meta,t,factor=4,sample_id_dedup) # larger factor, less number of clusters
{
  t0=t[,meta$group_ID]
  if(is.null(sample_id_dedup)){
    k0=round(dim(t0)[2]/factor)
    if(k0<1){k0=1}
    cluster_n=cutree(hclust(dist(t(t0))),k=k0)
  }else{
    sample_id=sample_id_dedup[meta$group_ID]
    cluster_n=rep(0,ncol(t0))
    names(cluster_n)=meta$group_ID
    for(s in 1:length(unique(sample_id))){
      sample_id_tmp=unique(sample_id)[s]
      t0_tmp=t0[,sample_id==sample_id_tmp]
      k0=round(dim(t0_tmp)[2]/(factor*length(unique(sample_id))))
      if(k0<1){k0=1}
      cluster_n_tmp=cutree(hclust(dist(t(t0_tmp))),k=k0)
      cluster_n_tmp=cluster_n_tmp+max(cluster_n)
      cluster_n[sample_id==sample_id_tmp]=cluster_n_tmp
    }
  }
  cluster=rep("",length(cluster_n))
  b_dummy=rep(1,dim(t)[1])
  
  for (i in unique(cluster_n))
    {cluster[cluster_n==i]=find_center(t0[,cluster_n==i,drop=F],b_dummy)}
  cluster
}

# initialization
initialize<-function(t,cdr3,e,hyper_priors,sample_id,b,seed_num=123)
{
  # setting up
  set.seed(as.numeric(seed_num))
  lambda=hyper_priors$lambda;xi=hyper_priors$xi;g=hyper_priors$g
  tau=hyper_priors$tau;u=hyper_priors$u;v=hyper_priors$v;
  initialize_cluster_factor=hyper_priors$initialize_cluster_factor
  
  # preprocess t
  if(!is.null(sample_id)){
    cdr3=paste(cdr3,sample_id,sep=';')
    sample_id_dedup=sample_id[!duplicated(cdr3)]
    names(sample_id_dedup)=cdr3[!duplicated(cdr3)]
  }else{
    sample_id_dedup=NULL
  }
  t=t[,!duplicated(cdr3),drop=F]
  colnames(t)=cdr3[!duplicated(cdr3)]

  # meta object
  # group ID is just the CDR3 sequence, cluster ID is the CDR3 sequence of the centers
  meta=data.frame(barcode=colnames(e),group_ID=cdr3,cluster_number=NA,stringsAsFactors = F)
  meta$cluster_number=initialize_cluster(meta,t,initialize_cluster_factor,sample_id_dedup)
  
  ## initialize random/placeholder variables
  # simple ones
  # the initialization of b is tricky
  # note this is not sampling from its distribution, but it is ok for initialization
  if(is.null(b)){
    b=apply(t,1,var)/10 
  }
  #b[]=mean(b)
  K=length(unique(meta$cluster_number))
  
  # phi
  phi0=aggregate(meta$group_ID,by=list(meta$cluster_number),function(x) length(unique(x)))
  phi=phi0[,2]
  names(phi)=phi0[,1]
  
  # t0
  # this is not exactly "right", but not wrong, either
  t00=aggregate(t(t[,meta$group_ID]),by=list(meta$cluster_number),mean) 
  t0=sapply(1:dim(t00)[1],function(i) as.numeric(unlist(t00[i,-1])),simplify=F)
  names(t0)=t00$Group.1
  t0=t0[names(phi)]
  
  # de, dt
  tmp=as.matrix(dist(t(e)))
  tmp=aggregate(tmp,by=list(meta$group_ID),mean)
  rownames(tmp)=tmp[,1]
  tmp=as.matrix(tmp[,-1])
  tmp=aggregate(t(tmp),by=list(meta$group_ID),mean)
  rownames(tmp)=tmp[,1]
  master_dist_e=as.matrix(tmp[,-1])
  master_dist_e=master_dist_e[colnames(t),colnames(t)]
  
  dt=de=list()
  coefs=c()
  
  for(k in 1:K)
  {
    c=names(phi)[k]
    group=unique(meta$group_ID[meta$cluster_number==c])
    de[[c]]=named_c(NULL,master_dist_e[names(phi)[k],group],group)
    dt[[c]]=named_c(NULL,colSums((t[,group,drop=F]-t[,c])^2/b/2),group)
    coefs=c(coefs,coef(lm(de[[c]]~dt[[c]]))[2]) # for estimating good values of a and ak
  }
 
  # sigma, this is also tricky
  sigma=mean(sapply(1:K,function(k) sd(de[[k]])),na.rm=T)*5
  if (is.na(sigma)) {sigma=1}
  
  # a, ak
  # note this is not sampling from its distribution
  coefs=coefs[which(coefs>0)]
  a=ifelse(is.na(mean(coefs)),1,mean(coefs))
  ak=rnorm(K,a,ifelse(is.na(sd(coefs)),1,sd(coefs))) 
  names(ak)=names(phi)
  
  # return 
  return(list(t=t,meta=meta,a=a,b=b,sigma=sigma,K=K,
    ak=ak,phi=phi,t0=t0,master_dist_e=master_dist_e,de=de,dt=dt))
}
