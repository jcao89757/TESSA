Tessa<-function(e,cdr3,t,hyper_priors,max_iter,sample_id=NULL,save=NULL,b=NULL,seed_num=123)
{
  # initialization
  cat('\nInitialization\n')
  options(warn=2)
  if ((!is.null(save)) && (!file.exists(save))) {dir.create(save)}
  if (!is.null(sample_id)){
    if(length(sample_id)!=ncol(t)){
      print('Unmatched sample IDs!')
      break
    }
  }
  if(!is.null(b)){
    preset_b=TRUE
  }else{
    preset_b=FALSE
  }
  print(seed_num)
  initialized=initialize(t,cdr3,e,hyper_priors,sample_id,b,seed_num=seed_num)
  t=initialized$t;meta=initialized$meta;a=initialized$a;b=initialized$b;sigma=initialized$sigma;
  K=initialized$K;ak=initialized$ak;phi=initialized$phi;t0=initialized$t0;
  master_dist_e=initialized$master_dist_e;de=initialized$de;dt=initialized$dt
  
  lambda=hyper_priors$lambda;xi=hyper_priors$xi;g=hyper_priors$g
  tau=hyper_priors$tau;u=hyper_priors$u;v=hyper_priors$v;
  initialize_cluster_factor=hyper_priors$initialize_cluster_factor
  
  # some intermediate variables for computational efficiencies
  meta_dedup=meta[!duplicated(meta$group_ID),] # meta: cell level; meta_dedup: group level
  rownames(meta_dedup)=meta_dedup$group_ID
  meta_dedup=meta_dedup[colnames(t),] 
  updated_recent=c() # for recording acceptance rates of b
  mean_t=rowMeans(t)
  
  # MCMC
  for (iter in 1:max_iter)
  {
    # report
    cat(paste('\nIteration round:',iter,"\n"))
    cat(paste("  # clusters:",K,"\n"))
    cluster_rate=1-sum(sapply(dt,function(x) length(x)==1))/length(unlist(dt))
    cat(paste("  Clustering rate:",round(cluster_rate,d=3),"\n"))
    
    # Dirichlet process
    results=DP(meta_dedup,meta,t0,dt,de,ak,phi,t,lambda,g,sigma,b,master_dist_e,K,
      a,xi,mean_t,sample_id)
    t0=results$t0;phi=results$phi;de=results$de;dt=results$dt;ak=results$ak;
    c=results$c;K=results$K;meta_dedup$cluster_number=results$cluster_number
    
    # other parameters
    t0=update_t0(t,meta_dedup,K,lambda,b,phi,t0)
    ak=update_ak(K,dt,de,sigma,a,g,phi,ak)
    regression_loss=sapply(1:K,function(k) sum((de[[k]]-ak[k]*dt[[k]])^2))
    sigma=update_sigma(u,v,K,phi,g,ak,a,de,dt,regression_loss)
    a=update_a(tau,K,g,sigma,ak)
    b_result=update_b(u,v,phi,t,t0,meta_dedup,K,dt,de,ak,sigma,regression_loss,b,preset_b)
    b=b_result$b;dt=b_result$dt
    updated_recent=c(b_result$updated,updated_recent)
    if (length(updated_recent)>200) {updated_recent=updated_recent[1:200]}
    cat(paste("  Recent b acceptance rate:",round(sum(updated_recent)/length(updated_recent),d=2),"\n"))
  }
  
  # wrap up
  options(warn=0)
  meta$cluster_number=meta_dedup[meta$group_ID,"cluster_number"]
  tessa_results=list(b=b,meta=meta,master_dist_e=master_dist_e,a=a,ak=ak,sigma=sigma,dt=dt,de=de,
    t=t,meta_dedup=meta_dedup,phi=phi,K=K)
  save(tessa_results,file=paste(save,"/tessa_final.RData",sep=""))
  write.csv(meta,file=paste(save,"/result_meta.csv",sep=""),quote=F,row.names=FALSE)
  return(tessa_results)
}
