update_t0=function(t,meta_dedup,K,lambda,b,phi,t0)
{
  tmp=aggregate(t(t),by=list(meta_dedup$cluster_number),sum)
  rownames(tmp)=tmp[,1]
  tmp=t(as.matrix(tmp[,-1]))
  tmp=tmp[,names(phi),drop=F]
  for (k in 1:K) # this may be further optimized
    {t0[[k]]=rnorm(dim(t)[1],tmp[,k]/(b/lambda+phi[k]),sqrt(1/(1/lambda+phi[k]/b)))}
  t0
}

update_ak=function(K,dt,de,sigma,a,g,phi,ak)
{
  AK=sapply(1:K,function(k) sum(dt[[k]]^2))+1/g
  BK=sapply(1:K,function(k) sum(dt[[k]]*de[[k]]))+a/g
  ak_new=rnorm(K,BK/AK,sigma/sqrt(AK))
  ak[phi>1]=ak_new[phi>1]
  ak
}

update_sigma=function(u,v,K,phi,g,ak,a,de,dt,regression_loss)
{
  C=u+1+(K+sum(phi))/2
  D=v+sum((ak-a)^2)/2/g+sum(regression_loss)/2
  sqrt(rinvgamma(1, shape=C-1, scale=D))
}

update_a=function(tau,K,g,sigma,ak)
{
  E=1/tau+K/(g*sigma^2)
  D=sum(ak)/(g*sigma^2)
  a=rtnorm(mu=D/E,sd=sqrt(1/E),a=0,b=Inf)
  a
}

update_b=function(u,v,phi,t,t0,meta_dedup,K,dt,de,ak,sigma,regression_loss,b,preset_b)
{
  if(!preset_b){
    # new beta
    MH_alpha=u+sum(phi)/2
    MH_beta=v+rowSums((t-unlist(t0[meta_dedup[colnames(t),"cluster_number"]]))^2)/2
    bnew=sapply(1:length(b),function(q) rinvgamma(1,shape=MH_alpha,scale=MH_beta[q]))
    
    # new dt
    dtnew=dt
    for(k in 1:K)
    {
      c=names(phi)[k]
      group=names(dt[[k]])
      dtnew[[k]]=colSums((t[,group,drop=F]-t[,c])^2/bnew/2)  
    }
    
    # F
    regression_loss_new=sapply(1:K,function(k) sum((de[[k]]-ak[k]*dtnew[[k]])^2))
    F=sum(regression_loss_new-regression_loss)/2/sigma^2
    
    updated=0
    if (runif(1,0,1)<min(1,exp(-F)))
    {
      updated=1
      b=bnew
      dt=dtnew
    }
  }else{
    b=b;dt=dt;updated=0
  }
  return(list(b=b,dt=dt,updated=updated))
}

calulate_DPprob<-function(phi,xi,group_ID,b,t0_new,ak_new,master_dist_e,sigma,t,ak,t0)
{
  t_group_ID=t[,group_ID]
  part1=c(phi,xi)
  part2=c(sapply(t0,function(x) sum(-(x-t_group_ID)^2/2/b)),
    sum(-(t0_new-t_group_ID)^2/2/b))
  de_test=c(master_dist_e[group_ID,names(phi)],master_dist_e[group_ID,group_ID])
  dt_test=c(colSums((t_group_ID-t[,names(phi),drop=F])^2/b/2),0)  
  part3=(de_test-dt_test*c(ak,ak_new))^2/(-2*sigma^2)
  tmp=part2+part3
  tmp=tmp-max(tmp)
  part1*exp(tmp)
}

DP<-function(meta_dedup,meta,t0,dt,de,ak,phi,t,lambda,g,sigma,b,master_dist_e,K,
  a,xi,mean_t,sample_id)
{
  if(!is.null(sample_id)){
    sample_id_dedup=sample_id[!duplicated(meta$group_ID)]
    names(sample_id_dedup)=meta_dedup$group_ID
  }else{
    sample_id_dedup=NULL
  }
  for (group_ID in 1:dim(meta_dedup)[1])
  {
    # remove group from old cluster
    group_to_operate=meta_dedup$group_ID[group_ID]
    kth_cluster_ind=which(names(ak)==meta_dedup$cluster_number[group_ID])
    
    if (phi[kth_cluster_ind]==1)
    {
      K=K-1
      phi=phi[-kth_cluster_ind]
      t0=t0[-kth_cluster_ind]
      dt=dt[-kth_cluster_ind]
      de=de[-kth_cluster_ind]
      ak=ak[-kth_cluster_ind]
    }else
    {
      phi[kth_cluster_ind]=phi[kth_cluster_ind]-1
    }
    
    # create a new cluster
    t0_new=rnorm(n=dim(t)[1],mean_t,sqrt(lambda))
    ak_new=rnorm(1,a,sigma*sqrt(g))
    
    # the assignment
    prob=calulate_DPprob(phi,xi,group_ID,b,t0_new,ak_new,master_dist_e,sigma,t,ak,t0)
    if(!is.null(sample_id)){
      sample_id_cluster=sapply(names(prob)[-length(prob)],function(name) strsplit(name,split = ';')[[1]][2])
      prob2rm=sample_id_cluster!=sample_id_dedup[group_ID]
      prob2rm=c(prob2rm,FALSE)
      if(sum(prob[!prob2rm])==0){
        prob[!prob2rm]=1
        print(paste('Random selection:',sample_id_dedup[group_ID]))
      }
      prob[prob2rm]=0
    }
    prob2sc=prob/sum(prob)
    new_cluster=sample(1:length(prob),1,prob=prob2sc)
    old_cluster_name=meta_dedup$cluster_number[group_ID]
    
    # move the group to operate from the old cluster into the new cluster
    if (new_cluster==length(prob)) # the newly created cluster
    {
      meta_dedup$cluster_number[group_ID]=group_to_operate # this CDR3 becomes the center naturally
      K=K+1
      # modify the old cluster
      keep=meta_dedup$cluster_number==old_cluster_name & meta_dedup$group_ID!=group_to_operate
      go_back=F
    }else # one of the old clusters
    {
      new_cluster_name=names(prob)[new_cluster]
      meta_dedup$cluster_number[group_ID]=new_cluster_name
      phi[new_cluster_name]=phi[new_cluster_name]+1
      # modify the old cluster
      keep=meta_dedup$cluster_number==old_cluster_name
      go_back=new_cluster_name==old_cluster_name # the picked cluster is the original cluster
    }
    
    # modify the old cluster
    if (sum(keep)>0 && (!go_back)) 
    {
      old_cluster_name1=find_center(t[,keep,drop=F],b)
      meta_dedup$cluster_number[keep]=old_cluster_name1
      which_to_update=which(names(phi)==old_cluster_name)
      names(phi)[which_to_update]=old_cluster_name1
      names(t0)=names(ak)=names(dt)=names(de)=names(phi)
      old_cluster_name=old_cluster_name1
      tmp=names(de[[old_cluster_name]])
      groups_to_operate=tmp[tmp!=group_to_operate]
      de_old_cluster_name=master_dist_e[groups_to_operate,old_cluster_name]
      dt_old_cluster_name=colSums((t[,groups_to_operate,drop=F]-t[,old_cluster_name])^2/b/2)
      de[[old_cluster_name]]=named_c(NULL,de_old_cluster_name,groups_to_operate)
      dt[[old_cluster_name]]=named_c(NULL,dt_old_cluster_name,groups_to_operate)
    }
      
    # modify the new cluster
    if (new_cluster==length(prob)) # new cluster
    {
      ak=named_c(ak,ak_new,group_to_operate)
      phi=named_c(phi,1,group_to_operate)
      t0[[group_to_operate]]=t0_new
      tmp=master_dist_e[group_to_operate,group_to_operate]
      de[[group_to_operate]]=named_c(NULL,tmp,group_to_operate)
      dt[[group_to_operate]]=named_c(NULL,0,group_to_operate)
    }else if (!go_back) # one of the old cluster, and also a different one
    {
      keep=meta_dedup$cluster_number==new_cluster_name
      new_cluster_name1=find_center(t[,keep,drop=F],b)
      meta_dedup$cluster_number[keep]=new_cluster_name1
      which_to_update=which(names(phi)==new_cluster_name)
      names(phi)[which_to_update]=new_cluster_name1
      names(t0)=names(ak)=names(dt)=names(de)=names(phi)
      new_cluster_name=new_cluster_name1
      groups_to_operate=c(names(de[[new_cluster_name]]),group_to_operate)
      de_new_cluster_name=master_dist_e[groups_to_operate,new_cluster_name]
      dt_new_cluster_name=colSums((t[,groups_to_operate,drop=F]-t[,new_cluster_name])^2/b/2)
      de[[new_cluster_name]]=named_c(NULL,de_new_cluster_name,groups_to_operate)
      dt[[new_cluster_name]]=named_c(NULL,dt_new_cluster_name,groups_to_operate)
    }
  }
  
  return(list(cluster_number=meta_dedup$cluster_number,phi=phi,de=de,K=K,
    dt=dt,ak=ak,t0=t0))
}
