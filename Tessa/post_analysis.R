getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
plot_tessa<-function(tessa_results,folder,labels=NA)
{
  library(data.table)
  # load data
  meta_dedup=tessa_results$meta_dedup
  meta=tessa_results$meta
  master_dist_e=tessa_results$master_dist_e
  b=tessa_results$b
  t=tessa_results$t
  K=tessa_results$K
  dt=tessa_results$dt
  de=tessa_results$de
  ak=tessa_results$ak
  phi=tessa_results$phi
  if (!file.exists(folder)) {dir.create(folder)}
  # are cluster centers more expanded than other members of the clusters?
  pdf(paste(folder,"/clone_size.pdf",sep=""),width=4,height=4)
  tmp=aggregate(meta$barcode,by=list(meta$group_ID),length)
  meta_dedup$n=NA
  meta_dedup[tmp$Group.1,"n"]=tmp$x
  if(!is.na(labels)){
    tmp2=aggregate(labels,by=list(meta$group_ID),getmode)
    tmp_keep=data.table(cbind(meta$group_ID,labels))
    tmp_table=table(tmp_keep)/tmp$x
    tmp_table=tmp_table[match(tmp2$Group.1,row.names(tmp_table)),]
    keep=sapply(1:nrow(tmp_table),function(i) tmp_table[i,tmp2$x[i]])
    meta_dedup[tmp2$Group.1,"labels"]=tmp2$x
    meta_dedup[tmp2$Group.1[keep<0.6],"labels"]=999
  }
  meta_dedup0=meta_dedup[meta_dedup$cluster_number %in% names(phi[phi>2]),]
  keep=meta_dedup0$group_ID!=meta_dedup0$cluster_number
  non_center_n=aggregate(meta_dedup0[keep,"n"],by=list(meta_dedup0$cluster_number[keep]),median)
  if(!is.na(labels)){
    meta_dedup0$col=c('salmon','blue','grey')[factor(meta_dedup0$labels,levels = c(1,0,999))]
    #salmon: CD8; blue: CD4; 'grey: mixed center (CD8 and CD4 cells -> same TCRs) 
    #or the center expresses same levels of CD8 and CD4'
    plot(non_center_n$x,meta_dedup0[non_center_n$Group.1,"n"],xlab="Clone size of non-center TCR groups",
         ylab="Clone size of center TCR groups",pch=19,col=meta_dedup0[non_center_n$Group.1,"col"])
    abline(0,1,col="red")
  }else{
    plot(non_center_n$x,meta_dedup0[non_center_n$Group.1,"n"],xlab="Clone size of non-center TCR groups",
       ylab="Clone size of center TCR groups",pch=19)
    abline(0,1,col="red")
  }
  dev.off()
  mean(non_center_n[,2])
  mean(meta_dedup0[non_center_n$Group.1,"n"])
  
  # expression-TCR distance plot 
  pdf(paste(folder,"/exp_TCR_pair_plot.pdf",sep=""),width=6,height=6)
  par(mfrow=c(3,2))
  for (k in 1:K) 
  {
    plot(de[[k]],dt[[k]],main=paste("Cluster:",names(dt)[k]),pch=19,
      col=1+as.numeric(as.factor(names(de[[k]]))),xlab="Expression dist.",
      ylab="TCR dist.")
    segments(x0=0,y0=0,y1=max(dt[[k]]),x1=ak[k]*max(dt[[k]]),lwd=2,lty=2)
  }
  dev.off()
  
  # density of TCR distances
  pdf(paste(folder,"/TCR_dist_density.pdf",sep=""),width=4,height=3)
  plot(density(as.matrix(dist(t(t/sqrt(b))))^2/2),xlab="TCR distances",ylab="",lwd=3)
  lines(density(unlist(dt)),lwd=3,col="red")
  dev.off()
  
  # exploratory plot at the TCR level
  pdf(paste(folder,"/TCR_explore.pdf",sep=""),width=20,height=20)
  pca_t=prcomp(t(t/sqrt(b)),scale.=F)$x
  plot(pca_t[,1],pca_t[,2],type="n",xlab="PC1",ylab="PC2")
  tmp=as.numeric(as.factor(meta_dedup$cluster_number))
  names(tmp)=meta_dedup$group_ID
  text(pca_t[,1],pca_t[,2],label=tmp,col=tmp,cex=0.5+0.5*(rownames(pca_t) %in% names(dt)))
  for (k in 1:K)
  {
    for (group in names(dt[[k]]))
    {
      segments(x0=pca_t[names(dt)[k],1],y0=pca_t[names(dt)[k],2],
        x1=pca_t[group,1],y1=pca_t[group,2],col=tmp[group])
    }
  }
  dev.off()
}

predict_tessa<-function(tessa_results,t_new,cutoff=NA)
{
  b=tessa_results$b
  if(is.na(cutoff)){
    cutofflist=sapply(tessa_results$dt,function(x) quantile(x,0.5))
    cutoff=quantile(cutofflist[cutofflist!=0],0.4)
  }
  cluster_new=cutree(hclust(dist(t(t_new/sqrt(b)),method='manhattan'),method='single'),h=cutoff)
}

plot_tessaClsuters=function(tessa_results,save_path){
  meta=tessa_results$meta
  meta_dedup=tessa_results$meta_dedup
  b=tessa_results$b
  tb=table(meta$group_ID)
  meta_dedup$clonalsize=tb[match(meta_dedup$group_ID,names(tb))]
  t=tessa_results$t
  phi=tessa_results$phi
  keep=names(phi)[phi>quantile(phi,0.75)]
  t=t[,meta_dedup$cluster_number%in%keep]
  meta_dedup=meta_dedup[meta_dedup$cluster_number%in%keep,]
  set.seed(123)
  data2plot=data.frame(Rtsne(t(t/sqrt(b)),dims = 2)$Y)
  colnames(data2plot)=c('tSNE_1','tSNE_2')
  row.names(data2plot)=meta_dedup$group_ID
  data2plot$clonalsize=as.numeric(meta_dedup$clonalsize)
  colorbase=colorRampPalette(rainbow(8))(length(keep))
  g=ggplot()+theme_bw(base_size = 12)+guides(size=FALSE)+
    xlim(range(data2plot$tSNE_1))+ylim(range(data2plot$tSNE_2))
  for(i in 1:length(unique(meta_dedup$cluster_number))){
    select_cluster=unique(meta_dedup$cluster_number)[i]
    meta_dedup_tmp=meta_dedup[meta_dedup$cluster_number==select_cluster,]
    data_tmp=data2plot[meta_dedup$cluster_number==select_cluster,]
    center=which(meta_dedup_tmp$group_ID==meta_dedup_tmp$cluster_number)
    data_tmp$ct1=data_tmp$tSNE_1[center]
    data_tmp$ct2=data_tmp$tSNE_2[center]
    g=g+geom_segment(data = data_tmp,color='grey80',aes(x=tSNE_1,y=tSNE_2,xend=ct1,yend=ct2),size=0.1)+
      geom_point(data = data_tmp,color=colorbase[i],alpha=0.8,aes(x=tSNE_1,y=tSNE_2,size=sqrt(clonalsize)))
  }
  ggsave(filename =paste(save_path,'/cluster_figure.pdf',sep = ''),plot = g,height=8,width=8)
}

  
