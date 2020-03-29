# truncated normal 
rtnorm<-function(mu,sd=1,a=-Inf,b=Inf)
{
  #truncated normal distribution sampling function.
  F<-runif(n=length(mu))
  Fa<-pnorm((a-mu)/sd,0,sd=1)
  Fa[a==-Inf] <-0
  Fb<-pnorm((b-mu)/sd,0,sd=1)
  Fb[b==Inf]<-1
  y<-mu+sd*qnorm(F*(Fb-Fa)+Fa)
  y
}

named_c<-function(vector,new_element,new_name)
{
  names(new_element)=new_name
  c(vector,new_element)
}