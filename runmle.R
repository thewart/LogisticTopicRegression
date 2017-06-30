library(gtools)
library(rstan)
readyparallel(4)
topetho <- stan_model("~/code/topetho/categorical_topetho.stan")

#get Y from R2julia
dat <- list(n=nrow(Y),K=10,B=ncol(Y)-ncovcols,Bs=sapply(Y[,-c(1:ncovcols),with=F],max),Y=as.matrix(Y[,-c(1:ncovcols),with=F]),alpha_p=1,alpha_t=1)

fit <- foreach(1:42) %dopar% { library(gtools); library(rstan)
  init <- list(pi=rdirichlet(1,alpha = rep(1,dat$K)) %>% as.vector(),
               theta_raw=sapply(Y[,-(1:ncovcols),with=F],function(x) table(x) %>% prop.table) %>% unlist() %>% matrix(nrow=dat$K,ncol=sum(dat$Bs),byrow = T))
  init$theta_raw <- init$theta_raw * pmax(1-rnorm(length(init$theta_raw),sd=0.5),0.01)
  moo <- optimizing(topetho,dat,verbose=T,init=init,as_vector=F,iter=600)
  return(moo)
}