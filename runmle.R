library(gtools)
library(rstan)
readyparallel()
topetho <- stan_model("~/code/topetho/categorical_topetho.stan")

n <- 4000
K <- 10
p <- rdirichlet(1,rep(2,K)) %>% as.vector()
nk <- rmultinom(1,size = n,prob = p) %>% as.vector()

pik <- matrix((1-0.9)/K,K,K)
diag(pik) <- 0.9

Y <- matrix(nrow=0,ncol=2)
for (i in 1:K) Y <- rbind(Y, cbind( (rmultinom(nk[i],1,pik[,i]) == 1) %>% apply(MARGIN = 2,FUN = which),
       (rmultinom(nk[i],1,pik[,K-i+1]) == 1) %>% apply(MARGIN = 2,FUN = which)))

dat <- list(n=n,K=K,B=2,Bs=rep(K,2),Y=Y,alpha_p=1,alpha_t=1)
moo <- optimizing(topetho,dat,verbose=T,vector)

qplot(y=moo$par[str_detect(names(moo$par),"^r\\[")] %>% matrix(nrow=n,ncol=K) %>% apply(MARGIN = 1,which.max),x=1:n)

dat <- list(n=nrow(Y),K=6,B=ncol(Y)-ncovcols,Bs=sapply(Y[,-c(1:ncovcols),with=F],max),Y=as.matrix(Y[,-c(1:ncovcols),with=F]),alpha_p=1,alpha_t=1)

init <- list(pi=rdirichlet(1,alpha = rep(1,dat$K)) %>% as.vector(),
             theta_raw=sapply(Y[,-(1:ncovcols),with=F],function(x) table(x) %>% prop.table) %>% unlist() %>% matrix(nrow=dat$K,ncol=sum(dat$Bs),byrow = T))
init$theta_raw <- init$theta_raw * pmax(1-rnorm(length(init$theta_raw),sd=0.5),0.01)
moo <- optimizing(topetho,dat,verbose=T,init=init,as_vector=F)

pm <- moo$value + nrow(moo$hessian)/2 * log(2*pi) - 0.5*(determinant(-moo$hessian,logarithm=T) %>% (function(x) c(x$sign*x$modulus)))


foo3 <- foreach(1:42) %dopar% { library(gtools); library(rstan)
  init <- list(pi=rdirichlet(1,alpha = rep(1,dat$K)) %>% as.vector(),
               theta_raw=sapply(Y[,-(1:ncovcols),with=F],function(x) table(x) %>% prop.table) %>% unlist() %>% matrix(nrow=dat$K,ncol=sum(dat$Bs),byrow = T))
  #init$theta_raw <- init$theta_raw * pmax(1-rnorm(length(init$theta_raw),sd=0.5),0.01)
  moo <- optimizing(topetho,dat,verbose=T,init=init,as_vector=F,iter=600)
  return(moo)
}


guh <- data.table(FocalID=Y$FocalID,z=apply(foo2[[30]]$par$r,1,which.max))
juh <- merge(guh[,.(y=sum(z==5),baseline=sum(z==7)),by=FocalID],Xdf,by="FocalID")
juh$FocalID2 <- juh$FocalID
glm(cbind(y,baseline) ~ sex*poly(age,2) + group,data=juh,family=binomial) %>% logLik()
cheetah.mm(cbind(y,baseline) ~ sex*poly(age,2) + group + (1|FocalID),data=juh,family=binomial) %>% logLik()
cheetah.mm(cbind(y,baseline) ~ sex*poly(age,2) + group + (1|FocalID) + (1|FocalID2),pedigree = list(FocalID=A),data=juh,family=binomial) %>% logLik()
