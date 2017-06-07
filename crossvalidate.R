#source('~/code/LogisticTopicRegression/R2julia.R')
library(mvtnorm)
library(rstan)
readyparallel(4)
rstan_options(auto_write = TRUE)
topetho <- stan_model("~/code/topetho/categorical_topetho.stan")

## divide folds
nfold <- 10
assfold <- vector("numeric",nrow(Y))
ID <- Y[,unique(FocalID)]
for (i in 1:length(ID)) {
  IDi <- Y[,which(FocalID==ID[i])]
  n <- length(IDi)
  blksz <- n %/%nfold
  afi <- c(rep(1:nfold,blksz),sample(1:nfold,n %% nfold))
  assfold[IDi] <- sample(afi,n)
}

## alternatively, load folds
assfold <- read.csv("~/analysis/logtopreg/crossvalidation/folds_5k.csv",header = F)[,1]
nfold <- max(assfold)

## factor analyis

fnum <- c(5,10,15)
lppd <- matrix(nrow=nfold,ncol=length(fnum))
for (i in 1:nfold) {
  Y1 <- Y[assfold!=i,-(1:ncovcols),with=F]
  Y2 <- Y[assfold==i,-(1:ncovcols),with=F]

  Ysd <- sapply(Y1,sd)
  Ymu <- sapply(Y1,mean)
  for (j in 1:length(fnum)) {
    fa <- factanal(Y1,factors = fnum[j],control=list(nstart=50,lower=1e-8))
    rhohat <- fa$loadings %*% t(fa$loadings) + diag(fa$uniquenesses)
    Sigma <- diag(Ysd) %*% rhohat %*% diag(Ysd)
    
    lppd[i,j] <- apply(as.matrix(Y2),MARGIN = 1,dmvnorm,mean = Ymu,sigma = Sigma,log = T) %>% sum
  }
  #lppd[i,length(fnum)+1] <- apply(as.matrix(Y2),MARGIN = 1,dmvnorm,mean = Ymu,sigma = cov(Y1),log = T) %>% sum
  
}

## get starting points for topic model
r <- list()
for (i in 1:nfold) {
  Y1 <- Y[assfold!=i,-(1:ncovcols),with=F]
  dat <- list(n=nrow(Y1),K=10,B=ncol(Y1),Bs=sapply(Y1,max),Y=Y1,alpha_p=1,alpha_t=1)
  initout <- foreach(1:12) %dopar% { library(gtools); library(rstan)
    init <- list(pi=rdirichlet(1,alpha = rep(1,dat$K)) %>% as.vector(),
                 theta_raw=sapply(Y[,-(1:ncovcols),with=F],function(x) table(x) %>% prop.table) %>% unlist() %>% matrix(nrow=dat$K,ncol=sum(dat$Bs),byrow = T))
    init$theta_raw <- init$theta_raw * pmax(1-rnorm(length(init$theta_raw),sd=0.5),0.01)
    moo <- optimizing(topetho,dat,verbose=T,init=init,as_vector=F,iter=200)
    return(moo)
  }
  
  bstft <- sapply(initout,function(x) x$value) %>% which.max
  r[[i]] <- initout[[bstft]]$par$r
  r[[i]] <- cbind(r[[i]],i)
}
r <- do.call(rbind,r)
