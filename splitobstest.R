obsdat <- readfocfiles(fpath,c("F","F","HH","R"))
obs <- obsdat[,unique(Observation)]
cutoff <- 310
splitobsdat <- list()
for (i in 1:length(obs)) splitobsdat[[i]] <- splitstate(obsdat[Observation==obs[i]],cutoff) #function at bottom
splitobsdat <- do.call(rbind,splitobsdat)
splitobsdat[RelativeEventTime>=cutoff,Observation:=paste0(Observation,".2")]
splitobsdat[RelativeEventTime<cutoff,Observation:=paste0(Observation,".1")]
splitobsdat[RelativeEventTime>=cutoff,RelativeEventTime:=RelativeEventTime-cutoff]
ptetho <- defaultpoint2()[type!="misc" & !(behavior%in%c("Vigilnce","PsCnTerm","GrmTerm","GrmPrsnt"))]
stetho <- defaultstate2()[type!="misc" & state!="Corral"]
Y <- collectfocal(splitobsdat,ptetho,stetho,state.maxsec = 320)

#resume R2julia here

filt <- Y[,lapply(.SD,function(x) mean(x>0)),.SD=eventslices(names(Y),ptetho)] > 0.005
filt <- colnames(filt)[!filt] %>% c(stetho[baseline==T,behavior]) %>% unique()
filt <- c(filt,"ScanProx","ScanProxInd")
filt <- filt[-c(1,5)]
Y <- Y[,-filt,with=F]


dat <- list(n=nrow(Y),K=10,B=ncol(Y)-ncovcols,Bs=sapply(Y[,-c(1:ncovcols),with=F],max),Y=as.matrix(Y[,-c(1:ncovcols),with=F]),alpha_p=1,alpha_t=1)

foo3 <- foreach(1:8) %dopar% { library(gtools); library(rstan)
  init <- list(pi=gtools::rdirichlet(1,alpha = rep(1,dat$K)) %>% as.vector(),
               theta_raw=sapply(Y[,-(1:ncovcols),with=F],function(x) table(x) %>% prop.table) %>% unlist() %>% matrix(nrow=dat$K,ncol=sum(dat$Bs),byrow = T))
  init$theta_raw <- init$theta_raw * pmax(1-rnorm(length(init$theta_raw),sd=0.5),0.01)
  moo <- optimizing(topetho,dat,verbose=T,init=init,as_vector=F,iter=500)
  return(moo)
}


splitstate <- function(obsdat,cutoff=330) {
  target <- obsdat[,RelativeEventTime < cutoff & (RelativeEventTime+Duration)>cutoff & Duration>0]
  repacts <- obsdat[target==T]
  repacts[,Duration:=Duration-cutoff+RelativeEventTime]
  repacts[,RelativeEventTime:=cutoff]
  obsdat[target,Duration:=cutoff-RelativeEventTime]
  newdat <- rbind(obsdat,repacts)
  setkey(newdat,"RelativeEventTime")
  return(newdat)
}
