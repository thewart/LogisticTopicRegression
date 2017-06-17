library(mvtnorm)
source('~/code/LogisticTopicRegression/R2julia.R')
ID <- Y[,unique(FocalID)]


fa1 <- factanal(Y[,-(1:ncovcols)],factors = 10,scores = "regression",control=list(nstart=50,lower=1e-8))
scores <- cbind(Y[,c(1,2,4)],fa1$scores) %>% melt()

#repeatability -- obs level
repid <- scores[,length(unique(Year)),by=FocalID]
repscores <- scores[FocalID %in% repid[V1==2,FocalID],mean(value),by=.(FocalID,Year,variable)]
repscores <- dcast(repscores,FocalID+variable~ Year)
repscores <- repscores[,cor(`2012`,`2013`),by=variable]
repscores$model <- "Factor model 1"

#repeatability -- org level
Y2 <- melt(Y[,-c(3,5)])[,mean(value),by=.(FocalID,Year,variable)] %>% dcast(formula=FocalID+Year~variable)
fa2 <- factanal(Y2[,-(1:2)],factors = 10,scores = "regression",control=list(nstart=50,lower=1e-8))
repscores2 <- cbind(Y2[,1:2],fa2$scores)[FocalID %in% repid[V1==2,FocalID]]
repscores2 <- melt(repscores2) %>% dcast(formula = FocalID+variable ~ Year)
repscores2 <- repscores2[,cor(`2012`,`2013`),by=variable]
repscores2$model <- "Factor model 2"

#get repeatability from topic model
Xdf <- data.table(FocalID=Y$FocalID,sex=getsex(Y$FocalID),age=getage(Y$FocalID,Y$Year),group=Y$Group,year=Y$Year)
Xdf[,age:=mean(unique(age)),by=FocalID]
Xdf <- merge(Xdf,drank,by.x = "FocalID",by.y="ID")
Xdf <- unique(Xdf)
path <- "/home/seth/analysis/logtopreg/fit_repeatability/"
d <- scan(paste0(path,"dims.csv")) %>% as.list
names(d) <- c("n","K","p","b","nsave")
etho <- rbind(ptetho[,c(1,3),with=F],stetho[,c(1,4),with=F])
etho <- rbind(etho,data.table(behavior=c("ScanProx","ScanProxInd"),type=rep("Affiliative",2)))
etadat <- data.table(FocalID=Xdf$FocalID %>% rep(.,rep(d$K,d$n)),
                     year=Xdf$year %>% rep(.,rep(d$K,d$n)),
                     topic=paste0("S",1:10),
                     iter=rep(1:d$nsave,rep(d$n*d$K,d$nsave)),
                     eta=scan(paste0(path,"eta.csv")))
etadat[,prob:=exp(eta)/sum(exp(eta)),by=.(FocalID,year,iter)]
etasumm <- etadat[iter>100,.(prob=mean(prob),eta=mean(eta)),by=.(FocalID,year,topic)]
etarep <- etasumm[FocalID %in% repid[V1==2,FocalID]] %>% dcast(FocalID + topic ~ year,value.var="prob")
reptopic <- etarep[,cor(`2012`,`2013`),by=.(topic)]
reptopic$model <- "State model"

ggplot(rbind(reptopic[,-1],repscores[,-1],repscores2[,-1]),aes(x=V1,fill=model)) + geom_histogram(position="dodge",bins=10) + theme_classic() + xlab("Correlation between 2012 and 2013 phenoypes") + ylab("# of States/Factors")

#### heritability
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
lmma <- stan_model("~/code/multivarlmm/lmm_animal.stan")
LA <- t(chol(A))
standat <- list(N=nrow(X),P=ncol(X),X=X,L_A=LA,R=nrow(A))

mscores <- scores[,mean(value),by=.(variable,FocalID)]
facts <- mscores[,unique(variable)]
fit <- list()
for (i in 1:10){
  standat$Y=mscores[variable==facts[i],V1]
  fit[[i]] <- sampling(lmma,data=standat,chains=2,iter=2000,warmup=500,thin=2,
                       pars=c("alpha","beta","sigma_eps","sigma_g","h2"))
}



Y2 <- melt(Y[,-c(3,5)])[,mean(value),by=.(FocalID,variable)] %>% dcast(formula=FocalID~variable)
fa2 <- factanal(Y2[,-(1)],factors = 10,scores = "regression",control=list(nstart=50,lower=1e-8))
scores2 <- data.table(FocalID=Y2$FocalID,fa2$scores)
fit2 <- list()
for (i in 1:10) {
  standat$Y <- scores2[,i+1,with=F] %>% unlist()
  fit2[[i]] <- sampling(lmma,data=standat,chains=2,iter=2000,warmup=500,thin=2,
                pars=c("alpha","beta","sigma_eps","sigma_g","h2"))
}


##### simulate fa data
Ysd <- sapply(Y[,-(1:ncovcols)],sd)
Ymu <- sapply(Y[,-(1:ncovcols)],mean)
Sigma <- diag(Ysd) %*% (fa1$loadings %*% t(fa1$loadings) + diag(fa1$uniquenesses)) %*% diag(Ysd)
#Sigma <- cov(Y[,-(1:ncovcols)])
Yrep_fa <- rmvnorm(1e6,Ymu,Sigma) %>% as.data.table()
for (i in 1:length(Yrep_fa)) {
  maxi <- max(Y[,-(1:ncovcols)][[i]])
  Yrep_fa[[i]][Yrep_fa[[i]]<1] <- 1
  Yrep_fa[[i]][Yrep_fa[[i]] > maxi] <- maxi
  Yrep_fa[[i]] <- round(Yrep_fa[[i]])
}
names(Yrep_fa) <- behname
Yrep_fa <- Yrep_fa[,.(Travel,Feed,`NonContactAgg(receive)`)]
Yrep_fa[,`NonContactAgg(receive)`:=`NonContactAgg(receive)` > 1]
Yrep_fa$Model <- "FA1"

#get real data
Ytmp <- Y[,-(1:ncovcols)]
names(Ytmp) <- behname
Ytmp <- Ytmp[,.(Travel,Feed,`NonContactAgg(receive)`)]
Ytmp[,`NonContactAgg(receive)`:=`NonContactAgg(receive)` > 1]
Ytmp$Model <- "Cayo data"

#simulate topic model
pl <- Y[,-(1:ncovcols)] %>% sapply(max)
thetadat <- data.table(value=scan("~/analysis/logtopreg/fitFKKR_A_fixed_2/topicparams.csv"),
                       behavior=rep(behname,pl),
                       level=lapply(pl,function(x) 1:x) %>% unlist() %>% as.factor(),
                       topic=rep(topicord_eta,each=sum(pl)),
                       iter=rep(1:d$nsave,each=d$K*sum(pl)))

thetamu <- thetadat[,.(value=mean(value)),by=.(topic,behavior,level)]
#thetamu <- thetamu[behavior %in% c("Scratch","Approach(give)","Groom(give)")]
thetamu <- thetamu[behavior %in% c("Travel","Feed","NonContactAgg(receive)")]

nsamp <- 125
Yrep_top <- data.table(Feed=vector(),`NonContactAgg(receive)`=vector(),Travel=vector())
for (j in 1:length(ID)) {
  ksamp <- rmultinom(n=1,
                     size=Y[,length(unique(Observation)),by=FocalID]$V1[j] * nsamp,
                     prob=etasumm[FocalID==ID[j],prob])
  
  for (i in 1:d$K) { 
    if (!ksamp[i]) next
    tmp <- thetamu[topic==topicord_eta[i],(rmultinom(ksamp[i],1,value)==1) %>% apply(2,which),by=.(behavior)]
    tmp[,samp:=1:ksamp[i]]
    Yrep_top <- rbind(Yrep_top,dcast(tmp,samp ~ behavior,value.var = "V1")[,-1])
  }
}
Yrep_top[,`NonContactAgg(receive)`:=`NonContactAgg(receive)` > 1]
setcolorder(Yrep_top,c(3,1,2))
Yrep_top$Model <- "State model"
repdat <- rbind(Yrep_top,Ytmp,Yrep_fa)[,.(Travel=mean(Travel),n=length(Travel),ste=2*sd(Travel)/sqrt(length(Travel))),by=.(Feed,`NonContactAgg(receive)`,Model)]
repdat[Model!="Cayo data",ste:=0]
repdat[,Model:=ordered(Model,levels=unique(Model)[c(2,1,3)])]
repdat[,`NonContactAgg(receive)`:=ordered(`NonContactAgg(receive)`,levels=c(T,F),labels=c("High","Low"))]
repdat <- repdat[!(Feed==5 & `NonContactAgg(receive)`=="High")]
ggplot(repdat,aes(x=Feed,y=Travel,color=`NonContactAgg(receive)`)) + geom_point(size=1,position=position_dodge(width=0.25)) +
  geom_line(position=position_dodge(width=0.25)) + geom_errorbar(aes(ymin=Travel-ste,ymax=Travel+ste),width=0,position=position_dodge(width=0.25)) +
  facet_wrap(~Model,scale="free") + scale_x_continuous(minor_breaks = NULL) + scale_color_discrete(name="Aggression\nreceived")
