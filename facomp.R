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
#etasumm <- etadat[iter>100,.(prob=mean(prob),eta=mean(eta)),by=.(FocalID,year,topic)]
etarep <- etadat[FocalID %in% repid[V1==2,FocalID]] %>% dcast(FocalID + topic + iter ~ year)
reptopic <- etarep[,cor(`2012`,`2013`),by=.(iter,topic)]
reptopic$model <- "State model"

ggplot(rbind(reptopic[,-1],repscores[,-1],repscores2[,-1]),aes(x=V1,fill=model)) + geom_histogram(position="dodge",bins=10) + theme_classic() + xlab("Correlation between 2012 and 2013 phenoypes") + ylab("# of States/Factors")
