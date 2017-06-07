mar <- Xdf[sex=="m",range(age)]
far <- Xdf[sex=="f",range(age)]
simdf <- data.table(sex=rep(c("f","m"),c(200,200)),
                    age=c(seq(far[1],far[2],length.out = 200),seq(mar[1],mar[2],length.out = 200)),
                    rank=Xdf[,mean(rank)],group="F")
tpdat <- betadat[iter>250,.(FocalID=Xdf$FocalID,tphat=as.vector(X %*% beta)),by=.(iter,topic)]
tpsimhat <- tpdat[,.(ID=1:400,age=simdf$age,sex=simdf$sex,etahat=lm(tphat ~ group +sex*poly(age,2) + sex*poly(rank,2),dat=cbind(Xdf,tphat)) %>% predict(newdata=simdf)),
                  ,by=.(iter,topic)]
tpsimhat <- merge(tpsimhat,mudatref)
tpsimhat[,etahat:=etahat+mu]; tpsimhat[,mu:=NULL]
tpsimhat[,prob:=exp(etahat)/sum(exp(etahat)),by=.(iter,ID)]
tphatsumm <- tpsimhat[,.(age=age[1],sex=sex[1],mu=mean(prob),lb=quantile(prob,0.05),ub=quantile(prob,0.95)),by=.(ID,topic)]
ggplot(tphatsumm,aes(x=age,y=mu,color=sex)) + geom_line() + geom_ribbon(aes(ymin=lb,ymax=ub,fill=sex),alpha=0.25,linetype=0) +
  facet_wrap(~topic,nrow = 2,scales = "free_y") + theme_light() + scale_x_continuous(name="Age",minor_breaks = NULL) + scale_y_continuous(name="Probability",minor_breaks = NULL)


simdf <- data.table(sex=rep(c("f","m"),c(3,3)),
                    rank=Xdf[,mean(rank)],
                    age=Xdf[,mean(age)],group=c("F","HH","R"))
tpsimhat <- tpdat[,.(ID=1:6,group=simdf$group,sex=simdf$sex,etahat=lm(tphat ~ group +sex*poly(age,2) + sex*poly(rank,2),dat=cbind(Xdf,tphat)) %>% predict(newdata=simdf)),
                  ,by=.(iter,topic)]
tpsimhat <- merge(tpsimhat,mudat)
tpsimhat[,etahat:=etahat+mu]; tpsimhat[,mu:=NULL]
tpsimhat[,prob:=exp(etahat)/sum(exp(etahat)),by=.(iter,ID)]
tphatsumm <- tpsimhat[,.(group=group[1],sex=sex[1],mu=mean(prob),lb=quantile(prob,0.05),ub=quantile(prob,0.95)),by=.(ID,topic)]
ggplot(tphatsumm,aes(x=group,y=mu,color=sex)) + 
  geom_point(position=position_dodge(width=0.5)) + geom_errorbar(aes(ymin=lb,ymax=ub),width=0,position = position_dodge(width = 0.5)) + theme_light() + facet_wrap(~topic,scale="free_y",nrow = 2) +
  scale_y_continuous(name = "Probability",minor_breaks = NULL) + xlab("Social Group") + scale_color_discrete("",labels=c("Female","Male"))

