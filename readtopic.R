<<<<<<< HEAD
path <- "/home/seth/analysis/logtopreg/fitFKKR_forpaper/"
=======
library(ggrepel)
path <- "/home/seth/analysis/logtopreg/fitFKKR_A_fixed_2/"
>>>>>>> raneffcov
d <- scan(paste0(path,"dims.csv")) %>% as.list
names(d) <- c("n","K","p","b","nsave")
etho <- rbind(ptetho[,c(1,3),with=F],stetho[,c(1,4),with=F])
etho <- rbind(etho,data.table(behavior=c("ScanProx","ScanProxInd"),type=rep("Affiliative",2)))

mudat <- fread(paste0(path,"mu.csv"),header = F)
mudat$iter <-1:d$nsave
mudat <- melt(mudat,id.vars="iter",variable.name = "topic",value.name = "mu")
mudat[,prob:=exp(mu)/ sum(exp(mu)),by=iter]
topicord <- mudat[iter>100,mean(prob),by=topic][,rank(-V1)] %>% ordered()
levels(topicord) <- paste0("S",levels(topicord))
mudat[,topic:=rep(topicord,rep(d$nsave,d$K))]

nzmu <- mudat[topic %in% mudat[iter>101,mean(prob),by=topic][,topic[V1>0.001]]]

#monkey-specific topic probs
etadat <- data.table(FocalID=Xdf$FocalID %>% rep(.,rep(d$K,d$n)),
                  topic=topicord,
                  iter=rep(1:d$nsave,rep(d$n*d$K,d$nsave)),
                  eta=scan(paste0(path,"eta.csv")))
etadat[,prob:=exp(eta)/sum(exp(eta)),by=.(FocalID,iter)]
etamu <- etadat[iter>100,.(prob=mean(prob),std=sd(prob)),by=.(iter,topic)][,.(prob=mean(prob),std=mean(std)),by=topic]
topicord_eta <- etamu[,rank(-prob)] %>% ordered()
levels(topicord_eta) <- paste0("S",levels(topicord_eta))
etadat[,topic:=topicord_eta]
etamu[,topic:=topicord_eta]
etasumm <- etadat[iter>100,.(prob=mean(prob),eta=mean(eta)),by=.(FocalID,topic)]
etadat[,prob:=exp(eta)/sum(exp(eta)),by=.(FocalID,iter)]

udat <- data.table(FocalID=Xdf$FocalID,
                  topic=rep(topicord_eta,rep(d$n,d$K)),
                  iter=rep(1:d$nsave,rep(d$n*d$K,d$nsave)),
                  u=scan(paste0(path,"u.csv")))


#population topic weight probabilities
ggplot(merge(etasumm,Xdf,by="FocalID"),aes(x=topic,y=prob)) + geom_boxplot() +
  ylab("Probability") + xlab("Behavioral state") + scale_y_continuous(limits = c(0,0.655),expand = c(0,0)) +
  theme_light() + theme(panel.grid.major.x = element_blank())

#potpulation topic weight variabilities
ggplot(etadat[,.(var=sd(prob)),by=.(topic,iter)][,.(mean(var),quantile(var,0.025),quantile(var,0.975)),by=topic],aes(x=topic,y=V1)) + geom_point() + geom_errorbar(aes(ymin=V2,ymax=V3),width=0)

#topic weights against covariates
ggplot(merge(etasumm[topic=="S5"],Xdf,by="FocalID"),aes(y=prob,x=rank,color=sex)) + 
  geom_point() + geom_smooth(method="lm",formula=y ~ poly(x,2),level=0.999) +
  xlab("Age") + ylab("State 3 Probability") + theme_light() + scale_color_discrete(labels=c("Female","Male"),guide=guide_legend(title=""))

ggplot(merge(etasumm[topic=="S4"],Xdf,by="FocalID"),aes(y=prob,x=group)) + 
  geom_boxplot() + xlab("Social Group") + ylab("State 10 Probability") + theme_light() 

ggplot(etamu,aes(x=topic,y=prob)) + geom_point() + geom_errorbar(ymin=prob-std,ymax=prob+std)


#topic contents
topic <- data.table(behav=names(Y)[-(1:ncovcols)], #see bottom for renamed behaviors
                    topic=rep(topicord_eta,rep(d$b,d$K)),
                    iter=rep(1:d$nsave,rep(d$K*d$b,d$nsave)),
                    value=scan(paste0(path,"topicmean.csv")),
                    std=scan(paste0(path,"topicvar.csv"))%>%sqrt)
topic[,type:=etho[str_detect(.BY[[1]],paste0("^",etho$behavior)),type],by=behav]
topic[,behav:=behname]
topic[,type:=factor(type)]

#topic zscores
topicz <- topic[iter>100,.(topic,value=(value-mean(value))/sd(value)),by=.(iter,type,behav)]

topicmeans <- topic[iter>100,.(value=mean(value)),by=.(behav,topic)] %>% 
  dcast(topic ~ behav) %>% 
  merge(nzmu[iter>100,.(prob=mean(prob)),by=topic],by="topic")

#visualize entire topics
topicsummary <- topic[iter > 100 & topic %in% nzmu$topic,
                      .(value=mean(value)),by=.(topic,type,behav)] %>% .[value > 1.05,.(topic,type,value=(value-1)/max(value-1)),by=behav]
#remove redundantish behaviors to declutter
#topicsummary <- topicsummary[!str_detect(behav,"^Groom[PET]")]

ggplot(topicsummary[topic %in% paste0("S",1:10),consolidate(behav,value,0.33),by=c("topic","type")],aes(y=value,x="",fill=type,label=behav)) + geom_point() +
  geom_label_repel(fontface="bold",size=4,force = 12) + coord_cartesian(ylim = c(0,1)) + theme_light(base_size = 16) + 
  scale_fill_brewer(drop=F,palette = "Accent", guide=guide_legend(title="")) + scale_y_continuous(name="Behavior (relative rate)",minor_breaks = NULL) + 
  facet_wrap(~topic,nrow=2) + scale_x_discrete(name="",breaks=NULL)
ggplot(topicsummary[topic %in% c("S10")],aes(y=value,x=topic,fill=type,label=behav)) + geom_point() +
  geom_label_repel(fontface="bold",size=4,force=5) + coord_cartesian(ylim = c(0,1)) + theme_classic() + 
  scale_fill_brewer(drop=F,palette = "Accent", guide=guide_legend(title="")) + 
  ylab("Behavior (relative rate)") + xlab("") #save as 8x5
ggplot(topicsummary[topic %in% c("S5")],aes(y=value,x=topic,fill=type,label=behav)) + geom_point() +
  geom_label_repel() #fontface="bold",size=4,force=3) + coord_cartesian(ylim = c(0,1)) + theme_classic() + 
  #scale_fill_brewer(drop=F,palette = "Accent", guide=guide_legend(title="")) + 
  #ylab("Behavior (relative rate)") + xlab("") #save as 7x5 

#behavior means/vars
topic[iter>100 & str_detect(behav,"displace"),.(mean(value),quantile(value,0.01),quantile(value,0.99)),by=.(topic,behav)] %>%
  ggplot(aes(x=topic,y=V1,color=behav)) + geom_point() + geom_errorbar(aes(ymin=V2,ymax=V3),width=0)
topic[iter>100 & behav=="ScanProx",.(mean(std),quantile(std,0.01),quantile(std,0.99)),by=.(topic,behav)] %>%
  ggplot(aes(x=topic,y=V1)) + geom_point() + geom_errorbar(aes(ymin=V2,ymax=V3),width=0)
topic[iter>100 & behav=="ScanProxInd",.(mean(value),mean(std)),by=.(topic,behav)] %>%
  ggplot(aes(x=topic,y=V1)) + geom_point() + geom_errorbar(aes(ymin=V1-V2,ymax=V1+V2),width=0)

topicz[iter>500,.(mean(value),quantile(value,0.025),quantile(value,0.975)),by=.(type,topic,behav)][
  topic %in% c("S4","S5") & type !="Self-Directed"] %>%
  ggplot(aes(x=V1,y=behav,color=topic)) + geom_point(position = position_dodgev(height=0.2)) + 
  geom_errorbarh(aes(xmin=V2,xmax=V3),width=0,position = position_dodgev(0.2)) +
  facet_grid(type~.,scales="free_y") + scale_x_continuous(minor_breaks = NULL)

#heritability
bl = "S2"
sigdat <- data.table(topic=rep(topicord_eta,rep(d$nsave,d$K)),
                 iter=rep(1:d$nsave,d$K),
                 sigu=fread(paste0(path,"sigma.csv"))[,(d$K+2):(d$K*2+1),with=FALSE] %>% unlist,
                 sigeta=fread(paste0(path,"sigma.csv"))[,2:(d$K+1),with=FALSE] %>% unlist)
sigdat[,h2:=sigdat[topic==bl,sigeta*sigu] + sigeta*sigu,by=topic]
sigdat[,h2:=h2/(h2+sigdat[topic==bl,sigeta] + sigeta),by=topic]
ggplot(sigdat[iter>250 & topic!=bl,.(Heritability=median(h2),lb=quantile(h2,0.17),llb=quantile(h2,0.025),ub=quantile(h2,0.83),uub=quantile(h2,0.975)),by=topic],
       aes(y=Heritability,ymax=ub,ymin=lb,x=topic)) + geom_errorbar(width=0,size=2) + 
  geom_errorbar(aes(ymin=llb,ymax=uub),width=0,size=0.25) + xlab("Behavioral state") +
  geom_point(shape = 21, colour = "black", fill = "white", size = 2, stroke = 2) +
  theme_light() + theme(panel.grid.major.x = element_blank())

#regression coefficients
betadat <- data.table(topic=rep(topicord_eta,rep(d$p,d$K)),
                   coeff=rep(ordered(colnames(X),levels=colnames(X)),d$K*d$nsave),
                   iter=rep(1:d$nsave,rep(d$p*d$K,d$nsave)),
                   beta=scan(paste0(path,"beta.csv")))
#betadat[,beta:=beta-betadat[topic==bl,beta],by=topic]
#plot means and cis 
ggplot(betadat[iter>250 & topic %in% nzmu$topic,
            .(mu=mean(beta),lb=quantile(beta,0.025),ub=quantile(beta,0.975)),by=.(coeff,topic)],
       aes(x=topic,y=mu,ymin=lb,ymax=ub,color=coeff)) + geom_point(position = position_dodge(0.5)) +
  geom_errorbar(width=0,position = position_dodge(0.5)) + geom_hline(yintercept = 0)
ggplot(betadat[iter>250 & topic=="S10",
               .(mu=mean(beta),lb=quantile(beta,0.025),ub=quantile(beta,0.975)),by=.(coeff,topic)],
       aes(x=mu,y=coeff,xmin=lb,xmax=ub)) + geom_point(position = position_dodge(0.5)) +
  geom_errorbarh(width=0) + geom_vline(xintercept = 0) + xlab("Regression coefficient for S10") + ylab("Covariate")
ggplot(betadat[iter>250 & topic !=bl,.(mu=mean(beta),lb=quantile(beta,0.025),ub=quantile(beta,0.975)),by=.(coeff,topic)],
       aes(x=mu,y=coeff,xmin=lb,xmax=ub)) + geom_point() + facet_wrap(~topic) +
  geom_errorbarh(height=0) + geom_vline(xintercept = 0) + scale_y_discrete("Covariate",labels = coefname) + 
  scale_x_continuous("Regression coefficient", minor_breaks = NULL)

#plot predicted functions
#mudatref <- mudat[,.(iter,mu=mu-mudat[topic==bl,mu]),by=topic]

mar <- Xdf[sex=="m",range(rank)]
far <- Xdf[sex=="f",range(rank)]
simdf <- data.table(sex=rep(c("f","m"),c(200,200)),
                    rank=c(seq(far[1],far[2],length.out = 200),seq(mar[1],mar[2],length.out = 200)),
                    age=Xdf[,mean(age)],group="F")
tpdat <- betadat[iter>100,.(FocalID=Xdf$FocalID,tphat=as.vector(X %*% beta)),by=.(iter,topic)]
tpsimhat <- tpdat[,.(ID=1:400,rank=simdf$rank,sex=simdf$sex,etahat=lm(tphat ~ group +sex*poly(age,2) + sex*poly(rank,2),dat=cbind(Xdf,tphat)) %>% predict(newdata=simdf)),
      ,by=.(iter,topic)]
tpsimhat <- merge(tpsimhat,mudat)
tpsimhat[,etahat:=etahat+mu]; tpsimhat[,mu:=NULL]
tpsimhat[,prob:=exp(etahat)/sum(exp(etahat)),by=.(iter,ID)]
tphatsumm <- tpsimhat[,.(rank=rank[1],sex=sex[1],mu=mean(prob),lb=quantile(prob,0.05),ub=quantile(prob,0.95)),by=.(ID,topic)]
ggplot(tphatsumm,aes(x=rank,y=mu,color=sex)) + geom_line() + geom_ribbon(aes(ymin=lb,ymax=ub,fill=sex),alpha=0.25,linetype=0) +
  facet_wrap(~topic,nrow = 2,scales = "free_y") + theme_light() + scale_y_continuous(name="Probability",minor_breaks = NULL) +
  scale_x_continuous(name="Dominance Rank",minor_breaks = NULL,breaks=1:3,labels=c("L","M","H")) + 
  scale_color_discrete("",labels=c("Female","Male")) + scale_fill_discrete("",labels=c("Female","Male"))

#plot R2
#etadatref <- etadat[,.(FocalID,iter,eta=eta-etadat[topic==bl,eta],prob=prob),by=topic]
tpdat2 <- merge(tpdat,mudat[,-"prob"])

xbdat <- copy(tpdat2)
xbdat[,tphat:=tphat+mu]
xbdat[,prob:=exp(tphat)/sum(exp(tphat)),by=.(iter,FocalID)]
r2dat <- merge(xbdat[,-"mu"],etadat,by = c("iter","topic","FocalID"))

r2sum <- r2dat[,.(r2=1-var(prob.y-prob.x)/var(prob.y)),by=c("iter","topic")]
ggplot(r2sum[,.(R2=mean(r2),lb=quantile(r2,0.17),llb=quantile(r2,0.025),ub=quantile(r2,0.83),uub=quantile(r2,0.95)),by=topic],
       aes(y=R2,ymax=ub,ymin=lb,x=topic)) + geom_errorbar(width=0,size=2) + 
  geom_errorbar(aes(ymin=llb,ymax=uub),width=0,size=0.25) + xlab("Behavioral state") +
  geom_point(shape = 21, colour = "black", fill = "white", size = 2, stroke = 2) + 
  theme_light() + theme(panel.grid.major.x = element_blank()) + geom_hline(yintercept = 0) + ylab(expression(R^2))

h2dat <- merge(tpdat2,udat,by=c("iter","topic","FocalID"))
h2dat[,tphat:=tphat+mu]
h2dat[,probxb:=exp(tphat)/sum(exp(tphat)),.(iter,FocalID)]
h2dat[,tphat:=tphat+u]
h2dat[,probxbu:=exp(tphat)/sum(exp(tphat)),.(iter,FocalID)]
r2dat <- merge(h2dat[,-c("mu","u")],etadat,by = c("iter","topic","FocalID"))

r2sum <- r2dat[,.(r2xbu=1-var(prob-probxbu)/var(prob),r2xb=1-var(prob-probxb)/var(prob)),by=c("iter","topic")]
r2sum <- melt(r2sum,measure.vars = c(3,4),value.name = "r2")
ggplot(r2sum[,.(R2=mean(r2),lb=quantile(r2,0.17),llb=quantile(r2,0.025),ub=quantile(r2,0.83),uub=quantile(r2,0.95)),by=.(variable,topic)],
       aes(y=R2,ymax=ub,ymin=lb,x=topic,color=variable)) + geom_errorbar(width=0,size=2,position= position_dodge(0.4)) + 
  geom_errorbar(aes(ymin=llb,ymax=uub),width=0,size=0.25,position= position_dodge(0.4)) + xlab("Behavioral state") +
  geom_point(shape = 21, fill = "white", size = 2, stroke = 2, position= position_dodge(0.4)) + ylab(expression(R^2)) +
  theme_light() + scale_color_discrete("",labels=c("Covariates \nand genetics \n","Covariates only"))

r2sum <- r2dat[,.(h2=1-var(prob-probxbu)/var(prob-probxb)),by=c("iter","topic")]
ggplot(r2sum[,.(R2=mean(h2),lb=quantile(h2,0.17),llb=quantile(h2,0.025),ub=quantile(h2,0.83),uub=quantile(h2,0.95)),by=topic],
       aes(y=R2,ymax=ub,ymin=lb,x=topic)) + geom_errorbar(width=0,size=2) + 
  geom_errorbar(aes(ymin=llb,ymax=uub),width=0,size=0.25) + xlab("Behavioral state") +
  geom_point(shape = 21, colour = "black", fill = "white", size = 2, stroke = 2) + 
  theme_light() + theme(panel.grid.major.x = element_blank()) + geom_hline(yintercept = 0) + ylab(expression("Pseudo H"^2))


#demonstrate an interaction effect
ggplot(Yraw[Vigilnce<median(Vigilnce)],aes(x=`Approach:initiate(focal)`,y=`threat:direct'n(give)`)) + geom_jitter() + geom_smooth(method="lm",formula=y~x)
ggplot(Yraw[Vigilnce>median(Vigilnce)],aes(x=`Approach:initiate(focal)`,y=`threat:direct'n(give)`)) + geom_jitter() + geom_smooth(method="lm",formula=y~x)
Yraw[Vigilnce>median(Vigilnce),lm(`Approach:initiate(focal)`~`threat:direct'n(give)`) %>% summary]
Yraw[Vigilnce<median(Vigilnce),lm(`Approach:initiate(focal)`~`threat:direct'n(give)`) %>% summary]

interdemo_raw <- Yraw[,.(Vigilnce,threat=`threat:direct'n(give)` + `threat:direct'n(receive)`,approach=`Approach:initiate(focal)`+`Approach:initiate(partner)`)]
interdemo <- copy(interdemo_raw)
interdemo[,threat:=cut(threat,c(0,1,2,3,6)-0.5,ordered_result = T)]
interdemo[,Vigilnce:=findInterval(Vigilnce,quantile(Vigilnce,c(0.33,0.66))) %>% ordered(.,0:2)]
idplt <- interdemo[,.(mu=mean(approach),sde=sqrt(var(approach)/length(approach))),by=.(Vigilnce,threat)]
ggplot(idplt,aes(x=threat,y=mu,color=Vigilnce)) + geom_point(position = position_dodge(0.4)) + 
  geom_errorbar(aes(ymin=mu-sde,ymax=mu+sde),width=0,position = position_dodge(0.4)) + theme_light() +
  geom_smooth(data=interdemo,aes(x=as.numeric(threat),y=approach,color=Vigilnce),method="lm",formula = y~x,se=F) + 
  ylab("Approach") + xlab("Threat") + scale_color_discrete(labels=c("low","medium","high")) + scale_x_discrete(labels=c("0","1","2",">2"))

#show that we can capture it in topics
ggplot(topicmeans,aes(y=`Approach:initiate(focal)` + `Approach:initiate(partner)` - 2,x=`threat:direct'n(give)` + `threat:direct'n(give)` - 2,
                      color=findInterval(Vigilnce,quantile(Vigilnce,c(0.33,0.66))) %>% ordered(.,0:2))) + 
  geom_point(size=5) + theme_light() + ylab("Approach") + xlab("Threat") + scale_color_discrete(guide=guide_legend(title="Vigilnce"), labels=c("low","medium","high"))



#rename behaviors
behname <- c("Scratch",
             "SelfGroom",
             "Threat(give)",
             "Threat(receive)",
             "Avoid(give)",
             "Avoid(receive)",
             "Displacement(give)",
             "Displacement(receive)",
             "FearGrimace(give)",
             "FearGrimace(receive)",
             "Submit(give)",
             "Submit(receive)",
             "NonContactAgg(give)",
             "NonContactAgg(receive)",
             "ContactAggression(give)",
             "ContactAggression(receive)",
             "Approach(give)",
             "Approach(receive)",
             "Leave(displace)",
             "Leave(give)",
             "Leave(receive)",
           #  "PassContEnd(give)",
            # "PassContEnd(receive)",
             "AffilVocal(give)",
             "AffilVocal(receive)",
            # "GroomEnd(displace)",
             # "GroomEnd(give)",
             # "GroomEnd(receive)",
             # "GroomPresent(give)",
             # "GroomPresent(receive)",
             "Groom(give)",
             "Groom(receive)",
             "Feed",
             "Travel",
             #"InFeedCorral",
             "PassiveContact",
             "SocialProximity",
             "ProximityGroupSize")

coefname <- c("Group HH","Group R","Male","Age",expression(Age^2),"Rank",expression(Rank^2),"Male x Age",expression("Male x "*Age^2),"Male x Rank",expression("Male x "*Rank^2))
             
consolidate <- function(behav,value,thresh=0.2)
{
  distect <- str_detect(behav,"displace")
  if (any(distect)) {
    disb <- behav[distect]
    behav <- behav[!distect]
    disv <- value[distect]
    value <- value[!distect]
  }
  
  for (i in 1:length(behav))
  {
    stem <- str_split(behav[i],"[\\(\\)]",simplify = T)
    if ( (length(stem) == 1)) next
    
    bi <- str_detect(behav,paste0("^",stem[1],"\\(")) %>% which()
    if (length(bi)==1) next
    
    if (length(bi)>2) browser()
    
    if (abs(diff(value[bi])) < thresh) {
      value[i] <- mean(value[bi])
      value <- value[-bi[2]]
      behav[i] <- paste0(stem[1],"(give/rec)")
      behav <- behav[-bi[2]]
    }
  }
  
  if (any(distect)) {
    value <- c(value,disv)
    behav <- c(behav,disb)
  }
  return(list(value=value,behav=behav))
  
}

##sparse data demo
ptetho$modifier <- NA
Yfig <- collectfocal(fpath,ptetho,stetho,group = c("F","F","HH","R"))
Yfig[,Groom:=GroomGIVE+GroomGET]
Yfig <- Yfig[,-c(2:5,18:20,23:26)]
setnames(Yfig,c("Observation",
                "Scratch",
                "SelfGroom",
                "Threat",
                "Avoid",
                "Displace",
                "FearGrimace",
                "Submit",
                "Noncontact Aggression",
                "Contact Aggression",
                "Approach",
                "Leave",
                "Affiliative Vocalization",
                "Feed",
                "Travel",
                "Groom"))
ggplot(melt(Yfig),aes(x=value)) + geom_histogram() + facet_wrap(~variable,scale="free",nrow=3) +
  scale_x_continuous("Count or Duration",minor_breaks = NULL) + scale_y_continuous("# of Focal Observations",minor_breaks = NULL)
