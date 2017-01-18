path <- "/home/seth/analysis/logtopreg/fitFKKR_forpaper/"
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
etadat <- data.table(FocalID=Y[,unique(FocalID)] %>% rep(.,rep(d$K,d$n)),
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

#population topic weight probabilities
ggplot(merge(etasumm,Xdf,by="FocalID"),aes(x=topic,y=prob)) + geom_boxplot() +
  ylab("Probability") + xlab("Behavioral state") + theme_light() + scale_y_continuous(limits = c(0,0.625),expand = c(0,0))

#topic weights against covariates
ggplot(merge(etasumm[topic=="S5"],Xdf,by="FocalID"),aes(y=prob,x=rank,color=sex)) + 
  geom_point() + geom_smooth(method="lm",formula=y ~ poly(x,2),level=0.999) +
  xlab("Age") + ylab("State 3 Probability") + theme_light() + scale_color_discrete(labels=c("Female","Male"),guide=guide_legend(title=""))

ggplot(merge(etasumm[topic=="S10"],Xdf,by="FocalID"),aes(y=prob,x=group)) + 
  geom_boxplot() + xlab("Social Group") + ylab("State 10 Probability") + theme_light() 

ggplot(etamu,aes(x=topic,y=prob)) + geom_point() + geom_errorbar(ymin=prob-std,ymax=prob+std)


#topic contents
topic <- data.table(behav=names(Y)[-(1:5)], #see bottom for renamed behaviors
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
topicsummary <- topicsummary[!str_detect(behav,"^Groom[PET]")]

ggplot(topicsummary[topic %in% c("S1","S3")],aes(y=value,x=topic,fill=type,label=behav)) + geom_point() +
  geom_label_repel(fontface="bold",size=4,force=1) + coord_cartesian(ylim = c(0,1)) + theme_classic() + 
  scale_fill_brewer(drop=F,palette = "Accent", guide=guide_legend(title="")) +
  ylab("Behavior (relative rate)") + xlab("State") #7x5
ggplot(topicsummary[topic %in% c("S10")],aes(y=value,x=topic,fill=type,label=behav)) + geom_point() +
  geom_label_repel(fontface="bold",size=4,force=8) + coord_cartesian(ylim = c(0,1)) + theme_classic() + 
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

topicz[iter>100,.(mean(value),quantile(value,0.01),quantile(value,0.99)),by=.(type,topic,behav)][
  topic %in% c("S4","S5") & type=="Agonistic"] %>%
  ggplot(aes(x=behav,y=V1,color=topic)) + geom_point() + geom_errorbar(aes(ymin=V2,ymax=V3),width=0)

#heritability
h2dat <- data.table(topic=rep(topicord_eta,rep(d$nsave,d$K)),
                 iter=rep(1:d$nsave,rep(d$p*d$K,d$nsave)),
                 h2=fread("~/analysis/logtopreg/fitFKKR/sigma.csv")[,(d$K+2):(d$K*2+1),with=FALSE] %>% unlist)
h2dat[,h2:=h2/(h2+1)]
ggplot(h2dat[iter>250,.(Heritability=median(h2),lb=quantile(h2,0.25),llb=quantile(h2,0.1),ub=quantile(h2,0.75),uub=quantile(h2,0.9)),by=topic],
       aes(y=Heritability,ymax=ub,ymin=lb,x=topic)) + geom_errorbar(width=0,size=2) + 
  geom_errorbar(aes(ymin=llb,ymax=uub),width=0,size=0.25) + xlab("Behavioral state") +
  geom_point(shape = 21, colour = "black", fill = "white", size = 2, stroke = 2) + theme_light()

#regression coefficients
betadat <- data.table(topic=rep(topicord,rep(d$p,d$K)),
                   coeff=rep(1:d$p,d$K*d$nsave) %>% as.character(),
                   iter=rep(1:d$nsave,rep(d$p*d$K,d$nsave)),
                   beta=scan(paste0(path,"beta.csv")))

#plot means and cis 
ggplot(betadat[iter>100 & topic %in% nzmu$topic,
            .(mu=mean(beta),lb=quantile(beta,0.025),ub=quantile(beta,0.975)),by=.(coeff,topic)],
       aes(x=topic,y=mu,ymin=lb,ymax=ub,color=coeff)) + geom_point(position = position_dodge(0.5)) +
  geom_errorbar(width=0,position = position_dodge(0.5)) + geom_hline(yintercept = 0)

#get individual level estimates
tphat <- betadat[iter>100,.(FocalID=Xdf$FocalID,tphat=as.vector(X %*% beta)),by=.(iter,topic)]
tphat <- merge(tphat,mudat[,-"prob",with=FALSE],by=c("topic","iter")) %>% 
  .[,tphat:=tphat+mu] %>% .[,mu:=NULL]
tphat[,prob:=exp(tphat)/sum(exp(tphat)),by=.(iter,FocalID)]
tphatsumm <- tphat[,.(mu=mean(prob),lb=quantile(prob,0.025),ub=quantile(prob,0.975)),by=.(FocalID,topic)]
ggplot(merge(tphatsumm[topic=="T3"],Xdf,by="FocalID"),aes(x=age,y=mu,color=sex)) + 
  geom_line() + geom_ribbon(aes(ymin=lb,ymax=ub,fill=sex),alpha=0.5,linetype=0)


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
ggplot(idplt,aes(x=threat,y=mu,color=Vigilnce)) + geom_point(position = position_dodge(0.2)) + 
  geom_errorbar(aes(ymin=mu-sde,ymax=mu+sde),width=0,position = position_dodge(0.2)) + theme_light() +
  geom_smooth(data=interdemo,aes(x=as.numeric(threat),y=approach,color=Vigilnce),method="lm",formula = y~x,se=F) + 
  ylab("Approach") + xlab("Threat") + scale_color_discrete(labels=c("low","medium","high")) + scale_x_discrete(labels=c("0","1","2",">2"))

#show that we can capture it in topics
ggplot(topicmeans,aes(y=`Approach:initiate(focal)` + `Approach:initiate(partner)` - 2,x=`threat:direct'n(give)` + `threat:direct'n(give)` - 2,
                      color=findInterval(Vigilnce,quantile(Vigilnce,c(0.33,0.66))) %>% ordered(.,0:2))) + 
  geom_point(size=5) + theme_light() + ylab("Approach") + xlab("Threat") + scale_color_discrete(guide=guide_legend(title="Vigilnce"), labels=c("low","medium","high"))



#rename behaviors
behname <- c("Scratch",
             "SelfGroom",
             "Vigilance",
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
             "PassContEnd(give)",
             "PassContEnd(receive)",
             "AffilVocal(give)",
             "AffilVocal(receive)",
            # "GroomEnd(displace)",
             "GroomEnd(give)",
             "GroomEnd(receive)",
             "GroomPresent(give)",
             "GroomPresent(receive)",
             "Groom(give)",
             "Groom(receive)",
             "Feed",
             "Travel",
             "InFeedCorral",
             "PassiveContact",
             "SocialProximity",
             "ProximityGroupSize")
             
