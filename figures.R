##requires R2julia and readtopic

####rename behaviors####
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
             "AffilVocal(give)",
             "AffilVocal(receive)",
             "Groom(give)",
             "Groom(receive)",
             "Feed",
             "Travel",
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

####population topic weight probabilities ####
etasumm <- etadat[,.(prob=mean(prob)),by=.(FocalID,topic)]
ggplot(merge(etasumm,Xdf,by="FocalID"),aes(x=topic,y=prob)) + geom_boxplot() +
  ylab("Probability") + xlab("Behavioral state") + scale_y_continuous(limits = c(0,0.655),expand = c(0,0)) +
  theme_light() + theme(panel.grid.major.x = element_blank())


### topic contents ####
library(ggrepel)
topic[,behav:=behname]
topicsummary <- topic[,.(value=mean(value)),by=.(topic,type,behav)] %>% .[value > 1.075,.(topic,type,value=(value-1)/max(value-1)),by=behav]

ggplot(topicsummary[topic %in% paste0("S",1:10),consolidate(behav,value,0.33),by=c("topic","type")],aes(y=value,x="",fill=type,label=behav)) + geom_point() +
  geom_label_repel(fontface="bold",size=4,force = 12) + coord_cartesian(ylim = c(0,1)) + theme_light(base_size = 16) + 
  scale_fill_brewer(drop=F,palette = "Accent", guide=guide_legend(title="")) + scale_y_continuous(name="Behavior (relative rate)",minor_breaks = NULL) + 
  facet_wrap(~topic,nrow=2) + scale_x_discrete(name="",breaks=NULL)


### regression coefficients ####
bl <- "S2"
betadat_bl <- copy(betadat)
betadat_bl[,beta:=beta-betadat[topic==bl,beta],by=topic]

ggplot(betadat_bl[topic !=bl,.(mu=mean(beta),lb=quantile(beta,0.025),ub=quantile(beta,0.975)),by=.(coeff,topic)],
       aes(x=mu,y=coeff,xmin=lb,xmax=ub)) + geom_point() + facet_wrap(~topic) +
  geom_errorbarh(height=0) + geom_vline(xintercept = 0) + scale_y_discrete("Covariate",labels = coefname) + 
  scale_x_continuous("Regression coefficient", minor_breaks = NULL)

### functions of covariates ####
tpdat <- betadat[.(FocalID=Xdf$FocalID,tphat=as.vector(X %*% beta)),by=.(iter,topic)]
whatplt <- "rank"

if (whatplt=="rank"){
  mar <- Xdf[sex=="m",range(rank)]
  far <- Xdf[sex=="f",range(rank)]
  simdf <- data.table(sex=rep(c("f","m"),c(200,200)),
                      rank=c(seq(far[1],far[2],length.out = 200),seq(mar[1],mar[2],length.out = 200)),
                      age=Xdf[,mean(age)],group="F")
} else if (whatplt=="age") {
  mar <- Xdf[sex=="m",range(age)]
  far <- Xdf[sex=="f",range(age)]
  simdf <- data.table(sex=rep(c("f","m"),c(200,200)),
                      age=c(seq(far[1],far[2],length.out = 200),seq(mar[1],mar[2],length.out = 200)),
                      rank=Xdf[,mean(rank)],group="F")
} else if (whatplt=="group") {
  simdf <- data.table(sex=rep(c("f","m"),c(3,3)),
                      rank=Xdf[,mean(rank)],
                      age=Xdf[,mean(age)],group=c("F","HH","R"))
}

tpsimhat <- tpdat[,.(ID=1:nrow(simdf),rank=simdf$rank,sex=simdf$sex,etahat=lm(tphat ~ group +sex*poly(age,2) + sex*poly(rank,2),dat=cbind(Xdf,tphat)) %>% predict(newdata=simdf)),
                  ,by=.(iter,topic)]
tpsimhat <- merge(tpsimhat,mudat)
tpsimhat[,etahat:=etahat+mu]; tpsimhat[,mu:=NULL]
tpsimhat[,prob:=exp(etahat)/sum(exp(etahat)),by=.(iter,ID)]

if (whatplt=="rank"){
  tphatsumm <- tpsimhat[,.(rank=rank[1],sex=sex[1],mu=mean(prob),lb=quantile(prob,0.05),ub=quantile(prob,0.95)),by=.(ID,topic)]
} else if (whatplt=="age") {
  tphatsumm <- tpsimhat[,.(age=age[1],sex=sex[1],mu=mean(prob),lb=quantile(prob,0.05),ub=quantile(prob,0.95)),by=.(ID,topic)]
} else if (whatplt=="group") {
  tphatsumm <- tpsimhat[,.(group
                           
                           =group[1],sex=sex[1],mu=mean(prob),lb=quantile(prob,0.05),ub=quantile(prob,0.95)),by=.(ID,topic)]
}
ggplot(tphatsumm,aes(x=rank,y=mu,color=sex)) + geom_line() + geom_ribbon(aes(ymin=lb,ymax=ub,fill=sex),alpha=0.25,linetype=0) +
  facet_wrap(~topic,nrow = 2,scales = "free_y") + theme_light() + scale_y_continuous(name="Probability",minor_breaks = NULL) +
  scale_x_continuous(name="Dominance Rank",minor_breaks = NULL,breaks=1:3,labels=c("L","M","H")) + 
  scale_color_discrete("",labels=c("Female","Male")) + scale_fill_discrete("",labels=c("Female","Male"))

## pseudo-h2 ####
tpdat2 <- merge(tpdat,mudat[,-"prob"]) #see previous section for tpdat
h2dat <- merge(tpdat2,udat,by=c("iter","topic","FocalID"))
h2dat[,tphat:=tphat+mu]
h2dat[,probxb:=exp(tphat)/sum(exp(tphat)),.(iter,FocalID)]
h2dat[,tphat:=tphat+u]
h2dat[,probxbu:=exp(tphat)/sum(exp(tphat)),.(iter,FocalID)]
r2dat <- merge(h2dat[,-c("mu","u")],etadat,by = c("iter","topic","FocalID"))

r2sum <- r2dat[,.(h2=1-var(prob-probxbu)/var(prob-probxb)),by=c("iter","topic")]
ggplot(r2sum[,.(R2=mean(h2),lb=quantile(h2,0.17),llb=quantile(h2,0.025),ub=quantile(h2,0.83),uub=quantile(h2,0.95)),by=topic],
       aes(y=R2,ymax=ub,ymin=lb,x=topic)) + geom_errorbar(width=0,size=2) + 
  geom_errorbar(aes(ymin=llb,ymax=uub),width=0,size=0.25) + xlab("Behavioral state") +
  geom_point(shape = 21, colour = "black", fill = "white", size = 2, stroke = 2) + 
  theme_light() + theme(panel.grid.major.x = element_blank()) + geom_hline(yintercept = 0) + ylab(expression("Pseudo H"^2))
