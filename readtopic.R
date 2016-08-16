path <- "/home/seth/analysis/logtopreg/fit43/"
d <- scan(paste0(path,"dims.csv")) %>% as.list
names(d) <- c("n","K","p","b","nsave")
mu <- fread(paste0(path,"mu.csv"),header = F)
mu$iter <-1:d$nsave
mu <- melt(mu,id.vars="iter",variable.name = "topic")
mu[,value:=exp(value)/ sum(exp(value)),by=iter]
nzmu <- mu[topic %in% mu[iter>101,mean(value),by=topic][,topic[V1>0.001]]]

topic <- data.table(behav=names(Y)[-(1:3)],
                    topic=paste0("V",1:d$K) %>% rep(rep(d$b,d$K)) %>% ordered(.,unique(.)),
                    iter=rep(1:d$nsave,rep(d$K*d$b,d$nsave)),
                    value=scan(paste0(path,"topicmean.csv")))

topicmeans <- topic[iter>100,.(value=mean(value)),by=.(behav,topic)] %>% 
  dcast(topic ~ behav) %>% 
  merge(nzmu[iter>100,.(value=mean(value)),by=topic],by="topic")

#visualize bivariate relationships between behaviors
ggplot(topicmeans,aes(x=Feed,y=Vigilnce,color=topic,size=value)) + geom_point()
ggplot(topicmeans,aes(y=`noncontactAgg:direct'n(give)`,x=`Approach:initiate(focal)`,color=topic,size=value)) + geom_point()

#visualize entire topics
topicsummary <- topic[iter > 100 & topic %in% nzmu$topic,.(topic,value=(value-1)/max(value-1)),by=.(behav,iter)] %>% 
  .[,.(value=mean(value)),by=.(topic,behav)]

#monkey-specific topic probs
eta <- data.table(FocalID=Y[,unique(FocalID)] %>% rep(.,rep(d$K,d$n)),
                  topic=paste0("V",1:d$K) %>% ordered(.,unique(.)),
                  iter=rep(1:d$nsave,rep(d$n*d$K,d$nsave)),
                  value=scan(paste0(path,"eta.csv")))
eta[,value:=exp(value)/sum(exp(value)),by=.(FocalID,iter)]
eta <- merge(eta,Xdf[,.(FocalID,sex)],by="FocalID")
ggplot(eta[iter>100,.(value=mean(value),sex=unique(sex)),by=.(FocalID,topic)],aes(x=topic,y=value,color=sex)) + geom_jitter(width=0.25,alpha=0.5)

#residuals
resid <-eta[iter>100,mean(value),by=.(FocalID,topic)][,.(FocalID,topic,value=V1-mu[iter>100,mean(value),by=topic]$V1)]

#regression coefficients
beta <- data.table(topic=paste0("V",1:d$K) %>% rep(.,rep(d$p,d$K)),
                   coeff=rep(1:d$p,rep(d$K,d$p)),
                   iter=rep(1:d$nsave,rep(d$p*d$K,d$nsave)),
                   beta=scan(paste0(path,"topicmean.csv")))
tphat <- beta[,.(FocalID=Xdf$FocalID,tphat=X %*% beta),by=.(topic,iter)]
tphat <- merge(tphat,mu,by="topic") %>% .[,tphat:=tphat+mu] %>% .[,mu:=NULL]