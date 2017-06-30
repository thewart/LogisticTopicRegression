path <- "/home/seth/analysis/logtopreg/fitFKKR_A_fixed_2/"
d <- scan(paste0(path,"dims.csv")) %>% as.list
names(d) <- c("n","K","p","b","nsave")
etho <- rbind(ptetho[,c(1,3),with=F],stetho[,c(1,4),with=F])
etho <- rbind(etho,data.table(behavior=c("ScanProx","ScanProxInd"),type=rep("Affiliative",2)))
warmup <- 100

##read in phenotypes -- need for topic order####
etadat <- data.table(FocalID=Xdf$FocalID %>% rep(.,rep(d$K,d$n)),
                     topic=1:10,
                     iter=rep(1:d$nsave,rep(d$n*d$K,d$nsave)),
                     eta=scan(paste0(path,"eta.csv")))
etadat <- etadat[iter>warmup]
etadat[,prob:=exp(eta)/sum(exp(eta)),by=.(FocalID,iter)]
etamu <- etadat[,.(prob=mean(prob),std=sd(prob)),by=.(iter,topic)][,.(prob=mean(prob),std=mean(std)),by=topic]
topicord <- etamu[,rank(-prob)] %>% ordered()
levels(topicord) <- paste0("S",levels(topicord))
etadat[,topic:=NULL]
etadat[,topic:=topicord]

## population mean phenotype ####
mudat <- fread(paste0(path,"mu.csv"),header = F)
mudat$iter <-1:d$nsave
mudat <- melt(mudat,id.vars="iter",variable.name = "topic",value.name = "mu")
mudat[,prob:=exp(mu)/ sum(exp(mu)),by=iter]
mudat[,topic:=rep(topicord,rep(d$nsave,d$K))]
mudat <- mudat[iter>warmup]

##  genetic effects ####
udat <- data.table(FocalID=Xdf$FocalID,
                  topic=rep(topicord,rep(d$n,d$K)),
                  iter=rep(1:d$nsave,rep(d$n*d$K,d$nsave)),
                  u=scan(paste0(path,"u.csv")))
udat <- udat[iter>warmup]

## topic contents -- means and vars ####
topic <- data.table(behav=names(Y)[-(1:ncovcols)],
                    topic=rep(topicord,rep(d$b,d$K)),
                    iter=rep(1:d$nsave,rep(d$K*d$b,d$nsave)),
                    value=scan(paste0(path,"topicmean.csv")),
                    std=scan(paste0(path,"topicvar.csv"))%>%sqrt)
topic[,type:=etho[str_detect(.BY[[1]],paste0("^",etho$behavior)),type],by=behav]
topic[,type:=factor(type)]
topic <- topic[iter>warmup]

## variance components -- remember all are scaled by sigeta in model! ####
sigdat <- data.table(topic=rep(topicord,rep(d$nsave,d$K)),
                 iter=rep(1:d$nsave,d$K),
                 sigu=fread(paste0(path,"sigma.csv"))[,(d$K+2):(d$K*2+1),with=FALSE] %>% unlist,
                 sigeta=fread(paste0(path,"sigma.csv"))[,2:(d$K+1),with=FALSE] %>% unlist)
sigdat <- sigdat[iter>warmup]

## regression coefficients -- remember must choose baseline to be interpretable! ####
betadat <- data.table(topic=rep(topicord,rep(d$p,d$K)),
                   coeff=rep(ordered(colnames(X),levels=colnames(X)),d$K*d$nsave),
                   iter=rep(1:d$nsave,rep(d$p*d$K,d$nsave)),
                   beta=scan(paste0(path,"beta.csv")))
betadat <- betadat[iter>warmup]