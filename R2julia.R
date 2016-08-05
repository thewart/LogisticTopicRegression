source("/home/seth/code/LogisticTopicRegression/parsefocaldata.R")
source("/home/seth/Dropbox/monkeybris/rscript/getAge.R")
basepath <- "~/Dropbox/focaldata_processed/"
fpath <- paste0(basepath,c("F2013/Txtexports_all_processed.csv"))
ptetho <- defaultpoint2()
#ptetho <- c("AffVoc","Approach","FearGrm","Leave","Submit","Vigilnce","avoid","contactAgg","displace","noncontactAgg","threat")
stetho <- defaultstate()
Y <- collectfocal(fpath,ptetho,stetho)
Y <- Y[FocalID %in% Y[,length(Observation),by=FocalID][V1>10,FocalID]] 
filt <- Y[,lapply(.SD,function(x) mean(x>0)),.SD=eventslices(Y,ptetho)] > 0.005
filt <- colnames(filt)[filt]
Y <- cbind(Y[,.(FocalID,Observation,Year)],Y[,lapply(.SD,function(x) 
  cut(x,c(0,quantile(x,c(0.05,0.2,0.4,0.6,0.8,0.95,1))+0.5) %>% unique,include.lowest = T) %>% as.numeric),
  .SD=filt])
X <- data.table(FocalID=Y$FocalID,sex=getsex(Y$FocalID),age=getage(Y$FocalID,Y$Year[1]),group=getgroup(Y$FocalID))

crushed <- do.call("paste0",X)
obsgroup <- lapply(crushed[!duplicated(X)],function(x) crushed==x)
X <- unique(X)
X <- model.matrix( ~ sex*poly(age,2),X)[,-1]
X <- apply(X,2,function(x) (x-mean(x))/sd(x))

