source("/home/seth/code/LogisticTopicRegression/parsefocaldata.R")
source("/home/seth/Dropbox/monkeybris/rscript/getAge.R")
basepath <- "~/Dropbox/focaldata_processed/"
fpath <- paste0(basepath,c("F2013/Txtexports_all_processed.csv"))
ptetho <- defaultpoint2()[type!="misc"]
#ptetho <- c("AffVoc","Approach","FearGrm","Leave","Submit","Vigilnce","avoid","contactAgg","displace","noncontactAgg","threat")
stetho <- defaultstate2()[type!="misc"]
Y <- collectfocal(fpath,ptetho,stetho)
Y <- Y[FocalID %in% Y[,length(Observation),by=FocalID][V1>10,FocalID]] 

#remove behaviors that happen too infrequenlty and baseline state behaviors
filt <- Y[,lapply(.SD,function(x) mean(x>0)),.SD=eventslices(names(Y),ptetho)] > 0.005
filt <- colnames(filt)[!filt] %>% c(stetho[baseline==T,behavior]) %>% unique()
Y <- Y[,-filt,with=F]
Yraw <- copy(Y)

#discritize!!
Y <- cbind(Y[,.(FocalID,Observation,Year)],Y[,lapply(.SD,function(x) 
  cut(x,c(0,quantile(x,c(0.05,0.2,0.4,0.6,0.8,0.95,1))+0.5) %>% unique,include.lowest = T) %>% as.numeric),
  .SD=-(1:4)])

#collect covariates
Xdf <- data.table(FocalID=Y$FocalID,sex=getsex(Y$FocalID),age=getage(Y$FocalID,Y$Year[1]),group=getgroup(Y$FocalID))

#generate covariate design matrix
crushed <- do.call("paste0",Xdf)
obsgroup <- lapply(crushed[!duplicated(Xdf)],function(x) crushed==x)
Xdf <- unique(Xdf)
X <- model.matrix( ~ sex*poly(age,2),Xdf)[,-1]
X <- apply(X,2,function(x) (x-mean(x))/sd(x))

