source("/home/seth/code/LogisticTopicRegression/parsefocaldata.R")
source("/home/seth/Dropbox/monkeybris/rscript/getAge.R")
basepath <- "~/Dropbox/focaldata_processed/"
fpath <- paste0(basepath,c("F2013/Txtexports_all_processed.csv"))
ptetho <- defaultpoint()[-12]
stetho <- defaultstate()[-4]
Y <- collectfocal(fpath,ptetho,stetho)
Y <- Y[FocalID %in% Y[,length(Observation),by=FocalID][V1>10,FocalID]]
X <- data.table(FocalID=Y$FocalID,sex=getsex(Y$FocalID),age=getage(Y$FocalID,Y$Year[1]),group=getgroup(Y$FocalID))

crushed <- do.call("paste0",X)
obsgroup <- lapply(crushed[!duplicated(X)],function(x) crushed==x)
X <- unique(X)
X <- model.matrix( ~ sex*poly(age,2),X)[,-1]
X <- apply(X,2,function(x) (x-mean(x))/sd(x))
