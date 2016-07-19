derepeat <- function(behav) {
  n <- length(behav)
  x <- !logical(n)
  if (n>1) {
    for (i in 2:n) 
      if (behav[i-1]==behav[i]) x[i] <- F
  }
  return(x)
}

statvect <- function(time,behav,states) {
  out <- array(0,length(states))
  for (i in 1:length(behav))
    out[match(behav[i],states)] <- time[i]
  return(out)
}

statprep <- function(stat,dat,maxsec=630,nsec=10,ctmc=F) {
  sdat <- dat[Behavior %in% stat,.(Observation,FocalID,Behavior,RelativeEventTime,Duration)]
  sdat <- sdat[RelativeEventTime<maxsec]
  zilch <- dat[!(Observation %in% sdat$Observation),unique(Observation)]
  sdat <- sdat[,if ((length(Behavior)>1) & (RelativeEventTime[2] - RelativeEventTime[1]) < 2.6)
    .(Behavior[-1],RelativeEventTime[-1],Duration[-1])
    else .(Behavior,RelativeEventTime,Duration),by=c("Observation")]
  
  if (!ctmc) {
    sdat[,V3:=pmin(V3,maxsec-V2)]
    sdat <- sdat[,.(time=sum(round(V3/nsec))),by=c("Observation","V1")]
    sdat <- sdat[,statvect(time,V1,stat),by=c("Observation")]
    sdat$behavior <- stat
    
    if (length(zilch) > 0) {
      zilch <- expand.grid(stat,0,zilch)[,3:1] %>% as.data.table()
      setnames(zilch,names(sdat))
      sdat <- rbind(sdat,zilch)
    }
    
    out <- dcast(melt(sdat), Observation ~ behavior)
    
  } else {
    sdat <- sdat[sdat[,derepeat(V1),by=c("Observation")]$V1]
    sdat[,V2:=V2/(max(V2)+1)]
    setnames(sdat,c(3,4),c("Behavior","StartTime"))
    sdat[,V3:=NULL]
    out <- sdat
  }
  return(out)
}

countprep <- function(behav,dat) {
  
  pt <- dat[,length(EventName),by=c("Observation","Behavior")]
  pt <- dcast(pt,Observation ~ Behavior,fill=0)
  pt <- pt[,c(1,which(names(pt) %in% behav)),with=F]
  return(pt)
}

defaultstate <- function() list(act=c("Rest","GroomGIVE","GroomGET","Feed","Travel"),
                                corr=c("OutCorral","InCorral"),
                                passetho=c("nopascon","passcont"),
                                inf=c("NoGrmInf","GromInf"))

defaultpoint <- function() c("Scratch","SelfGrm","Vigilnce","threat","avoid","FearGrm","Submit","noncontactAgg","contactAgg","Approach","Leave","InsptInf",
                             "AffVoc","GrmTerm","GrmPrsnt","ChinThrst")

collectfocal <- function(files,ptetho=NULL,stetho=NULL,nsec=60)
{
  if (is.null(stetho))
    stetho <- defaultstate()
  
  if (is.null(ptetho))   
    ptetho <- defaultpoint()
  
  bdat <- list()
  for (i in 1:length(files)) {
    bdat[[i]] <- fread(files[i])
    bdat[[i]] <- copy(bdat[[i]][!overtime & !BadObs])
  }
  bdat <- do.call(rbind,bdat)
  
  #construct state data
  ss <- sapply(stetho,length)
  sdat <- lapply(stetho,statprep,dat=bdat,nsec=nsec)
  Xs <- sdat[[1]]
  if (length(sdat)>1) for (i in 2:length(sdat)) Xs <- merge(Xs,sdat[[i]])
  #sdat <- unlist(lapply(sdat,as.vector)) %>% matrix(ncol=sum(sapply(statetho,length)))
  #colnames(sdat) <- unlist(statetho)
  
  #construct count data
  Xc <- countprep(ptetho,bdat)
  X <- merge(Xc,Xs)
  l <- length(X)
  X <- merge(X,unique(bdat[,cbind(Observation,FocalID,Observer,Year)]),by="Observation")
  setcolorder(X,c(1,(l+1):(l+3),2:l))
  setcolorder(X,c(names(X)[1:4],ptetho,unlist(stetho)))
  setkey(X,"FocalID")
  return(X)
}
