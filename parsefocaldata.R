defaultstate <- function() list(act=c("Rest","GroomGIVE","GroomGET","Feed","Travel"),
                                corr=c("OutCorral","InCorral"),
                                passetho=c("nopascon","passcont"),
                                inf=c("NoGrmInf","GromInf"))

defaultstate2 <- function() {
  self <- "Self-Directed"
  aff <- "affiliative"
  naff <- "non-affiliative"
  misc <- "misc"
  data.table(behavior=c("Rest","GroomGIVE","GroomGET","Feed","Travel",
                                                  "OutCorral","InCorral","nopascon","passcont","NoGrmInf","GromInf"),
                                       state=c(rep("Activity",5),rep("Corral",2),rep("PassiveContact",2),rep("GroomInfant",2)),
                                       baseline=c(T,F,F,F,F,T,F,T,F,T,F),
                                       type=c(self,aff,aff,self,self,self,self,naff,aff,misc,misc))
}

defaultpoint <- function() data.table(behavior=c("Scratch","SelfGrm","Vigilnce","threat","avoid","displace","FearGrm","Submit",
                                                 "noncontactAgg","contactAgg","Approach","Leave","PsCnTerm","InsptInf",
                                                 "AffVoc","GrmTerm","GrmPrsnt","ChinThrst","SexPresent","Mount"),
                                      modifier=c(NA,NA,NA,"Direction","Winner","Winner","Direction","Direction",
                                                 "Direction","Direction","Initiate","Initiate","Initiate",NA,
                                                 "Direction","Initiate","Direction","Direction",NA,NA))

defaultpoint2 <- function() {
  direct <- "direct'n\\(\\w+\\)"
  winner <- "winner\\?\\(\\w+\\)"
  init <- "initiate\\(\\w+\\)"
  self <- "Self-Directed"
  aff <- "affiliative"
  naff <- "non-affiliative"
  misc <- "misc"
  data.table(behavior=c("Scratch","SelfGrm","Vigilnce","threat","avoid","displace","FearGrm","Submit",
             "noncontactAgg","contactAgg","Approach","Leave","PsCnTerm","InsptInf",
             "AffVoc","GrmTerm","GrmPrsnt","ChinThrst","SexPresent","Mount"),
             modifier=c(NA,NA,NA,direct,winner,winner,direct,direct,direct,direct,
                        init,init,init,NA,direct,init,direct,direct,NA,NA),
             type=c(self,self,self,naff,naff,naff,naff,naff,naff,naff,aff,aff,aff,misc,aff,aff,aff,misc,misc,misc))
               
  # return(list(Scratch=NA, SelfGrm=NA, Vigilnce=NA, threat=direct, avoid=winner, displace=winner, FearGrm=direct, 
  #             Submit=direct, noncontactAgg=direct, contactAgg=direct, Approach=init, Leave=init, PsCnTerm=init, 
  #             InsptInf=NA, AffVoc=direct, GrmTerm=init, GrmPrsnt=direct, ChinThrst=direct, SexPresent=NA, Mount= NA))
}

#identify behaviors that are derivatives of top-level ethogram behaviors
eventslices <- function(X,ptetho) {
  if (is.list(ptetho)) ptetho <- ptetho$behavior
  return(names(X)[sapply(ptetho,function(x) str_detect(names(X),paste0("^",x)) %>% which) %>% unlist])
}

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
    
    out <- melt(sdat)[,variable:=NULL]
    
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
  pt <- pt[,c(1,sapply(names(pt),function(x) str_detect(x,pattern = behav) %>% any) %>% which),with=F]
  setcolorder(pt,c("Observation",eventslices(pt,behav)))
  return(pt)
}

#divide behaviors in ptetho based on matching to modifiers
eventsplit <- function(behav,ptetho,modifier,n) {
  if (!(behav %in% ptetho$behavior) || (ptetho[behavior==behav,is.na(modifier)]))
    return(rep(behav,n))
  else { 
    guh <- str_replace_all(modifier," ","") %>% str_extract(pattern=ptetho[behavior==behav,modifier])
    guh[is.na(guh)] <- rep(behav,sum(is.na(guh)))
    guh[!is.na(guh)] <- paste0(behav,":",guh[!is.na(guh)])
    return(guh)
  }
}
                             
collectfocal <- function(files,ptetho=NULL,stetho=NULL,nsec=10)
{
  if (is.null(stetho))
    stetho <- defaultstate()
  
  if (is.null(ptetho))   
    ptetho <- defaultpoint2()
  
  bdat <- list()
  for (i in 1:length(files)) {
    bdat[[i]] <- fread(files[i])
    bdat[[i]] <- copy(bdat[[i]][!overtime & !BadObs])
  }
  bdat <- do.call(rbind,bdat)
  
  #construct state data
  ss <- sapply(stetho,length)
  Xs <- stetho[,statprep(behavior,dat = bdat,nsec = nsec),by=state] %>% 
    dcast(Observation ~ behavior,value.var="value")

  #construct count data
  
  if (is.list(ptetho)) {
    bdat[,Behavior:=eventsplit(.BY[[1]],ptetho,BehaviorModifier,.N),by=Behavior]
    ptetho <- ptetho$behavior
  }
  
  Xc <- countprep(ptetho,bdat)
  X <- merge(Xc,Xs)
  l <- length(X)
  X <- merge(X,unique(bdat[,cbind(Observation,FocalID,Observer,Year)]),by="Observation")
  setcolorder(X,c(1,(l+1):(l+3),2:l))
  setkey(X,"FocalID")
  return(X)
}