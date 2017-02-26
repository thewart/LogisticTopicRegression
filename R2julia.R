source("/home/seth/code/LogisticTopicRegression/parsefocaldata.R")
source("/home/seth/Dropbox/monkeybris/rscript/getAge.R")
basepath <- "~/Dropbox/focaldata_processed/"
fpath <- paste0(basepath,c("F2013/Txtexports_all_processed.csv",
                          "F2012/Txtexports_all_processed.csv",
                           "HH2014/Txtexports_all_processed.csv",
                           "R2015/Txtexports_all_processed.csv"))
ptetho <- defaultpoint2()[type!="misc" & !(behavior%in%c("Vigilnce","PsCnTerm","GrmTerm","GrmPrsnt"))]
#ptetho <- c("AffVoc","Approach","FearGrm","Leave","Submit","Vigilnce","avoid","contactAgg","displace","noncontactAgg","threat")
stetho <- defaultstate2()[type!="misc" & state!="Corral"]
Y <- collectfocal(fpath,ptetho,stetho,group = c("F","F","HH","R"))
Y <- Y[FocalID %in% Y[,length(Observation),by=FocalID][V1>10,FocalID]]

ncovcols <- 5

#check for animals that switched groups and pick one group with more obs (so far it is always F)
redund <- Y[,length(unique(Group)),by=FocalID][V1>1,FocalID]
Y <- Y[!( (FocalID %in% redund) & (Group!="F"))]

#remove behaviors that happen too infrequenlty and baseline state behaviors
filt <- Y[,lapply(.SD,function(x) mean(x>0)),.SD=eventslices(names(Y),ptetho)] > 0.005
filt <- colnames(filt)[!filt] %>% c(stetho[baseline==T,behavior]) %>% unique()
Y <- Y[,-filt,with=F]
Yraw <- copy(Y)

#discritize!!
Y <- cbind(Y[,1:ncovcols],Y[,lapply(.SD,function(x) 
  cut(x,c(0,quantile(x,c(0.01,0.2,0.4,0.6,0.8,0.99,1))+0.5) %>% unique,include.lowest = T)),
  .SD=-(1:ncovcols)])
Y <- cbind(Y[,1:ncovcols],Y[,lapply(.SD,as.numeric),.SD=-(1:ncovcols)])

#collect covariates
library(xlsx)
drank <- read.xlsx("~/Dropbox/Subjects_attributes, dominance, etc/Dominance Hierarchies/DOMINANCE_ALLSUBJECTS_LONGLIST.xlsx",1) %>% as.data.table()
drank[,ID:=toupper(ID) %>% str_replace_all(.,"_","")]
drank <- drank[YEAR %in% c("2012","2013","2014","2015")]
drank[,ORD_RANK:=ordered(ORD_RANK,levels=c("L","M","H"))]
drank <- drank[,.(rank=as.numeric(ORD_RANK) %>% mean()),by=ID]

Xdf <- data.table(FocalID=Y$FocalID,sex=getsex(Y$FocalID),age=getage(Y$FocalID,Y$Year),group=Y$Group,year=Y$Year)
Xdf[,age:=mean(unique(age)),by=FocalID]
Xdf <- merge(Xdf,drank,by.x = "FocalID",by.y="ID")
#generate covariate design matrix

#pedigree!!!!
load("~/Dropbox/Pedigree and Life-History Data/pedigreeKW2016.RData")
mia <- Xdf[!(FocalID %in% bigped$id),FocalID]
Xdf <- Xdf[!(FocalID %in% mia),]
Y <- Y[!(FocalID %in% mia)]

crushed <- do.call("paste0",Xdf)
obsgroup <- lapply(crushed[!duplicated(Xdf)],function(x) crushed==x)
Xdf <- unique(Xdf)
#Xdf[group %in% c("BB","KK","S","V"),group:="O"]
X <- model.matrix( ~ group +sex*poly(age,2) + sex*poly(rank,2),Xdf)[,-1]
#+ group*poly(rank,2) + group*poly(age,2)
X[,-(1:3)] <- apply(X[,-(1:3)],2,function(x) (x-mean(x))/sd(x))

#A <- 2*kinship2::kinship(bigped$id,dadid=bigped$sire,momid=bigped$dam)
#A <- A[match(Xdf$FocalID,rownames(A)),match(Xdf$FocalID,rownames(A))]
#Z <- chol(A)
