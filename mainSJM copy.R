rm(list=ls())
timestart<-Sys.time()
options(warn=-1)
setwd("/Users/xuanwang/Dropbox (Harvard University)/Xuan/Surrogate_Longitudinal/simu")
# setwd("/home/xw127/")
source('funsSJM.R')
library(survival);library(MASS);#library("survminer")
library("JMbayes2");#library(rstanarm)
library(dplyr)

#### example:linear 
load('dat.example.rda')
X=dat$X;delta=dat$delta;T=dat$T;obsT=dat$obsT;Y=dat$Y

n=length(dat$X)
dN=(obsT>0); 
tn=max(rowSums(dN));
m=rowSums(dN)

array=as.vector(obsT)
sortIndex = sort(array, index.return=TRUE) 
I=sortIndex$ix[c(which(sortIndex$x>0)-1,n*tn)]
T1=array[I]
J=sapply(1:length(array),function(k){min(which(T1==array[k]))})

n=nrow(obsT);tn=ncol(obsT)
Y_=matrix(0,n,length(T1));dN_=Y_;
for (i in 1:n){
  tryCatch(
    { ID=J[seq(i,(i+n*(m[i]-1)),by=n)];
    Y_[i,ID]=Y[i,1:m[i]];
    dN_[i,ID]=1; }
    , error = function(e){})
}

out=est.linear(dat,re=5)


#### example: nonlinear
load('dat.example.nonlinear.rda')
X=dat$X;delta=dat$delta;T=dat$T;obsT=dat$obsT;Y=dat$Y

n=length(dat$X)
dN=(obsT>0); 
tn=max(rowSums(dN));
m=rowSums(dN)

array=as.vector(obsT)
sortIndex = sort(array, index.return=TRUE) 
I=sortIndex$ix[c(which(sortIndex$x>0)-1,n*tn)]
T1=array[I]
J=sapply(1:length(array),function(k){min(which(T1==array[k]))})

n=nrow(obsT);tn=ncol(obsT)
Y_=matrix(0,n,length(T1));dN_=Y_;
for (i in 1:n){
  tryCatch(
    { ID=J[seq(i,(i+n*(m[i]-1)),by=n)];
    Y_[i,ID]=Y[i,1:m[i]];
    dN_[i,ID]=1; }
    , error = function(e){})
}

gap=.1; nn=length(seq(0,15,gap))
out=est.nonlinear(dat,re=5)
 

timeend<-Sys.time()
runningtime<-timeend-timestart
print(runningtime) 
