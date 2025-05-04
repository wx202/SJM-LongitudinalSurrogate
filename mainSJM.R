#load libraries and functions
set.seed(100)
timestart<-Sys.time()
options(warn=-1)
source('funsSJM.R')
library(survival);library(MASS);
library("JMbayes2");
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

#run the function, computationally intensive
#re is the number of replications for the bootstrap
#re is set to 5 here for illustration but should be larger in practice
#outputs stacked vector of estimates and standard errors, out.nice is nice display

out=est.linear(dat,re=5)
out.nice = as.data.frame(matrix(out, nrow = 2,byrow = TRUE))
rownames(out.nice) = c("Estimate","SE")
out.nice

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

#run the function, computationally intensive
#re is the number of replications for the bootstrap
#re is set to 5 here for illustration but should be larger in practice
#out.nice is nice display

out=est.nonlinear(dat,re=5)
out.nice = data.frame(as.numeric(out$delta.s.hat),out$delta.s.se, out$delta.s.low,out$delta.s.up)
colnames(out.nice) = c("Estimate","SE","Lower 95%", "Upper 95%")
out.nice

jtimeend<-Sys.time()
runningtime<-timeend-timestart
print(runningtime) 
