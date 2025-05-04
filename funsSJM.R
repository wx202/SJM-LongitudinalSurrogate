VTM<-function(vc, dm){
  matrix(vc, ncol=length(vc), nrow=dm, byrow=T)
}


gen.bootstrap.weights=function( n, num.perturb=500){
  sapply(1:num.perturb,function(x) sample(1:n,n,replace=TRUE))
}


simulate <- function(n,beta0,psy1,psy2){
## this function is used to generate the simulation data, where X,Wd and
## Wr are  covariates for the longitudinal model,the survival model and
## the observation process model, resp.
## obsT is a n*max(m_i,i=1,...n) matrix denoting the observation time points;
## T1 is the ordered diferent time points in obsT;
  
tau=15;
gamma0=0;
eta0=0.5;

X=rbinom(n, 1, 0.5)
Wd=X;Wr=X;
epsilon <- log(-log(runif(n)))
D=10*exp(epsilon-X*eta0);  ##with \Lambda_0(t)=t/4;
C <- pmin(runif(n, 5, 25), tau)
T=pmin(D,C); delta=as.numeric(D<=C);

v1=psy1*exp(epsilon);
v2 <- runif(n, 1.9, 2.1)*exp(-psy2*epsilon/2);
ends <- c(0, tau)
Lambda <- matrix(0, n, length(ends) - 1)
for (j in 1:(length(ends)-1)){ Lambda[,j]=v2*exp(X*gamma0) }

out=poisson2(T,Lambda=matrix(20,nrow=n,ncol=1),ends,X*gamma0,v2);
obsT=out$obsT; dN=out$dN; m=out$m
tn=max(m);
array=as.vector(obsT)
sortIndex = sort(array, index.return=TRUE) 
I=sortIndex$ix[c(which(sortIndex$x>0)-1,n*tn)]
T1=array[I]
J=sapply(1:length(array),function(k){min(which(T1==array[k]))})

epsilon1<-matrix(rnorm(n*tn, mean = rep(rnorm(n, mean = 0, sd = 1),tn), sd = rep(0.2*(1:tn),each=n) ),n,tn)
if (set==1){
  Y=t(VTM(X*beta0,tn))+t(VTM(X*beta0,tn))*obsT+epsilon1;
}else if (set==2){
  Y=t(VTM(v1,tn))+t(VTM(X*beta0,tn))+t(VTM(X*beta0,tn))*obsT+epsilon1;
}else if (set==3){
  Y=.2*t(VTM(v1,tn))*obsT+t(VTM(X*beta0,tn))+t(VTM(X*beta0,tn))*obsT+epsilon1;
}else if (set==4){
  Y= t(VTM(rexp(n,1),tn))+.2*t(VTM(v1,tn))*obsT+t(VTM(X*beta0,tn))+t(VTM(X*beta0,tn))*obsT+epsilon1;
}else if (set==44){
  Y= t(VTM(rexp(n,1),tn))+.2*t(VTM(v1,tn))*obsT+t(VTM(X*beta0,tn))+t(VTM(X*beta0,tn))*10*log(obsT+1)+epsilon1;
}else if (set==5){
  b0 = matrix(pmax(rnorm(n*tn,mean = 50, sd = 16), rep(15,n*tn)),n,tn)
  b1=rnorm(n*tn,mean = -2, sd = 2.75)-0*v1
  epsilon1=matrix(rnorm(n*tn,mean = 0, sd = sqrt(0.667*abs(b0 + b1* obsT/4 + +t(VTM(X*beta0*8,tn))*obsT/4))),n, tn)
  Y= b0+b1*obsT/4+t(VTM(X*beta0*8,tn))*obsT/4+epsilon1;
 }else if (set==6){
  b0 = matrix(pmax(rnorm(n*tn,mean = 50, sd = 16), rep(15,n*tn)),n,tn)
  b1=rnorm(n*tn,mean = -2, sd = 2.75)-5*v1
  epsilon1=matrix(rnorm(n*tn,mean = 0, sd = sqrt(0.667*abs(b0 + b1* obsT/4 + +t(VTM(X*beta0*8,tn))*obsT/4))),n, tn)
  Y= b0+b1*obsT/4+t(VTM(X*beta0*8,tn))*obsT/4+epsilon1;
}


out=list('X'=X,'Wd'=Wd,'Wr'=Wr,'obsT'=obsT,'delta'=delta,'T1'=T1,'dN'=dN,'Y'=Y,'I'=I,'J'=J,'m'=m,'T'=T,
         'gamma0'=gamma0,'eta0'=eta0)
}


poisson=function(C=T,Lambda,ends,Xg=X*gamma0,v2){
     n=nrow(Lambda);k=ncol(Lambda);tn=1000;
    obsT=matrix(0,n,tn);
    for (i in 1:n){
    id=1; st=0;  J=1;  delta=1; 
    while(J<=k&&id<=tn){     
      if(delta==1){ u <- runif(1);x=-log(u)/Lambda[i,J]; }
      if(st+x>ends[J+1]){
        if(J+1<=k){
          x=(x+st-ends(J+1))*Lambda(i,J)/Lambda(i,(J+1)); st=ends(J+1); J=J+1;delta=0;
        } else {break}
      }else{ st=st+x 
      if (st>C[i]){ break }
      u <- runif(1)
      if(u<=Lambda[i,J]/Lambda[i,J]){     
      obsT[i,id]=st; id=id+1; }
      delta=1; 
      }
    }
    }
    
    dN=(obsT>0); 
    tn=max(rowSums(dN));
    obsT=obsT[,1:tn];
    dN=dN[,1:tn];
    m=rowSums(dN)
  
    out=list('obsT'=obsT,'dN'=dN,'m'=m)
}


poisson2=function(C=T,Lambda=matrix(10,nrow=n,ncol=1),ends,Xg=X*gamma0,v2){
  n=nrow(Lambda);k=ncol(Lambda);tn=1000;
  obsT=matrix(0,n,tn);
  for (i in 1:n){
    id=1; st=id-1;  J=1;  delta=1; 
    while(J<=k&&id<=tn){     
      if(delta==1){ u <- runif(1);x=-log(u)/Lambda[i,J]; }
      if(st+x>ends[J+1]){
        if(J+1<=k){
          x=(x+st-ends(J+1))*Lambda(i,J)/Lambda(i,(J+1)); st=ends(J+1); J=J+1;delta=0;
          } else {break}
        }else{
          st=st+x 
          if (st>C[i]){ break }
          obsT[i,id]=st; id=id+1; st=id-1
        }
    }
  }
  
  dN=(obsT>0); 
  tn=max(rowSums(dN));
  obsT=obsT[,1:tn];
  dN=dN[,1:tn];
  m=rowSums(dN)
  
  out=list('obsT'=obsT,'dN'=dN,'m'=m)
}


resam<- function(v,dat){
  ## Cox
  dat=data.frame('time'=T,'status'=delta,'X'=X)
  fit=coxph(Surv(time, status) ~ X,data = dat,weights = v) 
  eta_hat=fit$coefficients
  aa=survfit(fit, newdata=data.frame('X'=0) )
  cumhaz.tt=data.frame('time'=aa$time,'cumhaz'=aa$cumhaz) 
  tt=T1
  ind=sapply(1:length(tt), function(k){which.min((tt[k]-cumhaz.tt$time)[(tt[k]-cumhaz.tt$time)>=0])})
  Lambda0_T1=(cumhaz.tt$cumhaz[as.numeric(ind)]); 
  Lambda0_T1[is.na(Lambda0_T1)]=0
  ind=sapply(1:n, function(k){which.min((T[k]-cumhaz.tt$time)[(T[k]-cumhaz.tt$time)>=0])})
  Lambda0_T=(cumhaz.tt$cumhaz[ind])
  
  ## Estimate beta
  dYbar=matrix(NA,n,length(tt))
  Xbar=matrix(NA,n,length(tt))
  dGbar=matrix(NA,n,length(tt))
  Xbar2=matrix(NA,n,length(tt))
  dGbar2=matrix(NA,n,length(tt))
  for (i in 1:n){
    for (k in 1:length(tt)){
      if (T1[k] <= T[i]) {
        psai=(log(Lambda0_T)+eta_hat*X>=log(Lambda0_T1)[k]+eta_hat*X[i])*(eta_hat*X[i]>=eta_hat*X)
        dYbar[i,k]=sum(as.numeric(v)*Y_[,k]*dN_[,k]*psai)/sum(as.numeric(v)*psai)
        dGbar[i,k]=sum(as.numeric(v)*X*dN_[,k]*psai)/sum(as.numeric(v)*psai)
        Xbar[i,k]=sum(as.numeric(v)*X*psai)/sum(as.numeric(v)*psai)
        dGbar2[i,k]=sum(as.numeric(v)*X*dN_[,k]*psai*T1[k])/sum(as.numeric(v)*psai)
        Xbar2[i,k]=sum(as.numeric(v)*X*psai*T1[k])/sum(as.numeric(v)*psai)
      }
    }
  }
  
  nu=rep(NA,length(tt))
  nu2=rep(NA,length(tt))
  sigma=matrix(0,2,2)
  for (k in 1:length(tt)){
    nu[k]=sum(as.numeric(v)*(X-Xbar[,k])*(T>=tt[k])*(Y_[,k]*dN_[,k]-dYbar[,k]),na.rm = T)
    nu2[k]=sum(as.numeric(v)*(X*tt[k]-Xbar2[,k])*(T>=tt[k])*(Y_[,k]*dN_[,k]-dYbar[,k]),na.rm = T)
    tmp1=rbind(as.numeric(v)*as.numeric((X-Xbar[,k])*(T>=tt[k])),
               as.numeric(v)*as.numeric((X*tt[k]-Xbar2[,k])*(T>=tt[k]))); tmp1[is.na(tmp1)]=0
    tmp2=cbind((X*dN_[,k]-dGbar[,k]),(X*tt[k]*dN_[,k]-dGbar2[,k])); tmp2[is.na(tmp2)]=0
    sigma=sigma+tmp1%*%tmp2
  }
  
  est=ginv(sigma)%*%c(sum(nu),sum(nu2))   
  
  out=c((eta_hat),est[1],est[2],ginv(sigma[2,2])%*%c(sum(nu2)))
  
}


resam.nonlinear<- function(v,dat){
  ## Cox
  dat=data.frame('time'=T,'status'=delta,'X'=X)
  fit=coxph(Surv(time, status) ~ X,data = dat,weights = v) 
  eta_hat=fit$coefficients
  aa=survfit(fit, newdata=data.frame('X'=0) )
  cumhaz.tt=data.frame('time'=aa$time,'cumhaz'=aa$cumhaz) 
  tt=T1
  ind=sapply(1:length(tt), function(k){which.min((tt[k]-cumhaz.tt$time)[(tt[k]-cumhaz.tt$time)>=0])})
  Lambda0_T1=(cumhaz.tt$cumhaz[as.numeric(ind)]); 
  Lambda0_T1[is.na(Lambda0_T1)]=0
  ind=sapply(1:n, function(k){which.min((T[k]-cumhaz.tt$time)[(T[k]-cumhaz.tt$time)>=0])})
  Lambda0_T=(cumhaz.tt$cumhaz[ind])
  
  ## Estimate beta
  dYbar=matrix(NA,n,length(tt))
  Xbar=matrix(NA,n,length(tt))
  dGbar=matrix(NA,n,length(tt))
  Xbar2=matrix(NA,n,length(tt))
  dGbar2=matrix(NA,n,length(tt))
  for (i in 1:n){
    for (k in 1:length(tt)){
      if (T1[k] <= T[i]) {
        psai=(log(Lambda0_T)+eta_hat*X>=log(Lambda0_T1)[k]+eta_hat*X[i])*(eta_hat*X[i]>=eta_hat*X)
        dYbar[i,k]=sum(as.numeric(v)*Y_[,k]*dN_[,k]*psai)/sum(as.numeric(v)*psai)
        dGbar[i,k]=sum(as.numeric(v)*X*dN_[,k]*psai)/sum(as.numeric(v)*psai)
        Xbar[i,k]=sum(as.numeric(v)*X*psai)/sum(as.numeric(v)*psai)
        dGbar2[i,k]=sum(as.numeric(v)*X*dN_[,k]*psai*T1[k])/sum(as.numeric(v)*psai)
        Xbar2[i,k]=sum(as.numeric(v)*X*psai*T1[k])/sum(as.numeric(v)*psai)
      }
    }
  }
  
  nu=rep(NA,length(tt))
  nu2=rep(NA,length(tt))
  sigma=matrix(0,2,2)
  for (k in 1:length(tt)){
    nu[k]=sum(as.numeric(v)*(X-Xbar[,k])*(T>=tt[k])*(Y_[,k]*dN_[,k]-dYbar[,k]),na.rm = T)
    nu2[k]=sum(as.numeric(v)*(X*tt[k]-Xbar2[,k])*(T>=tt[k])*(Y_[,k]*dN_[,k]-dYbar[,k]),na.rm = T)
    tmp1=rbind(as.numeric(v)*as.numeric((X-Xbar[,k])*(T>=tt[k])),
               as.numeric(v)*as.numeric((X*tt[k]-Xbar2[,k])*(T>=tt[k]))); tmp1[is.na(tmp1)]=0
    tmp2=cbind((X*dN_[,k]-dGbar[,k]),(X*tt[k]*dN_[,k]-dGbar2[,k])); tmp2[is.na(tmp2)]=0
    sigma=sigma+tmp1%*%tmp2
  }
  
  est=ginv(sigma)%*%c(sum(nu),sum(nu2)) 
  
  ## nonlinear: splines

  phi = data.frame(bs(tt,df=5))
  dYbar=matrix(NA,n,length(tt))
  Xbar=matrix(NA,n,length(tt));dGbar=matrix(NA,n,length(tt))
  Xbar1=matrix(NA,n,length(tt));dGbar1=matrix(NA,n,length(tt))
  Xbar2=matrix(NA,n,length(tt));dGbar2=matrix(NA,n,length(tt))
  Xbar3=matrix(NA,n,length(tt));dGbar3=matrix(NA,n,length(tt))
  Xbar4=matrix(NA,n,length(tt));dGbar4=matrix(NA,n,length(tt))
  Xbar5=matrix(NA,n,length(tt));dGbar5=matrix(NA,n,length(tt))
  for (i in 1:n){
    for (k in 1:length(tt)){
      if (T1[k] <= T[i]) {
        psai=(log(Lambda0_T)+eta_hat*X>=log(Lambda0_T1)[k]+eta_hat*X[i])*(eta_hat*X[i]>=eta_hat*X)
        dYbar[i,k]=sum(as.numeric(v)*Y_[,k]*dN_[,k]*psai)/sum(as.numeric(v)*psai)
        dGbar[i,k]=sum(as.numeric(v)*X*dN_[,k]*psai)/sum(as.numeric(v)*psai)
        Xbar[i,k]=sum(as.numeric(v)*X*psai)/sum(as.numeric(v)*psai)
        dGbar1[i,k]=sum(as.numeric(v)*X*dN_[,k]*psai*phi[k,1])/sum(as.numeric(v)*psai)
        Xbar1[i,k]=sum(as.numeric(v)*X*psai*phi[k,1])/sum(as.numeric(v)*psai)
        dGbar2[i,k]=sum(as.numeric(v)*X*dN_[,k]*psai*phi[k,2])/sum(as.numeric(v)*psai)
        Xbar2[i,k]=sum(as.numeric(v)*X*psai*phi[k,2])/sum(as.numeric(v)*psai)
        dGbar3[i,k]=sum(as.numeric(v)*X*dN_[,k]*psai*phi[k,3])/sum(as.numeric(v)*psai)
        Xbar3[i,k]=sum(as.numeric(v)*X*psai*phi[k,3])/sum(as.numeric(v)*psai)
        dGbar4[i,k]=sum(as.numeric(v)*X*dN_[,k]*psai*phi[k,4])/sum(as.numeric(v)*psai)
        Xbar4[i,k]=sum(as.numeric(v)*X*psai*phi[k,4])/sum(as.numeric(v)*psai)
        dGbar5[i,k]=sum(as.numeric(v)*X*dN_[,k]*psai*phi[k,5])/sum(as.numeric(v)*psai)
        Xbar5[i,k]=sum(as.numeric(v)*X*psai*phi[k,5])/sum(as.numeric(v)*psai)
      }
    }
  }

  nu=rep(NA,length(tt));nu1=rep(NA,length(tt))
  nu2=rep(NA,length(tt));nu3=rep(NA,length(tt))
  nu4=rep(NA,length(tt));nu5=rep(NA,length(tt))
  sigma=matrix(0,6,6)
  for (k in 1:length(tt)){
    nu[k]=sum(as.numeric(v)*(X-Xbar[,k])*(T>=tt[k])*(Y_[,k]*dN_[,k]-dYbar[,k]),na.rm = T)
    nu1[k]=sum(as.numeric(v)*(X*phi[k,1]-Xbar1[,k])*(T>=tt[k])*(Y_[,k]*dN_[,k]-dYbar[,k]),na.rm = T)
    nu2[k]=sum(as.numeric(v)*(X*phi[k,2]-Xbar2[,k])*(T>=tt[k])*(Y_[,k]*dN_[,k]-dYbar[,k]),na.rm = T)
    nu3[k]=sum(as.numeric(v)*(X*phi[k,3]-Xbar3[,k])*(T>=tt[k])*(Y_[,k]*dN_[,k]-dYbar[,k]),na.rm = T)
    nu4[k]=sum(as.numeric(v)*(X*phi[k,4]-Xbar4[,k])*(T>=tt[k])*(Y_[,k]*dN_[,k]-dYbar[,k]),na.rm = T)
    nu5[k]=sum(as.numeric(v)*(X*phi[k,5]-Xbar5[,k])*(T>=tt[k])*(Y_[,k]*dN_[,k]-dYbar[,k]),na.rm = T)
    tmp1=rbind(as.numeric(v)*(X-Xbar[,k])*(T>=tt[k]),as.numeric(v)*(X*phi[k,1]-Xbar1[,k])*(T>=tt[k]),
               as.numeric(v)*(X*phi[k,2]-Xbar2[,k])*(T>=tt[k]),
               as.numeric(v)*(X*phi[k,3]-Xbar3[,k])*(T>=tt[k]),as.numeric(v)*(X*phi[k,4]-Xbar4[,k])*(T>=tt[k]),
               as.numeric(v)*(X*phi[k,5]-Xbar5[,k])*(T>=tt[k]))
    tmp1[is.na(tmp1)]=0
    tmp2=cbind((X*dN_[,k]-dGbar[,k]),(X*phi[k,1]*dN_[,k]-dGbar1[,k]),(X*phi[k,2]*dN_[,k]-dGbar2[,k]),
               (X*phi[k,3]*dN_[,k]-dGbar3[,k]),(X*phi[k,4]*dN_[,k]-dGbar4[,k]),(X*phi[k,5]*dN_[,k]-dGbar5[,k]))
    tmp2[is.na(tmp2)]=0
    sigma=sigma+tmp1%*%tmp2
  }
  est.n=ginv(sigma)%*%c(sum(nu),sum(nu1),sum(nu2),sum(nu3),sum(nu4),sum(nu5))

  ind=sapply(seq(0,quantile(apply(obsT, 1, max),.75),gap), function(k){which.min(abs(tt-k))})
  delta.s.true=(10*log(tt+1)/tt)[ind]
  delta.s=(as.matrix(phi)%*%as.matrix(est.n[-1],ncol=1)/tt)[ind]
  delta.s.hat=matrix(NA,1,length(seq(0,15,gap)))
  delta.s.hat[1:length(delta.s)]=delta.s

  delta.gfr.true=(1+10*log(tt+1))[ind]
  delta.gfr=(est.n[1]+as.matrix(phi)%*%as.matrix(est.n[-1],ncol=1))[ind]
  delta.gfr.hat=matrix(NA,1,length(seq(0,15,gap)))
  delta.gfr.hat[ 1:length(delta.gfr)]=delta.gfr
 
  
  out=c((eta_hat),est[1],est[2],delta.s.hat,delta.gfr.hat)
  
}


resam.nonlinear.boot<- function(index,dat){
  T=T[index];
  delta=delta[index];
  X=X[index]
  obsT=obsT[index,];
  dN=dN[index,];
  Y=Y[index,]
  m=m[index]

  array=as.vector(obsT)
  sortIndex = sort(array, index.return=TRUE)
  I=sortIndex$ix[c(which(sortIndex$x>0)-1,n*tn)]
  T1=array[I]
  # T1=unique(T1)
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
  
  ## Cox
  dat=data.frame('time'=T,'status'=delta,'X'=X)
  fit=coxph(Surv(time, status) ~ X,data = dat) 
  eta_hat=fit$coefficients
  aa=survfit(fit, newdata=data.frame('X'=0) )
  cumhaz.tt=data.frame('time'=aa$time,'cumhaz'=aa$cumhaz) 
  tt=T1
  ind=sapply(1:length(tt), function(k){which.min((tt[k]-cumhaz.tt$time)[(tt[k]-cumhaz.tt$time)>=0])})
  Lambda0_T1=(cumhaz.tt$cumhaz[as.numeric(ind)]); 
  Lambda0_T1[is.na(Lambda0_T1)]=0
  ind=sapply(1:n, function(k){which.min((T[k]-cumhaz.tt$time)[(T[k]-cumhaz.tt$time)>=0])})
  Lambda0_T=(cumhaz.tt$cumhaz[ind])
  
  ## Estimate beta
  dYbar=matrix(NA,n,length(tt))
  Xbar=matrix(NA,n,length(tt))
  dGbar=matrix(NA,n,length(tt))
  Xbar2=matrix(NA,n,length(tt))
  dGbar2=matrix(NA,n,length(tt))
  for (i in 1:n){
    for (k in 1:length(tt)){
      if (T1[k] <= T[i]) {
        psai=(log(Lambda0_T)+eta_hat*X>=log(Lambda0_T1)[k]+eta_hat*X[i])*(eta_hat*X[i]>=eta_hat*X)
        dYbar[i,k]=sum(Y_[,k]*dN_[,k]*psai)/sum(psai)
        dGbar[i,k]=sum(X*dN_[,k]*psai)/sum(psai)
        Xbar[i,k]=sum(X*psai)/sum(psai)
        dGbar2[i,k]=sum(X*dN_[,k]*psai*T1[k])/sum(psai)
        Xbar2[i,k]=sum(X*psai*T1[k])/sum(psai)
      }
    }
  }
  
  nu=rep(NA,length(tt))
  nu2=rep(NA,length(tt))
  sigma=matrix(0,2,2)
  for (k in 1:length(tt)){
    nu[k]=sum((X-Xbar[,k])*(T>=tt[k])*(Y_[,k]*dN_[,k]-dYbar[,k]),na.rm = T)
    nu2[k]=sum((X*tt[k]-Xbar2[,k])*(T>=tt[k])*(Y_[,k]*dN_[,k]-dYbar[,k]),na.rm = T)
    tmp1=rbind(as.numeric((X-Xbar[,k])*(T>=tt[k])),as.numeric((X*tt[k]-Xbar2[,k])*(T>=tt[k]))); tmp1[is.na(tmp1)]=0
    tmp2=cbind((X*dN_[,k]-dGbar[,k]),(X*tt[k]*dN_[,k]-dGbar2[,k])); tmp2[is.na(tmp2)]=0
    sigma=sigma+tmp1%*%tmp2
  }
  
  est=ginv(sigma)%*%c(sum(nu),sum(nu2)) 
  
  ## nonlinear: splines
  phi = data.frame(bs(tt,df=5))
  dYbar=matrix(NA,n,length(tt))
  Xbar=matrix(NA,n,length(tt));dGbar=matrix(NA,n,length(tt))
  Xbar1=matrix(NA,n,length(tt));dGbar1=matrix(NA,n,length(tt))
  Xbar2=matrix(NA,n,length(tt));dGbar2=matrix(NA,n,length(tt))
  Xbar3=matrix(NA,n,length(tt));dGbar3=matrix(NA,n,length(tt))
  Xbar4=matrix(NA,n,length(tt));dGbar4=matrix(NA,n,length(tt))
  Xbar5=matrix(NA,n,length(tt));dGbar5=matrix(NA,n,length(tt))
  for (i in 1:n){
    for (k in 1:length(tt)){
      if (T1[k] <= T[i]) {
        psai=(log(Lambda0_T)+eta_hat*X>=log(Lambda0_T1)[k]+eta_hat*X[i])*(eta_hat*X[i]>=eta_hat*X)
        dYbar[i,k]=sum(Y_[,k]*dN_[,k]*psai)/sum(psai)
        dGbar[i,k]=sum(X*dN_[,k]*psai)/sum(psai)
        Xbar[i,k]=sum(X*psai)/sum(psai)
        dGbar1[i,k]=sum(X*dN_[,k]*psai*phi[k,1])/sum(psai)
        Xbar1[i,k]=sum(X*psai*phi[k,1])/sum(psai)
        dGbar2[i,k]=sum(X*dN_[,k]*psai*phi[k,2])/sum(psai)
        Xbar2[i,k]=sum(X*psai*phi[k,2])/sum(psai)
        dGbar3[i,k]=sum(X*dN_[,k]*psai*phi[k,3])/sum(psai)
        Xbar3[i,k]=sum(X*psai*phi[k,3])/sum(psai)
        dGbar4[i,k]=sum(X*dN_[,k]*psai*phi[k,4])/sum(psai)
        Xbar4[i,k]=sum(X*psai*phi[k,4])/sum(psai)
        dGbar5[i,k]=sum(X*dN_[,k]*psai*phi[k,5])/sum(psai)
        Xbar5[i,k]=sum(X*psai*phi[k,5])/sum(psai)
      }
    }
  }
  
  nu=rep(NA,length(tt));nu1=rep(NA,length(tt))
  nu2=rep(NA,length(tt));nu3=rep(NA,length(tt))
  nu4=rep(NA,length(tt));nu5=rep(NA,length(tt))
  sigma=matrix(0,6,6)
  for (k in 1:length(tt)){
    nu[k]=sum((X-Xbar[,k])*(T>=tt[k])*(Y_[,k]*dN_[,k]-dYbar[,k]),na.rm = T)
    nu1[k]=sum((X*phi[k,1]-Xbar1[,k])*(T>=tt[k])*(Y_[,k]*dN_[,k]-dYbar[,k]),na.rm = T)
    nu2[k]=sum((X*phi[k,2]-Xbar2[,k])*(T>=tt[k])*(Y_[,k]*dN_[,k]-dYbar[,k]),na.rm = T)
    nu3[k]=sum((X*phi[k,3]-Xbar3[,k])*(T>=tt[k])*(Y_[,k]*dN_[,k]-dYbar[,k]),na.rm = T)
    nu4[k]=sum((X*phi[k,4]-Xbar4[,k])*(T>=tt[k])*(Y_[,k]*dN_[,k]-dYbar[,k]),na.rm = T)
    nu5[k]=sum((X*phi[k,5]-Xbar5[,k])*(T>=tt[k])*(Y_[,k]*dN_[,k]-dYbar[,k]),na.rm = T)
    tmp1=rbind((X-Xbar[,k])*(T>=tt[k]),(X*phi[k,1]-Xbar1[,k])*(T>=tt[k]),(X*phi[k,2]-Xbar2[,k])*(T>=tt[k]),
               (X*phi[k,3]-Xbar3[,k])*(T>=tt[k]),(X*phi[k,4]-Xbar4[,k])*(T>=tt[k]),(X*phi[k,5]-Xbar5[,k])*(T>=tt[k]))
    tmp1[is.na(tmp1)]=0
    tmp2=cbind((X*dN_[,k]-dGbar[,k]),(X*phi[k,1]*dN_[,k]-dGbar1[,k]),(X*phi[k,2]*dN_[,k]-dGbar2[,k]),
               (X*phi[k,3]*dN_[,k]-dGbar3[,k]),(X*phi[k,4]*dN_[,k]-dGbar4[,k]),(X*phi[k,5]*dN_[,k]-dGbar5[,k]))
    tmp2[is.na(tmp2)]=0
    sigma=sigma+tmp1%*%tmp2
  }
  est.n=ginv(sigma)%*%c(sum(nu),sum(nu1),sum(nu2),sum(nu3),sum(nu4),sum(nu5)) 
  ind=sapply(seq(0,quantile(apply(obsT, 1, max),.75),gap), function(k){which.min(abs(tt-k))})
  delta.s.true=(10*log(tt+1)/tt)[ind]
  delta.s=(as.matrix(phi)%*%as.matrix(est.n[-1],ncol=1)/tt)[ind]
  delta.s.hat=matrix(NA,1,length(seq(0,15,gap)))
  delta.s.hat[1:length(delta.s)]=delta.s
  delta.gfr.true=(1+10*log(tt+1))[ind]
  delta.gfr=(est.n[1]+as.matrix(phi)%*%as.matrix(est.n[-1],ncol=1))[ind]
  delta.gfr.hat=matrix(NA,1,length(seq(0,15,gap)))
  delta.gfr.hat[1:length(delta.gfr)]=delta.gfr
  
  out=c((eta_hat),est[1],est[2],delta.s.hat,delta.gfr.hat)
  
}


resam.mc<- function(index,dat){
  
T=T[index]
delta=delta[index]
X=X[index]

dN_=dN_[index,]
Y_=Y_[index,]

## Cox
dat=data.frame('time'=T,'status'=delta,'X'=X)
fit=coxph(Surv(time, status) ~ X,data = dat) 
eta_hat=fit$coefficients
aa=survfit(fit, newdata=data.frame('X'=0) )
cumhaz.tt=data.frame('time'=aa$time,'cumhaz'=aa$cumhaz) 
tt=T1
ind=sapply(1:length(tt), function(k){which.min((tt[k]-cumhaz.tt$time)[(tt[k]-cumhaz.tt$time)>=0])})
Lambda0_T1=(cumhaz.tt$cumhaz[as.numeric(ind)]); 
Lambda0_T1[is.na(Lambda0_T1)]=0
ind=sapply(1:n, function(k){which.min((T[k]-cumhaz.tt$time)[(T[k]-cumhaz.tt$time)>=0])})
Lambda0_T=(cumhaz.tt$cumhaz[ind])

## Estimate beta
dYbar=matrix(NA,n,length(tt))
Xbar=matrix(NA,n,length(tt))
dGbar=matrix(NA,n,length(tt))
Xbar2=matrix(NA,n,length(tt))
dGbar2=matrix(NA,n,length(tt))
for (i in 1:n){
  for (k in 1:length(tt)){
    if (T1[k] <= T[i]) {
      psai=(log(Lambda0_T)+eta_hat*X>=log(Lambda0_T1)[k]+eta_hat*X[i])*(eta_hat*X[i]>=eta_hat*X)
      dYbar[i,k]=sum(Y_[,k]*dN_[,k]*psai)/sum(psai)
      dGbar[i,k]=sum(X*dN_[,k]*psai)/sum(psai)
      Xbar[i,k]=sum(X*psai)/sum(psai)
      dGbar2[i,k]=sum(X*dN_[,k]*psai*T1[k])/sum(psai)
      Xbar2[i,k]=sum(X*psai*T1[k])/sum(psai)
    }
  }
}

nu=rep(NA,length(tt))
nu2=rep(NA,length(tt))
sigma=matrix(0,2,2)
for (k in 1:length(tt)){
  nu[k]=sum((X-Xbar[,k])*(T>=tt[k])*(Y_[,k]*dN_[,k]-dYbar[,k]),na.rm = T)
  nu2[k]=sum((X*tt[k]-Xbar2[,k])*(T>=tt[k])*(Y_[,k]*dN_[,k]-dYbar[,k]),na.rm = T)
  tmp1=rbind((X-Xbar[,k])*(T>=tt[k]),(X*tt[k]-Xbar2[,k])*(T>=tt[k])); tmp1[is.na(tmp1)]=0
  tmp2=cbind((X*dN_[,k]-dGbar[,k]),(X*tt[k]*dN_[,k]-dGbar2[,k])); tmp2[is.na(tmp2)]=0
  sigma=sigma+tmp1%*%tmp2
}

est=ginv(sigma)%*%c(sum(nu),sum(nu2)) 

dM=matrix(NA,n,length(tt))
for (k in 1:length(tt)){
  tmp=(T>=tt[k])*(Y_[,k]*dN_[,k]-dYbar[,k]); tmp[is.na(tmp)]=0
  tmp2=(T>=tt[k])*cbind((X*dN_[,k]-dGbar[,k]),(X*tt[k]*dN_[,k]-dGbar2[,k])); tmp2[is.na(tmp2)]=0
  dM[,k]=tmp-tmp2%*%est
}
l_x0=apply( apply(t(VTM((X<=0),length(tt)))*dM,1,cumsum), 2, mean ) 
l_x1=apply( apply(t(VTM((X<=1),length(tt)))*dM,1,cumsum), 2, mean ) 
l_sup=max(abs(l_x0),abs(l_x1)) 

out=c(l_sup)
}


est.linear=function(dat,re){
  ## Cox
  dat=data.frame('time'=T,'status'=delta,'X'=X)
  fit=coxph(Surv(time, status) ~ X,data = dat) 
  eta_hat=fit$coefficients
  aa=survfit(fit, newdata=data.frame('X'=0) )
  cumhaz.tt=data.frame('time'=aa$time,'cumhaz'=aa$cumhaz) 
  tt=T1
  ind=sapply(1:length(tt), function(k){which.min((tt[k]-cumhaz.tt$time)[(tt[k]-cumhaz.tt$time)>=0])})
  Lambda0_T1=(cumhaz.tt$cumhaz[as.numeric(ind)]); 
  Lambda0_T1[is.na(Lambda0_T1)]=0
  ind=sapply(1:n, function(k){which.min((T[k]-cumhaz.tt$time)[(T[k]-cumhaz.tt$time)>=0])})
  Lambda0_T=(cumhaz.tt$cumhaz[ind])
  
  ## Estimate beta
  dYbar=matrix(NA,n,length(tt))
  Xbar=matrix(NA,n,length(tt))
  dGbar=matrix(NA,n,length(tt))
  Xbar2=matrix(NA,n,length(tt))
  dGbar2=matrix(NA,n,length(tt))
  for (i in 1:n){
    for (k in 1:length(tt)){
      if (T1[k] <= T[i]) {
        psai=(log(Lambda0_T)+eta_hat*X>=log(Lambda0_T1)[k]+eta_hat*X[i])*(eta_hat*X[i]>=eta_hat*X)
        dYbar[i,k]=sum(Y_[,k]*dN_[,k]*psai)/sum(psai)
        dGbar[i,k]=sum(X*dN_[,k]*psai)/sum(psai)
        Xbar[i,k]=sum(X*psai)/sum(psai)
        dGbar2[i,k]=sum(X*dN_[,k]*psai*T1[k])/sum(psai)
        Xbar2[i,k]=sum(X*psai*T1[k])/sum(psai)
      }
    }
  }
  
  nu=rep(NA,length(tt))
  nu2=rep(NA,length(tt))
  sigma=matrix(0,2,2)
  for (k in 1:length(tt)){
    nu[k]=sum((X-Xbar[,k])*(T>=tt[k])*(Y_[,k]*dN_[,k]-dYbar[,k]),na.rm = T)
    nu2[k]=sum((X*tt[k]-Xbar2[,k])*(T>=tt[k])*(Y_[,k]*dN_[,k]-dYbar[,k]),na.rm = T)
    tmp1=rbind(as.numeric((X-Xbar[,k])*(T>=tt[k])),as.numeric((X*tt[k]-Xbar2[,k])*(T>=tt[k]))); tmp1[is.na(tmp1)]=0
    tmp2=cbind((X*dN_[,k]-dGbar[,k]),(X*tt[k]*dN_[,k]-dGbar2[,k])); tmp2[is.na(tmp2)]=0
    sigma=sigma+tmp1%*%tmp2
  }
  
  est=c((eta_hat),ginv(sigma)%*%c(sum(nu),sum(nu2)))
  esttom=c((eta_hat),ginv(sigma[2,2])%*%c(sum(nu2)))
  
  #### variance estimation
  v=matrix(rexp(n*re),nrow=n)
  temp=apply(v,2,resam,dat)
  est.se=apply(temp,1,sd)[1:3]
  esttom.se=apply(temp,1,sd)[c(1,4)]
  
  out=c(est,est.se)
}

est.nonlinear=function(dat,re){
  dat=data.frame('time'=T,'status'=delta,'X'=X)
  fit=coxph(Surv(time, status) ~ X,data = dat) 
  eta_hat=fit$coefficients
  aa=survfit(fit, newdata=data.frame('X'=0) )
  cumhaz.tt=data.frame('time'=aa$time,'cumhaz'=aa$cumhaz) 
  tt=T1
  ind=sapply(1:length(tt), function(k){which.min((tt[k]-cumhaz.tt$time)[(tt[k]-cumhaz.tt$time)>=0])})
  Lambda0_T1=(cumhaz.tt$cumhaz[as.numeric(ind)]); 
  Lambda0_T1[is.na(Lambda0_T1)]=0
  ind=sapply(1:n, function(k){which.min((T[k]-cumhaz.tt$time)[(T[k]-cumhaz.tt$time)>=0])})
  Lambda0_T=(cumhaz.tt$cumhaz[ind])
  
  ## nonlinear: splines
  phi = data.frame(bs(tt,df=5))
  dat.spline=data.frame(tt=tt,phi=phi)
  dYbar=matrix(NA,n,length(tt))
  Xbar=matrix(NA,n,length(tt));dGbar=matrix(NA,n,length(tt))
  Xbar1=matrix(NA,n,length(tt));dGbar1=matrix(NA,n,length(tt))
  Xbar2=matrix(NA,n,length(tt));dGbar2=matrix(NA,n,length(tt))
  Xbar3=matrix(NA,n,length(tt));dGbar3=matrix(NA,n,length(tt))
  Xbar4=matrix(NA,n,length(tt));dGbar4=matrix(NA,n,length(tt))
  Xbar5=matrix(NA,n,length(tt));dGbar5=matrix(NA,n,length(tt))
  for (i in 1:n){
    for (k in 1:length(tt)){
      if (T1[k] <= T[i]) {
        psai=(log(Lambda0_T)+eta_hat*X>=log(Lambda0_T1)[k]+eta_hat*X[i])*(eta_hat*X[i]>=eta_hat*X)
        dYbar[i,k]=sum(Y_[,k]*dN_[,k]*psai)/sum(psai)
        dGbar[i,k]=sum(X*dN_[,k]*psai)/sum(psai)
        Xbar[i,k]=sum(X*psai)/sum(psai)
        dGbar1[i,k]=sum(X*dN_[,k]*psai*phi[k,1])/sum(psai)
        Xbar1[i,k]=sum(X*psai*phi[k,1])/sum(psai)
        dGbar2[i,k]=sum(X*dN_[,k]*psai*phi[k,2])/sum(psai)
        Xbar2[i,k]=sum(X*psai*phi[k,2])/sum(psai)
        dGbar3[i,k]=sum(X*dN_[,k]*psai*phi[k,3])/sum(psai)
        Xbar3[i,k]=sum(X*psai*phi[k,3])/sum(psai)
        dGbar4[i,k]=sum(X*dN_[,k]*psai*phi[k,4])/sum(psai)
        Xbar4[i,k]=sum(X*psai*phi[k,4])/sum(psai)
        dGbar5[i,k]=sum(X*dN_[,k]*psai*phi[k,5])/sum(psai)
        Xbar5[i,k]=sum(X*psai*phi[k,5])/sum(psai)
      }
    }
  }

  
  nu=rep(NA,length(tt));nu1=rep(NA,length(tt))
  nu2=rep(NA,length(tt));nu3=rep(NA,length(tt))
  nu4=rep(NA,length(tt));nu5=rep(NA,length(tt))
  sigma=matrix(0,6,6)
  for (k in 1:length(tt)){
    nu[k]=sum((X-Xbar[,k])*(T>=tt[k])*(Y_[,k]*dN_[,k]-dYbar[,k]),na.rm = T)
    nu1[k]=sum((X*phi[k,1]-Xbar1[,k])*(T>=tt[k])*(Y_[,k]*dN_[,k]-dYbar[,k]),na.rm = T)
    nu2[k]=sum((X*phi[k,2]-Xbar2[,k])*(T>=tt[k])*(Y_[,k]*dN_[,k]-dYbar[,k]),na.rm = T)
    nu3[k]=sum((X*phi[k,3]-Xbar3[,k])*(T>=tt[k])*(Y_[,k]*dN_[,k]-dYbar[,k]),na.rm = T)
    nu4[k]=sum((X*phi[k,4]-Xbar4[,k])*(T>=tt[k])*(Y_[,k]*dN_[,k]-dYbar[,k]),na.rm = T)
    nu5[k]=sum((X*phi[k,5]-Xbar5[,k])*(T>=tt[k])*(Y_[,k]*dN_[,k]-dYbar[,k]),na.rm = T)
    tmp1=rbind((X-Xbar[,k])*(T>=tt[k]),(X*phi[k,1]-Xbar1[,k])*(T>=tt[k]),(X*phi[k,2]-Xbar2[,k])*(T>=tt[k]),
               (X*phi[k,3]-Xbar3[,k])*(T>=tt[k]),(X*phi[k,4]-Xbar4[,k])*(T>=tt[k]),(X*phi[k,5]-Xbar5[,k])*(T>=tt[k]))
    tmp1[is.na(tmp1)]=0
    tmp2=cbind((X*dN_[,k]-dGbar[,k]),(X*phi[k,1]*dN_[,k]-dGbar1[,k]),(X*phi[k,2]*dN_[,k]-dGbar2[,k]),
               (X*phi[k,3]*dN_[,k]-dGbar3[,k]),(X*phi[k,4]*dN_[,k]-dGbar4[,k]),(X*phi[k,5]*dN_[,k]-dGbar5[,k]))
    tmp2[is.na(tmp2)]=0
    sigma=sigma+tmp1%*%tmp2
  }
  est.n=ginv(sigma)%*%c(sum(nu),sum(nu1),sum(nu2),sum(nu3),sum(nu4),sum(nu5))
  ind=sapply(seq(0,quantile(apply(obsT, 1, max),.75),gap), function(k){which.min(abs(tt-k))})
  delta.s.true=(10*log(tt+1)/tt)[ind]
  delta.s=(as.matrix(phi)%*%as.matrix(est.n[-1],ncol=1)/tt)[ind]
  delta.s.hat=matrix(NA,1,length(seq(0,15,gap)))
  delta.s.hat[1:length(delta.s)]=delta.s
  
  delta.gfr.true=(1+10*log(tt+1))[ind]
  delta.gfr=(est.n[1]+as.matrix(phi)%*%as.matrix(est.n[-1],ncol=1))[ind]
  delta.gfr.hat=matrix(NA,1,length(seq(0,15,gap)))
  delta.gfr.hat[1:length(delta.gfr)]=delta.gfr
  

  v=matrix(rexp(n*re),nrow=n)
  temp=apply(v,2,resam.nonlinear,dat)
  est.se=apply(temp,1,sd)[1:3]
  delta.s.se=apply(temp,1,sd,na.rm=TRUE)[4:(4+nn-1)]
  delta.s.up=apply(temp,1,quantile,0.975,na.rm=TRUE)[4:(4+nn-1)]
  delta.s.low=apply(temp,1,quantile,0.025,na.rm=TRUE)[4:(4+nn-1)]
  
  out=list('delta.s.hat'=delta.s.hat,'delta.s.se'=delta.s.se,
           'delta.s.up'=delta.s.up,'delta.s.low'=delta.s.low)
}
