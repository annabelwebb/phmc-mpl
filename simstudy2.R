#simulation study 1
#Y ~ Weib(1,3) giving h0(y) = 3y^2

#generate data
gen_data=function(n,a1,a2,ev,b_true,g_true,b_type,g_type,overlap){
  #n is the sample size
  #a1 and a2 are scalars defining interval censoring
  #ev is the event proportion
  #b_true is a vector of the true values of beta (logistic regression)
  #g_true is a vector of the true values of gamma (ph regression)
  #b_type and g_type are vectors specifying the distributions covariates should come from
  #possible values: "u" for uniform, "n" for normal, "b" for binomial
  
  #create covariates
  p=length(b_true)
  q=length(g_true)
  
  if(length(b_type)!=(length(b_true)-1) | length(g_true)!=length(g_type)){
    stop("Incorrect dimensions in covariate specifications.")
  }
  
  if(length(overlap)!=(p+q-1)){
    stop("Incorrect dimensions in covariate overlap specifications.")
  }
  
  #make z
  Z=matrix(0,nrow=n,ncol=(p-1))
  
  z_norms=which(b_type=="n")
  z_unifs=which(b_type=="u")
  z_binoms=which(b_type=="b")
  
  if(any(b_type=="n")){
    for(nc in 1:length(z_norms)){
      z=rnorm(n,60,6)
      ind=z_norms[nc]
      Z[,ind]=z
    }
  }
  if(any(b_type=="u")){
    for(uc in 1:length(z_unifs)){
      z=runif(n)
      ind=z_unifs[uc]
      Z[,ind]=z
    }
  }
  if(any(b_type=="b")){
    for(bc in 1:length(z_binoms)){
      z=rbinom(n,1,0.5)
      ind=z_binoms[bc]
      Z[,ind]=z
    }
  }
  
  #check for duplicates between z and x and then create X
  dups=duplicated(overlap)
  if(any(dups)){
    X=matrix(0,nrow=n,ncol=q)
    
    z.ov=overlap[1:(p-1)]
    x.ov=overlap[-c(1:(p-1))]
    x.dups=dups[-c(1:(p-1))]
    x.dups.ind=which(x.dups)
    
    for(d in 1:sum(x.dups)){
      ind=x.dups.ind[d]
      var=x.ov[ind]
      wh=which(z.ov==var)
      X[,ind]=Z[,wh]
    }
    
    if(sum(x.dups)<q){
      x.new.ind=which(!x.dups)
      for(nv in 1:length(x.new.ind)){
        ind=x.new.ind[nv]
        type=g_type[ind]
        if(type=="n"){
          x=rnorm(n,60,6)
          X[,ind]=x
        }else if(type=="b"){
          x=rbinom(n,1,0.5)
          X[,ind]=x
        }else if(type=="u"){
          x=runif(n)
          X[,ind]=x
        }
      }
    }
  }else{
    X=matrix(0,nrow=n,ncol=q)
    
    x_norms=which(g_type=="n")
    x_unifs=which(g_type=="u")
    x_binoms=which(g_type=="b")
    
    if(any(g_type=="n")){
      for(nc in 1:length(x_norms)){
        x=rnorm(n,60,6)
        ind=x_norms[nc]
        X[,ind]=x
      }
    }
    if(any(g_type=="u")){
      for(uc in 1:length(x_unifs)){
        x=runif(n)
        ind=x_unifs[uc]
        X[,ind]=x
      }
    }
    if(any(g_type=="b")){
      for(bc in 1:length(z_binoms)){
        x=rbinom(n,1,0.5)
        ind=x_binoms[bc]
        X[,ind]=x
      }
    }
  }
  
  mean_X=apply(X,2,mean)
  X_std=X-rep(mean_X,each=n)
  mean_Z=apply(Z,2,mean)
  Z_std=Z-rep(mean_Z,each=n)
  
  #create pi for each i via logistic regression
  pi_i=rep(0,n)
  for(i in 1:n){
    pi_i[i]=exp(b_true[1]+sum(b_true[-c(1)]*Z_std[i,]))/(1+exp(b_true[1]+sum(b_true[-c(1)]*Z_std[i,])))
  }
  
  U_C=runif(n)
  cure_i=as.numeric(U_C<pi_i)
  uncure.n=length(cure_i[cure_i])
  cure.n=length(cure_i[!cure_i])
  
  #create event times for uncured fraction
  U_Y = runif(uncure.n)
  Y_i.u = rep(0, uncure.n)
  X_std_uncure=X_std[cure_i==1,]
  for(i in 1:uncure.n){
    Y_i.u[i] = (-log(U_Y[i])/exp(sum(g_true*X_std_uncure[i,])))^(1/3) #F^-1 method to sample from weibull(1, 3)
  }
  U_E = runif(uncure.n)
  U_L = runif(uncure.n)
  U_R = runif(uncure.n)
  TL_i.u=TR_i.u=rep(0,uncure.n)
  
  #conditions
  cond1=as.numeric(U_E<ev)
  cond2=as.numeric(ev<=U_E & a1*U_L<= Y_i.u & Y_i.u<=a2*U_R)
  cond3=as.numeric(ev<=U_E & a2*U_R < Y_i.u)
  cond4=as.numeric(ev<=U_E & Y_i.u<a1*U_L)
  cond5=as.numeric(ev<=U_E & Y_i.u<a1*U_L)
  
  for(i in 1:uncure.n){
    TL_i.u[i]=(Y_i.u[i]^cond1[i])*((a1*U_L[i])^cond2[i])*((a2*U_R[i])^cond3[i])*(0^cond4[i])
    TR_i.u[i]=(Y_i.u[i]^cond1[i])*((a1*U_L[i])^cond5[i])*((a2*U_R[i])^cond2[i])*(Inf^cond3[i])
  }
  
  #generate censoring times for cured fraction
  TL_i.c = rexp(cure.n,1/2)
  TR_i.c = rep(Inf, cure.n)
  
  #combine populations
  TL_i=TR_i = rep(0,n)
  TL_i[cure_i==1]=TL_i.u
  TL_i[cure_i==0]=TL_i.c
  TR_i[cure_i==1]=TR_i.u
  TR_i[cure_i==0]=TR_i.c
  
  data=matrix(c(TL_i,TR_i,cure_i,Z,X),ncol=(2+p+q),byrow=FALSE)
  data
  
}

#[a1,a2] = [0.9,1.3]
#2 covariates in PH model (one normal, one binomial)
#gamma = [0.75, -0.5]
#2 covariates in logistic model, same as covariates in PH model i.e. X = Z
#three different levels of cure rate i.e. three different sets of beta
#1: beta = [1.3863, 1, -1] i.e. 80% non cured in baseline group
#2: beta = [0, 1, -1] i.e. 50% non cured in baseline group
#3: beta = [-1.3863, 1, 1] i.e. 20% non cured in baseline group
#event proportion in the noncured group given by pi_e = 0%, 25%, 50%
#sample size n = 100, 500, 2000 with respective no. of knots ??? 10 for all??

#generate (to start) 10 samples
library(survival)


beta_save=matrix(rep(0,4*100),ncol=4)
gamma_save=matrix(rep(0,4*100),ncol=4)
censor_save=matrix(rep(0,4*100),ncol=4) #event, right, left, interval

for(r in 1:100){
  sim.data=gen_data(500,0.9,1.3,ev=0.25,c(1,1),c(0.5),c("b"),c("b"),c("z1","x1"))
  for(i in 1:500){
    if(sim.data[i,1]==0 & sim.data[i,2]!=Inf){
      sim.data[i,1]=NA
    }
  }
  
  sim.data.uncure=sim.data[sim.data[,3]==1,]
  sim.surv1=Surv(time=sim.data.uncure[,1],time2=sim.data.uncure[,2],type="interval2")
  censor_save[r,1]=sum(sim.surv1[,3]==1)
  censor_save[r,2]=sum(sim.surv1[,3]==0)
  censor_save[r,3]=sum(sim.surv1[,3]==2)
  censor_save[r,4]=sum(sim.surv1[,3]==3)
  
  sim.surv=Surv(time=sim.data[,1],time2=sim.data[,2],type="interval2")
  test=phmc_mpl(sim.surv~sim.data[,5]+sim.data[,6],pi.formula = ~sim.data[,4],data=as.data.frame(sim.data),phmc_mpl.control(n.knots=c(6,0),maxIter = c(1,1000,1001)))
  
  beta_save[r,1]=test$beta[1]
  beta_save[r,2]=test$beta[2]
  beta_save[r,3]=test$se$se_H[1]
  beta_save[r,4]=test$se$se_H[2]
  
  gamma_save[r,1]=test$gamma[1]
  gamma_save[r,2]=test$se$se_H[3]
  
  
  print(c(r,test$se$se_H[1:3]))
  
}
beta_save
gamma_save
censor_save

mean(beta_save[!is.na(beta_save[,3]),1])
mean(beta_save[!is.na(beta_save[,3]),2])
mean(beta_save[!is.na(beta_save[,3]),3])
mean(beta_save[!is.na(beta_save[,3]),4])
sd(beta_save[!is.na(beta_save[,3]),1])
sd(beta_save[!is.na(beta_save[,3]),2])

mean(gamma_save[!is.na(gamma_save[,3]),1])
mean(gamma_save[!is.na(gamma_save[,3]),2])
mean(gamma_save[!is.na(gamma_save[,3]),3])
mean(gamma_save[!is.na(gamma_save[,3]),4])
sd(gamma_save[!is.na(gamma_save[,3]),1])
sd(gamma_save[!is.na(gamma_save[,3]),2])

mean(censor_save[,1])/362.64
mean(censor_save[,2])/362.64
mean(censor_save[,3])/362.64
mean(censor_save[,4])/362.64


mean(apply(censor_save,1,sum))




while(!is.na(test$se$se_H[1])){
  sim.data=gen_data(500,0.9,1.3,ev=0.25,c(1,1),c(-0.75,0.5),c("b"),c("n","b"),c("z1","x1","z1"))
  for(i in 1:500){
    if(sim.data[i,1]==0 & sim.data[i,2]!=Inf){
      sim.data[i,1]=NA
    }
  }
  
  sim.data.uncure=sim.data[sim.data[,3]==1,]
  sim.surv1=Surv(time=sim.data.uncure[,1],time2=sim.data.uncure[,2],type="interval2")
  
  sim.surv=Surv(time=sim.data[,1],time2=sim.data[,2],type="interval2")
  test=phmc_mpl(sim.surv~sim.data[,5]+sim.data[,6],pi.formula = ~sim.data[,4],data=as.data.frame(sim.data),phmc_mpl.control(n.knots=c(6,0),maxIter = c(1,1000,1001)))
  print(test$se$se_H[1])
}








