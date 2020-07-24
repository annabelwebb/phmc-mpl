#simulation study 1
#Y ~ Weib(1,3) giving h0(y) = 3y^2

#generate data ------
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
    Y_i.u[i] = (-log(U_Y[i])/exp(sum(g_true*X_std_uncure[i])))^(1/3) #F^-1 method to sample from weibull(1, 3)
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
  
  data=data.frame(TL_i,TR_i,cure_i,Z,X)
  colnames(data)=c("TL","TR","cure","Z","X")
  data
  
}









# 1. n=2000,event=0.5,a=(0.9,1.3),beta=[1,1] -----
beta_save1=matrix(rep(0,4*1000),ncol=4)
gamma_save1=matrix(rep(0,2*1000),ncol=2)
censor_save1=matrix(rep(0,4*1000),ncol=4) #event, right, left, interval
converge_save1=matrix(rep(0,2*1000),ncol=2) #converged yes/no, number of iterations
theta_save1=matrix(rep(0,11*1000),ncol=11)
knots_save1=matrix(rep(0,10*1000),ncol=10)

for(r in 914:1000){
  sim.data=gen_data(2000,0.9,1.3,ev=0.5,c(1,1),c(0.5),c("b"),c("b"),c("z1","x1"))
  for(i in 1:2000){
    if(sim.data[i,1]==0 & sim.data[i,2]!=Inf){
      sim.data[i,1]=NA
    }
  }
  
  sim.data.uncure=sim.data[sim.data[,3]==1,]
  sim.surv1=Surv(time=sim.data.uncure[,1],time2=sim.data.uncure[,2],type="interval2")
  censor_save1[r,1]=sum(sim.surv1[,3]==1)
  censor_save1[r,2]=sum(sim.surv1[,3]==0)
  censor_save1[r,3]=sum(sim.surv1[,3]==2)
  censor_save1[r,4]=sum(sim.surv1[,3]==3)
  
  sim.surv=Surv(time=sim.data[,1],time2=sim.data[,2],type="interval2")
  test=phmc_mpl(sim.surv~sim.data$X,pi.formula = ~sim.data$Z,data=sim.data,phmc_mpl.control(n.knots=c(8,0),maxIter = c(1,5000,5001)))
  
  beta_save1[r,1]=test$beta[1]
  beta_save1[r,2]=test$beta[2]
  beta_save1[r,3]=test$se$se_H[1]
  beta_save1[r,4]=test$se$se_H[2]
  
  gamma_save1[r,1]=test$gamma[1]
  gamma_save1[r,2]=test$se$se_H[3]
  
  converge_save1[r,1]=test$convergence[1]
  converge_save1[r,2]=test$convergence[3]
  
  theta_save1[r,]=test$theta
  
  knots_save1[r,]=test$knots$Alpha
  
  print(c(r,test$se$se_H[1] ,test$convergence[1]))
  
}







# 2. n=2000,event=0.5,a=(0.9,1.3),beta=[0,1] ------

beta_save2=matrix(rep(0,4*1000),ncol=4)
gamma_save2=matrix(rep(0,2*1000),ncol=2)
censor_save2=matrix(rep(0,4*1000),ncol=4) #event, right, left, interval
converge_save2=matrix(rep(0,2*1000),ncol=2) #converged yes/no, number of iterations
theta_save2=matrix(rep(0,7*1000),ncol=7)
knots_save2=matrix(rep(0,6*1000),ncol=6)


for(r in 855:1000){
  sim.data=gen_data(2000,0.9,1.3,ev=0.5,c(0,1),c(0.5),c("b"),c("b"),c("z1","x1"))
  for(i in 1:2000){
    if(sim.data[i,1]==0 & sim.data[i,2]!=Inf){
      sim.data[i,1]=NA
    }
  }
  
  sim.data.uncure=sim.data[sim.data[,3]==1,]
  sim.surv1=Surv(time=sim.data.uncure[,1],time2=sim.data.uncure[,2],type="interval2")
  censor_save2[r,1]=sum(sim.surv1[,3]==1)
  censor_save2[r,2]=sum(sim.surv1[,3]==0)
  censor_save2[r,3]=sum(sim.surv1[,3]==2)
  censor_save2[r,4]=sum(sim.surv1[,3]==3)
  
  sim.surv=Surv(time=sim.data[,1],time2=sim.data[,2],type="interval2")
  test=phmc_mpl(sim.surv~sim.data$X,pi.formula = ~sim.data$Z,data=sim.data,phmc_mpl.control(n.knots=c(4,0),maxIter = c(1,10000,10001)))
  
  beta_save2[r,1]=test$beta[1]
  beta_save2[r,2]=test$beta[2]
  beta_save2[r,3]=test$se$se_H[1]
  beta_save2[r,4]=test$se$se_H[2]
  
  gamma_save2[r,1]=test$gamma[1]
  gamma_save2[r,2]=test$se$se_H[3]
  
  converge_save2[r,1]=test$convergence[1]
  converge_save2[r,2]=test$convergence[3]
  
  theta_save2[r,]=test$theta
  
  knots_save2[r,]=test$knots$Alpha
  
  print(c(r,test$se$se_H[1] ,test$convergence[1]))
  
}




# 3. n=2000,event=0.5,a=(0.9,1.3),beta=[-1,1]------

# 4. n=500,event=0.5,a=(0.9,1.3),beta=[1,1]--------


beta_save4=matrix(rep(0,4*100),ncol=4)
gamma_save4=matrix(rep(0,2*100),ncol=2)
censor_save4=matrix(rep(0,4*100),ncol=4) #event, right, left, interval
converge_save4=matrix(rep(0,2*100),ncol=2) #converged yes/no, number of iterations
theta_save4=matrix(rep(0,9*100),ncol=9)
knots_save4=matrix(rep(0,8*100),ncol=8)

for(r in 1:100){
  sim.data=gen_data(500,0.9,1.3,ev=0.5,c(1,1),c(0.5),c("b"),c("b"),c("z1","x1"))
  for(i in 1:500){
    if(sim.data[i,1]==0 & sim.data[i,2]!=Inf){
      sim.data[i,1]=NA
    }
  }
  
  sim.data.uncure=sim.data[sim.data[,3]==1,]
  sim.surv1=Surv(time=sim.data.uncure[,1],time2=sim.data.uncure[,2],type="interval2")
  censor_save4[r,1]=sum(sim.surv1[,3]==1)
  censor_save4[r,2]=sum(sim.surv1[,3]==0)
  censor_save4[r,3]=sum(sim.surv1[,3]==2)
  censor_save4[r,4]=sum(sim.surv1[,3]==3)
  
  sim.surv=Surv(time=sim.data[,1],time2=sim.data[,2],type="interval2")
  test=phmc_mpl(sim.surv~sim.data$X,pi.formula = ~sim.data$Z,data=sim.data,phmc_mpl.control(n.knots=c(6,0),maxIter = c(1,10000,10001)))
  
  beta_save4[r,1]=test$beta[1]
  beta_save4[r,2]=test$beta[2]
  beta_save4[r,3]=test$se$se_H[1]
  beta_save4[r,4]=test$se$se_H[2]
  
  gamma_save4[r,1]=test$gamma[1]
  gamma_save4[r,2]=test$se$se_H[3]
  
  converge_save4[r,1]=test$convergence[1]
  converge_save4[r,2]=test$convergence[3]
  
  theta_save4[r,]=test$theta
  
  knots_save4[r,]=test$knots$Alpha
  
  print(c(4,r,test$se$se_H[1] ,test$convergence[1]))
  
}

dfbeta_save4=data.frame(beta_save4)
dfgamma_save4=data.frame(gamma_save4)
dfcensor_save4=data.frame(censor_save4) #event, right, left, interval
dfconverge_save4=data.frame(converge_save4) #converged yes/no, number of iterations
dftheta_save4=data.frame(theta_save4)
dfknots_save4=data.frame(knots_save4)

write.table(dfbeta_save4,"beta4.csv",sep=",",dec = ".",row.names = F)
write.table(dfgamma_save4,"gamma4.csv",sep = ",",dec=".",row.names = F)
write.table(dfcensor_save4,"censor4.csv",sep = ",",dec=".",row.names = F)
write.table(dfconverge_save4,"converge4.csv",sep=",",dec=".",row.names = F)
write.table(dftheta_save4,"theta4.csv",sep = ",",dec=".",row.names = F)
write.table(dfknots_save4,"knots4.csv",sep = ",",dec=".",row.names = F)




# 5. n==500,event=0.5,a=(0.9,1.3),beta=[0,1]------


beta_save5=matrix(rep(0,4*100),ncol=4)
gamma_save5=matrix(rep(0,2*100),ncol=2)
censor_save5=matrix(rep(0,4*100),ncol=4) #event, right, left, interval
converge_save5=matrix(rep(0,2*100),ncol=2) #converged yes/no, number of iterations
theta_save5=matrix(rep(0,6*100),ncol=6)
knots_save5=matrix(rep(0,5*100),ncol=5)

for(r in 1:100){
  sim.data=gen_data(500,0.9,1.3,ev=0.5,c(0,1),c(0.5),c("b"),c("b"),c("z1","x1"))
  for(i in 1:500){
    if(sim.data[i,1]==0 & sim.data[i,2]!=Inf){
      sim.data[i,1]=NA
    }
  }
  
  sim.data.uncure=sim.data[sim.data[,3]==1,]
  sim.surv1=Surv(time=sim.data.uncure[,1],time2=sim.data.uncure[,2],type="interval2")
  censor_save5[r,1]=sum(sim.surv1[,3]==1)
  censor_save5[r,2]=sum(sim.surv1[,3]==0)
  censor_save5[r,3]=sum(sim.surv1[,3]==2)
  censor_save5[r,4]=sum(sim.surv1[,3]==3)
  
  sim.surv=Surv(time=sim.data[,1],time2=sim.data[,2],type="interval2")
  test=phmc_mpl(sim.surv~sim.data$X,pi.formula = ~sim.data$Z,data=sim.data,phmc_mpl.control(n.knots=c(3,0),maxIter = c(1,10000,10001)))
  
  beta_save5[r,1]=test$beta[1]
  beta_save5[r,2]=test$beta[2]
  beta_save5[r,3]=test$se$se_H[1]
  beta_save5[r,4]=test$se$se_H[2]
  
  gamma_save5[r,1]=test$gamma[1]
  gamma_save5[r,2]=test$se$se_H[3]
  
  converge_save5[r,1]=test$convergence[1]
  converge_save5[r,2]=test$convergence[3]
  
  theta_save5[r,]=test$theta
  
  knots_save5[r,]=test$knots$Alpha
  
  print(c(5, r,test$se$se_H[1] ,test$convergence[1]))
  
}

dfbeta_save5=data.frame(beta_save5)
dfgamma_save5=data.frame(gamma_save5)
dfcensor_save5=data.frame(censor_save5) #event, right, left, interval
dfconverge_save5=data.frame(converge_save5) #converged yes/no, number of iterations
dftheta_save5=data.frame(theta_save5)
dfknots_save5=data.frame(knots_save5)

write.table(dfbeta_save5,"beta5.csv",sep=",",dec = ".",row.names = F)
write.table(dfgamma_save5,"gamma5.csv",sep = ",",dec=".",row.names = F)
write.table(dfcensor_save5,"censor5.csv",sep = ",",dec=".",row.names = F)
write.table(dfconverge_save5,"converge5.csv",sep=",",dec=".",row.names = F)
write.table(dftheta_save5,"theta5.csv",sep = ",",dec=".",row.names = F)
write.table(dfknots_save5,"knots5.csv",sep = ",",dec=".",row.names = F)


# 6. n=500,event=0.5,a=(0.9,1.3),beta=[-1,1]------



beta_save6=matrix(rep(0,4*100),ncol=4)
gamma_save6=matrix(rep(0,2*100),ncol=2)
censor_save6=matrix(rep(0,4*100),ncol=4) #event, right, left, interval
converge_save6=matrix(rep(0,2*100),ncol=2) #converged yes/no, number of iterations
theta_save6=matrix(rep(0,4*100),ncol=4)
knots_save6=matrix(rep(0,3*100),ncol=3)

for(r in 74:100){
  sim.data=gen_data(500,0.9,1.3,ev=0.5,c(-1,1),c(0.5),c("b"),c("b"),c("z1","x1"))
  for(i in 1:500){
    if(sim.data[i,1]==0 & sim.data[i,2]!=Inf){
      sim.data[i,1]=NA
    }
  }
  
  sim.data.uncure=sim.data[sim.data[,3]==1,]
  sim.surv1=Surv(time=sim.data.uncure[,1],time2=sim.data.uncure[,2],type="interval2")
  censor_save6[r,1]=sum(sim.surv1[,3]==1)
  censor_save6[r,2]=sum(sim.surv1[,3]==0)
  censor_save6[r,3]=sum(sim.surv1[,3]==2)
  censor_save6[r,4]=sum(sim.surv1[,3]==3)
  
  sim.surv=Surv(time=sim.data[,1],time2=sim.data[,2],type="interval2")
  test=phmc_mpl(sim.surv~sim.data$X,pi.formula = ~sim.data$Z,data=sim.data,phmc_mpl.control(n.knots=c(1,0),maxIter = c(1,10000,10001)))
  
  beta_save6[r,1]=test$beta[1]
  beta_save6[r,2]=test$beta[2]
  beta_save6[r,3]=test$se$se_H[1]
  beta_save6[r,4]=test$se$se_H[2]
  
  gamma_save6[r,1]=test$gamma[1]
  gamma_save6[r,2]=test$se$se_H[3]
  
  converge_save6[r,1]=test$convergence[1]
  converge_save6[r,2]=test$convergence[3]
  
  theta_save6[r,]=test$theta
  
  knots_save6[r,]=test$knots$Alpha
  
  print(c(6, r,test$se$se_H[1] ,test$convergence[1]))
  
}

dfbeta_save6=data.frame(beta_save6)
dfgamma_save6=data.frame(gamma_save6)
dfcensor_save6=data.frame(censor_save6) #event, right, left, interval
dfconverge_save6=data.frame(converge_save6) #converged yes/no, number of iterations
dftheta_save6=data.frame(theta_save6)
dfknots_save6=data.frame(knots_save6)

write.table(dfbeta_save6,"beta6.csv",sep=",",dec = ".",row.names = F)
write.table(dfgamma_save6,"gamma6.csv",sep = ",",dec=".",row.names = F)
write.table(dfcensor_save6,"censor6.csv",sep = ",",dec=".",row.names = F)
write.table(dfconverge_save6,"converge6.csv",sep=",",dec=".",row.names = F)
write.table(dftheta_save6,"theta6.csv",sep = ",",dec=".",row.names = F)
write.table(dfknots_save6,"knots6.csv",sep = ",",dec=".",row.names = F)
 
# 13. n=500,event=0.25,a=(0.9,1.3),beta=[1,1]-----


beta_save13=matrix(rep(0,4*100),ncol=4)
gamma_save13=matrix(rep(0,2*100),ncol=2)
censor_save13=matrix(rep(0,4*100),ncol=4) #event, right, left, interval
converge_save13=matrix(rep(0,2*100),ncol=2) #converged yes/no, number of iterations
theta_save13=matrix(rep(0,6*100),ncol=6)
knots_save13=matrix(rep(0,5*100),ncol=5)

for(r in 1:100){
  sim.data=gen_data(500,0.9,1.3,ev=0.25,c(1,1),c(0.5),c("b"),c("b"),c("z1","x1"))
  for(i in 1:500){
    if(sim.data[i,1]==0 & sim.data[i,2]!=Inf){
      sim.data[i,1]=NA
    }
  }
  
  sim.data.uncure=sim.data[sim.data[,3]==1,]
  sim.surv1=Surv(time=sim.data.uncure[,1],time2=sim.data.uncure[,2],type="interval2")
  censor_save13[r,1]=sum(sim.surv1[,3]==1)
  censor_save13[r,2]=sum(sim.surv1[,3]==0)
  censor_save13[r,3]=sum(sim.surv1[,3]==2)
  censor_save13[r,4]=sum(sim.surv1[,3]==3)
  
  sim.surv=Surv(time=sim.data[,1],time2=sim.data[,2],type="interval2")
  test=phmc_mpl(sim.surv~sim.data$X,pi.formula = ~sim.data$Z,data=sim.data,phmc_mpl.control(n.knots=c(3,0),maxIter = c(1,10000,10001)))
  
  beta_save13[r,1]=test$beta[1]
  beta_save13[r,2]=test$beta[2]
  beta_save13[r,3]=test$se$se_H[1]
  beta_save13[r,4]=test$se$se_H[2]
  
  gamma_save13[r,1]=test$gamma[1]
  gamma_save13[r,2]=test$se$se_H[3]
  
  converge_save13[r,1]=test$convergence[1]
  converge_save13[r,2]=test$convergence[3]
  
  theta_save13[r,]=test$theta
  
  knots_save13[r,]=test$knots$Alpha
  
  print(c(13,r,test$se$se_H[1] ,test$convergence[1]))
  
}

dfbeta_save13=data.frame(beta_save13)
dfgamma_save13=data.frame(gamma_save13)
dfcensor_save13=data.frame(censor_save13) #event, right, left, interval
dfconverge_save13=data.frame(converge_save13) #converged yes/no, number of iterations
dftheta_save13=data.frame(theta_save13)
dfknots_save13=data.frame(knots_save13)

write.table(dfbeta_save13,"beta13.csv",sep=",",dec = ".",row.names = F)
write.table(dfgamma_save13,"gamma13.csv",sep = ",",dec=".",row.names = F)
write.table(dfcensor_save13,"censor13.csv",sep = ",",dec=".",row.names = F)
write.table(dfconverge_save13,"converge13.csv",sep=",",dec=".",row.names = F)
write.table(dftheta_save13,"theta13.csv",sep = ",",dec=".",row.names = F)
write.table(dfknots_save13,"knots13.csv",sep = ",",dec=".",row.names = F)


# 14. n=500,event=0.25,a=(0.9,1.3),beta=[0,1]-------


beta_save14=matrix(rep(0,4*100),ncol=4)
gamma_save14=matrix(rep(0,2*100),ncol=2)
censor_save14=matrix(rep(0,4*100),ncol=4) #event, right, left, interval
converge_save14=matrix(rep(0,2*100),ncol=2) #converged yes/no, number of iterations
theta_save14=matrix(rep(0,5*100),ncol=5)
knots_save14=matrix(rep(0,4*100),ncol=4)

for(r in 1:100){
  sim.data=gen_data(500,0.9,1.3,ev=0.25,c(0,1),c(0.5),c("b"),c("b"),c("z1","x1"))
  for(i in 1:500){
    if(sim.data[i,1]==0 & sim.data[i,2]!=Inf){
      sim.data[i,1]=NA
    }
  }
  
  sim.data.uncure=sim.data[sim.data[,3]==1,]
  sim.surv1=Surv(time=sim.data.uncure[,1],time2=sim.data.uncure[,2],type="interval2")
  censor_save14[r,1]=sum(sim.surv1[,3]==1)
  censor_save14[r,2]=sum(sim.surv1[,3]==0)
  censor_save14[r,3]=sum(sim.surv1[,3]==2)
  censor_save14[r,4]=sum(sim.surv1[,3]==3)
  
  sim.surv=Surv(time=sim.data[,1],time2=sim.data[,2],type="interval2")
  test=phmc_mpl(sim.surv~sim.data$X,pi.formula = ~sim.data$Z,data=sim.data,phmc_mpl.control(n.knots=c(2,0),maxIter = c(1,10000,10001)))
  
  beta_save14[r,1]=test$beta[1]
  beta_save14[r,2]=test$beta[2]
  beta_save14[r,3]=test$se$se_H[1]
  beta_save14[r,4]=test$se$se_H[2]
  
  gamma_save14[r,1]=test$gamma[1]
  gamma_save14[r,2]=test$se$se_H[3]
  
  converge_save14[r,1]=test$convergence[1]
  converge_save14[r,2]=test$convergence[3]
  
  theta_save14[r,]=test$theta
  
  knots_save14[r,]=test$knots$Alpha
  
  print(c(14,r,test$se$se_H[1] ,test$convergence[1]))
  
}

dfbeta_save14=data.frame(beta_save14)
dfgamma_save14=data.frame(gamma_save14)
dfcensor_save14=data.frame(censor_save14) #event, right, left, interval
dfconverge_save14=data.frame(converge_save14) #converged yes/no, number of iterations
dftheta_save14=data.frame(theta_save14)
dfknots_save14=data.frame(knots_save14)

write.table(dfbeta_save14,"beta14.csv",sep=",",dec = ".",row.names = F)
write.table(dfgamma_save14,"gamma14.csv",sep = ",",dec=".",row.names = F)
write.table(dfcensor_save14,"censor14.csv",sep = ",",dec=".",row.names = F)
write.table(dfconverge_save14,"converge14.csv",sep=",",dec=".",row.names = F)
write.table(dftheta_save14,"theta14.csv",sep = ",",dec=".",row.names = F)
write.table(dfknots_save14,"knots14.csv",sep = ",",dec=".",row.names = F)


# 15. n=500,event=0.25,a=(0.9,1.3),beta=[-1,1]------


beta_save15=matrix(rep(0,4*100),ncol=4)
gamma_save15=matrix(rep(0,2*100),ncol=2)
censor_save15=matrix(rep(0,4*100),ncol=4) #event, right, left, interval
converge_save15=matrix(rep(0,2*100),ncol=2) #converged yes/no, number of iterations
theta_save15=matrix(rep(0,4*100),ncol=4)
knots_save15=matrix(rep(0,3*100),ncol=3)

for(r in 1:100){
  sim.data=gen_data(500,0.9,1.3,ev=0.25,c(-1,1),c(0.5),c("b"),c("b"),c("z1","x1"))
  for(i in 1:500){
    if(sim.data[i,1]==0 & sim.data[i,2]!=Inf){
      sim.data[i,1]=NA
    }
  }
  
  sim.data.uncure=sim.data[sim.data[,3]==1,]
  sim.surv1=Surv(time=sim.data.uncure[,1],time2=sim.data.uncure[,2],type="interval2")
  censor_save15[r,1]=sum(sim.surv1[,3]==1)
  censor_save15[r,2]=sum(sim.surv1[,3]==0)
  censor_save15[r,3]=sum(sim.surv1[,3]==2)
  censor_save15[r,4]=sum(sim.surv1[,3]==3)
  
  sim.surv=Surv(time=sim.data[,1],time2=sim.data[,2],type="interval2")
  test=phmc_mpl(sim.surv~sim.data$X,pi.formula = ~sim.data$Z,data=sim.data,phmc_mpl.control(n.knots=c(1,0),maxIter = c(1,10000,10001)))
  
  beta_save15[r,1]=test$beta[1]
  beta_save15[r,2]=test$beta[2]
  beta_save15[r,3]=test$se$se_H[1]
  beta_save15[r,4]=test$se$se_H[2]
  
  gamma_save15[r,1]=test$gamma[1]
  gamma_save15[r,2]=test$se$se_H[3]
  
  converge_save15[r,1]=test$convergence[1]
  converge_save15[r,2]=test$convergence[3]
  
  theta_save15[r,]=test$theta
  
  knots_save15[r,]=test$knots$Alpha
  
  print(c(15,r,test$se$se_H[1] ,test$convergence[1]))
  
}

dfbeta_save15=data.frame(beta_save15)
dfgamma_save15=data.frame(gamma_save15)
dfcensor_save15=data.frame(censor_save15) #event, right, left, interval
dfconverge_save15=data.frame(converge_save15) #converged yes/no, number of iterations
dftheta_save15=data.frame(theta_save15)
dfknots_save15=data.frame(knots_save15)

write.table(dfbeta_save15,"beta15.csv",sep=",",dec = ".",row.names = F)
write.table(dfgamma_save15,"gamma15.csv",sep = ",",dec=".",row.names = F)
write.table(dfcensor_save15,"censor15.csv",sep = ",",dec=".",row.names = F)
write.table(dfconverge_save15,"converge15.csv",sep=",",dec=".",row.names = F)
write.table(dftheta_save15,"theta15.csv",sep = ",",dec=".",row.names = F)
write.table(dfknots_save15,"knots15.csv",sep = ",",dec=".",row.names = F)





# 19. n=2000,event=0,a=(0.9,1.3),beta=[1,1]---------




beta_save19=matrix(rep(0,4*100),ncol=4)
gamma_save19=matrix(rep(0,2*100),ncol=2)
censor_save19=matrix(rep(0,4*100),ncol=4) #event, right, left, interval
converge_save19=matrix(rep(0,2*100),ncol=2) #converged yes/no, number of iterations
theta_save19=matrix(rep(0,4*100),ncol=4)
knots_save19=matrix(rep(0,3*100),ncol=3)

for(r in 1:100){
  sim.data=gen_data(2000,0.5,1.1,ev=0,c(1,1),c(0.5),c("b"),c("b"),c("z1","x1"))
  for(i in 1:2000){
    if(sim.data[i,1]==0 & sim.data[i,2]!=Inf){
      sim.data[i,1]=NA
    }
  }
  
  sim.data.uncure=sim.data[sim.data[,3]==1,]
  sim.surv1=Surv(time=sim.data.uncure[,1],time2=sim.data.uncure[,2],type="interval2")
  censor_save19[r,1]=sum(sim.surv1[,3]==1)
  censor_save19[r,2]=sum(sim.surv1[,3]==0)
  censor_save19[r,3]=sum(sim.surv1[,3]==2)
  censor_save19[r,4]=sum(sim.surv1[,3]==3)
  
  sim.surv=Surv(time=sim.data[,1],time2=sim.data[,2],type="interval2")
  test=phmc_mpl(sim.surv~sim.data$X,pi.formula = ~sim.data$Z,data=sim.data,phmc_mpl.control(n.knots=c(1,0),maxIter = c(1,1000,1001)))
  
  beta_save19[r,1]=test$beta[1]
  beta_save19[r,2]=test$beta[2]
  beta_save19[r,3]=test$se$se_H[1]
  beta_save19[r,4]=test$se$se_H[2]
  
  gamma_save19[r,1]=test$gamma[1]
  gamma_save19[r,2]=test$se$se_H[3]
  
  converge_save19[r,1]=test$convergence[1]
  converge_save19[r,2]=test$convergence[3]
  
  theta_save19[r,]=test$theta
  
  knots_save19[r,]=test$knots$Alpha
  
  print(c(19,r,test$beta[1],test$se$se_H[1] ,test$convergence[1]))
  
}

dfbeta_save19=data.frame(beta_save19)
dfgamma_save19=data.frame(gamma_save19)
dfcensor_save19=data.frame(censor_save19) #event, right, left, interval
dfconverge_save19=data.frame(converge_save19) #converged yes/no, number of iterations
dftheta_save19=data.frame(theta_save19)
dfknots_save19=data.frame(knots_save19)

write.table(dfbeta_save19,"beta19.csv",sep=",",dec = ".",row.names = F)
write.table(dfgamma_save19,"gamma19.csv",sep = ",",dec=".",row.names = F)
write.table(dfcensor_save19,"censor19.csv",sep = ",",dec=".",row.names = F)
write.table(dfconverge_save19,"converge19.csv",sep=",",dec=".",row.names = F)
write.table(dftheta_save19,"theta19.csv",sep = ",",dec=".",row.names = F)
write.table(dfknots_save19,"knots19.csv",sep = ",",dec=".",row.names = F)






# 22. n=500,event=0,a=(0.9,1.3),beta=[1,1] -------

beta_save22=matrix(rep(0,4*100),ncol=4)
gamma_save22=matrix(rep(0,2*100),ncol=2)
censor_save22=matrix(rep(0,4*100),ncol=4) #event, right, left, interval
converge_save22=matrix(rep(0,2*100),ncol=2) #converged yes/no, number of iterations
theta_save22=matrix(rep(0,4*100),ncol=4)
knots_save22=matrix(rep(0,3*100),ncol=3)

for(r in 1:100){
  sim.data=gen_data(500,0.9,1.3,ev=0,c(1,1),c(0.5),c("b"),c("b"),c("z1","x1"))
  for(i in 1:500){
    if(sim.data[i,1]==0 & sim.data[i,2]!=Inf){
      sim.data[i,1]=NA
    }
  }
  
  sim.data.uncure=sim.data[sim.data[,3]==1,]
  sim.surv1=Surv(time=sim.data.uncure[,1],time2=sim.data.uncure[,2],type="interval2")
  censor_save22[r,1]=sum(sim.surv1[,3]==1)
  censor_save22[r,2]=sum(sim.surv1[,3]==0)
  censor_save22[r,3]=sum(sim.surv1[,3]==2)
  censor_save22[r,4]=sum(sim.surv1[,3]==3)
  
  sim.surv=Surv(time=sim.data[,1],time2=sim.data[,2],type="interval2")
  test=phmc_mpl(sim.surv~sim.data$X,pi.formula = ~sim.data$Z,data=sim.data,phmc_mpl.control(n.knots=c(1,0),maxIter = c(1,20000,20001)))
  
  beta_save22[r,1]=test$beta[1]
  beta_save22[r,2]=test$beta[2]
  beta_save22[r,3]=test$se$se_H[1]
  beta_save22[r,4]=test$se$se_H[2]
  
  gamma_save22[r,1]=test$gamma[1]
  gamma_save22[r,2]=test$se$se_H[3]
  
  converge_save22[r,1]=test$convergence[1]
  converge_save22[r,2]=test$convergence[3]
  
  theta_save22[r,]=test$theta
  
  knots_save22[r,]=test$knots$Alpha
  
  print(c(22,r,test$se$se_H[1] ,test$convergence[1]))
  
}

dfbeta_save22=data.frame(beta_save22)
dfgamma_save22=data.frame(gamma_save22)
dfcensor_save22=data.frame(censor_save22) #event, right, left, interval
dfconverge_save22=data.frame(converge_save22) #converged yes/no, number of iterations
dftheta_save22=data.frame(theta_save22)
dfknots_save22=data.frame(knots_save22)

write.table(dfbeta_save22,"beta22.csv",sep=",",dec = ".",row.names = F)
write.table(dfgamma_save22,"gamma22.csv",sep = ",",dec=".",row.names = F)
write.table(dfcensor_save22,"censor22.csv",sep = ",",dec=".",row.names = F)
write.table(dfconverge_save22,"converge22.csv",sep=",",dec=".",row.names = F)
write.table(dftheta_save22,"theta22.csv",sep = ",",dec=".",row.names = F)
write.table(dfknots_save22,"knots22.csv",sep = ",",dec=".",row.names = F)




# 23. n=500,event=0,a=(0.9,1.3),beta=[0,1] -------

beta_save23=matrix(rep(0,4*100),ncol=4)
gamma_save23=matrix(rep(0,2*100),ncol=2)
censor_save23=matrix(rep(0,4*100),ncol=4) #event, right, left, interval
converge_save23=matrix(rep(0,2*100),ncol=2) #converged yes/no, number of iterations
theta_save23=matrix(rep(0,4*100),ncol=4)
knots_save23=matrix(rep(0,3*100),ncol=3)

for(r in 1:100){
  sim.data=gen_data(500,0.9,1.3,ev=0,c(0,1),c(0.5),c("b"),c("b"),c("z1","x1"))
  for(i in 1:500){
    if(sim.data[i,1]==0 & sim.data[i,2]!=Inf){
      sim.data[i,1]=NA
    }
  }
  
  sim.data.uncure=sim.data[sim.data[,3]==1,]
  sim.surv1=Surv(time=sim.data.uncure[,1],time2=sim.data.uncure[,2],type="interval2")
  censor_save23[r,1]=sum(sim.surv1[,3]==1)
  censor_save23[r,2]=sum(sim.surv1[,3]==0)
  censor_save23[r,3]=sum(sim.surv1[,3]==2)
  censor_save23[r,4]=sum(sim.surv1[,3]==3)
  
  sim.surv=Surv(time=sim.data[,1],time2=sim.data[,2],type="interval2")
  test=phmc_mpl(sim.surv~sim.data$X,pi.formula = ~sim.data$Z,data=sim.data,phmc_mpl.control(n.knots=c(1,0),maxIter = c(1,20000,20001)))
  
  beta_save23[r,1]=test$beta[1]
  beta_save23[r,2]=test$beta[2]
  beta_save23[r,3]=test$se$se_H[1]
  beta_save23[r,4]=test$se$se_H[2]
  
  gamma_save23[r,1]=test$gamma[1]
  gamma_save23[r,2]=test$se$se_H[3]
  
  converge_save23[r,1]=test$convergence[1]
  converge_save23[r,2]=test$convergence[3]
  
  theta_save23[r,]=test$theta
  
  knots_save23[r,]=test$knots$Alpha
  
  print(c(23,r,test$se$se_H[1] ,test$convergence[1]))
  
}

dfbeta_save23=data.frame(beta_save23)
dfgamma_save23=data.frame(gamma_save23)
dfcensor_save23=data.frame(censor_save23) #event, right, left, interval
dfconverge_save23=data.frame(converge_save23) #converged yes/no, number of iterations
dftheta_save23=data.frame(theta_save23)
dfknots_save23=data.frame(knots_save23)

write.table(dfbeta_save23,"beta23.csv",sep=",",dec = ".",row.names = F)
write.table(dfgamma_save23,"gamma23.csv",sep = ",",dec=".",row.names = F)
write.table(dfcensor_save23,"censor23.csv",sep = ",",dec=".",row.names = F)
write.table(dfconverge_save23,"converge23.csv",sep=",",dec=".",row.names = F)
write.table(dftheta_save23,"theta23.csv",sep = ",",dec=".",row.names = F)
write.table(dfknots_save23,"knots23.csv",sep = ",",dec=".",row.names = F)


# 24. n=500,event=0,a=(0.9,1.3),beta=[-1,1] -------

beta_save24=matrix(rep(0,4*100),ncol=4)
gamma_save24=matrix(rep(0,2*100),ncol=2)
censor_save24=matrix(rep(0,4*100),ncol=4) #event, right, left, interval
converge_save24=matrix(rep(0,2*100),ncol=2) #converged yes/no, number of iterations
theta_save24=matrix(rep(0,4*100),ncol=4)
knots_save24=matrix(rep(0,3*100),ncol=3)

for(r in 1:100){
  sim.data=gen_data(500,0.9,1.3,ev=0,c(-1,1),c(0.5),c("b"),c("b"),c("z1","x1"))
  for(i in 1:500){
    if(sim.data[i,1]==0 & sim.data[i,2]!=Inf){
      sim.data[i,1]=NA
    }
  }
  
  sim.data.uncure=sim.data[sim.data[,3]==1,]
  sim.surv1=Surv(time=sim.data.uncure[,1],time2=sim.data.uncure[,2],type="interval2")
  censor_save24[r,1]=sum(sim.surv1[,3]==1)
  censor_save24[r,2]=sum(sim.surv1[,3]==0)
  censor_save24[r,3]=sum(sim.surv1[,3]==2)
  censor_save24[r,4]=sum(sim.surv1[,3]==3)
  
  sim.surv=Surv(time=sim.data[,1],time2=sim.data[,2],type="interval2")
  test=phmc_mpl(sim.surv~sim.data$X,pi.formula = ~sim.data$Z,data=sim.data,phmc_mpl.control(n.knots=c(1,0),maxIter = c(1,20000,20001)))
  
  beta_save24[r,1]=test$beta[1]
  beta_save24[r,2]=test$beta[2]
  beta_save24[r,3]=test$se$se_H[1]
  beta_save24[r,4]=test$se$se_H[2]
  
  gamma_save24[r,1]=test$gamma[1]
  gamma_save24[r,2]=test$se$se_H[3]
  
  converge_save24[r,1]=test$convergence[1]
  converge_save24[r,2]=test$convergence[3]
  
  theta_save24[r,]=test$theta
  
  knots_save24[r,]=test$knots$Alpha
  
  print(c(24,r,test$se$se_H[1] ,test$convergence[1]))
  
}

dfbeta_save24=data.frame(beta_save24)
dfgamma_save24=data.frame(gamma_save24)
dfcensor_save24=data.frame(censor_save24) #event, right, left, interval
dfconverge_save24=data.frame(converge_save24) #converged yes/no, number of iterations
dftheta_save24=data.frame(theta_save24)
dfknots_save24=data.frame(knots_save24)

write.table(dfbeta_save24,"beta24.csv",sep=",",dec = ".",row.names = F)
write.table(dfgamma_save24,"gamma24.csv",sep = ",",dec=".",row.names = F)
write.table(dfcensor_save24,"censor24.csv",sep = ",",dec=".",row.names = F)
write.table(dfconverge_save24,"converge24.csv",sep=",",dec=".",row.names = F)
write.table(dftheta_save24,"theta24.csv",sep = ",",dec=".",row.names = F)
write.table(dfknots_save24,"knots24.csv",sep = ",",dec=".",row.names = F)





# 25. n=100,event=0,a=(0.9,1.3),beta=[1,1] -------





# 26. n=100,event=0,a=(0.9,1.3),beta=[0,1] -------
# 27. n=100,event=0,a=(0.9,1.3),beta=[-1,1] -------









# 31. n=500.event=0.8,a=(0.9,1.3),beta=[1,1]-------

beta_save31=matrix(rep(0,4*100),ncol=4)
gamma_save31=matrix(rep(0,2*100),ncol=2)
censor_save31=matrix(rep(0,4*100),ncol=4) #event, right, left, interval
converge_save31=matrix(rep(0,2*100),ncol=2) #converged yes/no, number of iterations
theta_save31=matrix(rep(0,9*100),ncol=9)
knots_save31=matrix(rep(0,8*100),ncol=8)

for(r in 1:100){
  sim.data=gen_data(500,0.9,1.3,ev=0.8,c(1,1),c(0.5),c("b"),c("b"),c("z1","x1"))
  for(i in 1:500){
    if(sim.data[i,1]==0 & sim.data[i,2]!=Inf){
      sim.data[i,1]=NA
    }
  }
  
  sim.data.uncure=sim.data[sim.data[,3]==1,]
  sim.surv1=Surv(time=sim.data.uncure[,1],time2=sim.data.uncure[,2],type="interval2")
  censor_save31[r,1]=sum(sim.surv1[,3]==1)
  censor_save31[r,2]=sum(sim.surv1[,3]==0)
  censor_save31[r,3]=sum(sim.surv1[,3]==2)
  censor_save31[r,4]=sum(sim.surv1[,3]==3)
  
  sim.surv=Surv(time=sim.data[,1],time2=sim.data[,2],type="interval2")
  test=phmc_mpl(sim.surv~sim.data$X,pi.formula = ~sim.data$Z,data=sim.data,phmc_mpl.control(n.knots=c(6,0),maxIter = c(1,10000,10001)))
  
  beta_save31[r,1]=test$beta[1]
  beta_save31[r,2]=test$beta[2]
  beta_save31[r,3]=test$se$se_H[1]
  beta_save31[r,4]=test$se$se_H[2]
  
  gamma_save31[r,1]=test$gamma[1]
  gamma_save31[r,2]=test$se$se_H[3]
  
  converge_save31[r,1]=test$convergence[1]
  converge_save31[r,2]=test$convergence[3]
  
  theta_save31[r,]=test$theta
  
  knots_save31[r,]=test$knots$Alpha
  
  print(c(31,r,test$se$se_H[1] ,test$convergence[1]))
  
}

dfbeta_save31=data.frame(beta_save31)
dfgamma_save31=data.frame(gamma_save31)
dfcensor_save31=data.frame(censor_save31) #event, right, left, interval
dfconverge_save31=data.frame(converge_save31) #converged yes/no, number of iterations
dftheta_save31=data.frame(theta_save31)
dfknots_save31=data.frame(knots_save31)

write.table(dfbeta_save31,"beta31.csv",sep=",",dec = ".",row.names = F)
write.table(dfgamma_save31,"gamma31.csv",sep = ",",dec=".",row.names = F)
write.table(dfcensor_save31,"censor31.csv",sep = ",",dec=".",row.names = F)
write.table(dfconverge_save31,"converge31.csv",sep=",",dec=".",row.names = F)
write.table(dftheta_save31,"theta31.csv",sep = ",",dec=".",row.names = F)
write.table(dfknots_save31,"knots31.csv",sep = ",",dec=".",row.names = F)

# 32. n=500,event=0.8,a=(0.9,1.3),beta=[0,1]-------


beta_save32=matrix(rep(0,4*100),ncol=4)
gamma_save32=matrix(rep(0,2*100),ncol=2)
censor_save32=matrix(rep(0,4*100),ncol=4) #event, right, left, interval
converge_save32=matrix(rep(0,2*100),ncol=2) #converged yes/no, number of iterations
theta_save32=matrix(rep(0,7*100),ncol=7)
knots_save32=matrix(rep(0,6*100),ncol=6)

for(r in 1:100){
  sim.data=gen_data(500,0.9,1.3,ev=0.8,c(0,1),c(0.5),c("b"),c("b"),c("z1","x1"))
  for(i in 1:500){
    if(sim.data[i,1]==0 & sim.data[i,2]!=Inf){
      sim.data[i,1]=NA
    }
  }
  
  sim.data.uncure=sim.data[sim.data[,3]==1,]
  sim.surv1=Surv(time=sim.data.uncure[,1],time2=sim.data.uncure[,2],type="interval2")
  censor_save32[r,1]=sum(sim.surv1[,3]==1)
  censor_save32[r,2]=sum(sim.surv1[,3]==0)
  censor_save32[r,3]=sum(sim.surv1[,3]==2)
  censor_save32[r,4]=sum(sim.surv1[,3]==3)
  
  sim.surv=Surv(time=sim.data[,1],time2=sim.data[,2],type="interval2")
  test=phmc_mpl(sim.surv~sim.data$X,pi.formula = ~sim.data$Z,data=sim.data,phmc_mpl.control(n.knots=c(4,0),maxIter = c(1,10000,10001)))
  
  beta_save32[r,1]=test$beta[1]
  beta_save32[r,2]=test$beta[2]
  beta_save32[r,3]=test$se$se_H[1]
  beta_save32[r,4]=test$se$se_H[2]
  
  gamma_save32[r,1]=test$gamma[1]
  gamma_save32[r,2]=test$se$se_H[3]
  
  converge_save32[r,1]=test$convergence[1]
  converge_save32[r,2]=test$convergence[3]
  
  theta_save32[r,]=test$theta
  
  knots_save32[r,]=test$knots$Alpha
  
  print(c(32,r,test$se$se_H[1] ,test$convergence[1]))
  
}

dfbeta_save32=data.frame(beta_save32)
dfgamma_save32=data.frame(gamma_save32)
dfcensor_save32=data.frame(censor_save32) #event, right, left, interval
dfconverge_save32=data.frame(converge_save32) #converged yes/no, number of iterations
dftheta_save32=data.frame(theta_save32)
dfknots_save32=data.frame(knots_save32)

write.table(dfbeta_save32,"beta32.csv",sep=",",dec = ".",row.names = F)
write.table(dfgamma_save32,"gamma32.csv",sep = ",",dec=".",row.names = F)
write.table(dfcensor_save32,"censor32.csv",sep = ",",dec=".",row.names = F)
write.table(dfconverge_save32,"converge32.csv",sep=",",dec=".",row.names = F)
write.table(dftheta_save32,"theta32.csv",sep = ",",dec=".",row.names = F)
write.table(dfknots_save32,"knots32.csv",sep = ",",dec=".",row.names = F)

# 33. n=500,event=0.8,a=(0.9,1.3),beta=[-1,1]------


beta_save33=matrix(rep(0,4*100),ncol=4)
gamma_save33=matrix(rep(0,2*100),ncol=2)
censor_save33=matrix(rep(0,4*100),ncol=4) #event, right, left, interval
converge_save33=matrix(rep(0,2*100),ncol=2) #converged yes/no, number of iterations
theta_save33=matrix(rep(0,6*100),ncol=6)
knots_save33=matrix(rep(0,5*100),ncol=5)
f
for(r in 1:100){
  sim.data=gen_data(500,0.8,1.2,ev=0.75,c(-1,1),c(0.5),c("b"),c("b"),c("z1","x1"))
  for(i in 1:500){
    if(sim.data[i,1]==0 & sim.data[i,2]!=Inf){
      sim.data[i,1]=NA
    }
  }
  
  sim.data.uncure=sim.data[sim.data[,3]==1,]
  sim.surv1=Surv(time=sim.data.uncure[,1],time2=sim.data.uncure[,2],type="interval2")
  censor_save33[r,1]=sum(sim.surv1[,3]==1)
  censor_save33[r,2]=sum(sim.surv1[,3]==0)
  censor_save33[r,3]=sum(sim.surv1[,3]==2)
  censor_save33[r,4]=sum(sim.surv1[,3]==3)
  
  sim.surv=Surv(time=sim.data[,1],time2=sim.data[,2],type="interval2")
  test=phmc_mpl(sim.surv~sim.data$X,pi.formula = ~sim.data$Z,data=sim.data,phmc_mpl.control(n.knots=c(1,0),maxIter = c(1,10000,10001)))
  
  beta_save33[r,1]=test$beta[1]
  beta_save33[r,2]=test$beta[2]
  beta_save33[r,3]=test$se$se_H[1]
  beta_save33[r,4]=test$se$se_H[2]
  
  gamma_save33[r,1]=test$gamma[1]
  gamma_save33[r,2]=test$se$se_H[3]
  
  #converge_save33[r,1]=test$convergence[1]
  #converge_save33[r,2]=test$convergence[3]
  
  #theta_save33[r,]=test$theta
  
  #knots_save33[r,]=test$knots$Alpha
  
  print(c(33,r,test$se$se_H[1] ,test$beta[1]))
  
}

dfbeta_save33=data.frame(beta_save33)
dfgamma_save33=data.frame(gamma_save33)
dfcensor_save33=data.frame(censor_save33) #event, right, left, interval
dfconverge_save33=data.frame(converge_save33) #converged yes/no, number of iterations
dftheta_save33=data.frame(theta_save33)
dfknots_save33=data.frame(knots_save33)

write.table(dfbeta_save33,"beta33.csv",sep=",",dec = ".",row.names = F)
write.table(dfgamma_save33,"gamma33.csv",sep = ",",dec=".",row.names = F)
write.table(dfcensor_save33,"censor33.csv",sep = ",",dec=".",row.names = F)
write.table(dfconverge_save33,"converge33.csv",sep=",",dec=".",row.names = F)
write.table(dftheta_save33,"theta33.csv",sep = ",",dec=".",row.names = F)
write.table(dfknots_save33,"knots33.csv",sep = ",",dec=".",row.names = F)





#old stuff -------

beta_save=matrix(rep(0,4*100),ncol=4)
gamma_save=matrix(rep(0,2*100),ncol=2)
censor_save=matrix(rep(0,4*100),ncol=4) #event, right, left, interval
converge_save=matrix(rep(0,2*100),ncol=2) #converged yes/no, number of iterations


for(r in 1:100){
  sim.data=gen_data(500,1,1,ev=0.5,c(1,1),c(0.5),c("b"),c("b"),c("z1","x1"))
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
  test=phmc_mpl(sim.surv~sim.data[,5],pi.formula = ~sim.data[,4],data=as.data.frame(sim.data),phmc_mpl.control(n.knots=c(6,0),maxIter = c(1,10000,10001)))
  
  beta_save[r,1]=test$beta[1]
  beta_save[r,2]=test$beta[2]
  beta_save[r,3]=test$se$se_H[1]
  beta_save[r,4]=test$se$se_H[2]
  
  gamma_save[r,1]=test$gamma[1]
  gamma_save[r,2]=test$se$se_H[3]
  
  converge_save[r,1]=test$convergence[1]
  converge_save[r,2]=test$convergence[3]
  
  print(c(r,test$se$se_H[1] ,test$convergence[1]))
  
}
beta_save
gamma_save
censor_save
converge_save 

mean(beta_save[!is.na(beta_save[,3]),1])
mean(beta_save[!is.na(beta_save[,3]),2])
mean(beta_save[!is.na(beta_save[,3]),3])
mean(beta_save[!is.na(beta_save[,3]),4])
sd(beta_save[!is.na(beta_save[,3]),1])
sd(beta_save[!is.na(beta_save[,3]),2])

mean(gamma_save[!is.na(gamma_save[,2]),1])
mean(gamma_save[!is.na(gamma_save[,2]),2])
sd(gamma_save[!is.na(gamma_save[,2]),1])


mean(censor_save[,1])/362.64
mean(censor_save[,2])/362.64
mean(censor_save[,3])/362.64
mean(censor_save[,4])/362.64


mean(apply(censor_save,1,sum))




while(!is.na(test$se$se_H[1])){
  sim.data=gen_data(500,0.9,1.3,ev=0.25,c(1,1),c(0.5),c("b"),c("b"),c("z1","x1"))
  for(i in 1:500){
    if(sim.data[i,1]==0 & sim.data[i,2]!=Inf){
      sim.data[i,1]=NA
    }
  }
  
  sim.data.uncure=sim.data[sim.data[,3]==1,]
  sim.surv1=Surv(time=sim.data.uncure[,1],time2=sim.data.uncure[,2],type="interval2")
  
  sim.surv=Surv(time=sim.data[,1],time2=sim.data[,2],type="interval2")
  test=phmc_mpl(sim.surv~sim.data[,5],pi.formula = ~sim.data[,4],data=as.data.frame(sim.data),phmc_mpl.control(n.knots=c(6,0),maxIter = c(1,1000,1001)))
  print(test$se$se_H[1])
}



sim.data=gen_data(100,0.9,1.3,ev=0.5,c(1,1),c(0.5),c("b"),c("b"),c("z1","x1"))
for(i in 1:100){
  if(sim.data[i,1]==0 & sim.data[i,2]!=Inf){
    sim.data[i,1]=NA
  }
}

sim.data.uncure=sim.data[sim.data[,3]==1,]
sim.surv1=Surv(time=sim.data.uncure[,1],time2=sim.data.uncure[,2],type="interval2")

sim.surv=Surv(time=sim.data[,1],time2=sim.data[,2],type="interval2")
test=phmc_mpl(sim.surv~sim.data[,5],pi.formula = ~sim.data[,4],data=as.data.frame(sim.data),phmc_mpl.control(n.knots=c(5,0),maxIter = c(1,10000,10001)))
summary(test)





