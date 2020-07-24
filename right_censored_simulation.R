#right censored data only for comparative study with smcure


#generate data
gen_data_right=function(n,C_mean,b_true,g_true,b_type,g_type,overlap){
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
  
  #generate censoring times for cured population
  T_i.c = rexp(cure.n,1/2)
  
  #generate event & censoring times for non-cured population
  U_Y = runif(uncure.n)
  Y_i.u=C_i.u=T_i.u=censor_i.u = rep(0, uncure.n)
  X_std_uncure=X_std[cure_i==1,]
  for(i in 1:uncure.n){
    Y_i.u[i] = (-log(U_Y[i])/exp(sum(g_true*X_std_uncure[i])))^(1/3) #F^-1 method to sample from weibull(1, 3)
  }
  C_i.u = rexp(uncure.n,1/C_mean)
  
  
  #test if event < censoring for non-cured population
  for(i in 1:uncure.n){
    if(Y_i.u[i]<C_i.u[i]){
      T_i.u[i]=Y_i.u[i]
      censor_i.u[i]=1
    }else{
      T_i.u[i]=C_i.u[i]
    }
  }
  
  
  #combine populations and record censoring information (0=censored, 1=event)
  T_i=censor=rep(0,n)
  T_i[cure_i==0]=T_i.c
  T_i[cure_i==1]=T_i.u
  censor[cure_i==1]=censor_i.u
  
  
  data=data.frame(T_i,censor,cure_i,Z,X)
  colnames(data)=c("time","event","cure","Z","X1")
  data
  
}



sim.data=gen_data_right(500,4.2,c(1,1),c(0.5),c("b"),c("b"),c("z1","x1"))

sim.data.uncure=sim.data[sim.data[,3]==1,]
sim.surv1=Surv(time=sim.data.uncure[,1],event=sim.data.uncure[,2])


sim.surv=Surv(time=sim.data$time,event =sim.data$event)
test=phmc_mpl(sim.surv~sim.data$X1,pi.formula = ~sim.data$Z,data=sim.data,phmc_mpl.control(conv_limit = 1e-4, n.knots=c(6,0),maxIter = c(1,10000,10001)))
summary(test)

#library("smcure")
smc1 = smcure(sim.surv~X1, cureform=~Z, data=sim.data, model="ph")

# 1. beta = [1,1], n=2000, censor 0.5 ---------

r1_censorsave=matrix(rep(0,2*100),ncol=2)#whole dataset event proportion, non-cured event proportion
r1_bsave_mpl=matrix(rep(0,4*100),ncol=4)
r1_gsave_mpl=matrix(rep(0,2*100),ncol=2)
r1_bsave_em=matrix(rep(0,4*100),ncol=4)
r1_gsave_em=matrix(rep(0,2*100),ncol=2)

for(r in 1:100){
  sim.data=gen_data_right(2000,0.5,c(1,1),c(0.5),c("b"),c("b"),c("z1","x1"))
  
  sim.data.uncure=sim.data[sim.data[,3]==1,]
  sim.surv1=Surv(time=sim.data.uncure[,1],event=sim.data.uncure[,2])
  sim.surv=Surv(time=sim.data$time,event =sim.data$event)
  r1_censorsave[r,1]=sum(sim.surv[,2])/2000
  r1_censorsave[r,2]=sum(sim.surv1[,2])/nrow(sim.surv1)
  
  try1=try(phmc_mpl(sim.surv~sim.data$X1,pi.formula = ~sim.data$Z,data=sim.data,phmc_mpl.control(conv_limit = 1e-4, n.knots=c(5,0),maxIter = c(1,10000,10001))))
  try2=try(smcure(sim.surv~X1, cureform=~Z, data=sim.data, model="ph"))
  if(class(try1)!="try-error" & class(try2)!="try-error"){
    test=try1
    smc1=try2
    
    r1_bsave_mpl[r,1]=test$beta[1]
    r1_bsave_mpl[r,2]=test$beta[2]
    r1_bsave_mpl[r,3]=test$se$se_H[1]
    r1_bsave_mpl[r,4]=test$se$se_H[2]
    
    r1_gsave_mpl[r,1]=test$gamma[1]
    r1_gsave_mpl[r,2]=test$se$se_H[3]
    
    r1_bsave_em[r,1]=smc1$b[1]
    r1_bsave_em[r,2]=smc1$b[2]
    r1_bsave_em[r,3]=smc1$b_sd[1]
    r1_bsave_em[r,4]=smc1$b_sd[2]
    
    r1_gsave_em[r,1]=smc1$beta[1]
    r1_gsave_em[r,2]=smc1$beta_sd[1]
    
    print(c(r,"both 1",test$se$se_H[1]))
  } else if(class(try1)!="try-error"){
    test=try1
    r1_bsave_mpl[r,1]=test$beta[1]
    r1_bsave_mpl[r,2]=test$beta[2]
    r1_bsave_mpl[r,3]=test$se$se_H[1]
    r1_bsave_mpl[r,4]=test$se$se_H[2]
    
    r1_gsave_mpl[r,1]=test$gamma[1]
    r1_gsave_mpl[r,2]=test$se$se_H[3]
    print(c(r, "mpl 1",test$se$se_H[1]))
    
  } else if(class(try2)!="try-error"){
    smc1=try2
    r1_bsave_em[r,1]=smc1$b[1]
    r1_bsave_em[r,2]=smc1$b[2]
    r1_bsave_em[r,3]=smc1$b_sd[1]
    r1_bsave_em[r,4]=smc1$b_sd[2]
    
    r1_gsave_em[r,1]=smc1$beta[1]
    r1_gsave_em[r,2]=smc1$beta_sd[1]
    print(c(r, "em 1"))
  } else{
    print(c(r, "neither 1"))
  }
}


mean(r1_censorsave[!is.na(r1_bsave_mpl[,3])],1)
mean(r1_censorsave[!is.na(r1_bsave_mpl[,3])],2)
r1_censorsave=data.frame(r1_censorsave)
write.table(r1_censorsave, "r1censor.csv",sep=",",dec=".")

colnames(r1_bsave_mpl)=c("b0","b1","seb0","seb1")
r1_bsave_mpl=data.frame(r1_bsave_mpl)
write.table(r1_bsave_mpl, "r1bsave_mpl.csv",sep=",",dec=".",row.names = F)
colnames(r1_bsave_em)=c("b0","b1","seb0","seb1")
r1_bsave_em=data.frame(r1_bsave_em)
write.table(r1_bsave_em, "r1bsave_em.csv",sep=",",dec=".",row.names = F)
colnames(r1_gsave_mpl)=c("g1","seg1")
r1_gsave_mpl=data.frame(r1_gsave_mpl)
write.table(r1_gsave_mpl, "r1gsave_mpl.csv",sep=",",dec=".",row.names = F)
colnames(r1_gsave_em)=c("g1","seg1")
r1_gsave_em=data.frame(r1_gsave_em)
write.table(r1_gsave_em, "r1gsave_em.csv",sep=",",dec=".",row.names = F)


# 2. beta=[0,1], n=2000, censor 0.5 ---------

r2_censorsave=matrix(rep(0,2*100),ncol=2)#whole dataset event proportion, non-cured event proportion
r2_bsave_mpl=matrix(rep(0,4*100),ncol=4)
r2_gsave_mpl=matrix(rep(0,2*100),ncol=2)
r2_bsave_em=matrix(rep(0,4*100),ncol=4)
r2_gsave_em=matrix(rep(0,2*100),ncol=2)

for(r in 1:100){
  sim.data=gen_data_right(2000,0.5,c(0,1),c(0.5),c("b"),c("b"),c("z1","x1"))
  
  sim.data.uncure=sim.data[sim.data[,3]==1,]
  sim.surv1=Surv(time=sim.data.uncure[,1],event=sim.data.uncure[,2])
  sim.surv=Surv(time=sim.data$time,event =sim.data$event)
  r2_censorsave[r,1]=sum(sim.surv[,2])/2000
  r2_censorsave[r,2]=sum(sim.surv1[,2])/nrow(sim.surv1)
  
  try1=try(phmc_mpl(sim.surv~sim.data$X1,pi.formula = ~sim.data$Z,data=sim.data,phmc_mpl.control(conv_limit = 1e-4, n.knots=c(4,0),maxIter = c(1,10000,10001))))
  try2=try(smcure(sim.surv~X1, cureform=~Z, data=sim.data, model="ph"))
  if(class(try1)!="try-error" & class(try2)!="try-error"){
    test=try1
    smc1=try2
    
    r2_bsave_mpl[r,1]=test$beta[1]
    r2_bsave_mpl[r,2]=test$beta[2]
    r2_bsave_mpl[r,3]=test$se$se_H[1]
    r2_bsave_mpl[r,4]=test$se$se_H[2]
    
    r2_gsave_mpl[r,1]=test$gamma[1]
    r2_gsave_mpl[r,2]=test$se$se_H[3]
    
    r2_bsave_em[r,1]=smc1$b[1]
    r2_bsave_em[r,2]=smc1$b[2]
    r2_bsave_em[r,3]=smc1$b_sd[1]
    r2_bsave_em[r,4]=smc1$b_sd[2]
    
    r2_gsave_em[r,1]=smc1$beta[1]
    r2_gsave_em[r,2]=smc1$beta_sd[1]
    
    print(c(r,"both 2",test$se$se_H[1]))
  } else if(class(try1)!="try-error"){
    test=try1
    r2_bsave_mpl[r,1]=test$beta[1]
    r2_bsave_mpl[r,2]=test$beta[2]
    r2_bsave_mpl[r,3]=test$se$se_H[1]
    r2_bsave_mpl[r,4]=test$se$se_H[2]
    
    r2_gsave_mpl[r,1]=test$gamma[1]
    r2_gsave_mpl[r,2]=test$se$se_H[3]
    print(c(r, "mpl 2",test$se$se_H[1]))
    
  } else if(class(try2)!="try-error"){
    smc1=try2
    r2_bsave_em[r,1]=smc1$b[1]
    r2_bsave_em[r,2]=smc1$b[2]
    r2_bsave_em[r,3]=smc1$b_sd[1]
    r2_bsave_em[r,4]=smc1$b_sd[2]
    
    r2_gsave_em[r,1]=smc1$beta[1]
    r2_gsave_em[r,2]=smc1$beta_sd[1]
    print(c(r, "em 2"))
  } else{
    print(c(r, "neither 2"))
  }
}


mean(r2_censorsave[!is.na(r2_bsave_mpl[,3])],1)
mean(r2_censorsave[!is.na(r2_bsave_mpl[,3])],2)
r2_censorsave=data.frame(r2_censorsave)
write.table(r2_censorsave, "r2censor.csv",sep=",",dec=".")

colnames(r2_bsave_mpl)=c("b0","b1","seb0","seb1")
r2_bsave_mpl=data.frame(r2_bsave_mpl)
write.table(r2_bsave_mpl, "r2bsave_mpl.csv",sep=",",dec=".",row.names = F)
colnames(r2_bsave_em)=c("b0","b1","seb0","seb1")
r2_bsave_em=data.frame(r2_bsave_em)
write.table(r2_bsave_em, "r2bsave_em.csv",sep=",",dec=".",row.names = F)
colnames(r2_gsave_mpl)=c("g1","seg1")
r2_gsave_mpl=data.frame(r2_gsave_mpl)
write.table(r2_gsave_mpl, "r2gsave_mpl.csv",sep=",",dec=".",row.names = F)
colnames(r2_gsave_em)=c("g1","seg1")
r2_gsave_em=data.frame(r2_gsave_em)
write.table(r2_gsave_em, "r2gsave_em.csv",sep=",",dec=".",row.names = F)


# 3. beta = [-1,1], n=2000, censor 0.5 --------

r3_censorsave=matrix(rep(0,2*100),ncol=2)#whole dataset event proportion, non-cured event proportion
r3_bsave_mpl=matrix(rep(0,4*100),ncol=4)
r3_gsave_mpl=matrix(rep(0,2*100),ncol=2)
r3_bsave_em=matrix(rep(0,4*100),ncol=4)
r3_gsave_em=matrix(rep(0,2*100),ncol=2)

for(r in 1:100){
  sim.data=gen_data_right(2000,0.5,c(-1,1),c(0.5),c("b"),c("b"),c("z1","x1"))
  
  sim.data.uncure=sim.data[sim.data[,3]==1,]
  sim.surv1=Surv(time=sim.data.uncure[,1],event=sim.data.uncure[,2])
  sim.surv=Surv(time=sim.data$time,event =sim.data$event)
  r3_censorsave[r,1]=sum(sim.surv[,2])/2000
  r3_censorsave[r,2]=sum(sim.surv1[,2])/nrow(sim.surv1)
  
  try1=try(phmc_mpl(sim.surv~sim.data$X1,pi.formula = ~sim.data$Z,data=sim.data,phmc_mpl.control(conv_limit = 1e-4, n.knots=c(3,0),maxIter = c(1,10000,10001))))
  try2=try(smcure(sim.surv~X1, cureform=~Z, data=sim.data, model="ph"))
  if(class(try1)!="try-error" & class(try2)!="try-error"){
    test=try1
    smc1=try2
    
    r3_bsave_mpl[r,1]=test$beta[1]
    r3_bsave_mpl[r,2]=test$beta[2]
    r3_bsave_mpl[r,3]=test$se$se_H[1]
    r3_bsave_mpl[r,4]=test$se$se_H[2]
    
    r3_gsave_mpl[r,1]=test$gamma[1]
    r3_gsave_mpl[r,2]=test$se$se_H[3]
    
    r3_bsave_em[r,1]=smc1$b[1]
    r3_bsave_em[r,2]=smc1$b[2]
    r3_bsave_em[r,3]=smc1$b_sd[1]
    r3_bsave_em[r,4]=smc1$b_sd[2]
    
    r3_gsave_em[r,1]=smc1$beta[1]
    r3_gsave_em[r,2]=smc1$beta_sd[1]
    
    print(c(r,"both 3",test$se$se_H[1]))
  } else if(class(try1)!="try-error"){
    test=try1
    r3_bsave_mpl[r,1]=test$beta[1]
    r3_bsave_mpl[r,2]=test$beta[2]
    r3_bsave_mpl[r,3]=test$se$se_H[1]
    r3_bsave_mpl[r,4]=test$se$se_H[2]
    
    r3_gsave_mpl[r,1]=test$gamma[1]
    r3_gsave_mpl[r,2]=test$se$se_H[3]
    print(c(r, "mpl 3",test$se$se_H[1]))
    
  } else if(class(try2)!="try-error"){
    smc1=try2
    r3_bsave_em[r,1]=smc1$b[1]
    r3_bsave_em[r,2]=smc1$b[2]
    r3_bsave_em[r,3]=smc1$b_sd[1]
    r3_bsave_em[r,4]=smc1$b_sd[2]
    
    r3_gsave_em[r,1]=smc1$beta[1]
    r3_gsave_em[r,2]=smc1$beta_sd[1]
    print(c(r, "em 3"))
  } else{
    print(c(r, "neither 3"))
  }
}


mean(r3_censorsave[!is.na(r3_bsave_mpl[,3])],1)
mean(r3_censorsave[!is.na(r3_bsave_mpl[,3])],2)
r3_censorsave=data.frame(r3_censorsave)
write.table(r3_censorsave, "r3censor.csv",sep=",",dec=".")

colnames(r3_bsave_mpl)=c("b0","b1","seb0","seb1")
r3_bsave_mpl=data.frame(r3_bsave_mpl)
write.table(r3_bsave_mpl, "r3bsave_mpl.csv",sep=",",dec=".",row.names = F)
colnames(r3_bsave_em)=c("b0","b1","seb0","seb1")
r3_bsave_em=data.frame(r3_bsave_em)
write.table(r3_bsave_em, "r3bsave_em.csv",sep=",",dec=".",row.names = F)
colnames(r3_gsave_mpl)=c("g1","seg1")
r3_gsave_mpl=data.frame(r3_gsave_mpl)
write.table(r3_gsave_mpl, "r3gsave_mpl.csv",sep=",",dec=".",row.names = F)
colnames(r3_gsave_em)=c("g1","seg1")
r3_gsave_em=data.frame(r3_gsave_em)
write.table(r3_gsave_em, "r3gsave_em.csv",sep=",",dec=".",row.names = F)





# 4. beta = [1,1], n=500, censor 0.5 --------

r4_censorsave=matrix(rep(0,2*100),ncol=2)#whole dataset event proportion, non-cured event proportion
r4_bsave_mpl=matrix(rep(0,4*100),ncol=4)
r4_gsave_mpl=matrix(rep(0,2*100),ncol=2)
r4_bsave_em=matrix(rep(0,4*100),ncol=4)
r4_gsave_em=matrix(rep(0,2*100),ncol=2)


for(r in 1:100){
  sim.data=gen_data_right(500,0.5,c(1,1),c(0.5),c("b"),c("b"),c("z1","x1"))
  
  sim.data.uncure=sim.data[sim.data[,3]==1,]
  sim.surv1=Surv(time=sim.data.uncure[,1],event=sim.data.uncure[,2])
  sim.surv=Surv(time=sim.data$time,event =sim.data$event)
  r4_censorsave[r,1]=sum(sim.surv[,2])/500
  r4_censorsave[r,2]=sum(sim.surv1[,2])/nrow(sim.surv1)
  
  try1=try(phmc_mpl(sim.surv~sim.data$X1,pi.formula = ~sim.data$Z,data=sim.data,phmc_mpl.control(conv_limit = 1e-4, n.knots=c(3,0),maxIter = c(1,10000,10001))))
  try2=try(smcure(sim.surv~X1, cureform=~Z, data=sim.data, model="ph"))
  if(class(try1)!="try-error" & class(try2)!="try-error"){
    test=try1
    smc1=try2
    
    r4_bsave_mpl[r,1]=test$beta[1]
    r4_bsave_mpl[r,2]=test$beta[2]
    r4_bsave_mpl[r,3]=test$se$se_H[1]
    r4_bsave_mpl[r,4]=test$se$se_H[2]
    
    r4_gsave_mpl[r,1]=test$gamma[1]
    r4_gsave_mpl[r,2]=test$se$se_H[3]
    
    r4_bsave_em[r,1]=smc1$b[1]
    r4_bsave_em[r,2]=smc1$b[2]
    r4_bsave_em[r,3]=smc1$b_sd[1]
    r4_bsave_em[r,4]=smc1$b_sd[2]
    
    r4_gsave_em[r,1]=smc1$beta[1]
    r4_gsave_em[r,2]=smc1$beta_sd[1]
    
    print(c(r,"both 4",test$se$se_H[1]))
  } else if(class(try1)!="try-error"){
    test=try1
    r4_bsave_mpl[r,1]=test$beta[1]
    r4_bsave_mpl[r,2]=test$beta[2]
    r4_bsave_mpl[r,3]=test$se$se_H[1]
    r4_bsave_mpl[r,4]=test$se$se_H[2]
    
    r4_gsave_mpl[r,1]=test$gamma[1]
    r4_gsave_mpl[r,2]=test$se$se_H[3]
    print(c(r, "mpl 4",test$se$se_H[1]))
    
  } else if(class(try2)!="try-error"){
    smc1=try2
    r4_bsave_em[r,1]=smc1$b[1]
    r4_bsave_em[r,2]=smc1$b[2]
    r4_bsave_em[r,3]=smc1$b_sd[1]
    r4_bsave_em[r,4]=smc1$b_sd[2]
    
    r4_gsave_em[r,1]=smc1$beta[1]
    r4_gsave_em[r,2]=smc1$beta_sd[1]
    print(c(r, "em 4"))
  } else{
    print(c(r, "neither 4"))
  }
}


mean(r4_censorsave[!is.na(r4_bsave_mpl[,3])],1)
mean(r4_censorsave[!is.na(r4_bsave_mpl[,3])],2)
r4_censorsave=data.frame(r4_censorsave)
write.table(r4_censorsave, "r4censor.csv",sep=",",dec=".")



colnames(r4_bsave_mpl)=c("b0","b1","seb0","seb1")
r4_bsave_mpl=data.frame(r4_bsave_mpl)
write.table(r4_bsave_mpl, "r4bsave_mpl.csv",sep=",",dec=".",row.names = F)
colnames(r4_bsave_em)=c("b0","b1","seb0","seb1")
r4_bsave_em=data.frame(r4_bsave_em)
write.table(r4_bsave_em, "r4bsave_em.csv",sep=",",dec=".",row.names = F)
colnames(r4_gsave_mpl)=c("g1","seg1")
r4_gsave_mpl=data.frame(r4_gsave_mpl)
write.table(r4_gsave_mpl, "r4gsave_mpl.csv",sep=",",dec=".",row.names = F)
colnames(r4_gsave_em)=c("g1","seg1")
r4_gsave_em=data.frame(r4_gsave_em)
write.table(r4_gsave_em, "r4gsave_em.csv",sep=",",dec=".",row.names = F)

#mpl
mean(r4_bsave_mpl[!is.na(r4_bsave_mpl[,3]),1])
mean(r4_bsave_mpl[!is.na(r4_bsave_mpl[,3]),2])
mean(r4_bsave_mpl[!is.na(r4_bsave_mpl[,3]),3])
sd(r4_bsave_mpl[!is.na(r4_bsave_mpl[,3]),1])
mean(r4_bsave_mpl[!is.na(r4_bsave_mpl[,3]),4])
sd(r4_bsave_mpl[!is.na(r4_bsave_mpl[,3]),2])

mean(r4_gsave_mpl[!is.na(r4_gsave_mpl[,2]),1])
mean(r4_gsave_mpl[!is.na(r4_gsave_mpl[,2]),2])
sd(r4_gsave_mpl[!is.na(r4_gsave_mpl[,2]),1])

r4_bsave_mpl$b0_LL=r4_bsave_mpl$b0_UL=r4_bsave_mpl$b1_LL=r4_bsave_mpl$b1_UL=rep(0,100)
r4_gsave_mpl$g1_LL=r4_gsave_mpl$g1_UL=rep(0,100)

r4_bsave_mpl$b0_LL[!is.na(r4_bsave_mpl$seb0)]=r4_bsave_mpl$b0[!is.na(r4_bsave_mpl$seb0)]-1.96*r4_bsave_mpl$seb0[!is.na(r4_bsave_mpl$seb0)]
r4_bsave_mpl$b0_UL[!is.na(r4_bsave_mpl$seb0)]=r4_bsave_mpl$b0[!is.na(r4_bsave_mpl$seb0)]+1.96*r4_bsave_mpl$seb0[!is.na(r4_bsave_mpl$seb0)]

r4_bsave_mpl$b1_LL[!is.na(r4_bsave_mpl$seb0)]=r4_bsave_mpl$b1[!is.na(r4_bsave_mpl$seb0)]-1.96*r4_bsave_mpl$seb1[!is.na(r4_bsave_mpl$seb0)]
r4_bsave_mpl$b1_UL[!is.na(r4_bsave_mpl$seb0)]=r4_bsave_mpl$b1[!is.na(r4_bsave_mpl$seb0)]+1.96*r4_bsave_mpl$seb1[!is.na(r4_bsave_mpl$seb0)]

r4_bsave_mpl$b0_cov=r4_bsave_mpl$b1_cov=rep(0,100)
r4_bsave_mpl$b0_cov[!is.na(r4_bsave_mpl$seb0)]=as.numeric(r4_bsave_mpl$b0_LL[!is.na(r4_bsave_mpl$seb0)]<=1 & 1<=r4_bsave_mpl$b0_UL[!is.na(r4_bsave_mpl$seb0)])
r4_bsave_mpl$b1_cov[!is.na(r4_bsave_mpl$seb0)]=as.numeric(r4_bsave_mpl$b1_LL[!is.na(r4_bsave_mpl$seb0)]<=1 & 1<=r4_bsave_mpl$b1_UL[!is.na(r4_bsave_mpl$seb0)])

r4_gsave_mpl$g1_LL[!is.na(r4_gsave_mpl$seg1)]=r4_gsave_mpl$g1[!is.na(r4_gsave_mpl$seg1)]-1.96*r4_gsave_mpl$seg1[!is.na(r4_gsave_mpl$seg1)]
r4_gsave_mpl$g1_UL[!is.na(r4_gsave_mpl$seg1)]=r4_gsave_mpl$g1[!is.na(r4_gsave_mpl$seg1)]+1.96*r4_gsave_mpl$seg1[!is.na(r4_gsave_mpl$seg1)]
r4_gsave_mpl$g1_cov=rep(0,100)
r4_gsave_mpl$g1_cov[!is.na(r4_gsave_mpl$seg1)]=as.numeric(r4_gsave_mpl$g1_LL[!is.na(r4_gsave_mpl$seg1)]<=0.5 & 0.5<=r4_gsave_mpl$g1_UL[!is.na(r4_gsave_mpl$seg1)])

#em

mean(r4_bsave_em[,1])
mean(r4_bsave_em[,2])
mean(r4_bsave_em[,3])
sd(r4_bsave_em[,1])
mean(r4_bsave_em[,4])
sd(r4_bsave_em[,2])

mean(r4_gsave_em[,1])
mean(r4_gsave_em[,2])
sd(r4_gsave_em[,1])

r4_bsave_em$b0_LL=r4_bsave_em$b0_UL=r4_bsave_em$b1_LL=r4_bsave_em$b1_UL=rep(0,100)
r4_gsave_em$g1_LL=r4_gsave_em$g1_UL=rep(0,100)

r4_bsave_em$b0_LL=r4_bsave_em$b0-1.96*r4_bsave_em$seb0
r4_bsave_em$b0_UL=r4_bsave_em$b0+1.96*r4_bsave_em$seb0

r4_bsave_em$b1_LL=r4_bsave_em$b1-1.96*r4_bsave_em$seb1
r4_bsave_em$b1_UL=r4_bsave_em$b1+1.96*r4_bsave_em$seb1

r4_bsave_em$b0_cov=r4_bsave_em$b1_cov=rep(0,100)
r4_bsave_em$b0_cov=as.numeric(r4_bsave_em$b0_LL<=1 & 1<=r4_bsave_em$b0_UL)
r4_bsave_em$b1_cov=as.numeric(r4_bsave_em$b1_LL<=1 & 1<=r4_bsave_em$b1_UL)

r4_gsave_em$g1_LL=r4_gsave_em$g1-1.96*r4_gsave_em$seg1
r4_gsave_em$g1_UL=r4_gsave_em$g1+1.96*r4_gsave_em$seg1
r4_gsave_em$g1_cov=rep(0,100)
r4_gsave_em$g1_cov=as.numeric(r4_gsave_em$g1_LL<=0.5 & 0.5<=r4_gsave_em$g1_UL)



# 5. beta = [0,1], n=500, censor 0.5 ---------


r5_censorsave=matrix(rep(0,2*100),ncol=2)#whole dataset event proportion, non-cured event proportion
r5_bsave_mpl=matrix(rep(0,4*100),ncol=4)
r5_gsave_mpl=matrix(rep(0,2*100),ncol=2)
r5_bsave_em=matrix(rep(0,4*100),ncol=4)
r5_gsave_em=matrix(rep(0,2*100),ncol=2)

for(r in 1:100){
  sim.data=gen_data_right(500,0.5,c(0,1),c(0.5),c("b"),c("b"),c("z1","x1"))
  
  sim.data.uncure=sim.data[sim.data[,3]==1,]
  sim.surv1=Surv(time=sim.data.uncure[,1],event=sim.data.uncure[,2])
  sim.surv=Surv(time=sim.data$time,event =sim.data$event)
  r5_censorsave[r,1]=sum(sim.surv[,2])/500
  r5_censorsave[r,2]=sum(sim.surv1[,2])/nrow(sim.surv1)
  
  try1=try(phmc_mpl(sim.surv~sim.data$X1,pi.formula = ~sim.data$Z,data=sim.data,phmc_mpl.control(conv_limit = 1e-4, n.knots=c(2,0),maxIter = c(1,10000,10001))))
  try2=try(smcure(sim.surv~X1, cureform=~Z, data=sim.data, model="ph"))
  if(class(try1)!="try-error" & class(try2)!="try-error"){
    test=try1
    smc1=try2
    
    r5_bsave_mpl[r,1]=test$beta[1]
    r5_bsave_mpl[r,2]=test$beta[2]
    r5_bsave_mpl[r,3]=test$se$se_H[1]
    r5_bsave_mpl[r,4]=test$se$se_H[2]
    
    r5_gsave_mpl[r,1]=test$gamma[1]
    r5_gsave_mpl[r,2]=test$se$se_H[3]
    
    r5_bsave_em[r,1]=smc1$b[1]
    r5_bsave_em[r,2]=smc1$b[2]
    r5_bsave_em[r,3]=smc1$b_sd[1]
    r5_bsave_em[r,4]=smc1$b_sd[2]
    
    r5_gsave_em[r,1]=smc1$beta[1]
    r5_gsave_em[r,2]=smc1$beta_sd[1]
    
    print(c(r,"both 5",test$se$se_H[1]))
  } else if(class(try1)!="try-error"){
    test=try1
    r5_bsave_mpl[r,1]=test$beta[1]
    r5_bsave_mpl[r,2]=test$beta[2]
    r5_bsave_mpl[r,3]=test$se$se_H[1]
    r5_bsave_mpl[r,4]=test$se$se_H[2]
    
    r5_gsave_mpl[r,1]=test$gamma[1]
    r5_gsave_mpl[r,2]=test$se$se_H[3]
    print(c(r, "mpl 5",test$se$se_H[1]))
    
  } else if(class(try2)!="try-error"){
    smc1=try2
    r5_bsave_em[r,1]=smc1$b[1]
    r5_bsave_em[r,2]=smc1$b[2]
    r5_bsave_em[r,3]=smc1$b_sd[1]
    r5_bsave_em[r,4]=smc1$b_sd[2]
    
    r5_gsave_em[r,1]=smc1$beta[1]
    r5_gsave_em[r,2]=smc1$beta_sd[1]
    print(c(r, "em 5"))
  } else{
    print(c(r, "neither 5"))
  }
}


mean(r5_censorsave[!is.na(r5_bsave_mpl[,3])],1)
mean(r5_censorsave[!is.na(r5_bsave_mpl[,3])],2)
r5_censorsave=data.frame(r5_censorsave)
write.table(r5_censorsave, "r5censor.csv",sep=",",dec=".")

colnames(r5_bsave_mpl)=c("b0","b1","seb0","seb1")
r5_bsave_mpl=data.frame(r5_bsave_mpl)
write.table(r5_bsave_mpl, "r5bsave_mpl.csv",sep=",",dec=".",row.names = F)
colnames(r5_bsave_em)=c("b0","b1","seb0","seb1")
r5_bsave_em=data.frame(r5_bsave_em)
write.table(r5_bsave_em, "r5bsave_em.csv",sep=",",dec=".",row.names = F)
colnames(r5_gsave_mpl)=c("g1","seg1")
r5_gsave_mpl=data.frame(r5_gsave_mpl)
write.table(r5_gsave_mpl, "r5gsave_mpl.csv",sep=",",dec=".",row.names = F)
colnames(r5_gsave_em)=c("g1","seg1")
r5_gsave_em=data.frame(r5_gsave_em)
write.table(r5_gsave_em, "r5gsave_em.csv",sep=",",dec=".",row.names = F)



#mpl
mean(r5_bsave_mpl[!is.na(r5_bsave_mpl[,3]),1])
mean(r5_bsave_mpl[!is.na(r5_bsave_mpl[,3]),2])
mean(r5_bsave_mpl[!is.na(r5_bsave_mpl[,3]),3])
sd(r5_bsave_mpl[!is.na(r5_bsave_mpl[,3]),1])
mean(r5_bsave_mpl[!is.na(r5_bsave_mpl[,3]),4])
sd(r5_bsave_mpl[!is.na(r5_bsave_mpl[,3]),2])

mean(r5_gsave_mpl[!is.na(r5_gsave_mpl[,2]),1])
mean(r5_gsave_mpl[!is.na(r5_gsave_mpl[,2]),2])
sd(r5_gsave_mpl[!is.na(r5_gsave_mpl[,2]),1])

r5_bsave_mpl$b0_LL=r5_bsave_mpl$b0_UL=r5_bsave_mpl$b1_LL=r5_bsave_mpl$b1_UL=rep(0,100)
r5_gsave_mpl$g1_LL=r5_gsave_mpl$g1_UL=rep(0,100)

r5_bsave_mpl$b0_LL[!is.na(r5_bsave_mpl$seb0)]=r5_bsave_mpl$b0[!is.na(r5_bsave_mpl$seb0)]-1.96*r5_bsave_mpl$seb0[!is.na(r5_bsave_mpl$seb0)]
r5_bsave_mpl$b0_UL[!is.na(r5_bsave_mpl$seb0)]=r5_bsave_mpl$b0[!is.na(r5_bsave_mpl$seb0)]+1.96*r5_bsave_mpl$seb0[!is.na(r5_bsave_mpl$seb0)]

r5_bsave_mpl$b1_LL[!is.na(r5_bsave_mpl$seb0)]=r5_bsave_mpl$b1[!is.na(r5_bsave_mpl$seb0)]-1.96*r5_bsave_mpl$seb1[!is.na(r5_bsave_mpl$seb0)]
r5_bsave_mpl$b1_UL[!is.na(r5_bsave_mpl$seb0)]=r5_bsave_mpl$b1[!is.na(r5_bsave_mpl$seb0)]+1.96*r5_bsave_mpl$seb1[!is.na(r5_bsave_mpl$seb0)]

r5_bsave_mpl$b0_cov=r5_bsave_mpl$b1_cov=rep(0,100)
r5_bsave_mpl$b0_cov[!is.na(r5_bsave_mpl$seb0)]=as.numeric(r5_bsave_mpl$b0_LL[!is.na(r5_bsave_mpl$seb0)]<=0 & 0<=r5_bsave_mpl$b0_UL[!is.na(r5_bsave_mpl$seb0)])
r5_bsave_mpl$b1_cov[!is.na(r5_bsave_mpl$seb0)]=as.numeric(r5_bsave_mpl$b1_LL[!is.na(r5_bsave_mpl$seb0)]<=1 & 1<=r5_bsave_mpl$b1_UL[!is.na(r5_bsave_mpl$seb0)])

r5_gsave_mpl$g1_LL[!is.na(r5_gsave_mpl$seg1)]=r5_gsave_mpl$g1[!is.na(r5_gsave_mpl$seg1)]-1.96*r5_gsave_mpl$seg1[!is.na(r5_gsave_mpl$seg1)]
r5_gsave_mpl$g1_UL[!is.na(r5_gsave_mpl$seg1)]=r5_gsave_mpl$g1[!is.na(r5_gsave_mpl$seg1)]+1.96*r5_gsave_mpl$seg1[!is.na(r5_gsave_mpl$seg1)]
r5_gsave_mpl$g1_cov=rep(0,100)
r5_gsave_mpl$g1_cov[!is.na(r5_gsave_mpl$seg1)]=as.numeric(r5_gsave_mpl$g1_LL[!is.na(r5_gsave_mpl$seg1)]<=0.5 & 0.5<=r5_gsave_mpl$g1_UL[!is.na(r5_gsave_mpl$seg1)])


#em

mean(r5_bsave_em[,1])
mean(r5_bsave_em[,2])
mean(r5_bsave_em[,3])
sd(r5_bsave_em[,1])
mean(r5_bsave_em[,4])
sd(r5_bsave_em[,2])

mean(r5_gsave_em[,1])
mean(r5_gsave_em[,2])
sd(r5_gsave_em[,1])

r5_bsave_em$b0_LL=r5_bsave_em$b0_UL=r5_bsave_em$b1_LL=r5_bsave_em$b1_UL=rep(0,100)
r5_gsave_em$g1_LL=r5_gsave_em$g1_UL=rep(0,100)

r5_bsave_em$b0_LL=r5_bsave_em$b0-1.96*r5_bsave_em$seb0
r5_bsave_em$b0_UL=r5_bsave_em$b0+1.96*r5_bsave_em$seb0

r5_bsave_em$b1_LL=r5_bsave_em$b1-1.96*r5_bsave_em$seb1
r5_bsave_em$b1_UL=r5_bsave_em$b1+1.96*r5_bsave_em$seb1

r5_bsave_em$b0_cov=r5_bsave_em$b1_cov=rep(0,100)
r5_bsave_em$b0_cov=as.numeric(r5_bsave_em$b0_LL<=0 & 0<=r5_bsave_em$b0_UL)
r5_bsave_em$b1_cov=as.numeric(r5_bsave_em$b1_LL<=1 & 1<=r5_bsave_em$b1_UL)

r5_gsave_em$g1_LL=r5_gsave_em$g1-1.96*r5_gsave_em$seg1
r5_gsave_em$g1_UL=r5_gsave_em$g1+1.96*r5_gsave_em$seg1
r5_gsave_em$g1_cov=rep(0,100)
r5_gsave_em$g1_cov=as.numeric(r5_gsave_em$g1_LL<=0.5 & 0.5<=r5_gsave_em$g1_UL)

# 6. beta = [-1,1], n=500, censor 0.5 --------


r6_censorsave=matrix(rep(0,2*100),ncol=2)#whole dataset event proportion, non-cured event proportion
r6_bsave_mpl=matrix(rep(0,4*100),ncol=4)
r6_gsave_mpl=matrix(rep(0,2*100),ncol=2)
r6_bsave_em=matrix(rep(0,4*100),ncol=4)
r6_gsave_em=matrix(rep(0,2*100),ncol=2)


for(r in 1:100){
  sim.data=gen_data_right(500,0.5,c(-1,1),c(0.5),c("b"),c("b"),c("z1","x1"))
  
  sim.data.uncure=sim.data[sim.data[,3]==1,]
  sim.surv1=Surv(time=sim.data.uncure[,1],event=sim.data.uncure[,2])
  sim.surv=Surv(time=sim.data$time,event =sim.data$event)
  r6_censorsave[r,1]=sum(sim.surv[,2])/500
  r6_censorsave[r,2]=sum(sim.surv1[,2])/nrow(sim.surv1)
  
  try1=try(phmc_mpl(sim.surv~sim.data$X1,pi.formula = ~sim.data$Z,data=sim.data,phmc_mpl.control(conv_limit = 1e-4, n.knots=c(1,0),maxIter = c(1,10000,10001))))
  try2=try(smcure(sim.surv~X1, cureform=~Z, data=sim.data, model="ph"))
  if(class(try1)!="try-error" & class(try2)!="try-error"){
    test=try1
    smc1=try2
    
    r6_bsave_mpl[r,1]=test$beta[1]
    r6_bsave_mpl[r,2]=test$beta[2]
    r6_bsave_mpl[r,3]=test$se$se_H[1]
    r6_bsave_mpl[r,4]=test$se$se_H[2]
    
    r6_gsave_mpl[r,1]=test$gamma[1]
    r6_gsave_mpl[r,2]=test$se$se_H[3]
    
    r6_bsave_em[r,1]=smc1$b[1]
    r6_bsave_em[r,2]=smc1$b[2]
    r6_bsave_em[r,3]=smc1$b_sd[1]
    r6_bsave_em[r,4]=smc1$b_sd[2]
    
    r6_gsave_em[r,1]=smc1$beta[1]
    r6_gsave_em[r,2]=smc1$beta_sd[1]
    
    print(c(r,"both 6",test$se$se_H[1]))
  } else if(class(try1)!="try-error"){
    test=try1
    r6_bsave_mpl[r,1]=test$beta[1]
    r6_bsave_mpl[r,2]=test$beta[2]
    r6_bsave_mpl[r,3]=test$se$se_H[1]
    r6_bsave_mpl[r,4]=test$se$se_H[2]
    
    r6_gsave_mpl[r,1]=test$gamma[1]
    r6_gsave_mpl[r,2]=test$se$se_H[3]
    print(c(r, "mpl 6",test$se$se_H[1]))
    
  } else if(class(try2)!="try-error"){
    smc1=try2
    r6_bsave_em[r,1]=smc1$b[1]
    r6_bsave_em[r,2]=smc1$b[2]
    r6_bsave_em[r,3]=smc1$b_sd[1]
    r6_bsave_em[r,4]=smc1$b_sd[2]
    
    r6_gsave_em[r,1]=smc1$beta[1]
    r6_gsave_em[r,2]=smc1$beta_sd[1]
    print(c(r, "em 6"))
  } else{
    print(c(r, "neither 6"))
  }
}


mean(r6_censorsave[!is.na(r6_bsave_mpl[,3])],1)
mean(r6_censorsave[!is.na(r6_bsave_mpl[,3])],2)
r6_censorsave=data.frame(r6_censorsave)
write.table(r6_censorsave, "r6censor.csv",sep=",",dec=".")

colnames(r6_bsave_mpl)=c("b0","b1","seb0","seb1")
r6_bsave_mpl=data.frame(r6_bsave_mpl)
write.table(r6_bsave_mpl, "r6bsave_mpl.csv",sep=",",dec=".",row.names = F)
colnames(r6_bsave_em)=c("b0","b1","seb0","seb1")
r6_bsave_em=data.frame(r6_bsave_em)
write.table(r6_bsave_em, "r6bsave_em.csv",sep=",",dec=".",row.names = F)
colnames(r6_gsave_mpl)=c("g1","seg1")
r6_gsave_mpl=data.frame(r6_gsave_mpl)
write.table(r6_gsave_mpl, "r6gsave_mpl.csv",sep=",",dec=".",row.names = F)
colnames(r6_gsave_em)=c("g1","seg1")
r6_gsave_em=data.frame(r6_gsave_em)
write.table(r6_gsave_em, "r6gsave_em.csv",sep=",",dec=".",row.names = F)


#mpl
mean(r6_bsave_mpl[!is.na(r6_bsave_mpl[,3]),1])
mean(r6_bsave_mpl[!is.na(r6_bsave_mpl[,3]),2])
mean(r6_bsave_mpl[!is.na(r6_bsave_mpl[,3]),3])
sd(r6_bsave_mpl[!is.na(r6_bsave_mpl[,3]),1])
mean(r6_bsave_mpl[!is.na(r6_bsave_mpl[,3]),4])
sd(r6_bsave_mpl[!is.na(r6_bsave_mpl[,3]),2])

mean(r6_gsave_mpl[!is.na(r6_gsave_mpl[,2]),1])
mean(r6_gsave_mpl[!is.na(r6_gsave_mpl[,2]),2])
sd(r6_gsave_mpl[!is.na(r6_gsave_mpl[,2]),1])

r6_bsave_mpl$b0_LL=r6_bsave_mpl$b0_UL=r6_bsave_mpl$b1_LL=r6_bsave_mpl$b1_UL=rep(0,100)
r6_gsave_mpl$g1_LL=r6_gsave_mpl$g1_UL=rep(0,100)

r6_bsave_mpl$b0_LL[!is.na(r6_bsave_mpl$seb0)]=r6_bsave_mpl$b0[!is.na(r6_bsave_mpl$seb0)]-1.96*r6_bsave_mpl$seb0[!is.na(r6_bsave_mpl$seb0)]
r6_bsave_mpl$b0_UL[!is.na(r6_bsave_mpl$seb0)]=r6_bsave_mpl$b0[!is.na(r6_bsave_mpl$seb0)]+1.96*r6_bsave_mpl$seb0[!is.na(r6_bsave_mpl$seb0)]

r6_bsave_mpl$b1_LL[!is.na(r6_bsave_mpl$seb0)]=r6_bsave_mpl$b1[!is.na(r6_bsave_mpl$seb0)]-1.96*r6_bsave_mpl$seb1[!is.na(r6_bsave_mpl$seb0)]
r6_bsave_mpl$b1_UL[!is.na(r6_bsave_mpl$seb0)]=r6_bsave_mpl$b1[!is.na(r6_bsave_mpl$seb0)]+1.96*r6_bsave_mpl$seb1[!is.na(r6_bsave_mpl$seb0)]

r6_bsave_mpl$b0_cov=r6_bsave_mpl$b1_cov=rep(0,100)
r6_bsave_mpl$b0_cov[!is.na(r6_bsave_mpl$seb0)]=as.numeric(r6_bsave_mpl$b0_LL[!is.na(r6_bsave_mpl$seb0)]<=-1 & -1<=r6_bsave_mpl$b0_UL[!is.na(r6_bsave_mpl$seb0)])
r6_bsave_mpl$b1_cov[!is.na(r6_bsave_mpl$seb0)]=as.numeric(r6_bsave_mpl$b1_LL[!is.na(r6_bsave_mpl$seb0)]<=1 & 1<=r6_bsave_mpl$b1_UL[!is.na(r6_bsave_mpl$seb0)])

r6_gsave_mpl$g1_LL[!is.na(r6_gsave_mpl$seg1)]=r6_gsave_mpl$g1[!is.na(r6_gsave_mpl$seg1)]-1.96*r6_gsave_mpl$seg1[!is.na(r6_gsave_mpl$seg1)]
r6_gsave_mpl$g1_UL[!is.na(r6_gsave_mpl$seg1)]=r6_gsave_mpl$g1[!is.na(r6_gsave_mpl$seg1)]+1.96*r6_gsave_mpl$seg1[!is.na(r6_gsave_mpl$seg1)]
r6_gsave_mpl$g1_cov=rep(0,100)
r6_gsave_mpl$g1_cov[!is.na(r6_gsave_mpl$seg1)]=as.numeric(r6_gsave_mpl$g1_LL[!is.na(r6_gsave_mpl$seg1)]<=0.5 & 0.5<=r6_gsave_mpl$g1_UL[!is.na(r6_gsave_mpl$seg1)])


#em

mean(r6_bsave_em[,1])
mean(r6_bsave_em[,2])
mean(r6_bsave_em[,3])
sd(r6_bsave_em[,1])
mean(r6_bsave_em[,4])
sd(r6_bsave_em[,2])

mean(r6_gsave_em[,1])
mean(r6_gsave_em[,2])
sd(r6_gsave_em[,1])

r6_bsave_em$b0_LL=r6_bsave_em$b0_UL=r6_bsave_em$b1_LL=r6_bsave_em$b1_UL=rep(0,100)
r6_gsave_em$g1_LL=r6_gsave_em$g1_UL=rep(0,100)

r6_bsave_em$b0_LL=r6_bsave_em$b0-1.96*r6_bsave_em$seb0
r6_bsave_em$b0_UL=r6_bsave_em$b0+1.96*r6_bsave_em$seb0

r6_bsave_em$b1_LL=r6_bsave_em$b1-1.96*r6_bsave_em$seb1
r6_bsave_em$b1_UL=r6_bsave_em$b1+1.96*r6_bsave_em$seb1

r6_bsave_em$b0_cov=r6_bsave_em$b1_cov=rep(0,100)
r6_bsave_em$b0_cov=as.numeric(r6_bsave_em$b0_LL<=-1 & -1<=r6_bsave_em$b0_UL)
r6_bsave_em$b1_cov=as.numeric(r6_bsave_em$b1_LL<=1 & 1<=r6_bsave_em$b1_UL)

r6_gsave_em$g1_LL=r6_gsave_em$g1-1.96*r6_gsave_em$seg1
r6_gsave_em$g1_UL=r6_gsave_em$g1+1.96*r6_gsave_em$seg1
r6_gsave_em$g1_cov=rep(0,100)
r6_gsave_em$g1_cov=as.numeric(r6_gsave_em$g1_LL<=0.5 & 0.5<=r6_gsave_em$g1_UL)


# 7. beta = [1,1], n=100, censor 0.5 -------


r7_censorsave=matrix(rep(0,2*100),ncol=2)#whole dataset event proportion, non-cured event proportion
r7_bsave_mpl=matrix(rep(0,4*100),ncol=4)
r7_gsave_mpl=matrix(rep(0,2*100),ncol=2)
r7_bsave_em=matrix(rep(0,4*100),ncol=4)
r7_gsave_em=matrix(rep(0,2*100),ncol=2)


for(r in 1:100){
  sim.data=gen_data_right(100,0.5,c(1,1),c(0.5),c("b"),c("b"),c("z1","x1"))
  
  sim.data.uncure=sim.data[sim.data[,3]==1,]
  sim.surv1=Surv(time=sim.data.uncure[,1],event=sim.data.uncure[,2])
  sim.surv=Surv(time=sim.data$time,event =sim.data$event)
  r7_censorsave[r,1]=sum(sim.surv[,2])/100
  r7_censorsave[r,2]=sum(sim.surv1[,2])/nrow(sim.surv1)
  
  try1=try(phmc_mpl(sim.surv~sim.data$X1,pi.formula = ~sim.data$Z,data=sim.data,phmc_mpl.control(conv_limit = 1e-4, n.knots=c(1,0),maxIter = c(1,10000,10001))))
  try2=try(smcure(sim.surv~X1, cureform=~Z, data=sim.data, model="ph"))
  if(class(try1)!="try-error" & class(try2)!="try-error"){
    test=try1
    smc1=try2
    
    r7_bsave_mpl[r,1]=test$beta[1]
    r7_bsave_mpl[r,2]=test$beta[2]
    r7_bsave_mpl[r,3]=test$se$se_H[1]
    r7_bsave_mpl[r,4]=test$se$se_H[2]
    
    r7_gsave_mpl[r,1]=test$gamma[1]
    r7_gsave_mpl[r,2]=test$se$se_H[3]

    r7_bsave_em[r,1]=smc1$b[1]
    r7_bsave_em[r,2]=smc1$b[2]
    r7_bsave_em[r,3]=smc1$b_sd[1]
    r7_bsave_em[r,4]=smc1$b_sd[2]
    
    r7_gsave_em[r,1]=smc1$beta[1]
    r7_gsave_em[r,2]=smc1$beta_sd[1]
    
    print(c(r,"both",test$se$se_H[1]))
  } else if(class(try1)!="try-error"){
    test=try1
    r7_bsave_mpl[r,1]=test$beta[1]
    r7_bsave_mpl[r,2]=test$beta[2]
    r7_bsave_mpl[r,3]=test$se$se_H[1]
    r7_bsave_mpl[r,4]=test$se$se_H[2]
    
    r7_gsave_mpl[r,1]=test$gamma[1]
    r7_gsave_mpl[r,2]=test$se$se_H[3]
    print(c(r, "mpl",test$se$se_H[1]))
    
  } else if(class(try2)!="try-error"){
    smc1=try2
    r7_bsave_em[r,1]=smc1$b[1]
    r7_bsave_em[r,2]=smc1$b[2]
    r7_bsave_em[r,3]=smc1$b_sd[1]
    r7_bsave_em[r,4]=smc1$b_sd[2]
    
    r7_gsave_em[r,1]=smc1$beta[1]
    r7_gsave_em[r,2]=smc1$beta_sd[1]
    print(c(r, "em"))
  } else{
    print(c(r, "neither"))
  }
}


mean(r7_censorsave[!is.na(r7_bsave_mpl[,3])],1)
mean(r7_censorsave[!is.na(r7_bsave_mpl[,3])],2)
r7_censorsave=data.frame(r7_censorsave)
write.table(r7_censorsave, "r7censor.csv",sep=",",dec=".")

colnames(r7_bsave_mpl)=c("b0","b1","seb0","seb1")
r7_bsave_mpl=data.frame(r7_bsave_mpl)
write.table(r7_bsave_mpl, "r7bsave_mpl.csv",sep=",",dec=".",row.names = F)
colnames(r7_bsave_em)=c("b0","b1","seb0","seb1")
r7_bsave_em=data.frame(r7_bsave_em)
write.table(r7_bsave_em, "r7bsave_em.csv",sep=",",dec=".",row.names = F)
colnames(r7_gsave_mpl)=c("g1","seg1")
r7_gsave_mpl=data.frame(r7_gsave_mpl)
write.table(r7_gsave_mpl, "r7gsave_mpl.csv",sep=",",dec=".",row.names = F)
colnames(r7_gsave_em)=c("g1","seg1")
r7_gsave_em=data.frame(r7_gsave_em)
write.table(r7_gsave_em, "r7gsave_em.csv",sep=",",dec=".",row.names = F)


#mpl
mean(r7_bsave_mpl[!is.na(r7_bsave_mpl[,3]),1])
mean(r7_bsave_mpl[!is.na(r7_bsave_mpl[,3]),2])
mean(r7_bsave_mpl[!is.na(r7_bsave_mpl[,3]),3])
sd(r7_bsave_mpl[!is.na(r7_bsave_mpl[,3]),1])
mean(r7_bsave_mpl[!is.na(r7_bsave_mpl[,3]),4])
sd(r7_bsave_mpl[!is.na(r7_bsave_mpl[,3]),2])

mean(r7_gsave_mpl[!is.na(r7_gsave_mpl[,2]),1])
mean(r7_gsave_mpl[!is.na(r7_gsave_mpl[,2]),2])
sd(r7_gsave_mpl[!is.na(r7_gsave_mpl[,2]),1])

r7_bsave_mpl$b0_LL=r7_bsave_mpl$b0_UL=r7_bsave_mpl$b1_LL=r7_bsave_mpl$b1_UL=rep(0,100)
r7_gsave_mpl$g1_LL=r7_gsave_mpl$g1_UL=rep(0,100)

r7_bsave_mpl$b0_LL[!is.na(r7_bsave_mpl$seb0)]=r7_bsave_mpl$b0[!is.na(r7_bsave_mpl$seb0)]-1.96*r7_bsave_mpl$seb0[!is.na(r7_bsave_mpl$seb0)]
r7_bsave_mpl$b0_UL[!is.na(r7_bsave_mpl$seb0)]=r7_bsave_mpl$b0[!is.na(r7_bsave_mpl$seb0)]+1.96*r7_bsave_mpl$seb0[!is.na(r7_bsave_mpl$seb0)]

r7_bsave_mpl$b1_LL[!is.na(r7_bsave_mpl$seb0)]=r7_bsave_mpl$b1[!is.na(r7_bsave_mpl$seb0)]-1.96*r7_bsave_mpl$seb1[!is.na(r7_bsave_mpl$seb0)]
r7_bsave_mpl$b1_UL[!is.na(r7_bsave_mpl$seb0)]=r7_bsave_mpl$b1[!is.na(r7_bsave_mpl$seb0)]+1.96*r7_bsave_mpl$seb1[!is.na(r7_bsave_mpl$seb0)]

r7_bsave_mpl$b0_cov=r7_bsave_mpl$b1_cov=rep(0,100)
r7_bsave_mpl$b0_cov[!is.na(r7_bsave_mpl$seb0)]=as.numeric(r7_bsave_mpl$b0_LL[!is.na(r7_bsave_mpl$seb0)]<=1 & 1<=r7_bsave_mpl$b0_UL[!is.na(r7_bsave_mpl$seb0)])
r7_bsave_mpl$b1_cov[!is.na(r7_bsave_mpl$seb0)]=as.numeric(r7_bsave_mpl$b1_LL[!is.na(r7_bsave_mpl$seb0)]<=1 & 1<=r7_bsave_mpl$b1_UL[!is.na(r7_bsave_mpl$seb0)])

r7_gsave_mpl$g1_LL[!is.na(r7_gsave_mpl$seg1)]=r7_gsave_mpl$g1[!is.na(r7_gsave_mpl$seg1)]-1.96*r7_gsave_mpl$seg1[!is.na(r7_gsave_mpl$seg1)]
r7_gsave_mpl$g1_UL[!is.na(r7_gsave_mpl$seg1)]=r7_gsave_mpl$g1[!is.na(r7_gsave_mpl$seg1)]+1.96*r7_gsave_mpl$seg1[!is.na(r7_gsave_mpl$seg1)]
r7_gsave_mpl$g1_cov=rep(0,100)
r7_gsave_mpl$g1_cov[!is.na(r7_gsave_mpl$seg1)]=as.numeric(r7_gsave_mpl$g1_LL[!is.na(r7_gsave_mpl$seg1)]<=0.5 & 0.5<=r7_gsave_mpl$g1_UL[!is.na(r7_gsave_mpl$seg1)])


#em

mean(r7_bsave_em[,1])
mean(r7_bsave_em[,2])
mean(r7_bsave_em[,3])
sd(r7_bsave_em[,1])
mean(r7_bsave_em[,4])
sd(r7_bsave_em[,2])

mean(r7_gsave_em[,1])
mean(r7_gsave_em[,2])
sd(r7_gsave_em[,1])

r7_bsave_em$b0_LL=r7_bsave_em$b0_UL=r7_bsave_em$b1_LL=r7_bsave_em$b1_UL=rep(0,100)
r7_gsave_em$g1_LL=r7_gsave_em$g1_UL=rep(0,100)

r7_bsave_em$b0_LL=r7_bsave_em$b0-1.96*r7_bsave_em$seb0
r7_bsave_em$b0_UL=r7_bsave_em$b0+1.96*r7_bsave_em$seb0

r7_bsave_em$b1_LL=r7_bsave_em$b1-1.96*r7_bsave_em$seb1
r7_bsave_em$b1_UL=r7_bsave_em$b1+1.96*r7_bsave_em$seb1

r7_bsave_em$b0_cov=r7_bsave_em$b1_cov=rep(0,100)
r7_bsave_em$b0_cov=as.numeric(r7_bsave_em$b0_LL<=1 & 1<=r7_bsave_em$b0_UL)
r7_bsave_em$b1_cov=as.numeric(r7_bsave_em$b1_LL<=1 & 1<=r7_bsave_em$b1_UL)

r7_gsave_em$g1_LL=r7_gsave_em$g1-1.96*r7_gsave_em$seg1
r7_gsave_em$g1_UL=r7_gsave_em$g1+1.96*r7_gsave_em$seg1
r7_gsave_em$g1_cov=rep(0,100)
r7_gsave_em$g1_cov=as.numeric(r7_gsave_em$g1_LL<=0.5 & 0.5<=r7_gsave_em$g1_UL)


# 8. beta = [0,1], n=100, censor 0.5 ------

 
r8_censorsave=matrix(rep(0,2*100),ncol=2)#whole dataset event proportion, non-cured event proportion
r8_bsave_mpl=matrix(rep(0,4*100),ncol=4)
r8_gsave_mpl=matrix(rep(0,2*100),ncol=2)
r8_bsave_em=matrix(rep(0,4*100),ncol=4)
r8_gsave_em=matrix(rep(0,2*100),ncol=2)


for(r in 1:100){
  sim.data=gen_data_right(100,0.5,c(0,1),c(0.5),c("b"),c("b"),c("z1","x1"))
  
  sim.data.uncure=sim.data[sim.data[,3]==1,]
  sim.surv1=Surv(time=sim.data.uncure[,1],event=sim.data.uncure[,2])
  sim.surv=Surv(time=sim.data$time,event =sim.data$event)
  r8_censorsave[r,1]=sum(sim.surv[,2])/100
  r8_censorsave[r,2]=sum(sim.surv1[,2])/nrow(sim.surv1)
  
  try1=try(phmc_mpl(sim.surv~sim.data$X1,pi.formula = ~sim.data$Z,data=sim.data,phmc_mpl.control(conv_limit = 1e-4, n.knots=c(1,0),maxIter = c(1,10000,10001))))
  try2=try(smcure(sim.surv~X1, cureform=~Z, data=sim.data, model="ph"))
  if(class(try1)!="try-error" & class(try2)!="try-error"){
    test=try1
    smc1=try2
    
    r8_bsave_mpl[r,1]=test$beta[1]
    r8_bsave_mpl[r,2]=test$beta[2]
    r8_bsave_mpl[r,3]=test$se$se_H[1]
    r8_bsave_mpl[r,4]=test$se$se_H[2]
    
    r8_gsave_mpl[r,1]=test$gamma[1]
    r8_gsave_mpl[r,2]=test$se$se_H[3]
    
    r8_bsave_em[r,1]=smc1$b[1]
    r8_bsave_em[r,2]=smc1$b[2]
    r8_bsave_em[r,3]=smc1$b_sd[1]
    r8_bsave_em[r,4]=smc1$b_sd[2]
    
    r8_gsave_em[r,1]=smc1$beta[1]
    r8_gsave_em[r,2]=smc1$beta_sd[1]
    
    print(c(r,"both",test$se$se_H[1]))
  } else if(class(try1)!="try-error"){
    test=try1
    r8_bsave_mpl[r,1]=test$beta[1]
    r8_bsave_mpl[r,2]=test$beta[2]
    r8_bsave_mpl[r,3]=test$se$se_H[1]
    r8_bsave_mpl[r,4]=test$se$se_H[2]
    
    r8_gsave_mpl[r,1]=test$gamma[1]
    r8_gsave_mpl[r,2]=test$se$se_H[3]
    print(c(r, "mpl",test$se$se_H[1]))
    
  } else if(class(try2)!="try-error"){
    smc1=try2
    r8_bsave_em[r,1]=smc1$b[1]
    r8_bsave_em[r,2]=smc1$b[2]
    r8_bsave_em[r,3]=smc1$b_sd[1]
    r8_bsave_em[r,4]=smc1$b_sd[2]
    
    r8_gsave_em[r,1]=smc1$beta[1]
    r8_gsave_em[r,2]=smc1$beta_sd[1]
    print(c(r, "em"))
  } else{
    print(c(r, "neither"))
  }
}


mean(r8_censorsave[!is.na(r8_bsave_mpl[,3])],1)
mean(r8_censorsave[!is.na(r8_bsave_mpl[,3])],2)
r8_censorsave=data.frame(r8_censorsave)
write.table(r8_censorsave, "r8censor.csv",sep=",",dec=".")

colnames(r8_bsave_mpl)=c("b0","b1","seb0","seb1")
r8_bsave_mpl=data.frame(r8_bsave_mpl)
write.table(r8_bsave_mpl, "r8bsave_mpl.csv",sep=",",dec=".",row.names = F)
colnames(r8_bsave_em)=c("b0","b1","seb0","seb1")
r8_bsave_em=data.frame(r8_bsave_em)
write.table(r8_bsave_em, "r8bsave_em.csv",sep=",",dec=".",row.names = F)
colnames(r8_gsave_mpl)=c("g1","seg1")
r8_gsave_mpl=data.frame(r8_gsave_mpl)
write.table(r8_gsave_mpl, "r8gsave_mpl.csv",sep=",",dec=".",row.names = F)
colnames(r8_gsave_em)=c("g1","seg1")
r8_gsave_em=data.frame(r8_gsave_em)
write.table(r8_gsave_em, "r8gsave_em.csv",sep=",",dec=".",row.names = F)


#mpl
mean(r8_bsave_mpl[!is.na(r8_bsave_mpl[,3]),1])
mean(r8_bsave_mpl[!is.na(r8_bsave_mpl[,3]),2])
mean(r8_bsave_mpl[!is.na(r8_bsave_mpl[,3]),3])
sd(r8_bsave_mpl[!is.na(r8_bsave_mpl[,3]),1])
mean(r8_bsave_mpl[!is.na(r8_bsave_mpl[,3]),4])
sd(r8_bsave_mpl[!is.na(r8_bsave_mpl[,3]),2])

mean(r8_gsave_mpl[!is.na(r8_gsave_mpl[,2]),1])
mean(r8_gsave_mpl[!is.na(r8_gsave_mpl[,2]),2])
sd(r8_gsave_mpl[!is.na(r8_gsave_mpl[,2]),1])

r8_bsave_mpl$b0_LL=r8_bsave_mpl$b0_UL=r8_bsave_mpl$b1_LL=r8_bsave_mpl$b1_UL=rep(0,100)
r8_gsave_mpl$g1_LL=r8_gsave_mpl$g1_UL=rep(0,100)

r8_bsave_mpl$b0_LL[!is.na(r8_bsave_mpl$seb0)]=r8_bsave_mpl$b0[!is.na(r8_bsave_mpl$seb0)]-1.96*r8_bsave_mpl$seb0[!is.na(r8_bsave_mpl$seb0)]
r8_bsave_mpl$b0_UL[!is.na(r8_bsave_mpl$seb0)]=r8_bsave_mpl$b0[!is.na(r8_bsave_mpl$seb0)]+1.96*r8_bsave_mpl$seb0[!is.na(r8_bsave_mpl$seb0)]

r8_bsave_mpl$b1_LL[!is.na(r8_bsave_mpl$seb0)]=r8_bsave_mpl$b1[!is.na(r8_bsave_mpl$seb0)]-1.96*r8_bsave_mpl$seb1[!is.na(r8_bsave_mpl$seb0)]
r8_bsave_mpl$b1_UL[!is.na(r8_bsave_mpl$seb0)]=r8_bsave_mpl$b1[!is.na(r8_bsave_mpl$seb0)]+1.96*r8_bsave_mpl$seb1[!is.na(r8_bsave_mpl$seb0)]

r8_bsave_mpl$b0_cov=r8_bsave_mpl$b1_cov=rep(0,100)
r8_bsave_mpl$b0_cov[!is.na(r8_bsave_mpl$seb0)]=as.numeric(r8_bsave_mpl$b0_LL[!is.na(r8_bsave_mpl$seb0)]<=0 & 0<=r8_bsave_mpl$b0_UL[!is.na(r8_bsave_mpl$seb0)])
r8_bsave_mpl$b1_cov[!is.na(r8_bsave_mpl$seb0)]=as.numeric(r8_bsave_mpl$b1_LL[!is.na(r8_bsave_mpl$seb0)]<=1 & 1<=r8_bsave_mpl$b1_UL[!is.na(r8_bsave_mpl$seb0)])

r8_gsave_mpl$g1_LL[!is.na(r8_gsave_mpl$seg1)]=r8_gsave_mpl$g1[!is.na(r8_gsave_mpl$seg1)]-1.96*r8_gsave_mpl$seg1[!is.na(r8_gsave_mpl$seg1)]
r8_gsave_mpl$g1_UL[!is.na(r8_gsave_mpl$seg1)]=r8_gsave_mpl$g1[!is.na(r8_gsave_mpl$seg1)]+1.96*r8_gsave_mpl$seg1[!is.na(r8_gsave_mpl$seg1)]
r8_gsave_mpl$g1_cov=rep(0,100)
r8_gsave_mpl$g1_cov[!is.na(r8_gsave_mpl$seg1)]=as.numeric(r8_gsave_mpl$g1_LL[!is.na(r8_gsave_mpl$seg1)]<=0.5 & 0.5<=r8_gsave_mpl$g1_UL[!is.na(r8_gsave_mpl$seg1)])


#em

mean(r8_bsave_em[,1])
mean(r8_bsave_em[,2])
mean(r8_bsave_em[,3])
sd(r8_bsave_em[,1])
mean(r8_bsave_em[,4])
sd(r8_bsave_em[,2])

mean(r8_gsave_em[,1])
mean(r8_gsave_em[,2])
sd(r8_gsave_em[,1])

r8_bsave_em$b0_LL=r8_bsave_em$b0_UL=r8_bsave_em$b1_LL=r8_bsave_em$b1_UL=rep(0,100)
r8_gsave_em$g1_LL=r8_gsave_em$g1_UL=rep(0,100)

r8_bsave_em$b0_LL=r8_bsave_em$b0-1.96*r8_bsave_em$seb0
r8_bsave_em$b0_UL=r8_bsave_em$b0+1.96*r8_bsave_em$seb0

r8_bsave_em$b1_LL=r8_bsave_em$b1-1.96*r8_bsave_em$seb1
r8_bsave_em$b1_UL=r8_bsave_em$b1+1.96*r8_bsave_em$seb1

r8_bsave_em$b0_cov=r8_bsave_em$b1_cov=rep(0,100)
r8_bsave_em$b0_cov=as.numeric(r8_bsave_em$b0_LL<=0 & 0<=r8_bsave_em$b0_UL)
r8_bsave_em$b1_cov=as.numeric(r8_bsave_em$b1_LL<=1 & 1<=r8_bsave_em$b1_UL)

r8_gsave_em$g1_LL=r8_gsave_em$g1-1.96*r8_gsave_em$seg1
r8_gsave_em$g1_UL=r8_gsave_em$g1+1.96*r8_gsave_em$seg1
r8_gsave_em$g1_cov=rep(0,100)
r8_gsave_em$g1_cov=as.numeric(r8_gsave_em$g1_LL<=0.5 & 0.5<=r8_gsave_em$g1_UL)


# 9. beta = [-1,1], n=100, censor 0.5 -------


r9_censorsave=matrix(rep(0,2*100),ncol=2)#whole dataset event proportion, non-cured event proportion
r9_bsave_mpl=matrix(rep(0,4*100),ncol=4)
r9_gsave_mpl=matrix(rep(0,2*100),ncol=2)
r9_bsave_em=matrix(rep(0,4*100),ncol=4)
r9_gsave_em=matrix(rep(0,2*100),ncol=2)


for(r in 23:100){
  sim.data=gen_data_right(100,0.5,c(-1,1),c(0.5),c("b"),c("b"),c("z1","x1"))
  
  sim.data.uncure=sim.data[sim.data[,3]==1,]
  sim.surv1=Surv(time=sim.data.uncure[,1],event=sim.data.uncure[,2])
  sim.surv=Surv(time=sim.data$time,event =sim.data$event)
  r9_censorsave[r,1]=sum(sim.surv[,2])/100
  r9_censorsave[r,2]=sum(sim.surv1[,2])/nrow(sim.surv1)
  
  try1=try(phmc_mpl(sim.surv~sim.data$X1,pi.formula = ~sim.data$Z,data=sim.data,phmc_mpl.control(conv_limit = 1e-4, n.knots=c(1,0),maxIter = c(1,30000,30001))))
  try2=try(smcure(sim.surv~X1, cureform=~Z, data=sim.data, model="ph"))
  if(class(try1)!="try-error" & class(try2)!="try-error"){
    test=try1
    smc1=try2
    
    r9_bsave_mpl[r,1]=test$beta[1]
    r9_bsave_mpl[r,2]=test$beta[2]
    r9_bsave_mpl[r,3]=test$se$se_H[1]
    r9_bsave_mpl[r,4]=test$se$se_H[2]
    
    r9_gsave_mpl[r,1]=test$gamma[1]
    r9_gsave_mpl[r,2]=test$se$se_H[3]
    
    r9_bsave_em[r,1]=smc1$b[1]
    r9_bsave_em[r,2]=smc1$b[2]
    r9_bsave_em[r,3]=smc1$b_sd[1]
    r9_bsave_em[r,4]=smc1$b_sd[2]
    
    r9_gsave_em[r,1]=smc1$beta[1]
    r9_gsave_em[r,2]=smc1$beta_sd[1]
    
    print(c(r,"both y",test$se$se_H[1]))
  } else if(class(try1)!="try-error"){
    test=try1
    r9_bsave_mpl[r,1]=test$beta[1]
    r9_bsave_mpl[r,2]=test$beta[2]
    r9_bsave_mpl[r,3]=test$se$se_H[1]
    r9_bsave_mpl[r,4]=test$se$se_H[2]
    
    r9_gsave_mpl[r,1]=test$gamma[1]
    r9_gsave_mpl[r,2]=test$se$se_H[3]
    print(c(r, "mpl y",test$se$se_H[1]))
    
  } else if(class(try2)!="try-error"){
    smc1=try2
    r9_bsave_em[r,1]=smc1$b[1]
    r9_bsave_em[r,2]=smc1$b[2]
    r9_bsave_em[r,3]=smc1$b_sd[1]
    r9_bsave_em[r,4]=smc1$b_sd[2]
    
    r9_gsave_em[r,1]=smc1$beta[1]
    r9_gsave_em[r,2]=smc1$beta_sd[1]
    print(c(r, "em y"))
  } else{
    print(c(r, "neither y"))
  }
}


mean(r9_censorsave[!is.na(r9_bsave_mpl[,3])],1)
mean(r9_censorsave[!is.na(r9_bsave_mpl[,3])],2)
r9_censorsave=data.frame(r9_censorsave)
write.table(r9_censorsave, "r9censor.csv",sep=",",dec=".")

colnames(r9_bsave_mpl)=c("b0","b1","seb0","seb1")
r9_bsave_mpl=data.frame(r9_bsave_mpl)
write.table(r9_bsave_mpl, "r9bsave_mpl.csv",sep=",",dec=".",row.names = F)
colnames(r9_bsave_em)=c("b0","b1","seb0","seb1")
r9_bsave_em=data.frame(r9_bsave_em)
write.table(r9_bsave_em, "r9bsave_em.csv",sep=",",dec=".",row.names = F)
colnames(r9_gsave_mpl)=c("g1","seg1")
r9_gsave_mpl=data.frame(r9_gsave_mpl)
write.table(r9_gsave_mpl, "r9gsave_mpl.csv",sep=",",dec=".",row.names = F)
colnames(r9_gsave_em)=c("g1","seg1")
r9_gsave_em=data.frame(r9_gsave_em)
write.table(r9_gsave_em, "r9gsave_em.csv",sep=",",dec=".",row.names = F)


#mpl
mean(r9_bsave_mpl[!is.na(r9_bsave_mpl[,3]),1])
mean(r9_bsave_mpl[!is.na(r9_bsave_mpl[,3]),2])
mean(r9_bsave_mpl[!is.na(r9_bsave_mpl[,3]),3])
sd(r9_bsave_mpl[!is.na(r9_bsave_mpl[,3]),1])
mean(r9_bsave_mpl[!is.na(r9_bsave_mpl[,3]),4])
sd(r9_bsave_mpl[!is.na(r9_bsave_mpl[,3]),2])

mean(r9_gsave_mpl[!is.na(r9_gsave_mpl[,2]),1])
mean(r9_gsave_mpl[!is.na(r9_gsave_mpl[,2]),2])
sd(r9_gsave_mpl[!is.na(r9_gsave_mpl[,2]),1])

r9_bsave_mpl$b0_LL=r9_bsave_mpl$b0_UL=r9_bsave_mpl$b1_LL=r9_bsave_mpl$b1_UL=rep(0,100)
r9_gsave_mpl$g1_LL=r9_gsave_mpl$g1_UL=rep(0,100)

r9_bsave_mpl$b0_LL[!is.na(r9_bsave_mpl$seb0)]=r9_bsave_mpl$b0[!is.na(r9_bsave_mpl$seb0)]-1.96*r9_bsave_mpl$seb0[!is.na(r9_bsave_mpl$seb0)]
r9_bsave_mpl$b0_UL[!is.na(r9_bsave_mpl$seb0)]=r9_bsave_mpl$b0[!is.na(r9_bsave_mpl$seb0)]+1.96*r9_bsave_mpl$seb0[!is.na(r9_bsave_mpl$seb0)]

r9_bsave_mpl$b1_LL[!is.na(r9_bsave_mpl$seb0)]=r9_bsave_mpl$b1[!is.na(r9_bsave_mpl$seb0)]-1.96*r9_bsave_mpl$seb1[!is.na(r9_bsave_mpl$seb0)]
r9_bsave_mpl$b1_UL[!is.na(r9_bsave_mpl$seb0)]=r9_bsave_mpl$b1[!is.na(r9_bsave_mpl$seb0)]+1.96*r9_bsave_mpl$seb1[!is.na(r9_bsave_mpl$seb0)]

r9_bsave_mpl$b0_cov=r9_bsave_mpl$b1_cov=rep(0,100)
r9_bsave_mpl$b0_cov[!is.na(r9_bsave_mpl$seb0)]=as.numeric(r9_bsave_mpl$b0_LL[!is.na(r9_bsave_mpl$seb0)]<=-1 & -1<=r9_bsave_mpl$b0_UL[!is.na(r9_bsave_mpl$seb0)])
r9_bsave_mpl$b1_cov[!is.na(r9_bsave_mpl$seb0)]=as.numeric(r9_bsave_mpl$b1_LL[!is.na(r9_bsave_mpl$seb0)]<=1 & 1<=r9_bsave_mpl$b1_UL[!is.na(r9_bsave_mpl$seb0)])

r9_gsave_mpl$g1_LL[!is.na(r9_gsave_mpl$seg1)]=r9_gsave_mpl$g1[!is.na(r9_gsave_mpl$seg1)]-1.96*r9_gsave_mpl$seg1[!is.na(r9_gsave_mpl$seg1)]
r9_gsave_mpl$g1_UL[!is.na(r9_gsave_mpl$seg1)]=r9_gsave_mpl$g1[!is.na(r9_gsave_mpl$seg1)]+1.96*r9_gsave_mpl$seg1[!is.na(r9_gsave_mpl$seg1)]
r9_gsave_mpl$g1_cov=rep(0,100)
r9_gsave_mpl$g1_cov[!is.na(r9_gsave_mpl$seg1)]=as.numeric(r9_gsave_mpl$g1_LL[!is.na(r9_gsave_mpl$seg1)]<=0.5 & 0.5<=r9_gsave_mpl$g1_UL[!is.na(r9_gsave_mpl$seg1)])


#em

mean(r9_bsave_em[,1])
mean(r9_bsave_em[,2])
mean(r9_bsave_em[,3])
sd(r9_bsave_em[,1])
mean(r9_bsave_em[,4])
sd(r9_bsave_em[,2])

mean(r9_gsave_em[,1])
mean(r9_gsave_em[,2])
sd(r9_gsave_em[,1])

r9_bsave_em$b0_LL=r9_bsave_em$b0_UL=r9_bsave_em$b1_LL=r9_bsave_em$b1_UL=rep(0,100)
r9_gsave_em$g1_LL=r9_gsave_em$g1_UL=rep(0,100)

r9_bsave_em$b0_LL=r9_bsave_em$b0-1.96*r9_bsave_em$seb0
r9_bsave_em$b0_UL=r9_bsave_em$b0+1.96*r9_bsave_em$seb0

r9_bsave_em$b1_LL=r9_bsave_em$b1-1.96*r9_bsave_em$seb1
r9_bsave_em$b1_UL=r9_bsave_em$b1+1.96*r9_bsave_em$seb1

r9_bsave_em$b0_cov=r9_bsave_em$b1_cov=rep(0,100)
r9_bsave_em$b0_cov=as.numeric(r9_bsave_em$b0_LL<=-1 & -1<=r9_bsave_em$b0_UL)
r9_bsave_em$b1_cov=as.numeric(r9_bsave_em$b1_LL<=1 & 1<=r9_bsave_em$b1_UL)

r9_gsave_em$g1_LL=r9_gsave_em$g1-1.96*r9_gsave_em$seg1
r9_gsave_em$g1_UL=r9_gsave_em$g1+1.96*r9_gsave_em$seg1
r9_gsave_em$g1_cov=rep(0,100)
r9_gsave_em$g1_cov=as.numeric(r9_gsave_em$g1_LL<=0.5 & 0.5<=r9_gsave_em$g1_UL)


# 10. beta = [1,1], n=2000, censor 1.25 --------


r10_censorsave=matrix(rep(0,2*100),ncol=2)#whole dataset event proportion, non-cured event proportion
r10_bsave_mpl=matrix(rep(0,4*100),ncol=4)
r10_gsave_mpl=matrix(rep(0,2*100),ncol=2)
r10_bsave_em=matrix(rep(0,4*100),ncol=4)
r10_gsave_em=matrix(rep(0,2*100),ncol=2)


for(r in 18:100){
  sim.data=gen_data_right(2000,1.25,c(1,1),c(0.5),c("b"),c("b"),c("z1","x1"))
  
  sim.data.uncure=sim.data[sim.data[,3]==1,]
  sim.surv1=Surv(time=sim.data.uncure[,1],event=sim.data.uncure[,2])
  sim.surv=Surv(time=sim.data$time,event =sim.data$event)
  
  r10_censorsave[r,1]=sum(sim.surv[,2])/2000
  r10_censorsave[r,2]=sum(sim.surv1[,2])/nrow(sim.surv1)
  
  test=phmc_mpl(sim.surv~sim.data$X1,pi.formula = ~sim.data$Z,data=sim.data,phmc_mpl.control(conv_limit = 1e-4, n.knots=c(8,0),maxIter = c(1,10000,10001)))
  
  r10_bsave_mpl[r,1]=test$beta[1]
  r10_bsave_mpl[r,2]=test$beta[2]
  r10_bsave_mpl[r,3]=test$se$se_H[1]
  r10_bsave_mpl[r,4]=test$se$se_H[2]
  
  r10_gsave_mpl[r,1]=test$gamma[1]
  r10_gsave_mpl[r,2]=test$se$se_H[3]
  
  smc1 = smcure(sim.surv~X1, cureform=~Z, data=sim.data, model="ph")
  
  r10_bsave_em[r,1]=smc1$b[1]
  r10_bsave_em[r,2]=smc1$b[2]
  r10_bsave_em[r,3]=smc1$b_sd[1]
  r10_bsave_em[r,4]=smc1$b_sd[2]
  
  r10_gsave_em[r,1]=smc1$beta[1]
  r10_gsave_em[r,2]=smc1$beta_sd[1]
  
  print(c(r,test$se$se_H[1]))
  
}


mean(r10_censorsave[!is.na(r10_bsave_mpl[,3])],1)
mean(r10_censorsave[!is.na(r10_bsave_mpl[,3])],2)
r10_censorsave=data.frame(r10_censorsave)
write.table(r10_censorsave, "r10censor.csv",sep=",",dec=".")

colnames(r10_bsave_mpl)=c("b0","b1","seb0","seb1")
r10_bsave_mpl=data.frame(r10_bsave_mpl)
write.table(r10_bsave_mpl, "r10bsave_mpl.csv",sep=",",dec=".",row.names = F)
colnames(r10_bsave_em)=c("b0","b1","seb0","seb1")
r10_bsave_em=data.frame(r10_bsave_em)
write.table(r10_bsave_em, "r10bsave_em.csv",sep=",",dec=".",row.names = F)
colnames(r10_gsave_mpl)=c("g1","seg1")
r10_gsave_mpl=data.frame(r10_gsave_mpl)
write.table(r10_gsave_mpl, "r10gsave_mpl.csv",sep=",",dec=".",row.names = F)
colnames(r10_gsave_em)=c("g1","seg1")
r10_gsave_em=data.frame(r10_gsave_em)
write.table(r10_gsave_em, "r10gsave_em.csv",sep=",",dec=".",row.names = F)


#mpl
mean(r10_bsave_mpl[!is.na(r10_bsave_mpl[,3]),1])
mean(r10_bsave_mpl[!is.na(r10_bsave_mpl[,3]),2])
mean(r10_bsave_mpl[!is.na(r10_bsave_mpl[,3]),3])
sd(r10_bsave_mpl[!is.na(r10_bsave_mpl[,3]),1])
mean(r10_bsave_mpl[!is.na(r10_bsave_mpl[,3]),4])
sd(r10_bsave_mpl[!is.na(r10_bsave_mpl[,3]),2])

mean(r10_gsave_mpl[!is.na(r10_gsave_mpl[,2]),1])
mean(r10_gsave_mpl[!is.na(r10_gsave_mpl[,2]),2])
sd(r10_gsave_mpl[!is.na(r10_gsave_mpl[,2]),1])

r10_bsave_mpl$b0_LL=r10_bsave_mpl$b0_UL=r10_bsave_mpl$b1_LL=r10_bsave_mpl$b1_UL=rep(0,100)
r10_gsave_mpl$g1_LL=r10_gsave_mpl$g1_UL=rep(0,100)

r10_bsave_mpl$b0_LL[!is.na(r10_bsave_mpl$seb0)]=r10_bsave_mpl$b0[!is.na(r10_bsave_mpl$seb0)]-1.96*r10_bsave_mpl$seb0[!is.na(r10_bsave_mpl$seb0)]
r10_bsave_mpl$b0_UL[!is.na(r10_bsave_mpl$seb0)]=r10_bsave_mpl$b0[!is.na(r10_bsave_mpl$seb0)]+1.96*r10_bsave_mpl$seb0[!is.na(r10_bsave_mpl$seb0)]

r10_bsave_mpl$b1_LL[!is.na(r10_bsave_mpl$seb0)]=r10_bsave_mpl$b1[!is.na(r10_bsave_mpl$seb0)]-1.96*r10_bsave_mpl$seb1[!is.na(r10_bsave_mpl$seb0)]
r10_bsave_mpl$b1_UL[!is.na(r10_bsave_mpl$seb0)]=r10_bsave_mpl$b1[!is.na(r10_bsave_mpl$seb0)]+1.96*r10_bsave_mpl$seb1[!is.na(r10_bsave_mpl$seb0)]

r10_bsave_mpl$b0_cov=r10_bsave_mpl$b1_cov=rep(0,100)
r10_bsave_mpl$b0_cov[!is.na(r10_bsave_mpl$seb0)]=as.numeric(r10_bsave_mpl$b0_LL[!is.na(r10_bsave_mpl$seb0)]<=-1 & -1<=r10_bsave_mpl$b0_UL[!is.na(r10_bsave_mpl$seb0)])
r10_bsave_mpl$b1_cov[!is.na(r10_bsave_mpl$seb0)]=as.numeric(r10_bsave_mpl$b1_LL[!is.na(r10_bsave_mpl$seb0)]<=1 & 1<=r10_bsave_mpl$b1_UL[!is.na(r10_bsave_mpl$seb0)])

r10_gsave_mpl$g1_LL[!is.na(r10_gsave_mpl$seg1)]=r10_gsave_mpl$g1[!is.na(r10_gsave_mpl$seg1)]-1.96*r10_gsave_mpl$seg1[!is.na(r10_gsave_mpl$seg1)]
r10_gsave_mpl$g1_UL[!is.na(r10_gsave_mpl$seg1)]=r10_gsave_mpl$g1[!is.na(r10_gsave_mpl$seg1)]+1.96*r10_gsave_mpl$seg1[!is.na(r10_gsave_mpl$seg1)]
r10_gsave_mpl$g1_cov=rep(0,100)
r10_gsave_mpl$g1_cov[!is.na(r10_gsave_mpl$seg1)]=as.numeric(r10_gsave_mpl$g1_LL[!is.na(r10_gsave_mpl$seg1)]<=0.5 & 0.5<=r10_gsave_mpl$g1_UL[!is.na(r10_gsave_mpl$seg1)])


#em

mean(r10_bsave_em[,1])
mean(r10_bsave_em[,2])
mean(r10_bsave_em[,3])
sd(r10_bsave_em[,1])
mean(r10_bsave_em[,4])
sd(r10_bsave_em[,2])

mean(r10_gsave_em[,1])
mean(r10_gsave_em[,2])
sd(r10_gsave_em[,1])

r10_bsave_em$b0_LL=r10_bsave_em$b0_UL=r10_bsave_em$b1_LL=r10_bsave_em$b1_UL=rep(0,100)
r10_gsave_em$g1_LL=r10_gsave_em$g1_UL=rep(0,100)

r10_bsave_em$b0_LL=r10_bsave_em$b0-1.96*r10_bsave_em$seb0
r10_bsave_em$b0_UL=r10_bsave_em$b0+1.96*r10_bsave_em$seb0

r10_bsave_em$b1_LL=r10_bsave_em$b1-1.96*r10_bsave_em$seb1
r10_bsave_em$b1_UL=r10_bsave_em$b1+1.96*r10_bsave_em$seb1

r10_bsave_em$b0_cov=r10_bsave_em$b1_cov=rep(0,100)
r10_bsave_em$b0_cov=as.numeric(r10_bsave_em$b0_LL<=-1 & -1<=r10_bsave_em$b0_UL)
r10_bsave_em$b1_cov=as.numeric(r10_bsave_em$b1_LL<=1 & 1<=r10_bsave_em$b1_UL)

r10_gsave_em$g1_LL=r10_gsave_em$g1-1.96*r10_gsave_em$seg1
r10_gsave_em$g1_UL=r10_gsave_em$g1+1.96*r10_gsave_em$seg1
r10_gsave_em$g1_cov=rep(0,100)
r10_gsave_em$g1_cov=as.numeric(r10_gsave_em$g1_LL<=0.5 & 0.5<=r10_gsave_em$g1_UL)



# 13. beta = [1,1], n=500, censor 1.25 -------

r1_censorsave=matrix(rep(0,2*100),ncol=2)#whole dataset event proportion, non-cured event proportion
r1_bsave_mpl=matrix(rep(0,4*100),ncol=4)
r1_gsave_mpl=matrix(rep(0,2*100),ncol=2)
r1_bsave_em=matrix(rep(0,4*100),ncol=4)
r1_gsave_em=matrix(rep(0,2*100),ncol=2)

for(r in 1:100){
  sim.data=gen_data_right(500,1.25,c(1,1),c(0.5),c("b"),c("b"),c("z1","x1"))
  
  sim.data.uncure=sim.data[sim.data[,3]==1,]
  sim.surv1=Surv(time=sim.data.uncure[,1],event=sim.data.uncure[,2])
  sim.surv=Surv(time=sim.data$time,event =sim.data$event)
  
  #r1_censorsave[r,1]=sum(sim.surv[,2])/500
  #r1_censorsave[r,2]=sum(sim.surv1[,2])/nrow(sim.surv1)
  
  #test=phmc_mpl(sim.surv~sim.data$X1,pi.formula = ~sim.data$Z,data=sim.data,phmc_mpl.control(conv_limit = 1e-4, n.knots=c(6,0),maxIter = c(1,10000,10001)))
  
  #r1_bsave_mpl[r,1]=test$beta[1]
  #r1_bsave_mpl[r,2]=test$beta[2]
  #r1_bsave_mpl[r,3]=test$se$se_H[1]
  #r1_bsave_mpl[r,4]=test$se$se_H[2]
  
  #r1_gsave_mpl[r,1]=test$gamma[1]
  #r1_gsave_mpl[r,2]=test$se$se_H[3]
  
  smc1 = smcure(sim.surv~X1, cureform=~Z, data=sim.data, model="ph")
  
  r1_bsave_em[r,1]=smc1$b[1]
  r1_bsave_em[r,2]=smc1$b[2]
  r1_bsave_em[r,3]=smc1$b_sd[1]
  r1_bsave_em[r,4]=smc1$b_sd[2]
  
  r1_gsave_em[r,1]=smc1$beta[1]
  r1_gsave_em[r,2]=smc1$beta_sd[1]
  
  #print(c(r,sum(sim.surv1[,2])/nrow(sim.surv1),test$beta[1]))
  print(c(r,smc1$b[2],r1_bsave_em[r,2]))
}

mean(r1_censorsave[!is.na(r1_bsave_mpl[,3]),1])
mean(r1_censorsave[!is.na(r1_bsave_mpl[,3]),2])
r1_censorsave=data.frame(r1_censorsave)
write.table(r1_censorsave, "r1censor.csv",sep=",",dec=".")


colnames(r1_bsave_mpl)=c("b0","b1","seb0","seb1")
r1_bsave_mpl=data.frame(r1_bsave_mpl)
colnames(r1_bsave_em)=c("b0","b1","seb0","seb1")
r1_bsave_em=data.frame(r1_bsave_em)
colnames(r1_gsave_mpl)=c("g1","seg1")
r1_gsave_mpl=data.frame(r1_gsave_mpl)
colnames(r1_gsave_em)=c("g1","seg1")
r1_gsave_em=data.frame(r1_gsave_em)




#mpl
mean(r1_bsave_mpl[!is.na(r1_bsave_mpl[,3]),1])
mean(r1_bsave_mpl[!is.na(r1_bsave_mpl[,3]),2])
mean(r1_bsave_mpl[!is.na(r1_bsave_mpl[,3]),3])
sd(r1_bsave_mpl[!is.na(r1_bsave_mpl[,3]),1])
mean(r1_bsave_mpl[!is.na(r1_bsave_mpl[,3]),4])
sd(r1_bsave_mpl[!is.na(r1_bsave_mpl[,3]),2])

mean(r1_gsave_mpl[!is.na(r1_gsave_mpl[,2]),1])
mean(r1_gsave_mpl[!is.na(r1_gsave_mpl[,2]),2])
sd(r1_gsave_mpl[!is.na(r1_gsave_mpl[,2]),1])

r1_bsave_mpl$b0_LL=r1_bsave_mpl$b0_UL=r1_bsave_mpl$b1_LL=r1_bsave_mpl$b1_UL=rep(0,100)
r1_gsave_mpl$g1_LL=r1_gsave_mpl$g1_UL=rep(0,100)

r1_bsave_mpl$b0_LL[!is.na(r1_bsave_mpl$seb0)]=r1_bsave_mpl$b0[!is.na(r1_bsave_mpl$seb0)]-1.96*r1_bsave_mpl$seb0[!is.na(r1_bsave_mpl$seb0)]
r1_bsave_mpl$b0_UL[!is.na(r1_bsave_mpl$seb0)]=r1_bsave_mpl$b0[!is.na(r1_bsave_mpl$seb0)]+1.96*r1_bsave_mpl$seb0[!is.na(r1_bsave_mpl$seb0)]

r1_bsave_mpl$b1_LL[!is.na(r1_bsave_mpl$seb0)]=r1_bsave_mpl$b1[!is.na(r1_bsave_mpl$seb0)]-1.96*r1_bsave_mpl$seb1[!is.na(r1_bsave_mpl$seb0)]
r1_bsave_mpl$b1_UL[!is.na(r1_bsave_mpl$seb0)]=r1_bsave_mpl$b1[!is.na(r1_bsave_mpl$seb0)]+1.96*r1_bsave_mpl$seb1[!is.na(r1_bsave_mpl$seb0)]

r1_bsave_mpl$b0_cov=r1_bsave_mpl$b1_cov=rep(0,100)
r1_bsave_mpl$b0_cov[!is.na(r1_bsave_mpl$seb0)]=as.numeric(r1_bsave_mpl$b0_LL[!is.na(r1_bsave_mpl$seb0)]<=1 & 1<=r1_bsave_mpl$b0_UL[!is.na(r1_bsave_mpl$seb0)])
r1_bsave_mpl$b1_cov[!is.na(r1_bsave_mpl$seb0)]=as.numeric(r1_bsave_mpl$b1_LL[!is.na(r1_bsave_mpl$seb0)]<=1 & 1<=r1_bsave_mpl$b1_UL[!is.na(r1_bsave_mpl$seb0)])


r1_gsave_mpl$g1_LL[!is.na(r1_gsave_mpl$seg1)]=r1_gsave_mpl$g1[!is.na(r1_gsave_mpl$seg1)]-1.96*r1_gsave_mpl$seg1[!is.na(r1_gsave_mpl$seg1)]
r1_gsave_mpl$g1_UL[!is.na(r1_gsave_mpl$seg1)]=r1_gsave_mpl$g1[!is.na(r1_gsave_mpl$seg1)]+1.96*r1_gsave_mpl$seg1[!is.na(r1_gsave_mpl$seg1)]
r1_gsave_mpl$g1_cov=rep(0,100)
r1_gsave_mpl$g1_cov[!is.na(r1_gsave_mpl$seg1)]=as.numeric(r1_gsave_mpl$g1_LL[!is.na(r1_gsave_mpl$seg1)]<=0.5 & 0.5<=r1_gsave_mpl$g1_UL[!is.na(r1_gsave_mpl$seg1)])


write.table(r1_bsave_mpl, "r1bsave_mpl.csv",sep=",",dec=".",row.names = FALSE)
write.table(r1_gsave_mpl, "r1gsave_mpl.csv",sep=",",dec=".",row.names = FALSE)

#em

mean(r1_bsave_em[,1])
mean(r1_bsave_em[,2])
mean(r1_bsave_em[,3])
sd(r1_bsave_em[,1])
mean(r1_bsave_em[,4])
sd(r1_bsave_em[,2])

mean(r1_gsave_em[,1])
mean(r1_gsave_em[,2])
sd(r1_gsave_em[,1])

r1_bsave_em$b0_LL=r1_bsave_em$b0_UL=r1_bsave_em$b1_LL=r1_bsave_em$b1_UL=rep(0,100)
r1_gsave_em$g1_LL=r1_gsave_em$g1_UL=rep(0,100)

r1_bsave_em$b0_LL=r1_bsave_em$b0-1.96*r1_bsave_em$seb0
r1_bsave_em$b0_UL=r1_bsave_em$b0+1.96*r1_bsave_em$seb0

r1_bsave_em$b1_LL=r1_bsave_em$b1-1.96*r1_bsave_em$seb1
r1_bsave_em$b1_UL=r1_bsave_em$b1+1.96*r1_bsave_em$seb1

r1_bsave_em$b0_cov=r1_bsave_em$b1_cov=rep(0,100)
r1_bsave_em$b0_cov=as.numeric(r1_bsave_em$b0_LL<=1 & 1<=r1_bsave_em$b0_UL)
r1_bsave_em$b1_cov=as.numeric(r1_bsave_em$b1_LL<=1 & 1<=r1_bsave_em$b1_UL)

r1_gsave_em$g1_LL=r1_gsave_em$g1-1.96*r1_gsave_em$seg1
r1_gsave_em$g1_UL=r1_gsave_em$g1+1.96*r1_gsave_em$seg1
r1_gsave_em$g1_cov=rep(0,100)
r1_gsave_em$g1_cov=as.numeric(r1_gsave_em$g1_LL<=0.5 & 0.5<=r1_gsave_em$g1_UL)

write.table(r1_bsave_em, "r1bsave_em.csv",sep=",",dec=".",row.names = FALSE)
write.table(r1_gsave_em, "r1gsave_em.csv",sep=",",dec=".",row.names = FALSE)


# 14. beta = [0,1], n=500, censor 1.25 -------

r2_censorsave=matrix(rep(0,2*100),ncol=2)#whole dataset event proportion, non-cured event proportion
r2_bsave_mpl=matrix(rep(0,4*100),ncol=4)
r2_gsave_mpl=matrix(rep(0,2*100),ncol=2)
r2_bsave_em=matrix(rep(0,4*100),ncol=4)
r2_gsave_em=matrix(rep(0,2*100),ncol=2)

for(r in 1:100){
  sim.data=gen_data_right(500,1.25,c(0,1),c(0.5),c("b"),c("b"),c("z1","x1"))
  
  sim.data.uncure=sim.data[sim.data[,3]==1,]
  sim.surv1=Surv(time=sim.data.uncure[,1],event=sim.data.uncure[,2])
  sim.surv=Surv(time=sim.data$time,event =sim.data$event)
  
  r2_censorsave[r,1]=sum(sim.surv[,2])/500
  r2_censorsave[r,2]=sum(sim.surv1[,2])/nrow(sim.surv1)
  
  test=phmc_mpl(sim.surv~sim.data$X1,pi.formula = ~sim.data$Z,data=sim.data,phmc_mpl.control(conv_limit = 1e-4, n.knots=c(3,0),maxIter = c(1,10000,10001)))
  
  r2_bsave_mpl[r,1]=test$beta[1]
  r2_bsave_mpl[r,2]=test$beta[2]
  r2_bsave_mpl[r,3]=test$se$se_H[1]
  r2_bsave_mpl[r,4]=test$se$se_H[2]

  r2_gsave_mpl[r,1]=test$gamma[1]
  r2_gsave_mpl[r,2]=test$se$se_H[3]
  
  print(c(r,sum(sim.surv1[,2])/nrow(sim.surv1),test$se$se_H[1]))
  
  #smc1 = smcure(sim.surv~X1, cureform=~Z, data=sim.data, model="ph")
  
  #r2_bsave_em[r,1]=smc1$b[1]
  #r2_bsave_em[r,2]=smc1$b[2]
  #r2_bsave_em[r,3]=smc1$b_sd[1]
  #r2_bsave_em[r,4]=smc1$b_sd[2]
  
  #r2_gsave_em[r,1]=smc1$beta[1]
  #r2_gsave_em[r,2]=smc1$beta_sd[1]
  
  
  
}


mean(r2_censorsave[!is.na(r2_bsave_mpl[,3])],1)
mean(r2_censorsave[!is.na(r2_bsave_mpl[,3])],2)
r2_censorsave=data.frame(r2_censorsave)
write.table(r2_censorsave, "r2censor.csv",sep=",",dec=".")

colnames(r2_bsave_mpl)=c("b0","b1","seb0","seb1")
r2_bsave_mpl=data.frame(r2_bsave_mpl)
write.table(r2_bsave_mpl, "r2bsave_mpl.csv",sep=",",dec=".",row.names = FALSE)
colnames(r2_bsave_em)=c("b0","b1","seb0","seb1")
r2_bsave_em=data.frame(r2_bsave_em)
write.table(r2_bsave_em, "r2bsave_em.csv",sep=",",dec=".",row.names = FALSE)
colnames(r2_gsave_mpl)=c("g1","seg1")
r2_gsave_mpl=data.frame(r2_gsave_mpl)
write.table(r2_gsave_mpl, "r2gsave_mpl.csv",sep=",",dec=".",row.names = FALSE)
colnames(r2_gsave_em)=c("g1","seg1")
r2_gsave_em=data.frame(r2_gsave_em)
write.table(r2_gsave_em, "r2gsave_em.csv",sep=",",dec=".",row.names = FALSE)

#mpl
mean(r2_bsave_mpl[!is.na(r2_bsave_mpl[,3]),1])
mean(r2_bsave_mpl[!is.na(r2_bsave_mpl[,3]),2])
mean(r2_bsave_mpl[!is.na(r2_bsave_mpl[,3]),3])
sd(r2_bsave_mpl[!is.na(r2_bsave_mpl[,3]),1])
mean(r2_bsave_mpl[!is.na(r2_bsave_mpl[,3]),4])
sd(r2_bsave_mpl[!is.na(r2_bsave_mpl[,3]),2])

mean(r2_gsave_mpl[!is.na(r2_gsave_mpl[,2]),1])
mean(r2_gsave_mpl[!is.na(r2_gsave_mpl[,2]),2])
sd(r2_gsave_mpl[!is.na(r2_gsave_mpl[,2]),1])

r2_bsave_mpl$b0_LL=r2_bsave_mpl$b0_UL=r2_bsave_mpl$b1_LL=r2_bsave_mpl$b1_UL=rep(0,100)
r2_gsave_mpl$g1_LL=r2_gsave_mpl$g1_UL=rep(0,100)

r2_bsave_mpl$b0_LL[!is.na(r2_bsave_mpl$seb0)]=r2_bsave_mpl$b0[!is.na(r2_bsave_mpl$seb0)]-1.96*r2_bsave_mpl$seb0[!is.na(r2_bsave_mpl$seb0)]
r2_bsave_mpl$b0_UL[!is.na(r2_bsave_mpl$seb0)]=r2_bsave_mpl$b0[!is.na(r2_bsave_mpl$seb0)]+1.96*r2_bsave_mpl$seb0[!is.na(r2_bsave_mpl$seb0)]

r2_bsave_mpl$b1_LL[!is.na(r2_bsave_mpl$seb0)]=r2_bsave_mpl$b1[!is.na(r2_bsave_mpl$seb0)]-1.96*r2_bsave_mpl$seb1[!is.na(r2_bsave_mpl$seb0)]
r2_bsave_mpl$b1_UL[!is.na(r2_bsave_mpl$seb0)]=r2_bsave_mpl$b1[!is.na(r2_bsave_mpl$seb0)]+1.96*r2_bsave_mpl$seb1[!is.na(r2_bsave_mpl$seb0)]

r2_bsave_mpl$b0_cov=r2_bsave_mpl$b1_cov=rep(0,100)
r2_bsave_mpl$b0_cov[!is.na(r2_bsave_mpl$seb0)]=as.numeric(r2_bsave_mpl$b0_LL[!is.na(r2_bsave_mpl$seb0)]<=0 & 0<=r2_bsave_mpl$b0_UL[!is.na(r2_bsave_mpl$seb0)])
r2_bsave_mpl$b1_cov[!is.na(r2_bsave_mpl$seb0)]=as.numeric(r2_bsave_mpl$b1_LL[!is.na(r2_bsave_mpl$seb0)]<=1 & 1<=r2_bsave_mpl$b1_UL[!is.na(r2_bsave_mpl$seb0)])

r2_gsave_mpl$g1_LL[!is.na(r2_gsave_mpl$seg1)]=r2_gsave_mpl$g1[!is.na(r2_gsave_mpl$seg1)]-1.96*r2_gsave_mpl$seg1[!is.na(r2_gsave_mpl$seg1)]
r2_gsave_mpl$g1_UL[!is.na(r2_gsave_mpl$seg1)]=r2_gsave_mpl$g1[!is.na(r2_gsave_mpl$seg1)]+1.96*r2_gsave_mpl$seg1[!is.na(r2_gsave_mpl$seg1)]
r2_gsave_mpl$g1_cov=rep(0,100)
r2_gsave_mpl$g1_cov[!is.na(r2_gsave_mpl$seg1)]=as.numeric(r2_gsave_mpl$g1_LL[!is.na(r2_gsave_mpl$seg1)]<=0.5 & 0.5<=r2_gsave_mpl$g1_UL[!is.na(r2_gsave_mpl$seg1)])





#em

mean(r2_bsave_em[,1])
mean(r2_bsave_em[,2])
mean(r2_bsave_em[,3])
sd(r2_bsave_em[,1])
mean(r2_bsave_em[,4])
sd(r2_bsave_em[,2])

mean(r2_gsave_em[,1])
mean(r2_gsave_em[,2])
sd(r2_gsave_em[,1])

r2_bsave_em$b0_LL=r2_bsave_em$b0_UL=r2_bsave_em$b1_LL=r2_bsave_em$b1_UL=rep(0,100)
r2_gsave_em$g1_LL=r2_gsave_em$g1_UL=rep(0,100)

r2_bsave_em$b0_LL=r2_bsave_em$b0-1.96*r2_bsave_em$seb0
r2_bsave_em$b0_UL=r2_bsave_em$b0+1.96*r2_bsave_em$seb0

r2_bsave_em$b1_LL=r2_bsave_em$b1-1.96*r2_bsave_em$seb1
r2_bsave_em$b1_UL=r2_bsave_em$b1+1.96*r2_bsave_em$seb1

r2_bsave_em$b0_cov=r2_bsave_em$b1_cov=rep(0,100)
r2_bsave_em$b0_cov=as.numeric(r2_bsave_em$b0_LL<=0 & 0<=r2_bsave_em$b0_UL)
r2_bsave_em$b1_cov=as.numeric(r2_bsave_em$b1_LL<=1 & 1<=r2_bsave_em$b1_UL)

r2_gsave_em$g1_LL=r2_gsave_em$g1-1.96*r2_gsave_em$seg1
r2_gsave_em$g1_UL=r2_gsave_em$g1+1.96*r2_gsave_em$seg1
r2_gsave_em$g1_cov=rep(0,100)
r2_gsave_em$g1_cov=as.numeric(r2_gsave_em$g1_LL<=0.5 & 0.5<=r2_gsave_em$g1_UL)








# 15. beta = [-1,1], n=500, censor 1.25 -------

r3_censorsave=matrix(rep(0,2*100),ncol=2)#whole dataset event proportion, non-cured event proportion
r3_bsave_mpl=matrix(rep(0,4*100),ncol=4)
r3_gsave_mpl=matrix(rep(0,2*100),ncol=2)
r3_bsave_em=matrix(rep(0,4*100),ncol=4)
r3_gsave_em=matrix(rep(0,2*100),ncol=2)

for(r in 1:100){
  sim.data=gen_data_right(500,1.25,c(-1,1),c(0.5),c("b"),c("b"),c("z1","x1"))
  
  #sim.data.uncure=sim.data[sim.data[,3]==1,]
  #sim.surv1=Surv(time=sim.data.uncure[,1],event=sim.data.uncure[,2])
  sim.surv=Surv(time=sim.data$time,event =sim.data$event)
  
  #r3_censorsave[r,1]=sum(sim.surv[,2])/500
  #r3_censorsave[r,2]=sum(sim.surv1[,2])/nrow(sim.surv1)
  
  #test=phmc_mpl(sim.surv~sim.data$X1,pi.formula = ~sim.data$Z,data=sim.data,phmc_mpl.control(conv_limit = 1e-4, n.knots=c(1,0),maxIter = c(1,10000,10001)))
  
  #r3_bsave_mpl[r,1]=test$beta[1]
  #r3_bsave_mpl[r,2]=test$beta[2]
  #r3_bsave_mpl[r,3]=test$se$se_H[1]
  #r3_bsave_mpl[r,4]=test$se$se_H[2]
  
  #r3_gsave_mpl[r,1]=test$gamma[1]
  #r3_gsave_mpl[r,2]=test$se$se_H[3]

  #print(c(r,sum(sim.surv1[,2])/nrow(sim.surv1),test$se$se_H[1],test$knots$m))
  
  smc1 = smcure(sim.surv~X1, cureform=~Z, data=sim.data, model="ph")
  
  r3_bsave_em[r,1]=smc1$b[1]
  r3_bsave_em[r,2]=smc1$b[2]
  r3_bsave_em[r,3]=smc1$b_sd[1]
  r3_bsave_em[r,4]=smc1$b_sd[2]
  
  r3_gsave_em[r,1]=smc1$beta[1]
  r3_gsave_em[r,2]=smc1$beta_sd[1]
  
  
  
}



mean(r3_censorsave[!is.na(r3_bsave_mpl[,3])],1)
mean(r3_censorsave[!is.na(r3_bsave_mpl[,3])],2)
r3_censorsave=data.frame(r3_censorsave)
write.table(r3_censorsave, "r3censor.csv",sep=",",dec=".")



colnames(r3_bsave_mpl)=c("b0","b1","seb0","seb1")
r3_bsave_mpl=data.frame(r3_bsave_mpl)
write.table(r3_bsave_mpl, "r3bsave_mpl.csv",sep=",",dec=".",row.names = F)
colnames(r3_bsave_em)=c("b0","b1","seb0","seb1")
r3_bsave_em=data.frame(r3_bsave_em)
write.table(r3_bsave_em, "r3bsave_em.csv",sep=",",dec=".",row.names = F)
colnames(r3_gsave_mpl)=c("g1","seg1")
r3_gsave_mpl=data.frame(r3_gsave_mpl)
write.table(r3_gsave_mpl, "r3gsave_mpl.csv",sep=",",dec=".",row.names = F)
colnames(r3_gsave_em)=c("g1","seg1")
r3_gsave_em=data.frame(r3_gsave_em)
write.table(r3_gsave_em, "r3gsave_em.csv",sep=",",dec=".",row.names = F)

#mpl
mean(r3_bsave_mpl[!is.na(r3_bsave_mpl[,3]),1])
mean(r3_bsave_mpl[!is.na(r3_bsave_mpl[,3]),2])
mean(r3_bsave_mpl[!is.na(r3_bsave_mpl[,3]),3])
sd(r3_bsave_mpl[!is.na(r3_bsave_mpl[,3]),1])
mean(r3_bsave_mpl[!is.na(r3_bsave_mpl[,3]),4])
sd(r3_bsave_mpl[!is.na(r3_bsave_mpl[,3]),2])

mean(r3_gsave_mpl[!is.na(r3_gsave_mpl[,2]),1])
mean(r3_gsave_mpl[!is.na(r3_gsave_mpl[,2]),2])
sd(r3_gsave_mpl[!is.na(r3_gsave_mpl[,2]),1])

r3_bsave_mpl$b0_LL=r3_bsave_mpl$b0_UL=r3_bsave_mpl$b1_LL=r3_bsave_mpl$b1_UL=rep(0,100)
r3_gsave_mpl$g1_LL=r3_gsave_mpl$g1_UL=rep(0,100)

r3_bsave_mpl$b0_LL[!is.na(r3_bsave_mpl$seb0)]=r3_bsave_mpl$b0[!is.na(r3_bsave_mpl$seb0)]-1.96*r3_bsave_mpl$seb0[!is.na(r3_bsave_mpl$seb0)]
r3_bsave_mpl$b0_UL[!is.na(r3_bsave_mpl$seb0)]=r3_bsave_mpl$b0[!is.na(r3_bsave_mpl$seb0)]+1.96*r3_bsave_mpl$seb0[!is.na(r3_bsave_mpl$seb0)]

r3_bsave_mpl$b1_LL[!is.na(r3_bsave_mpl$seb0)]=r3_bsave_mpl$b1[!is.na(r3_bsave_mpl$seb0)]-1.96*r3_bsave_mpl$seb1[!is.na(r3_bsave_mpl$seb0)]
r3_bsave_mpl$b1_UL[!is.na(r3_bsave_mpl$seb0)]=r3_bsave_mpl$b1[!is.na(r3_bsave_mpl$seb0)]+1.96*r3_bsave_mpl$seb1[!is.na(r3_bsave_mpl$seb0)]

r3_bsave_mpl$b0_cov=r3_bsave_mpl$b1_cov=rep(0,100)
r3_bsave_mpl$b0_cov[!is.na(r3_bsave_mpl$seb0)]=as.numeric(r3_bsave_mpl$b0_LL[!is.na(r3_bsave_mpl$seb0)]<=-1 & -1<=r3_bsave_mpl$b0_UL[!is.na(r3_bsave_mpl$seb0)])
r3_bsave_mpl$b1_cov[!is.na(r3_bsave_mpl$seb0)]=as.numeric(r3_bsave_mpl$b1_LL[!is.na(r3_bsave_mpl$seb0)]<=1 & 1<=r3_bsave_mpl$b1_UL[!is.na(r3_bsave_mpl$seb0)])

r3_gsave_mpl$g1_LL[!is.na(r3_gsave_mpl$seg1)]=r3_gsave_mpl$g1[!is.na(r3_gsave_mpl$seg1)]-1.96*r3_gsave_mpl$seg1[!is.na(r3_gsave_mpl$seg1)]
r3_gsave_mpl$g1_UL[!is.na(r3_gsave_mpl$seg1)]=r3_gsave_mpl$g1[!is.na(r3_gsave_mpl$seg1)]+1.96*r3_gsave_mpl$seg1[!is.na(r3_gsave_mpl$seg1)]
r3_gsave_mpl$g1_cov=rep(0,100)
r3_gsave_mpl$g1_cov[!is.na(r3_gsave_mpl$seg1)]=as.numeric(r3_gsave_mpl$g1_LL[!is.na(r3_gsave_mpl$seg1)]<=0.5 & 0.5<=r3_gsave_mpl$g1_UL[!is.na(r3_gsave_mpl$seg1)])

#em

mean(r3_bsave_em[,1])
mean(r3_bsave_em[,2])
mean(r3_bsave_em[,3])
sd(r3_bsave_em[,1])
mean(r3_bsave_em[,4])
sd(r3_bsave_em[,2])

mean(r3_gsave_em[,1])
mean(r3_gsave_em[,2])
sd(r3_gsave_em[,1])

r3_bsave_em$b0_LL=r3_bsave_em$b0_UL=r3_bsave_em$b1_LL=r3_bsave_em$b1_UL=rep(0,100)
r3_gsave_em$g1_LL=r3_gsave_em$g1_UL=rep(0,100)

r3_bsave_em$b0_LL=r3_bsave_em$b0-1.96*r3_bsave_em$seb0
r3_bsave_em$b0_UL=r3_bsave_em$b0+1.96*r3_bsave_em$seb0

r3_bsave_em$b1_LL=r3_bsave_em$b1-1.96*r3_bsave_em$seb1
r3_bsave_em$b1_UL=r3_bsave_em$b1+1.96*r3_bsave_em$seb1

r3_bsave_em$b0_cov=r3_bsave_em$b1_cov=rep(0,100)
r3_bsave_em$b0_cov=as.numeric(r3_bsave_em$b0_LL<=-1 & -1<=r3_bsave_em$b0_UL)
r3_bsave_em$b1_cov=as.numeric(r3_bsave_em$b1_LL<=1 & 1<=r3_bsave_em$b1_UL)

r3_gsave_em$g1_LL=r3_gsave_em$g1-1.96*r3_gsave_em$seg1
r3_gsave_em$g1_UL=r3_gsave_em$g1+1.96*r3_gsave_em$seg1
r3_gsave_em$g1_cov=rep(0,100)
r3_gsave_em$g1_cov=as.numeric(r3_gsave_em$g1_LL<=0.5 & 0.5<=r3_gsave_em$g1_UL)

# 16. beta = [1,1], n=100, censor 1.25 --------

r16_censorsave=matrix(rep(0,2*100),ncol=2)#whole dataset event proportion, non-cured event proportion
r16_bsave_mpl=matrix(rep(0,4*100),ncol=4)
r16_gsave_mpl=matrix(rep(0,2*100),ncol=2)
r16_bsave_em=matrix(rep(0,4*100),ncol=4)
r16_gsave_em=matrix(rep(0,2*100),ncol=2)



for(r in 1:100){
  sim.data=gen_data_right(100,1.25,c(1,1),c(0.5),c("b"),c("b"),c("z1","x1"))
  
  sim.data.uncure=sim.data[sim.data[,3]==1,]
  sim.surv1=Surv(time=sim.data.uncure[,1],event=sim.data.uncure[,2])
  sim.surv=Surv(time=sim.data$time,event =sim.data$event)
  r16_censorsave[r,1]=sum(sim.surv[,2])/100
  r16_censorsave[r,2]=sum(sim.surv1[,2])/nrow(sim.surv1)
  
  try1=try(phmc_mpl(sim.surv~sim.data$X1,pi.formula = ~sim.data$Z,data=sim.data,phmc_mpl.control(conv_limit = 1e-4, n.knots=c(1,0),maxIter = c(1,30000,30001))))
  try2=try(smcure(sim.surv~X1, cureform=~Z, data=sim.data, model="ph"))
  if(class(try1)!="try-error" & class(try2)!="try-error"){
    test=try1
    smc1=try2
    
    r16_bsave_mpl[r,1]=test$beta[1]
    r16_bsave_mpl[r,2]=test$beta[2]
    r16_bsave_mpl[r,3]=test$se$se_H[1]
    r16_bsave_mpl[r,4]=test$se$se_H[2]
    
    r16_gsave_mpl[r,1]=test$gamma[1]
    r16_gsave_mpl[r,2]=test$se$se_H[3]
    
    r16_bsave_em[r,1]=smc1$b[1]
    r16_bsave_em[r,2]=smc1$b[2]
    r16_bsave_em[r,3]=smc1$b_sd[1]
    r16_bsave_em[r,4]=smc1$b_sd[2]
    
    r16_gsave_em[r,1]=smc1$beta[1]
    r16_gsave_em[r,2]=smc1$beta_sd[1]
    
    print(c(r,"both 16",test$se$se_H[1]))
  } else if(class(try1)!="try-error"){
    test=try1
    r16_bsave_mpl[r,1]=test$beta[1]
    r16_bsave_mpl[r,2]=test$beta[2]
    r16_bsave_mpl[r,3]=test$se$se_H[1]
    r16_bsave_mpl[r,4]=test$se$se_H[2]
    
    r16_gsave_mpl[r,1]=test$gamma[1]
    r16_gsave_mpl[r,2]=test$se$se_H[3]
    print(c(r, "mpl 16",test$se$se_H[1]))
    
  } else if(class(try2)!="try-error"){
    smc1=try2
    r16_bsave_em[r,1]=smc1$b[1]
    r16_bsave_em[r,2]=smc1$b[2]
    r16_bsave_em[r,3]=smc1$b_sd[1]
    r16_bsave_em[r,4]=smc1$b_sd[2]
    
    r16_gsave_em[r,1]=smc1$beta[1]
    r16_gsave_em[r,2]=smc1$beta_sd[1]
    print(c(r, "em 16"))
  } else{
    print(c(r, "neither 16"))
  }
}




mean(r16_censorsave[!is.na(r16_bsave_mpl[,3])],1)
mean(r16_censorsave[!is.na(r16_bsave_mpl[,3])],2)
r16_censorsave=data.frame(r16_censorsave)
write.table(r16_censorsave, "r16censor.csv",sep=",",dec=".")

colnames(r16_bsave_mpl)=c("b0","b1","seb0","seb1")
r16_bsave_mpl=data.frame(r16_bsave_mpl)
write.table(r16_bsave_mpl, "r16bsave_mpl.csv",sep=",",dec=".",row.names = F)
colnames(r16_bsave_em)=c("b0","b1","seb0","seb1")
r16_bsave_em=data.frame(r16_bsave_em)
write.table(r16_bsave_em, "r16bsave_em.csv",sep=",",dec=".",row.names = F)
colnames(r16_gsave_mpl)=c("g1","seg1")
r16_gsave_mpl=data.frame(r16_gsave_mpl)
write.table(r16_gsave_mpl, "r16gsave_mpl.csv",sep=",",dec=".",row.names = F)
colnames(r16_gsave_em)=c("g1","seg1")
r16_gsave_em=data.frame(r16_gsave_em)
write.table(r16_gsave_em, "r16gsave_em.csv",sep=",",dec=".",row.names = F)


#mpl
mean(r16_bsave_mpl[!is.na(r16_bsave_mpl[,3]),1])
mean(r16_bsave_mpl[!is.na(r16_bsave_mpl[,3]),2])
mean(r16_bsave_mpl[!is.na(r16_bsave_mpl[,3]),3])
sd(r16_bsave_mpl[!is.na(r16_bsave_mpl[,3]),1])
mean(r16_bsave_mpl[!is.na(r16_bsave_mpl[,3]),4])
sd(r16_bsave_mpl[!is.na(r16_bsave_mpl[,3]),2])

mean(r16_gsave_mpl[!is.na(r16_gsave_mpl[,2]),1])
mean(r16_gsave_mpl[!is.na(r16_gsave_mpl[,2]),2])
sd(r16_gsave_mpl[!is.na(r16_gsave_mpl[,2]),1])

r16_bsave_mpl$b0_LL=r16_bsave_mpl$b0_UL=r16_bsave_mpl$b1_LL=r16_bsave_mpl$b1_UL=rep(0,100)
r16_gsave_mpl$g1_LL=r16_gsave_mpl$g1_UL=rep(0,100)

r16_bsave_mpl$b0_LL[!is.na(r16_bsave_mpl$seb0)]=r16_bsave_mpl$b0[!is.na(r16_bsave_mpl$seb0)]-1.96*r16_bsave_mpl$seb0[!is.na(r16_bsave_mpl$seb0)]
r16_bsave_mpl$b0_UL[!is.na(r16_bsave_mpl$seb0)]=r16_bsave_mpl$b0[!is.na(r16_bsave_mpl$seb0)]+1.96*r16_bsave_mpl$seb0[!is.na(r16_bsave_mpl$seb0)]

r16_bsave_mpl$b1_LL[!is.na(r16_bsave_mpl$seb0)]=r16_bsave_mpl$b1[!is.na(r16_bsave_mpl$seb0)]-1.96*r16_bsave_mpl$seb1[!is.na(r16_bsave_mpl$seb0)]
r16_bsave_mpl$b1_UL[!is.na(r16_bsave_mpl$seb0)]=r16_bsave_mpl$b1[!is.na(r16_bsave_mpl$seb0)]+1.96*r16_bsave_mpl$seb1[!is.na(r16_bsave_mpl$seb0)]

r16_bsave_mpl$b0_cov=r16_bsave_mpl$b1_cov=rep(0,100)
r16_bsave_mpl$b0_cov[!is.na(r16_bsave_mpl$seb0)]=as.numeric(r16_bsave_mpl$b0_LL[!is.na(r16_bsave_mpl$seb0)]<=1 & 1<=r16_bsave_mpl$b0_UL[!is.na(r16_bsave_mpl$seb0)])
r16_bsave_mpl$b1_cov[!is.na(r16_bsave_mpl$seb0)]=as.numeric(r16_bsave_mpl$b1_LL[!is.na(r16_bsave_mpl$seb0)]<=1 & 1<=r16_bsave_mpl$b1_UL[!is.na(r16_bsave_mpl$seb0)])

r16_gsave_mpl$g1_LL[!is.na(r16_gsave_mpl$seg1)]=r16_gsave_mpl$g1[!is.na(r16_gsave_mpl$seg1)]-1.96*r16_gsave_mpl$seg1[!is.na(r16_gsave_mpl$seg1)]
r16_gsave_mpl$g1_UL[!is.na(r16_gsave_mpl$seg1)]=r16_gsave_mpl$g1[!is.na(r16_gsave_mpl$seg1)]+1.96*r16_gsave_mpl$seg1[!is.na(r16_gsave_mpl$seg1)]
r16_gsave_mpl$g1_cov=rep(0,100)
r16_gsave_mpl$g1_cov[!is.na(r16_gsave_mpl$seg1)]=as.numeric(r16_gsave_mpl$g1_LL[!is.na(r16_gsave_mpl$seg1)]<=0.5 & 0.5<=r16_gsave_mpl$g1_UL[!is.na(r16_gsave_mpl$seg1)])


#em

mean(r16_bsave_em[,1])
mean(r16_bsave_em[,2])
mean(r16_bsave_em[,3])
sd(r16_bsave_em[,1])
mean(r16_bsave_em[,4])
sd(r16_bsave_em[,2])

mean(r16_gsave_em[,1])
mean(r16_gsave_em[,2])
sd(r16_gsave_em[,1])

r16_bsave_em$b0_LL=r16_bsave_em$b0_UL=r16_bsave_em$b1_LL=r16_bsave_em$b1_UL=rep(0,100)
r16_gsave_em$g1_LL=r16_gsave_em$g1_UL=rep(0,100)

r16_bsave_em$b0_LL=r16_bsave_em$b0-1.96*r16_bsave_em$seb0
r16_bsave_em$b0_UL=r16_bsave_em$b0+1.96*r16_bsave_em$seb0

r16_bsave_em$b1_LL=r16_bsave_em$b1-1.96*r16_bsave_em$seb1
r16_bsave_em$b1_UL=r16_bsave_em$b1+1.96*r16_bsave_em$seb1

r16_bsave_em$b0_cov=r16_bsave_em$b1_cov=rep(0,100)
r16_bsave_em$b0_cov=as.numeric(r16_bsave_em$b0_LL<=1 & 1<=r16_bsave_em$b0_UL)
r16_bsave_em$b1_cov=as.numeric(r16_bsave_em$b1_LL<=1 & 1<=r16_bsave_em$b1_UL)

r16_gsave_em$g1_LL=r16_gsave_em$g1-1.96*r16_gsave_em$seg1
r16_gsave_em$g1_UL=r16_gsave_em$g1+1.96*r16_gsave_em$seg1
r16_gsave_em$g1_cov=rep(0,100)
r16_gsave_em$g1_cov=as.numeric(r16_gsave_em$g1_LL<=0.5 & 0.5<=r16_gsave_em$g1_UL)


# 17. beta=[0,1], n=100, censor 1.25 ------


r17_censorsave=matrix(rep(0,2*100),ncol=2)#whole dataset event proportion, non-cured event proportion
r17_bsave_mpl=matrix(rep(0,4*100),ncol=4)
r17_gsave_mpl=matrix(rep(0,2*100),ncol=2)
r17_bsave_em=matrix(rep(0,4*100),ncol=4)
r17_gsave_em=matrix(rep(0,2*100),ncol=2)


for(r in 1:100){
  sim.data=gen_data_right(100,1.25,c(0,1),c(0.5),c("b"),c("b"),c("z1","x1"))
  
  sim.data.uncure=sim.data[sim.data[,3]==1,]
  sim.surv1=Surv(time=sim.data.uncure[,1],event=sim.data.uncure[,2])
  sim.surv=Surv(time=sim.data$time,event =sim.data$event)
  r17_censorsave[r,1]=sum(sim.surv[,2])/100
  r17_censorsave[r,2]=sum(sim.surv1[,2])/nrow(sim.surv1)
  
  try1=try(phmc_mpl(sim.surv~sim.data$X1,pi.formula = ~sim.data$Z,data=sim.data,phmc_mpl.control(conv_limit = 1e-4, n.knots=c(1,0),maxIter = c(1,30000,30001))))
  try2=try(smcure(sim.surv~X1, cureform=~Z, data=sim.data, model="ph"))
  if(class(try1)!="try-error" & class(try2)!="try-error"){
    test=try1
    smc1=try2
    
    r17_bsave_mpl[r,1]=test$beta[1]
    r17_bsave_mpl[r,2]=test$beta[2]
    r17_bsave_mpl[r,3]=test$se$se_H[1]
    r17_bsave_mpl[r,4]=test$se$se_H[2]
    
    r17_gsave_mpl[r,1]=test$gamma[1]
    r17_gsave_mpl[r,2]=test$se$se_H[3]
    
    r17_bsave_em[r,1]=smc1$b[1]
    r17_bsave_em[r,2]=smc1$b[2]
    r17_bsave_em[r,3]=smc1$b_sd[1]
    r17_bsave_em[r,4]=smc1$b_sd[2]
    
    r17_gsave_em[r,1]=smc1$beta[1]
    r17_gsave_em[r,2]=smc1$beta_sd[1]
    
    print(c(r,"both 17",test$se$se_H[1]))
  } else if(class(try1)!="try-error"){
    test=try1
    r17_bsave_mpl[r,1]=test$beta[1]
    r17_bsave_mpl[r,2]=test$beta[2]
    r17_bsave_mpl[r,3]=test$se$se_H[1]
    r17_bsave_mpl[r,4]=test$se$se_H[2]
    
    r17_gsave_mpl[r,1]=test$gamma[1]
    r17_gsave_mpl[r,2]=test$se$se_H[3]
    print(c(r, "mpl 17",test$se$se_H[1]))
    
  } else if(class(try2)!="try-error"){
    smc1=try2
    r17_bsave_em[r,1]=smc1$b[1]
    r17_bsave_em[r,2]=smc1$b[2]
    r17_bsave_em[r,3]=smc1$b_sd[1]
    r17_bsave_em[r,4]=smc1$b_sd[2]
    
    r17_gsave_em[r,1]=smc1$beta[1]
    r17_gsave_em[r,2]=smc1$beta_sd[1]
    print(c(r, "em 17"))
  } else{
    print(c(r, "neither 17"))
  }
}




mean(r17_censorsave[!is.na(r17_bsave_mpl[,3])],1)
mean(r17_censorsave[!is.na(r17_bsave_mpl[,3])],2)
r17_censorsave=data.frame(r17_censorsave)
write.table(r17_censorsave, "r17censor.csv",sep=",",dec=".")

colnames(r17_bsave_mpl)=c("b0","b1","seb0","seb1")
r17_bsave_mpl=data.frame(r17_bsave_mpl)
write.table(r17_bsave_mpl, "r17bsave_mpl.csv",sep=",",dec=".",row.names = F)
colnames(r17_bsave_em)=c("b0","b1","seb0","seb1")
r17_bsave_em=data.frame(r17_bsave_em)
write.table(r17_bsave_em, "r17bsave_em.csv",sep=",",dec=".",row.names = F)
colnames(r17_gsave_mpl)=c("g1","seg1")
r17_gsave_mpl=data.frame(r17_gsave_mpl)
write.table(r17_gsave_mpl, "r17gsave_mpl.csv",sep=",",dec=".",row.names = F)
colnames(r17_gsave_em)=c("g1","seg1")
r17_gsave_em=data.frame(r17_gsave_em)
write.table(r17_gsave_em, "r17gsave_em.csv",sep=",",dec=".",row.names = F)


#mpl
mean(r17_bsave_mpl[!is.na(r17_bsave_mpl[,3]),1])
mean(r17_bsave_mpl[!is.na(r17_bsave_mpl[,3]),2])
mean(r17_bsave_mpl[!is.na(r17_bsave_mpl[,3]),3])
sd(r17_bsave_mpl[!is.na(r17_bsave_mpl[,3]),1])
mean(r17_bsave_mpl[!is.na(r17_bsave_mpl[,3]),4])
sd(r17_bsave_mpl[!is.na(r17_bsave_mpl[,3]),2])

mean(r17_gsave_mpl[!is.na(r17_gsave_mpl[,2]),1])
mean(r17_gsave_mpl[!is.na(r17_gsave_mpl[,2]),2])
sd(r17_gsave_mpl[!is.na(r17_gsave_mpl[,2]),1])

r17_bsave_mpl$b0_LL=r17_bsave_mpl$b0_UL=r17_bsave_mpl$b1_LL=r17_bsave_mpl$b1_UL=rep(0,100)
r17_gsave_mpl$g1_LL=r17_gsave_mpl$g1_UL=rep(0,100)

r17_bsave_mpl$b0_LL[!is.na(r17_bsave_mpl$seb0)]=r17_bsave_mpl$b0[!is.na(r17_bsave_mpl$seb0)]-1.96*r17_bsave_mpl$seb0[!is.na(r17_bsave_mpl$seb0)]
r17_bsave_mpl$b0_UL[!is.na(r17_bsave_mpl$seb0)]=r17_bsave_mpl$b0[!is.na(r17_bsave_mpl$seb0)]+1.96*r17_bsave_mpl$seb0[!is.na(r17_bsave_mpl$seb0)]

r17_bsave_mpl$b1_LL[!is.na(r17_bsave_mpl$seb0)]=r17_bsave_mpl$b1[!is.na(r17_bsave_mpl$seb0)]-1.96*r17_bsave_mpl$seb1[!is.na(r17_bsave_mpl$seb0)]
r17_bsave_mpl$b1_UL[!is.na(r17_bsave_mpl$seb0)]=r17_bsave_mpl$b1[!is.na(r17_bsave_mpl$seb0)]+1.96*r17_bsave_mpl$seb1[!is.na(r17_bsave_mpl$seb0)]

r17_bsave_mpl$b0_cov=r17_bsave_mpl$b1_cov=rep(0,100)
r17_bsave_mpl$b0_cov[!is.na(r17_bsave_mpl$seb0)]=as.numeric(r17_bsave_mpl$b0_LL[!is.na(r17_bsave_mpl$seb0)]<=0 & 0<=r17_bsave_mpl$b0_UL[!is.na(r17_bsave_mpl$seb0)])
r17_bsave_mpl$b1_cov[!is.na(r17_bsave_mpl$seb0)]=as.numeric(r17_bsave_mpl$b1_LL[!is.na(r17_bsave_mpl$seb0)]<=1 & 1<=r17_bsave_mpl$b1_UL[!is.na(r17_bsave_mpl$seb0)])

r17_gsave_mpl$g1_LL[!is.na(r17_gsave_mpl$seg1)]=r17_gsave_mpl$g1[!is.na(r17_gsave_mpl$seg1)]-1.96*r17_gsave_mpl$seg1[!is.na(r17_gsave_mpl$seg1)]
r17_gsave_mpl$g1_UL[!is.na(r17_gsave_mpl$seg1)]=r17_gsave_mpl$g1[!is.na(r17_gsave_mpl$seg1)]+1.96*r17_gsave_mpl$seg1[!is.na(r17_gsave_mpl$seg1)]
r17_gsave_mpl$g1_cov=rep(0,100)
r17_gsave_mpl$g1_cov[!is.na(r17_gsave_mpl$seg1)]=as.numeric(r17_gsave_mpl$g1_LL[!is.na(r17_gsave_mpl$seg1)]<=0.5 & 0.5<=r17_gsave_mpl$g1_UL[!is.na(r17_gsave_mpl$seg1)])


#em

mean(r17_bsave_em[,1])
mean(r17_bsave_em[,2])
mean(r17_bsave_em[,3])
sd(r17_bsave_em[,1])
mean(r17_bsave_em[,4])
sd(r17_bsave_em[,2])

mean(r17_gsave_em[,1])
mean(r17_gsave_em[,2])
sd(r17_gsave_em[,1])

r17_bsave_em$b0_LL=r17_bsave_em$b0_UL=r17_bsave_em$b1_LL=r17_bsave_em$b1_UL=rep(0,100)
r17_gsave_em$g1_LL=r17_gsave_em$g1_UL=rep(0,100)

r17_bsave_em$b0_LL=r17_bsave_em$b0-1.96*r17_bsave_em$seb0
r17_bsave_em$b0_UL=r17_bsave_em$b0+1.96*r17_bsave_em$seb0

r17_bsave_em$b1_LL=r17_bsave_em$b1-1.96*r17_bsave_em$seb1
r17_bsave_em$b1_UL=r17_bsave_em$b1+1.96*r17_bsave_em$seb1

r17_bsave_em$b0_cov=r17_bsave_em$b1_cov=rep(0,100)
r17_bsave_em$b0_cov=as.numeric(r17_bsave_em$b0_LL<=0 & 0<=r17_bsave_em$b0_UL)
r17_bsave_em$b1_cov=as.numeric(r17_bsave_em$b1_LL<=1 & 1<=r17_bsave_em$b1_UL)

r17_gsave_em$g1_LL=r17_gsave_em$g1-1.96*r17_gsave_em$seg1
r17_gsave_em$g1_UL=r17_gsave_em$g1+1.96*r17_gsave_em$seg1
r17_gsave_em$g1_cov=rep(0,100)
r17_gsave_em$g1_cov=as.numeric(r17_gsave_em$g1_LL<=0.5 & 0.5<=r17_gsave_em$g1_UL)


# 18. beta = [-1,1], n=100, censor 1.25 -----




r18_censorsave=matrix(rep(0,2*100),ncol=2)#whole dataset event proportion, non-cured event proportion
r18_bsave_mpl=matrix(rep(0,4*100),ncol=4)
r18_gsave_mpl=matrix(rep(0,2*100),ncol=2)
r18_bsave_em=matrix(rep(0,4*100),ncol=4)
r18_gsave_em=matrix(rep(0,2*100),ncol=2)


for(r in 1:100){
  sim.data=gen_data_right(100,1.25,c(-1,1),c(0.5),c("b"),c("b"),c("z1","x1"))
  
  sim.data.uncure=sim.data[sim.data[,3]==1,]
  sim.surv1=Surv(time=sim.data.uncure[,1],event=sim.data.uncure[,2])
  sim.surv=Surv(time=sim.data$time,event =sim.data$event)
  r18_censorsave[r,1]=sum(sim.surv[,2])/100
  r18_censorsave[r,2]=sum(sim.surv1[,2])/nrow(sim.surv1)
  
  try1=try(phmc_mpl(sim.surv~sim.data$X1,pi.formula = ~sim.data$Z,data=sim.data,phmc_mpl.control(conv_limit = 1e-4, n.knots=c(1,0),maxIter = c(1,30000,30001))))
  try2=try(smcure(sim.surv~X1, cureform=~Z, data=sim.data, model="ph"))
  if(class(try1)!="try-error" & class(try2)!="try-error"){
    test=try1
    smc1=try2
    
    r18_bsave_mpl[r,1]=test$beta[1]
    r18_bsave_mpl[r,2]=test$beta[2]
    r18_bsave_mpl[r,3]=test$se$se_H[1]
    r18_bsave_mpl[r,4]=test$se$se_H[2]
    
    r18_gsave_mpl[r,1]=test$gamma[1]
    r18_gsave_mpl[r,2]=test$se$se_H[3]
    
    r18_bsave_em[r,1]=smc1$b[1]
    r18_bsave_em[r,2]=smc1$b[2]
    r18_bsave_em[r,3]=smc1$b_sd[1]
    r18_bsave_em[r,4]=smc1$b_sd[2]
    
    r18_gsave_em[r,1]=smc1$beta[1]
    r18_gsave_em[r,2]=smc1$beta_sd[1]
    
    print(c(r,"both 18",test$se$se_H[1]))
  } else if(class(try1)!="try-error"){
    test=try1
    r18_bsave_mpl[r,1]=test$beta[1]
    r18_bsave_mpl[r,2]=test$beta[2]
    r18_bsave_mpl[r,3]=test$se$se_H[1]
    r18_bsave_mpl[r,4]=test$se$se_H[2]
    
    r18_gsave_mpl[r,1]=test$gamma[1]
    r18_gsave_mpl[r,2]=test$se$se_H[3]
    print(c(r, "mpl 18",test$se$se_H[1]))
    
  } else if(class(try2)!="try-error"){
    smc1=try2
    r18_bsave_em[r,1]=smc1$b[1]
    r18_bsave_em[r,2]=smc1$b[2]
    r18_bsave_em[r,3]=smc1$b_sd[1]
    r18_bsave_em[r,4]=smc1$b_sd[2]
    
    r18_gsave_em[r,1]=smc1$beta[1]
    r18_gsave_em[r,2]=smc1$beta_sd[1]
    print(c(r, "em 18"))
  } else{
    print(c(r, "neither 18"))
  }
}




mean(r18_censorsave[!is.na(r18_bsave_mpl[,3])],1)
mean(r18_censorsave[!is.na(r18_bsave_mpl[,3])],2)
r18_censorsave=data.frame(r18_censorsave)
write.table(r18_censorsave, "r18censor.csv",sep=",",dec=".")

colnames(r18_bsave_mpl)=c("b0","b1","seb0","seb1")
r18_bsave_mpl=data.frame(r18_bsave_mpl)
write.table(r18_bsave_mpl, "r18bsave_mpl.csv",sep=",",dec=".",row.names = F)
colnames(r18_bsave_em)=c("b0","b1","seb0","seb1")
r18_bsave_em=data.frame(r18_bsave_em)
write.table(r18_bsave_em, "r18bsave_em.csv",sep=",",dec=".",row.names = F)
colnames(r18_gsave_mpl)=c("g1","seg1")
r18_gsave_mpl=data.frame(r18_gsave_mpl)
write.table(r18_gsave_mpl, "r18gsave_mpl.csv",sep=",",dec=".",row.names = F)
colnames(r18_gsave_em)=c("g1","seg1")
r18_gsave_em=data.frame(r18_gsave_em)
write.table(r18_gsave_em, "r18gsave_em.csv",sep=",",dec=".",row.names = F)


#mpl
mean(r18_bsave_mpl[!is.na(r18_bsave_mpl[,3]),1])
mean(r18_bsave_mpl[!is.na(r18_bsave_mpl[,3]),2])
mean(r18_bsave_mpl[!is.na(r18_bsave_mpl[,3]),3])
sd(r18_bsave_mpl[!is.na(r18_bsave_mpl[,3]),1])
mean(r18_bsave_mpl[!is.na(r18_bsave_mpl[,3]),4])
sd(r18_bsave_mpl[!is.na(r18_bsave_mpl[,3]),2])

mean(r18_gsave_mpl[!is.na(r18_gsave_mpl[,2]),1])
mean(r18_gsave_mpl[!is.na(r18_gsave_mpl[,2]),2])
sd(r18_gsave_mpl[!is.na(r18_gsave_mpl[,2]),1])

r18_bsave_mpl$b0_LL=r18_bsave_mpl$b0_UL=r18_bsave_mpl$b1_LL=r18_bsave_mpl$b1_UL=rep(0,100)
r18_gsave_mpl$g1_LL=r18_gsave_mpl$g1_UL=rep(0,100)

r18_bsave_mpl$b0_LL[!is.na(r18_bsave_mpl$seb0)]=r18_bsave_mpl$b0[!is.na(r18_bsave_mpl$seb0)]-1.96*r18_bsave_mpl$seb0[!is.na(r18_bsave_mpl$seb0)]
r18_bsave_mpl$b0_UL[!is.na(r18_bsave_mpl$seb0)]=r18_bsave_mpl$b0[!is.na(r18_bsave_mpl$seb0)]+1.96*r18_bsave_mpl$seb0[!is.na(r18_bsave_mpl$seb0)]

r18_bsave_mpl$b1_LL[!is.na(r18_bsave_mpl$seb0)]=r18_bsave_mpl$b1[!is.na(r18_bsave_mpl$seb0)]-1.96*r18_bsave_mpl$seb1[!is.na(r18_bsave_mpl$seb0)]
r18_bsave_mpl$b1_UL[!is.na(r18_bsave_mpl$seb0)]=r18_bsave_mpl$b1[!is.na(r18_bsave_mpl$seb0)]+1.96*r18_bsave_mpl$seb1[!is.na(r18_bsave_mpl$seb0)]

r18_bsave_mpl$b0_cov=r18_bsave_mpl$b1_cov=rep(0,100)
r18_bsave_mpl$b0_cov[!is.na(r18_bsave_mpl$seb0)]=as.numeric(r18_bsave_mpl$b0_LL[!is.na(r18_bsave_mpl$seb0)]<=-1 & -1<=r18_bsave_mpl$b0_UL[!is.na(r18_bsave_mpl$seb0)])
r18_bsave_mpl$b1_cov[!is.na(r18_bsave_mpl$seb0)]=as.numeric(r18_bsave_mpl$b1_LL[!is.na(r18_bsave_mpl$seb0)]<=1 & 1<=r18_bsave_mpl$b1_UL[!is.na(r18_bsave_mpl$seb0)])

r18_gsave_mpl$g1_LL[!is.na(r18_gsave_mpl$seg1)]=r18_gsave_mpl$g1[!is.na(r18_gsave_mpl$seg1)]-1.96*r18_gsave_mpl$seg1[!is.na(r18_gsave_mpl$seg1)]
r18_gsave_mpl$g1_UL[!is.na(r18_gsave_mpl$seg1)]=r18_gsave_mpl$g1[!is.na(r18_gsave_mpl$seg1)]+1.96*r18_gsave_mpl$seg1[!is.na(r18_gsave_mpl$seg1)]
r18_gsave_mpl$g1_cov=rep(0,100)
r18_gsave_mpl$g1_cov[!is.na(r18_gsave_mpl$seg1)]=as.numeric(r18_gsave_mpl$g1_LL[!is.na(r18_gsave_mpl$seg1)]<=0.5 & 0.5<=r18_gsave_mpl$g1_UL[!is.na(r18_gsave_mpl$seg1)])


#em

mean(r18_bsave_em[,1])
mean(r18_bsave_em[,2])
mean(r18_bsave_em[,3])
sd(r18_bsave_em[,1])
mean(r18_bsave_em[,4])
sd(r18_bsave_em[,2])

mean(r18_gsave_em[,1])
mean(r18_gsave_em[,2])
sd(r18_gsave_em[,1])

r18_bsave_em$b0_LL=r18_bsave_em$b0_UL=r18_bsave_em$b1_LL=r18_bsave_em$b1_UL=rep(0,100)
r18_gsave_em$g1_LL=r18_gsave_em$g1_UL=rep(0,100)

r18_bsave_em$b0_LL=r18_bsave_em$b0-1.96*r18_bsave_em$seb0
r18_bsave_em$b0_UL=r18_bsave_em$b0+1.96*r18_bsave_em$seb0

r18_bsave_em$b1_LL=r18_bsave_em$b1-1.96*r18_bsave_em$seb1
r18_bsave_em$b1_UL=r18_bsave_em$b1+1.96*r18_bsave_em$seb1

r18_bsave_em$b0_cov=r18_bsave_em$b1_cov=rep(0,100)
r18_bsave_em$b0_cov=as.numeric(r18_bsave_em$b0_LL<=-1 & -1<=r18_bsave_em$b0_UL)
r18_bsave_em$b1_cov=as.numeric(r18_bsave_em$b1_LL<=1 & 1<=r18_bsave_em$b1_UL)

r18_gsave_em$g1_LL=r18_gsave_em$g1-1.96*r18_gsave_em$seg1
r18_gsave_em$g1_UL=r18_gsave_em$g1+1.96*r18_gsave_em$seg1
r18_gsave_em$g1_cov=rep(0,100)
r18_gsave_em$g1_cov=as.numeric(r18_gsave_em$g1_LL<=0.5 & 0.5<=r18_gsave_em$g1_UL)

# 22. beta=[1,1], n=500, censor 4.2 -------



r22_censorsave=matrix(rep(0,2*100),ncol=2)#whole dataset event proportion, non-cured event proportion
r22_bsave_mpl=matrix(rep(0,4*100),ncol=4)
r22_gsave_mpl=matrix(rep(0,2*100),ncol=2)
r22_bsave_em=matrix(rep(0,4*100),ncol=4)
r22_gsave_em=matrix(rep(0,2*100),ncol=2)


for(r in 1:100){
  sim.data=gen_data_right(500,4.2,c(1,1),c(0.5),c("b"),c("b"),c("z1","x1"))
  
  sim.data.uncure=sim.data[sim.data[,3]==1,]
  sim.surv1=Surv(time=sim.data.uncure[,1],event=sim.data.uncure[,2])
  sim.surv=Surv(time=sim.data$time,event =sim.data$event)
  
  r22_censorsave[r,1]=sum(sim.surv[,2])/500
  r22_censorsave[r,2]=sum(sim.surv1[,2])/nrow(sim.surv1)
  
  test=phmc_mpl(sim.surv~sim.data$X1,pi.formula = ~sim.data$Z,data=sim.data,phmc_mpl.control(conv_limit = 1e-4, n.knots=c(4,0),maxIter = c(1,10000,10001)))
  
  r22_bsave_mpl[r,1]=test$beta[1]
  r22_bsave_mpl[r,2]=test$beta[2]
  r22_bsave_mpl[r,3]=test$se$se_H[1]
  r22_bsave_mpl[r,4]=test$se$se_H[2]
  
  r22_gsave_mpl[r,1]=test$gamma[1]
  r22_gsave_mpl[r,2]=test$se$se_H[3]
  
  smc1 = smcure(sim.surv~X1, cureform=~Z, data=sim.data, model="ph")
  
  r22_bsave_em[r,1]=smc1$b[1]
  r22_bsave_em[r,2]=smc1$b[2]
  r22_bsave_em[r,3]=smc1$b_sd[1]
  r22_bsave_em[r,4]=smc1$b_sd[2]
  
  r22_gsave_em[r,1]=smc1$beta[1]
  r22_gsave_em[r,2]=smc1$beta_sd[1]
  
  print(c("22",r,test$se$se_H[1]))
  
}

print(c("done :)")) 

mean(r22_censorsave[!is.na(r22_bsave_mpl[,3])],1)
mean(r22_censorsave[!is.na(r22_bsave_mpl[,3])],2)
r22_censorsave=data.frame(r22_censorsave)
write.table(r22_censorsave, "r22censor.csv",sep=",",dec=".")

colnames(r22_bsave_mpl)=c("b0","b1","seb0","seb1")
r22_bsave_mpl=data.frame(r22_bsave_mpl)
write.table(r22_bsave_mpl, "r22bsave_mpl.csv",sep=",",dec=".",row.names = F)
colnames(r22_bsave_em)=c("b0","b1","seb0","seb1")
r22_bsave_em=data.frame(r22_bsave_em)
write.table(r22_bsave_em, "r22bsave_em.csv",sep=",",dec=".",row.names = F)
colnames(r22_gsave_mpl)=c("g1","seg1")
r22_gsave_mpl=data.frame(r22_gsave_mpl)
write.table(r22_gsave_mpl, "r22gsave_mpl.csv",sep=",",dec=".",row.names = F)
colnames(r22_gsave_em)=c("g1","seg1")
r22_gsave_em=data.frame(r22_gsave_em)
write.table(r22_gsave_em, "r22gsave_em.csv",sep=",",dec=".",row.names = F)


#mpl
mean(r22_bsave_mpl[!is.na(r22_bsave_mpl[,3]),1])
mean(r22_bsave_mpl[!is.na(r22_bsave_mpl[,3]),2])
mean(r22_bsave_mpl[!is.na(r22_bsave_mpl[,3]),3])
sd(r22_bsave_mpl[!is.na(r22_bsave_mpl[,3]),1])
mean(r22_bsave_mpl[!is.na(r22_bsave_mpl[,3]),4])
sd(r22_bsave_mpl[!is.na(r22_bsave_mpl[,3]),2])

mean(r22_gsave_mpl[!is.na(r22_gsave_mpl[,2]),1])
mean(r22_gsave_mpl[!is.na(r22_gsave_mpl[,2]),2])
sd(r22_gsave_mpl[!is.na(r22_gsave_mpl[,2]),1])

r22_bsave_mpl$b0_LL=r22_bsave_mpl$b0_UL=r22_bsave_mpl$b1_LL=r22_bsave_mpl$b1_UL=rep(0,100)
r22_gsave_mpl$g1_LL=r22_gsave_mpl$g1_UL=rep(0,100)

r22_bsave_mpl$b0_LL[!is.na(r22_bsave_mpl$seb0)]=r22_bsave_mpl$b0[!is.na(r22_bsave_mpl$seb0)]-1.96*r22_bsave_mpl$seb0[!is.na(r22_bsave_mpl$seb0)]
r22_bsave_mpl$b0_UL[!is.na(r22_bsave_mpl$seb0)]=r22_bsave_mpl$b0[!is.na(r22_bsave_mpl$seb0)]+1.96*r22_bsave_mpl$seb0[!is.na(r22_bsave_mpl$seb0)]

r22_bsave_mpl$b1_LL[!is.na(r22_bsave_mpl$seb0)]=r22_bsave_mpl$b1[!is.na(r22_bsave_mpl$seb0)]-1.96*r22_bsave_mpl$seb1[!is.na(r22_bsave_mpl$seb0)]
r22_bsave_mpl$b1_UL[!is.na(r22_bsave_mpl$seb0)]=r22_bsave_mpl$b1[!is.na(r22_bsave_mpl$seb0)]+1.96*r22_bsave_mpl$seb1[!is.na(r22_bsave_mpl$seb0)]

r22_bsave_mpl$b0_cov=r22_bsave_mpl$b1_cov=rep(0,100)
r22_bsave_mpl$b0_cov[!is.na(r22_bsave_mpl$seb0)]=as.numeric(r22_bsave_mpl$b0_LL[!is.na(r22_bsave_mpl$seb0)]<=1 & 1<=r22_bsave_mpl$b0_UL[!is.na(r22_bsave_mpl$seb0)])
r22_bsave_mpl$b1_cov[!is.na(r22_bsave_mpl$seb0)]=as.numeric(r22_bsave_mpl$b1_LL[!is.na(r22_bsave_mpl$seb0)]<=1 & 1<=r22_bsave_mpl$b1_UL[!is.na(r22_bsave_mpl$seb0)])

r22_gsave_mpl$g1_LL[!is.na(r22_gsave_mpl$seg1)]=r22_gsave_mpl$g1[!is.na(r22_gsave_mpl$seg1)]-1.96*r22_gsave_mpl$seg1[!is.na(r22_gsave_mpl$seg1)]
r22_gsave_mpl$g1_UL[!is.na(r22_gsave_mpl$seg1)]=r22_gsave_mpl$g1[!is.na(r22_gsave_mpl$seg1)]+1.96*r22_gsave_mpl$seg1[!is.na(r22_gsave_mpl$seg1)]
r22_gsave_mpl$g1_cov=rep(0,100)
r22_gsave_mpl$g1_cov[!is.na(r22_gsave_mpl$seg1)]=as.numeric(r22_gsave_mpl$g1_LL[!is.na(r22_gsave_mpl$seg1)]<=0.5 & 0.5<=r22_gsave_mpl$g1_UL[!is.na(r22_gsave_mpl$seg1)])


#em

mean(r22_bsave_em[,1])
mean(r22_bsave_em[,2])
mean(r22_bsave_em[,3])
sd(r22_bsave_em[,1])
mean(r22_bsave_em[,4])
sd(r22_bsave_em[,2])

mean(r22_gsave_em[,1])
mean(r22_gsave_em[,2])
sd(r22_gsave_em[,1])

r22_bsave_em$b0_LL=r22_bsave_em$b0_UL=r22_bsave_em$b1_LL=r22_bsave_em$b1_UL=rep(0,100)
r22_gsave_em$g1_LL=r22_gsave_em$g1_UL=rep(0,100)

r22_bsave_em$b0_LL=r22_bsave_em$b0-1.96*r22_bsave_em$seb0
r22_bsave_em$b0_UL=r22_bsave_em$b0+1.96*r22_bsave_em$seb0

r22_bsave_em$b1_LL=r22_bsave_em$b1-1.96*r22_bsave_em$seb1
r22_bsave_em$b1_UL=r22_bsave_em$b1+1.96*r22_bsave_em$seb1

r22_bsave_em$b0_cov=r22_bsave_em$b1_cov=rep(0,100)
r22_bsave_em$b0_cov=as.numeric(r22_bsave_em$b0_LL<=1 & 1<=r22_bsave_em$b0_UL)
r22_bsave_em$b1_cov=as.numeric(r22_bsave_em$b1_LL<=1 & 1<=r22_bsave_em$b1_UL)

r22_gsave_em$g1_LL=r22_gsave_em$g1-1.96*r22_gsave_em$seg1
r22_gsave_em$g1_UL=r22_gsave_em$g1+1.96*r22_gsave_em$seg1
r22_gsave_em$g1_cov=rep(0,100)
r22_gsave_em$g1_cov=as.numeric(r22_gsave_em$g1_LL<=0.5 & 0.5<=r22_gsave_em$g1_UL)

# 23. beta=[0,1], n=500, censor 4.2 -------


r23_censorsave=matrix(rep(0,2*100),ncol=2)#whole dataset event proportion, non-cured event proportion
r23_bsave_mpl=matrix(rep(0,4*100),ncol=4)
r23_gsave_mpl=matrix(rep(0,2*100),ncol=2)
r23_bsave_em=matrix(rep(0,4*100),ncol=4)
r23_gsave_em=matrix(rep(0,2*100),ncol=2)


for(r in 1:100){
  sim.data=gen_data_right(500,4.2,c(0,1),c(0.5),c("b"),c("b"),c("z1","x1"))
  
  sim.data.uncure=sim.data[sim.data[,3]==1,]
  sim.surv1=Surv(time=sim.data.uncure[,1],event=sim.data.uncure[,2])
  sim.surv=Surv(time=sim.data$time,event =sim.data$event)
  
  r23_censorsave[r,1]=sum(sim.surv[,2])/500
  r23_censorsave[r,2]=sum(sim.surv1[,2])/nrow(sim.surv1)
  
  test=phmc_mpl(sim.surv~sim.data$X1,pi.formula = ~sim.data$Z,data=sim.data,phmc_mpl.control(conv_limit = 1e-4, n.knots=c(2,0),maxIter = c(1,10000,10001)))
  
  r23_bsave_mpl[r,1]=test$beta[1]
  r23_bsave_mpl[r,2]=test$beta[2]
  r23_bsave_mpl[r,3]=test$se$se_H[1]
  r23_bsave_mpl[r,4]=test$se$se_H[2]
  
  r23_gsave_mpl[r,1]=test$gamma[1]
  r23_gsave_mpl[r,2]=test$se$se_H[3]
  
  smc1 = smcure(sim.surv~X1, cureform=~Z, data=sim.data, model="ph")
  
  r23_bsave_em[r,1]=smc1$b[1]
  r23_bsave_em[r,2]=smc1$b[2]
  r23_bsave_em[r,3]=smc1$b_sd[1]
  r23_bsave_em[r,4]=smc1$b_sd[2]
  
  r23_gsave_em[r,1]=smc1$beta[1]
  r23_gsave_em[r,2]=smc1$beta_sd[1]
  
  print(c("23",r,test$se$se_H[1]))
  
}


mean(r23_censorsave[!is.na(r23_bsave_mpl[,3])],1)
mean(r23_censorsave[!is.na(r23_bsave_mpl[,3])],2)
r23_censorsave=data.frame(r23_censorsave)
write.table(r23_censorsave, "r23censor.csv",sep=",",dec=".")

colnames(r23_bsave_mpl)=c("b0","b1","seb0","seb1")
r23_bsave_mpl=data.frame(r23_bsave_mpl)
write.table(r23_bsave_mpl, "r23bsave_mpl.csv",sep=",",dec=".",row.names = F)
colnames(r23_bsave_em)=c("b0","b1","seb0","seb1")
r23_bsave_em=data.frame(r23_bsave_em)
write.table(r23_bsave_em, "r23bsave_em.csv",sep=",",dec=".",row.names = F)
colnames(r23_gsave_mpl)=c("g1","seg1")
r23_gsave_mpl=data.frame(r23_gsave_mpl)
write.table(r23_gsave_mpl, "r23gsave_mpl.csv",sep=",",dec=".",row.names = F)
colnames(r23_gsave_em)=c("g1","seg1")
r23_gsave_em=data.frame(r23_gsave_em)
write.table(r23_gsave_em, "r23gsave_em.csv",sep=",",dec=".",row.names = F)


#mpl
mean(r23_bsave_mpl[!is.na(r23_bsave_mpl[,3]),1])
mean(r23_bsave_mpl[!is.na(r23_bsave_mpl[,3]),2])
mean(r23_bsave_mpl[!is.na(r23_bsave_mpl[,3]),3])
sd(r23_bsave_mpl[!is.na(r23_bsave_mpl[,3]),1])
mean(r23_bsave_mpl[!is.na(r23_bsave_mpl[,3]),4])
sd(r23_bsave_mpl[!is.na(r23_bsave_mpl[,3]),2])

mean(r23_gsave_mpl[!is.na(r23_gsave_mpl[,2]),1])
mean(r23_gsave_mpl[!is.na(r23_gsave_mpl[,2]),2])
sd(r23_gsave_mpl[!is.na(r23_gsave_mpl[,2]),1])

r23_bsave_mpl$b0_LL=r23_bsave_mpl$b0_UL=r23_bsave_mpl$b1_LL=r23_bsave_mpl$b1_UL=rep(0,100)
r23_gsave_mpl$g1_LL=r23_gsave_mpl$g1_UL=rep(0,100)

r23_bsave_mpl$b0_LL[!is.na(r23_bsave_mpl$seb0)]=r23_bsave_mpl$b0[!is.na(r23_bsave_mpl$seb0)]-1.96*r23_bsave_mpl$seb0[!is.na(r23_bsave_mpl$seb0)]
r23_bsave_mpl$b0_UL[!is.na(r23_bsave_mpl$seb0)]=r23_bsave_mpl$b0[!is.na(r23_bsave_mpl$seb0)]+1.96*r23_bsave_mpl$seb0[!is.na(r23_bsave_mpl$seb0)]

r23_bsave_mpl$b1_LL[!is.na(r23_bsave_mpl$seb0)]=r23_bsave_mpl$b1[!is.na(r23_bsave_mpl$seb0)]-1.96*r23_bsave_mpl$seb1[!is.na(r23_bsave_mpl$seb0)]
r23_bsave_mpl$b1_UL[!is.na(r23_bsave_mpl$seb0)]=r23_bsave_mpl$b1[!is.na(r23_bsave_mpl$seb0)]+1.96*r23_bsave_mpl$seb1[!is.na(r23_bsave_mpl$seb0)]

r23_bsave_mpl$b0_cov=r23_bsave_mpl$b1_cov=rep(0,100)
r23_bsave_mpl$b0_cov[!is.na(r23_bsave_mpl$seb0)]=as.numeric(r23_bsave_mpl$b0_LL[!is.na(r23_bsave_mpl$seb0)]<=0 & 0<=r23_bsave_mpl$b0_UL[!is.na(r23_bsave_mpl$seb0)])
r23_bsave_mpl$b1_cov[!is.na(r23_bsave_mpl$seb0)]=as.numeric(r23_bsave_mpl$b1_LL[!is.na(r23_bsave_mpl$seb0)]<=1 & 1<=r23_bsave_mpl$b1_UL[!is.na(r23_bsave_mpl$seb0)])

r23_gsave_mpl$g1_LL[!is.na(r23_gsave_mpl$seg1)]=r23_gsave_mpl$g1[!is.na(r23_gsave_mpl$seg1)]-1.96*r23_gsave_mpl$seg1[!is.na(r23_gsave_mpl$seg1)]
r23_gsave_mpl$g1_UL[!is.na(r23_gsave_mpl$seg1)]=r23_gsave_mpl$g1[!is.na(r23_gsave_mpl$seg1)]+1.96*r23_gsave_mpl$seg1[!is.na(r23_gsave_mpl$seg1)]
r23_gsave_mpl$g1_cov=rep(0,100)
r23_gsave_mpl$g1_cov[!is.na(r23_gsave_mpl$seg1)]=as.numeric(r23_gsave_mpl$g1_LL[!is.na(r23_gsave_mpl$seg1)]<=0.5 & 0.5<=r23_gsave_mpl$g1_UL[!is.na(r23_gsave_mpl$seg1)])


#em

mean(r23_bsave_em[,1])
mean(r23_bsave_em[,2])
mean(r23_bsave_em[,3])
sd(r23_bsave_em[,1])
mean(r23_bsave_em[,4])
sd(r23_bsave_em[,2])

mean(r23_gsave_em[,1])
mean(r23_gsave_em[,2])
sd(r23_gsave_em[,1])

r23_bsave_em$b0_LL=r23_bsave_em$b0_UL=r23_bsave_em$b1_LL=r23_bsave_em$b1_UL=rep(0,100)
r23_gsave_em$g1_LL=r23_gsave_em$g1_UL=rep(0,100)

r23_bsave_em$b0_LL=r23_bsave_em$b0-1.96*r23_bsave_em$seb0
r23_bsave_em$b0_UL=r23_bsave_em$b0+1.96*r23_bsave_em$seb0

r23_bsave_em$b1_LL=r23_bsave_em$b1-1.96*r23_bsave_em$seb1
r23_bsave_em$b1_UL=r23_bsave_em$b1+1.96*r23_bsave_em$seb1

r23_bsave_em$b0_cov=r23_bsave_em$b1_cov=rep(0,100)
r23_bsave_em$b0_cov=as.numeric(r23_bsave_em$b0_LL<=0 & 0<=r23_bsave_em$b0_UL)
r23_bsave_em$b1_cov=as.numeric(r23_bsave_em$b1_LL<=1 & 1<=r23_bsave_em$b1_UL)

r23_gsave_em$g1_LL=r23_gsave_em$g1-1.96*r23_gsave_em$seg1
r23_gsave_em$g1_UL=r23_gsave_em$g1+1.96*r23_gsave_em$seg1
r23_gsave_em$g1_cov=rep(0,100)
r23_gsave_em$g1_cov=as.numeric(r23_gsave_em$g1_LL<=0.5 & 0.5<=r23_gsave_em$g1_UL)

# 24. beta=[-1,1], n=500, censor 4.2 ---------



r24_censorsave=matrix(rep(0,2*100),ncol=2)#whole dataset event proportion, non-cured event proportion
r24_bsave_mpl=matrix(rep(0,4*100),ncol=4)
r24_gsave_mpl=matrix(rep(0,2*100),ncol=2)
r24_bsave_em=matrix(rep(0,4*100),ncol=4)
r24_gsave_em=matrix(rep(0,2*100),ncol=2)

for(r in 1:100){
  sim.data=gen_data_right(500,4.2,c(-1,1),c(0.5),c("b"),c("b"),c("z1","x1"))
  
  sim.data.uncure=sim.data[sim.data[,3]==1,]
  sim.surv1=Surv(time=sim.data.uncure[,1],event=sim.data.uncure[,2])
  sim.surv=Surv(time=sim.data$time,event =sim.data$event)
  r24_censorsave[r,1]=sum(sim.surv[,2])/500
  r24_censorsave[r,2]=sum(sim.surv1[,2])/nrow(sim.surv1)
  
  try1=try(phmc_mpl(sim.surv~sim.data$X1,pi.formula = ~sim.data$Z,data=sim.data,phmc_mpl.control(conv_limit = 1e-4, n.knots=c(1,0),maxIter = c(1,20000,20001))))
  try2=try(smcure(sim.surv~X1, cureform=~Z, data=sim.data, model="ph"))
  if(class(try1)!="try-error" & class(try2)!="try-error"){
    test=try1
    smc1=try2
    
    r24_bsave_mpl[r,1]=test$beta[1]
    r24_bsave_mpl[r,2]=test$beta[2]
    r24_bsave_mpl[r,3]=test$se$se_H[1]
    r24_bsave_mpl[r,4]=test$se$se_H[2]
    
    r24_gsave_mpl[r,1]=test$gamma[1]
    r24_gsave_mpl[r,2]=test$se$se_H[3]
    
    r24_bsave_em[r,1]=smc1$b[1]
    r24_bsave_em[r,2]=smc1$b[2]
    r24_bsave_em[r,3]=smc1$b_sd[1]
    r24_bsave_em[r,4]=smc1$b_sd[2]
    
    r24_gsave_em[r,1]=smc1$beta[1]
    r24_gsave_em[r,2]=smc1$beta_sd[1]
    
    print(c(r,"both 24A",test$se$se_H[1]))
  } else if(class(try1)!="try-error"){
    test=try1
    r24_bsave_mpl[r,1]=test$beta[1]
    r24_bsave_mpl[r,2]=test$beta[2]
    r24_bsave_mpl[r,3]=test$se$se_H[1]
    r24_bsave_mpl[r,4]=test$se$se_H[2]
    
    r24_gsave_mpl[r,1]=test$gamma[1]
    r24_gsave_mpl[r,2]=test$se$se_H[3]
    print(c(r, "mpl 24A",test$se$se_H[1]))
    
  } else if(class(try2)!="try-error"){
    smc1=try2
    r24_bsave_em[r,1]=smc1$b[1]
    r24_bsave_em[r,2]=smc1$b[2]
    r24_bsave_em[r,3]=smc1$b_sd[1]
    r24_bsave_em[r,4]=smc1$b_sd[2]
    
    r24_gsave_em[r,1]=smc1$beta[1]
    r24_gsave_em[r,2]=smc1$beta_sd[1]
    print(c(r, "em 24A"))
  } else{
    print(c(r, "neither 24A"))
  }
}




mean(r24_censorsave[!is.na(r24_bsave_mpl[,3])],1)
mean(r24_censorsave[!is.na(r24_bsave_mpl[,3])],2)
r24_censorsave=data.frame(r24_censorsave)
write.table(r24_censorsave, "r24censor.csv",sep=",",dec=".")

colnames(r24_bsave_mpl)=c("b0","b1","seb0","seb1")
r24_bsave_mpl=data.frame(r24_bsave_mpl)
write.table(r24_bsave_mpl, "r24bsave_mpl.csv",sep=",",dec=".",row.names = F)
colnames(r24_bsave_em)=c("b0","b1","seb0","seb1")
r24_bsave_em=data.frame(r24_bsave_em)
write.table(r24_bsave_em, "r24bsave_em.csv",sep=",",dec=".",row.names = F)
colnames(r24_gsave_mpl)=c("g1","seg1")
r24_gsave_mpl=data.frame(r24_gsave_mpl)
write.table(r24_gsave_mpl, "r24gsave_mpl.csv",sep=",",dec=".",row.names = F)
colnames(r24_gsave_em)=c("g1","seg1")
r24_gsave_em=data.frame(r24_gsave_em)
write.table(r24_gsave_em, "r24gsave_em.csv",sep=",",dec=".",row.names = F)


#mpl
mean(r24_bsave_mpl[!is.na(r24_bsave_mpl[,3]),1])
mean(r24_bsave_mpl[!is.na(r24_bsave_mpl[,3]),2])
mean(r24_bsave_mpl[!is.na(r24_bsave_mpl[,3]),3])
sd(r24_bsave_mpl[!is.na(r24_bsave_mpl[,3]),1])
mean(r24_bsave_mpl[!is.na(r24_bsave_mpl[,3]),4])
sd(r24_bsave_mpl[!is.na(r24_bsave_mpl[,3]),2])

mean(r24_gsave_mpl[!is.na(r24_gsave_mpl[,2]),1])
mean(r24_gsave_mpl[!is.na(r24_gsave_mpl[,2]),2])
sd(r24_gsave_mpl[!is.na(r24_gsave_mpl[,2]),1])

r24_bsave_mpl$b0_LL=r24_bsave_mpl$b0_UL=r24_bsave_mpl$b1_LL=r24_bsave_mpl$b1_UL=rep(0,100)
r24_gsave_mpl$g1_LL=r24_gsave_mpl$g1_UL=rep(0,100)

r24_bsave_mpl$b0_LL[!is.na(r24_bsave_mpl$seb0)]=r24_bsave_mpl$b0[!is.na(r24_bsave_mpl$seb0)]-1.96*r24_bsave_mpl$seb0[!is.na(r24_bsave_mpl$seb0)]
r24_bsave_mpl$b0_UL[!is.na(r24_bsave_mpl$seb0)]=r24_bsave_mpl$b0[!is.na(r24_bsave_mpl$seb0)]+1.96*r24_bsave_mpl$seb0[!is.na(r24_bsave_mpl$seb0)]

r24_bsave_mpl$b1_LL[!is.na(r24_bsave_mpl$seb0)]=r24_bsave_mpl$b1[!is.na(r24_bsave_mpl$seb0)]-1.96*r24_bsave_mpl$seb1[!is.na(r24_bsave_mpl$seb0)]
r24_bsave_mpl$b1_UL[!is.na(r24_bsave_mpl$seb0)]=r24_bsave_mpl$b1[!is.na(r24_bsave_mpl$seb0)]+1.96*r24_bsave_mpl$seb1[!is.na(r24_bsave_mpl$seb0)]

r24_bsave_mpl$b0_cov=r24_bsave_mpl$b1_cov=rep(0,100)
r24_bsave_mpl$b0_cov[!is.na(r24_bsave_mpl$seb0)]=as.numeric(r24_bsave_mpl$b0_LL[!is.na(r24_bsave_mpl$seb0)]<=-1 & -1<=r24_bsave_mpl$b0_UL[!is.na(r24_bsave_mpl$seb0)])
r24_bsave_mpl$b1_cov[!is.na(r24_bsave_mpl$seb0)]=as.numeric(r24_bsave_mpl$b1_LL[!is.na(r24_bsave_mpl$seb0)]<=1 & 1<=r24_bsave_mpl$b1_UL[!is.na(r24_bsave_mpl$seb0)])

r24_gsave_mpl$g1_LL[!is.na(r24_gsave_mpl$seg1)]=r24_gsave_mpl$g1[!is.na(r24_gsave_mpl$seg1)]-1.96*r24_gsave_mpl$seg1[!is.na(r24_gsave_mpl$seg1)]
r24_gsave_mpl$g1_UL[!is.na(r24_gsave_mpl$seg1)]=r24_gsave_mpl$g1[!is.na(r24_gsave_mpl$seg1)]+1.96*r24_gsave_mpl$seg1[!is.na(r24_gsave_mpl$seg1)]
r24_gsave_mpl$g1_cov=rep(0,100)
r24_gsave_mpl$g1_cov[!is.na(r24_gsave_mpl$seg1)]=as.numeric(r24_gsave_mpl$g1_LL[!is.na(r24_gsave_mpl$seg1)]<=0.5 & 0.5<=r24_gsave_mpl$g1_UL[!is.na(r24_gsave_mpl$seg1)])


#em

mean(r24_bsave_em[,1])
mean(r24_bsave_em[,2])
mean(r24_bsave_em[,3])
sd(r24_bsave_em[,1])
mean(r24_bsave_em[,4])
sd(r24_bsave_em[,2])

mean(r24_gsave_em[,1])
mean(r24_gsave_em[,2])
sd(r24_gsave_em[,1])

r24_bsave_em$b0_LL=r24_bsave_em$b0_UL=r24_bsave_em$b1_LL=r24_bsave_em$b1_UL=rep(0,100)
r24_gsave_em$g1_LL=r24_gsave_em$g1_UL=rep(0,100)

r24_bsave_em$b0_LL=r24_bsave_em$b0-1.96*r24_bsave_em$seb0
r24_bsave_em$b0_UL=r24_bsave_em$b0+1.96*r24_bsave_em$seb0

r24_bsave_em$b1_LL=r24_bsave_em$b1-1.96*r24_bsave_em$seb1
r24_bsave_em$b1_UL=r24_bsave_em$b1+1.96*r24_bsave_em$seb1

r24_bsave_em$b0_cov=r24_bsave_em$b1_cov=rep(0,100)
r24_bsave_em$b0_cov=as.numeric(r24_bsave_em$b0_LL<=-1 & -1<=r24_bsave_em$b0_UL)
r24_bsave_em$b1_cov=as.numeric(r24_bsave_em$b1_LL<=1 & 1<=r24_bsave_em$b1_UL)

r24_gsave_em$g1_LL=r24_gsave_em$g1-1.96*r24_gsave_em$seg1
r24_gsave_em$g1_UL=r24_gsave_em$g1+1.96*r24_gsave_em$seg1
r24_gsave_em$g1_cov=rep(0,100)
r24_gsave_em$g1_cov=as.numeric(r24_gsave_em$g1_LL<=0.5 & 0.5<=r24_gsave_em$g1_UL)


# 25. beta=[1,1], n=100, censor 4.2 ------------


r25_censorsave=matrix(rep(0,2*100),ncol=2)#whole dataset event proportion, non-cured event proportion
r25_bsave_mpl=matrix(rep(0,4*100),ncol=4)
r25_gsave_mpl=matrix(rep(0,2*100),ncol=2)
r25_bsave_em=matrix(rep(0,4*100),ncol=4)
r25_gsave_em=matrix(rep(0,2*100),ncol=2)


for(r in 65:100){
  sim.data=gen_data_right(100,4.2,c(1,1),c(0.5),c("b"),c("b"),c("z1","x1"))
  
  sim.data.uncure=sim.data[sim.data[,3]==1,]
  sim.surv1=Surv(time=sim.data.uncure[,1],event=sim.data.uncure[,2])
  sim.surv=Surv(time=sim.data$time,event =sim.data$event)
  
  r25_censorsave[r,1]=sum(sim.surv[,2])/100
  r25_censorsave[r,2]=sum(sim.surv1[,2])/nrow(sim.surv1)
  
  test=phmc_mpl(sim.surv~sim.data$X1,pi.formula = ~sim.data$Z,data=sim.data,phmc_mpl.control(conv_limit = 1e-4, n.knots=c(3,0),maxIter = c(1,10000,10001)))
  
  r25_bsave_mpl[r,1]=test$beta[1]
  r25_bsave_mpl[r,2]=test$beta[2]
  r25_bsave_mpl[r,3]=test$se$se_H[1]
  r25_bsave_mpl[r,4]=test$se$se_H[2]
  
  r25_gsave_mpl[r,1]=test$gamma[1]
  r25_gsave_mpl[r,2]=test$se$se_H[3]
  
  smc1 = smcure(sim.surv~X1, cureform=~Z, data=sim.data, model="ph")
  
  r25_bsave_em[r,1]=smc1$b[1]
  r25_bsave_em[r,2]=smc1$b[2]
  r25_bsave_em[r,3]=smc1$b_sd[1]
  r25_bsave_em[r,4]=smc1$b_sd[2]
  
  r25_gsave_em[r,1]=smc1$beta[1]
  r25_gsave_em[r,2]=smc1$beta_sd[1]
  
  print(c(r,r25_bsave_mpl[r,1],test$se$se_H[1]))
  
}


mean(r25_censorsave[!is.na(r25_bsave_mpl[,3])],1)
mean(r25_censorsave[!is.na(r25_bsave_mpl[,3])],2)
r25_censorsave=data.frame(r25_censorsave)
write.table(r25_censorsave, "r25censor.csv",sep=",",dec=".")

colnames(r25_bsave_mpl)=c("b0","b1","seb0","seb1")
r25_bsave_mpl=data.frame(r25_bsave_mpl)
write.table(r25_bsave_mpl, "r25bsave_mpl.csv",sep=",",dec=".",row.names = F)
colnames(r25_bsave_em)=c("b0","b1","seb0","seb1")
r25_bsave_em=data.frame(r25_bsave_em)
write.table(r25_bsave_em, "r25bsave_em.csv",sep=",",dec=".",row.names = F)
colnames(r25_gsave_mpl)=c("g1","seg1")
r25_gsave_mpl=data.frame(r25_gsave_mpl)
write.table(r25_gsave_mpl, "r25gsave_mpl.csv",sep=",",dec=".",row.names = F)
colnames(r25_gsave_em)=c("g1","seg1")
r25_gsave_em=data.frame(r25_gsave_em)
write.table(r25_gsave_em, "r25gsave_em.csv",sep=",",dec=".",row.names = F)


#mpl
mean(r25_bsave_mpl[!is.na(r25_bsave_mpl[,3]),1])
mean(r25_bsave_mpl[!is.na(r25_bsave_mpl[,3]),2])
mean(r25_bsave_mpl[!is.na(r25_bsave_mpl[,3]),3])
sd(r25_bsave_mpl[!is.na(r25_bsave_mpl[,3]),1])
mean(r25_bsave_mpl[!is.na(r25_bsave_mpl[,3]),4])
sd(r25_bsave_mpl[!is.na(r25_bsave_mpl[,3]),2])

mean(r25_gsave_mpl[!is.na(r25_gsave_mpl[,2]),1])
mean(r25_gsave_mpl[!is.na(r25_gsave_mpl[,2]),2])
sd(r25_gsave_mpl[!is.na(r25_gsave_mpl[,2]),1])

r25_bsave_mpl$b0_LL=r25_bsave_mpl$b0_UL=r25_bsave_mpl$b1_LL=r25_bsave_mpl$b1_UL=rep(0,100)
r25_gsave_mpl$g1_LL=r25_gsave_mpl$g1_UL=rep(0,100)

r25_bsave_mpl$b0_LL[!is.na(r25_bsave_mpl$seb0)]=r25_bsave_mpl$b0[!is.na(r25_bsave_mpl$seb0)]-1.96*r25_bsave_mpl$seb0[!is.na(r25_bsave_mpl$seb0)]
r25_bsave_mpl$b0_UL[!is.na(r25_bsave_mpl$seb0)]=r25_bsave_mpl$b0[!is.na(r25_bsave_mpl$seb0)]+1.96*r25_bsave_mpl$seb0[!is.na(r25_bsave_mpl$seb0)]

r25_bsave_mpl$b1_LL[!is.na(r25_bsave_mpl$seb0)]=r25_bsave_mpl$b1[!is.na(r25_bsave_mpl$seb0)]-1.96*r25_bsave_mpl$seb1[!is.na(r25_bsave_mpl$seb0)]
r25_bsave_mpl$b1_UL[!is.na(r25_bsave_mpl$seb0)]=r25_bsave_mpl$b1[!is.na(r25_bsave_mpl$seb0)]+1.96*r25_bsave_mpl$seb1[!is.na(r25_bsave_mpl$seb0)]

r25_bsave_mpl$b0_cov=r25_bsave_mpl$b1_cov=rep(0,100)
r25_bsave_mpl$b0_cov[!is.na(r25_bsave_mpl$seb0)]=as.numeric(r25_bsave_mpl$b0_LL[!is.na(r25_bsave_mpl$seb0)]<=1 & 1<=r25_bsave_mpl$b0_UL[!is.na(r25_bsave_mpl$seb0)])
r25_bsave_mpl$b1_cov[!is.na(r25_bsave_mpl$seb0)]=as.numeric(r25_bsave_mpl$b1_LL[!is.na(r25_bsave_mpl$seb0)]<=1 & 1<=r25_bsave_mpl$b1_UL[!is.na(r25_bsave_mpl$seb0)])

r25_gsave_mpl$g1_LL[!is.na(r25_gsave_mpl$seg1)]=r25_gsave_mpl$g1[!is.na(r25_gsave_mpl$seg1)]-1.96*r25_gsave_mpl$seg1[!is.na(r25_gsave_mpl$seg1)]
r25_gsave_mpl$g1_UL[!is.na(r25_gsave_mpl$seg1)]=r25_gsave_mpl$g1[!is.na(r25_gsave_mpl$seg1)]+1.96*r25_gsave_mpl$seg1[!is.na(r25_gsave_mpl$seg1)]
r25_gsave_mpl$g1_cov=rep(0,100)
r25_gsave_mpl$g1_cov[!is.na(r25_gsave_mpl$seg1)]=as.numeric(r25_gsave_mpl$g1_LL[!is.na(r25_gsave_mpl$seg1)]<=0.5 & 0.5<=r25_gsave_mpl$g1_UL[!is.na(r25_gsave_mpl$seg1)])


#em

mean(r25_bsave_em[,1])
mean(r25_bsave_em[,2])
mean(r25_bsave_em[,3])
sd(r25_bsave_em[,1])
mean(r25_bsave_em[,4])
sd(r25_bsave_em[,2])

mean(r25_gsave_em[,1])
mean(r25_gsave_em[,2])
sd(r25_gsave_em[,1])

r25_bsave_em$b0_LL=r25_bsave_em$b0_UL=r25_bsave_em$b1_LL=r25_bsave_em$b1_UL=rep(0,100)
r25_gsave_em$g1_LL=r25_gsave_em$g1_UL=rep(0,100)

r25_bsave_em$b0_LL=r25_bsave_em$b0-1.96*r25_bsave_em$seb0
r25_bsave_em$b0_UL=r25_bsave_em$b0+1.96*r25_bsave_em$seb0

r25_bsave_em$b1_LL=r25_bsave_em$b1-1.96*r25_bsave_em$seb1
r25_bsave_em$b1_UL=r25_bsave_em$b1+1.96*r25_bsave_em$seb1

r25_bsave_em$b0_cov=r25_bsave_em$b1_cov=rep(0,100)
r25_bsave_em$b0_cov=as.numeric(r25_bsave_em$b0_LL<=1 & 1<=r25_bsave_em$b0_UL)
r25_bsave_em$b1_cov=as.numeric(r25_bsave_em$b1_LL<=1 & 1<=r25_bsave_em$b1_UL)

r25_gsave_em$g1_LL=r25_gsave_em$g1-1.96*r25_gsave_em$seg1
r25_gsave_em$g1_UL=r25_gsave_em$g1+1.96*r25_gsave_em$seg1
r25_gsave_em$g1_cov=rep(0,100)
r25_gsave_em$g1_cov=as.numeric(r25_gsave_em$g1_LL<=0.5 & 0.5<=r25_gsave_em$g1_UL)

# 26. beta=[0,1], n=100, censor 4.2 ------




r26_censorsave=matrix(rep(0,2*100),ncol=2)#whole dataset event proportion, non-cured event proportion
r26_bsave_mpl=matrix(rep(0,4*100),ncol=4)
r26_gsave_mpl=matrix(rep(0,2*100),ncol=2)
r26_bsave_em=matrix(rep(0,4*100),ncol=4)
r26_gsave_em=matrix(rep(0,2*100),ncol=2)

for(r in 1:100){
  sim.data=gen_data_right(100,4.2,c(0,1),c(0.5),c("b"),c("b"),c("z1","x1"))
  
  sim.data.uncure=sim.data[sim.data[,3]==1,]
  sim.surv1=Surv(time=sim.data.uncure[,1],event=sim.data.uncure[,2])
  sim.surv=Surv(time=sim.data$time,event =sim.data$event)
  r26_censorsave[r,1]=sum(sim.surv[,2])/100
  r26_censorsave[r,2]=sum(sim.surv1[,2])/nrow(sim.surv1)
  
  try1=try(phmc_mpl(sim.surv~sim.data$X1,pi.formula = ~sim.data$Z,data=sim.data,phmc_mpl.control(conv_limit = 1e-4, n.knots=c(1,0),maxIter = c(1,20000,20001))))
  try2=try(smcure(sim.surv~X1, cureform=~Z, data=sim.data, model="ph"))
  if(class(try1)!="try-error" & class(try2)!="try-error"){
    test=try1
    smc1=try2
    
    r26_bsave_mpl[r,1]=test$beta[1]
    r26_bsave_mpl[r,2]=test$beta[2]
    r26_bsave_mpl[r,3]=test$se$se_H[1]
    r26_bsave_mpl[r,4]=test$se$se_H[2]
    
    r26_gsave_mpl[r,1]=test$gamma[1]
    r26_gsave_mpl[r,2]=test$se$se_H[3]
    
    r26_bsave_em[r,1]=smc1$b[1]
    r26_bsave_em[r,2]=smc1$b[2]
    r26_bsave_em[r,3]=smc1$b_sd[1]
    r26_bsave_em[r,4]=smc1$b_sd[2]
    
    r26_gsave_em[r,1]=smc1$beta[1]
    r26_gsave_em[r,2]=smc1$beta_sd[1]
    
    print(c(r,"both 26A",test$se$se_H[1]))
  } else if(class(try1)!="try-error"){
    test=try1
    r26_bsave_mpl[r,1]=test$beta[1]
    r26_bsave_mpl[r,2]=test$beta[2]
    r26_bsave_mpl[r,3]=test$se$se_H[1]
    r26_bsave_mpl[r,4]=test$se$se_H[2]
    
    r26_gsave_mpl[r,1]=test$gamma[1]
    r26_gsave_mpl[r,2]=test$se$se_H[3]
    print(c(r, "mpl 26A",test$se$se_H[1]))
    
  } else if(class(try2)!="try-error"){
    smc1=try2
    r26_bsave_em[r,1]=smc1$b[1]
    r26_bsave_em[r,2]=smc1$b[2]
    r26_bsave_em[r,3]=smc1$b_sd[1]
    r26_bsave_em[r,4]=smc1$b_sd[2]
    
    r26_gsave_em[r,1]=smc1$beta[1]
    r26_gsave_em[r,2]=smc1$beta_sd[1]
    print(c(r, "em 26A"))
  } else{
    print(c(r, "neither 26A"))
  }
}



mean(r26_censorsave[!is.na(r26_bsave_mpl[,3])],1)
mean(r26_censorsave[!is.na(r26_bsave_mpl[,3])],2)
r26_censorsave=data.frame(r26_censorsave)
write.table(r26_censorsave, "r26censor.csv",sep=",",dec=".")

colnames(r26_bsave_mpl)=c("b0","b1","seb0","seb1")
r26_bsave_mpl=data.frame(r26_bsave_mpl)
write.table(r26_bsave_mpl, "r26bsave_mpl.csv",sep=",",dec=".",row.names = F)
colnames(r26_bsave_em)=c("b0","b1","seb0","seb1")
r26_bsave_em=data.frame(r26_bsave_em)
write.table(r26_bsave_em, "r26bsave_em.csv",sep=",",dec=".",row.names = F)
colnames(r26_gsave_mpl)=c("g1","seg1")
r26_gsave_mpl=data.frame(r26_gsave_mpl)
write.table(r26_gsave_mpl, "r26gsave_mpl.csv",sep=",",dec=".",row.names = F)
colnames(r26_gsave_em)=c("g1","seg1")
r26_gsave_em=data.frame(r26_gsave_em)
write.table(r26_gsave_em, "r26gsave_em.csv",sep=",",dec=".",row.names = F)


#mpl
mean(r26_bsave_mpl[!is.na(r26_bsave_mpl[,3]),1])
mean(r26_bsave_mpl[!is.na(r26_bsave_mpl[,3]),2])
mean(r26_bsave_mpl[!is.na(r26_bsave_mpl[,3]),3])
sd(r26_bsave_mpl[!is.na(r26_bsave_mpl[,3]),1])
mean(r26_bsave_mpl[!is.na(r26_bsave_mpl[,3]),4])
sd(r26_bsave_mpl[!is.na(r26_bsave_mpl[,3]),2])

mean(r26_gsave_mpl[!is.na(r26_gsave_mpl[,2]),1])
mean(r26_gsave_mpl[!is.na(r26_gsave_mpl[,2]),2])
sd(r26_gsave_mpl[!is.na(r26_gsave_mpl[,2]),1])

r26_bsave_mpl$b0_LL=r26_bsave_mpl$b0_UL=r26_bsave_mpl$b1_LL=r26_bsave_mpl$b1_UL=rep(0,100)
r26_gsave_mpl$g1_LL=r26_gsave_mpl$g1_UL=rep(0,100)

r26_bsave_mpl$b0_LL[!is.na(r26_bsave_mpl$seb0)]=r26_bsave_mpl$b0[!is.na(r26_bsave_mpl$seb0)]-1.96*r26_bsave_mpl$seb0[!is.na(r26_bsave_mpl$seb0)]
r26_bsave_mpl$b0_UL[!is.na(r26_bsave_mpl$seb0)]=r26_bsave_mpl$b0[!is.na(r26_bsave_mpl$seb0)]+1.96*r26_bsave_mpl$seb0[!is.na(r26_bsave_mpl$seb0)]

r26_bsave_mpl$b1_LL[!is.na(r26_bsave_mpl$seb0)]=r26_bsave_mpl$b1[!is.na(r26_bsave_mpl$seb0)]-1.96*r26_bsave_mpl$seb1[!is.na(r26_bsave_mpl$seb0)]
r26_bsave_mpl$b1_UL[!is.na(r26_bsave_mpl$seb0)]=r26_bsave_mpl$b1[!is.na(r26_bsave_mpl$seb0)]+1.96*r26_bsave_mpl$seb1[!is.na(r26_bsave_mpl$seb0)]

r26_bsave_mpl$b0_cov=r26_bsave_mpl$b1_cov=rep(0,100)
r26_bsave_mpl$b0_cov[!is.na(r26_bsave_mpl$seb0)]=as.numeric(r26_bsave_mpl$b0_LL[!is.na(r26_bsave_mpl$seb0)]<=1 & 1<=r26_bsave_mpl$b0_UL[!is.na(r26_bsave_mpl$seb0)])
r26_bsave_mpl$b1_cov[!is.na(r26_bsave_mpl$seb0)]=as.numeric(r26_bsave_mpl$b1_LL[!is.na(r26_bsave_mpl$seb0)]<=1 & 1<=r26_bsave_mpl$b1_UL[!is.na(r26_bsave_mpl$seb0)])

r26_gsave_mpl$g1_LL[!is.na(r26_gsave_mpl$seg1)]=r26_gsave_mpl$g1[!is.na(r26_gsave_mpl$seg1)]-1.96*r26_gsave_mpl$seg1[!is.na(r26_gsave_mpl$seg1)]
r26_gsave_mpl$g1_UL[!is.na(r26_gsave_mpl$seg1)]=r26_gsave_mpl$g1[!is.na(r26_gsave_mpl$seg1)]+1.96*r26_gsave_mpl$seg1[!is.na(r26_gsave_mpl$seg1)]
r26_gsave_mpl$g1_cov=rep(0,100)
r26_gsave_mpl$g1_cov[!is.na(r26_gsave_mpl$seg1)]=as.numeric(r26_gsave_mpl$g1_LL[!is.na(r26_gsave_mpl$seg1)]<=0.5 & 0.5<=r26_gsave_mpl$g1_UL[!is.na(r26_gsave_mpl$seg1)])


#em

mean(r26_bsave_em[,1])
mean(r26_bsave_em[,2])
mean(r26_bsave_em[,3])
sd(r26_bsave_em[,1])
mean(r26_bsave_em[,4])
sd(r26_bsave_em[,2])

mean(r26_gsave_em[,1])
mean(r26_gsave_em[,2])
sd(r26_gsave_em[,1])

r26_bsave_em$b0_LL=r26_bsave_em$b0_UL=r26_bsave_em$b1_LL=r26_bsave_em$b1_UL=rep(0,100)
r26_gsave_em$g1_LL=r26_gsave_em$g1_UL=rep(0,100)

r26_bsave_em$b0_LL=r26_bsave_em$b0-1.96*r26_bsave_em$seb0
r26_bsave_em$b0_UL=r26_bsave_em$b0+1.96*r26_bsave_em$seb0

r26_bsave_em$b1_LL=r26_bsave_em$b1-1.96*r26_bsave_em$seb1
r26_bsave_em$b1_UL=r26_bsave_em$b1+1.96*r26_bsave_em$seb1

r26_bsave_em$b0_cov=r26_bsave_em$b1_cov=rep(0,100)
r26_bsave_em$b0_cov=as.numeric(r26_bsave_em$b0_LL<=1 & 1<=r26_bsave_em$b0_UL)
r26_bsave_em$b1_cov=as.numeric(r26_bsave_em$b1_LL<=1 & 1<=r26_bsave_em$b1_UL)

r26_gsave_em$g1_LL=r26_gsave_em$g1-1.96*r26_gsave_em$seg1
r26_gsave_em$g1_UL=r26_gsave_em$g1+1.96*r26_gsave_em$seg1
r26_gsave_em$g1_cov=rep(0,100)
r26_gsave_em$g1_cov=as.numeric(r26_gsave_em$g1_LL<=0.5 & 0.5<=r26_gsave_em$g1_UL)

# 27. beta=[-1,1], n=100, censor 4.2 ---------



r27_censorsave=matrix(rep(0,2*100),ncol=2)#whole dataset event proportion, non-cured event proportion
r27_bsave_mpl=matrix(rep(0,4*100),ncol=4)
r27_gsave_mpl=matrix(rep(0,2*100),ncol=2)
r27_bsave_em=matrix(rep(0,4*100),ncol=4)
r27_gsave_em=matrix(rep(0,2*100),ncol=2)


for(r in 1:100){
  sim.data=gen_data_right(100,4.2,c(-1,1),c(0.5),c("b"),c("b"),c("z1","x1"))
  
  sim.data.uncure=sim.data[sim.data[,3]==1,]
  sim.surv1=Surv(time=sim.data.uncure[,1],event=sim.data.uncure[,2])
  sim.surv=Surv(time=sim.data$time,event =sim.data$event)
  r27_censorsave[r,1]=sum(sim.surv[,2])/100
  r27_censorsave[r,2]=sum(sim.surv1[,2])/nrow(sim.surv1)
  
  try1=try(phmc_mpl(sim.surv~sim.data$X1,pi.formula = ~sim.data$Z,data=sim.data,phmc_mpl.control(conv_limit = 1e-4, n.knots=c(1,0),maxIter = c(1,30000,30001))))
  try2=try(smcure(sim.surv~X1, cureform=~Z, data=sim.data, model="ph"))
  if(class(try1)!="try-error" & class(try2)!="try-error"){
    test=try1
    smc1=try2
    
    r27_bsave_mpl[r,1]=test$beta[1]
    r27_bsave_mpl[r,2]=test$beta[2]
    r27_bsave_mpl[r,3]=test$se$se_H[1]
    r27_bsave_mpl[r,4]=test$se$se_H[2]
    
    r27_gsave_mpl[r,1]=test$gamma[1]
    r27_gsave_mpl[r,2]=test$se$se_H[3]
    
    r27_bsave_em[r,1]=smc1$b[1]
    r27_bsave_em[r,2]=smc1$b[2]
    r27_bsave_em[r,3]=smc1$b_sd[1]
    r27_bsave_em[r,4]=smc1$b_sd[2]
    
    r27_gsave_em[r,1]=smc1$beta[1]
    r27_gsave_em[r,2]=smc1$beta_sd[1]
    
    print(c(r,"both 27",test$se$se_H[1]))
  } else if(class(try1)!="try-error"){
    test=try1
    r27_bsave_mpl[r,1]=test$beta[1]
    r27_bsave_mpl[r,2]=test$beta[2]
    r27_bsave_mpl[r,3]=test$se$se_H[1]
    r27_bsave_mpl[r,4]=test$se$se_H[2]
    
    r27_gsave_mpl[r,1]=test$gamma[1]
    r27_gsave_mpl[r,2]=test$se$se_H[3]
    print(c(r, "mpl 27",test$se$se_H[1]))
    
  } else if(class(try2)!="try-error"){
    smc1=try2
    r27_bsave_em[r,1]=smc1$b[1]
    r27_bsave_em[r,2]=smc1$b[2]
    r27_bsave_em[r,3]=smc1$b_sd[1]
    r27_bsave_em[r,4]=smc1$b_sd[2]
    
    r27_gsave_em[r,1]=smc1$beta[1]
    r27_gsave_em[r,2]=smc1$beta_sd[1]
    print(c(r, "em 27"))
  } else{
    print(c(r, "neither 27"))
  }
}


mean(r27_censorsave[!is.na(r27_bsave_mpl[,3])],1)
mean(r27_censorsave[!is.na(r27_bsave_mpl[,3])],2)
r27_censorsave=data.frame(r27_censorsave)
write.table(r27_censorsave, "r27censor.csv",sep=",",dec=".")

colnames(r27_bsave_mpl)=c("b0","b1","seb0","seb1")
r27_bsave_mpl=data.frame(r27_bsave_mpl)
write.table(r27_bsave_mpl, "r27bsave_mpl.csv",sep=",",dec=".",row.names = F)
colnames(r27_bsave_em)=c("b0","b1","seb0","seb1")
r27_bsave_em=data.frame(r27_bsave_em)
write.table(r27_bsave_em, "r27bsave_em.csv",sep=",",dec=".",row.names = F)
colnames(r27_gsave_mpl)=c("g1","seg1")
r27_gsave_mpl=data.frame(r27_gsave_mpl)
write.table(r27_gsave_mpl, "r27gsave_mpl.csv",sep=",",dec=".",row.names = F)
colnames(r27_gsave_em)=c("g1","seg1")
r27_gsave_em=data.frame(r27_gsave_em)
write.table(r27_gsave_em, "r27gsave_em.csv",sep=",",dec=".",row.names = F)


#mpl
mean(r27_bsave_mpl[!is.na(r27_bsave_mpl[,3]),1])
mean(r27_bsave_mpl[!is.na(r27_bsave_mpl[,3]),2])
mean(r27_bsave_mpl[!is.na(r27_bsave_mpl[,3]),3])
sd(r27_bsave_mpl[!is.na(r27_bsave_mpl[,3]),1])
mean(r27_bsave_mpl[!is.na(r27_bsave_mpl[,3]),4])
sd(r27_bsave_mpl[!is.na(r27_bsave_mpl[,3]),2])

mean(r27_gsave_mpl[!is.na(r27_gsave_mpl[,2]),1])
mean(r27_gsave_mpl[!is.na(r27_gsave_mpl[,2]),2])
sd(r27_gsave_mpl[!is.na(r27_gsave_mpl[,2]),1])

r27_bsave_mpl$b0_LL=r27_bsave_mpl$b0_UL=r27_bsave_mpl$b1_LL=r27_bsave_mpl$b1_UL=rep(0,100)
r27_gsave_mpl$g1_LL=r27_gsave_mpl$g1_UL=rep(0,100)

r27_bsave_mpl$b0_LL[!is.na(r27_bsave_mpl$seb0)]=r27_bsave_mpl$b0[!is.na(r27_bsave_mpl$seb0)]-1.96*r27_bsave_mpl$seb0[!is.na(r27_bsave_mpl$seb0)]
r27_bsave_mpl$b0_UL[!is.na(r27_bsave_mpl$seb0)]=r27_bsave_mpl$b0[!is.na(r27_bsave_mpl$seb0)]+1.96*r27_bsave_mpl$seb0[!is.na(r27_bsave_mpl$seb0)]

r27_bsave_mpl$b1_LL[!is.na(r27_bsave_mpl$seb0)]=r27_bsave_mpl$b1[!is.na(r27_bsave_mpl$seb0)]-1.96*r27_bsave_mpl$seb1[!is.na(r27_bsave_mpl$seb0)]
r27_bsave_mpl$b1_UL[!is.na(r27_bsave_mpl$seb0)]=r27_bsave_mpl$b1[!is.na(r27_bsave_mpl$seb0)]+1.96*r27_bsave_mpl$seb1[!is.na(r27_bsave_mpl$seb0)]

r27_bsave_mpl$b0_cov=r27_bsave_mpl$b1_cov=rep(0,100)
r27_bsave_mpl$b0_cov[!is.na(r27_bsave_mpl$seb0)]=as.numeric(r27_bsave_mpl$b0_LL[!is.na(r27_bsave_mpl$seb0)]<=1 & 1<=r27_bsave_mpl$b0_UL[!is.na(r27_bsave_mpl$seb0)])
r27_bsave_mpl$b1_cov[!is.na(r27_bsave_mpl$seb0)]=as.numeric(r27_bsave_mpl$b1_LL[!is.na(r27_bsave_mpl$seb0)]<=1 & 1<=r27_bsave_mpl$b1_UL[!is.na(r27_bsave_mpl$seb0)])

r27_gsave_mpl$g1_LL[!is.na(r27_gsave_mpl$seg1)]=r27_gsave_mpl$g1[!is.na(r27_gsave_mpl$seg1)]-1.96*r27_gsave_mpl$seg1[!is.na(r27_gsave_mpl$seg1)]
r27_gsave_mpl$g1_UL[!is.na(r27_gsave_mpl$seg1)]=r27_gsave_mpl$g1[!is.na(r27_gsave_mpl$seg1)]+1.96*r27_gsave_mpl$seg1[!is.na(r27_gsave_mpl$seg1)]
r27_gsave_mpl$g1_cov=rep(0,100)
r27_gsave_mpl$g1_cov[!is.na(r27_gsave_mpl$seg1)]=as.numeric(r27_gsave_mpl$g1_LL[!is.na(r27_gsave_mpl$seg1)]<=0.5 & 0.5<=r27_gsave_mpl$g1_UL[!is.na(r27_gsave_mpl$seg1)])


#em

mean(r27_bsave_em[,1])
mean(r27_bsave_em[,2])
mean(r27_bsave_em[,3])
sd(r27_bsave_em[,1])
mean(r27_bsave_em[,4])
sd(r27_bsave_em[,2])

mean(r27_gsave_em[,1])
mean(r27_gsave_em[,2])
sd(r27_gsave_em[,1])

r27_bsave_em$b0_LL=r27_bsave_em$b0_UL=r27_bsave_em$b1_LL=r27_bsave_em$b1_UL=rep(0,100)
r27_gsave_em$g1_LL=r27_gsave_em$g1_UL=rep(0,100)

r27_bsave_em$b0_LL=r27_bsave_em$b0-1.96*r27_bsave_em$seb0
r27_bsave_em$b0_UL=r27_bsave_em$b0+1.96*r27_bsave_em$seb0

r27_bsave_em$b1_LL=r27_bsave_em$b1-1.96*r27_bsave_em$seb1
r27_bsave_em$b1_UL=r27_bsave_em$b1+1.96*r27_bsave_em$seb1

r27_bsave_em$b0_cov=r27_bsave_em$b1_cov=rep(0,100)
r27_bsave_em$b0_cov=as.numeric(r27_bsave_em$b0_LL<=1 & 1<=r27_bsave_em$b0_UL)
r27_bsave_em$b1_cov=as.numeric(r27_bsave_em$b1_LL<=1 & 1<=r27_bsave_em$b1_UL)

r27_gsave_em$g1_LL=r27_gsave_em$g1-1.96*r27_gsave_em$seg1
r27_gsave_em$g1_UL=r27_gsave_em$g1+1.96*r27_gsave_em$seg1
r27_gsave_em$g1_cov=rep(0,100)
r27_gsave_em$g1_cov=as.numeric(r27_gsave_em$g1_LL<=0.5 & 0.5<=r27_gsave_em$g1_UL)





