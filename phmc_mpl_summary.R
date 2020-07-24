summary.phmc_mpl=function(obj,full=FALSE){
  col.names = c("Estimate", "Std. Error", "z-value", "Pr(>|z|)")
  p=obj$dimensions$p
  q=obj$dimensions$q
  m=obj$dimensions$m
  seB=obj$se$se_H[1:p]
  seG=obj$se$se_H[(p+1):(p+q)]
  seT=obj$se$se_H[(p+q+1):(p+q+m)]
  matxB = cbind(obj$beta,seB,obj$beta/seB,2*(1-pnorm(abs(obj$beta/seB))))
  colnames(matxB)=col.names
  matxG = cbind(obj$gamma,seG,obj$gamma/seG,2*(1-pnorm(abs(obj$gamma/seG))))
  colnames(matxG)=col.names
  matxT = cbind(obj$theta,seT,obj$theta/seT,2*(1-pnorm(abs(obj$theta/seT))))
  colnames(matxT)=col.names
  out=list(Beta=matxB,Gamma=matxG,Theta=matxT,inf=list(call=obj$call,full=full, ploglik=obj$ploglik,maxIter=obj$maxIter))
  class(out) = "summary.phmc_mpl"
  out

}

print.phmc_mpl=function(x){
  #cat("\n")
  #print(x$call)
  cat("\nLog-likelihood : ",x$ploglik[1],"\n",sep="")    
  cat("\nLogistic regression parameters :\n")
  vect=c(x$Beta)
  #names(vect)=dimnames(x$data$X)[[2]]
  print(vect)
  cat("\nProportional hazards regression parameters :\n")
  vectg=c(x$Gamma)
  #names(vectg)
  print(vectg)
  cat("\nBaseline hazard parameters : \n")
  print(x$Theta)
  cat("\n")
}

print.summary.phmc_mpl=function(x,...){
  inf = x$inf
  cat("\n")
  print(inf$call)
  cat("\n-----\n\n")

  cat("Mixture Cure Proportional Hazards Model Fitted Using MPL","\n\n\n")
  cat("Penalized log-likelihood  :  ",inf$ploglik,"\n\n",sep="")
  cat(ifelse(inf$maxIter[1]==1,
             "Fixed smoothing value     :  ",
             "Estimated smoothing value :  "),
      inf$smooth,"\n",sep="")
  
  cat("Data             : ",inf$data,"\n",sep="")  
  
  cat("Logistic regression parameters : \n",sep="")
  printCoefmat(x$Beta, P.values=TRUE, has.Pvalue=TRUE)
  
  cat("\nProportional hazards regression parameters : \n",sep="")
  printCoefmat(x$Gamma, P.values=TRUE, has.Pvalue=TRUE) 
  
  if(inf$full){
    cat("\nBaseline hazard parameter vector : \n",sep="")
    #print(x$Theta)
    printCoefmat(x$Theta, P.values=TRUE, has.Pvalue=TRUE)
    #print theta and inference using printCoefmat, print information about active constraints?
  }
  
}




