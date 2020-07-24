##KNOTS##
knots_phmc = function(events,basis="msplines",n.knots=c(8,2), range.quant=c(0.075,0.9),order=3){
  n.events=length(events)
  range=range(events)
  if(n.knots[2]==0){
    Alpha=quantile(events,seq(0,1,length.out=(n.knots[1]+2)))
  }else{
    Alpha1=quantile(events,seq(0,range.quant[2],length.out=(n.knots[1]+1)))
    Alpha2=seq(quantile(events,range.quant[2]),range(events)[2],length=n.knots[2]+2)
    Alpha=c(Alpha1,Alpha2[-1])
  }
  n.Alpha  = length(Alpha)
  if(basis=="msplines"){
    m = n.Alpha+order-2
    list(m=m, Alpha=Alpha, Delta=rep(1,m))
  }else{
    stop("Choose msplines basis.")
  }
}


