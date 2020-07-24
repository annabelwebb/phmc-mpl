#using aml data from survival package
head(aml)

aml.surv=Surv(aml$time, event=aml$status)

coxph(aml.surv~aml$x,data=aml)
coxph_mpl(aml.surv~as.vector(aml$x), data=aml)
mpl.s=predict(coxph_mpl(aml.surv~aml$x, data=aml),type="survival",time=aml$time)
plot(coxph_mpl(aml.surv~aml$x, data=aml))


flexsurvreg(aml.surv~aml$x, data=aml, dist="exp")
flexsurvreg(aml.surv~aml$x, data=aml, dist="gompertz")

summary(flexsurvreg(aml.surv~aml$x, data=aml, dist="exp"),type="survival")
summary(flexsurvreg(aml.surv~aml$x, data=aml, dist="gompertz"),type="survival")



#generate interval censored data (use in non-cure model just to see....)
sim.data=gen_data(500,0.9,1.3,ev=0.5,c(0,1),c(1),c("u"),c("u"),c("z1","z1"))
for(i in 1:500){
  if(sim.data[i,1]==0 & sim.data[i,2]!=Inf){
    sim.data[i,1]=NA
  }
}
sim.surv=Surv(time=sim.data[,1],time2=sim.data[,2],type="interval2")
test=phmc_mpl(sim.surv~sim.data$X,pi.formula = ~sim.data$Z,data=sim.data,phmc_mpl.control(n.knots=c(1,0),maxIter = c(1,10000,10001)))
summary(test)

flexsurvreg(Surv(time=sim.data[,1],time2=sim.data[,2],type="interval2")~sim.data$X,data=sim.data,dist="gompertz")
flexsurvcure(sim.surv~sim.data$X, data=sim.data, dist="gompertz")

