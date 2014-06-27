set.seed(1)
simu <- simulateOutbreak(R=1,neg=100/365,pi=0.5)
record<-inferTTree(ptreeFromFullTree(simu),mcmcIterations=10000,startPi=0.5,updatePi=FALSE,testing=TRUE)
#par(mfrow=c(2,2))
#plot(sapply(record,function(x) x$pTTree+x$pPTree),ylab='Posterior probability',xlab='MCMC iterations',type='l')
#hist(sapply(record,function(x) x$pi),xlab='Posterior of pi',main='Prior is Unif(0,1)')
#hist(sapply(record,function(x) x$neg),xlab='Posterior of Ne*g',main='Prior is Exp(1)')
#hist(sapply(record,function(x) x$R),xlab='Posterior of R',main='Prior is Exp(1)')

ncases=c()
for (i in 5000:10000) {
  ptree=ttreeFromFullTree(record[[i]]$tree)
  ncases=c(ncases,nrow(ptree))
}
par(mfrow=c(2,1))
hist(ncases,breaks=0.5:50.5,main='Inferred number of cases')
plot(dbinom(10,1:50,0.5),type='h',main='Theoretical distribution for number of cases')
#TODO: the inferred and theoretical distributions should match

#TODO: check that in the testing conditions, transmission events are distributed uniformly on the ptree