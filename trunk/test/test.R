set.seed(1)
simu <- simulateOutbreak(R=1,neg=100/365,pi=0.5)
#simu[,1]=simu[,1]*10
pi=0.5
record<-inferTTree(ptreeFromFullTree(simu),mcmcIterations=10000,startPi=pi,updatePi=FALSE,testing=TRUE)
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
ncases2=ncases;ncases2[which(ncases2>100)]=100
hist(ncases2,breaks=0.5:100.5,main='Inferred number of cases')
plot(dbinom(10,1:100,pi),type='h',main='Theoretical distribution for number of cases')
#The inferred and theoretical distributions should match!

#TODO: check that in the testing conditions, transmission events are distributed uniformly on the ptree