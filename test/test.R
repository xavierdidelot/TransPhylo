set.seed(0)
simu <- simulateOutbreak(neg=100/365,pi=0.5,ninf=20)

pi=0.5
start <- Sys.time()
record<-inferTTree(ptreeFromFullTree(simu),mcmcIterations=10000,startPi=pi,updatePi=FALSE,testing=TRUE)
print(Sys.time()-start)

#par(mfrow=c(2,2))
#plot(sapply(record,function(x) x$pTTree+x$pPTree),ylab='Posterior probability',xlab='MCMC iterations',type='l')
#hist(sapply(record,function(x) x$pi),xlab='Posterior of pi',main='Prior is Unif(0,1)')
#hist(sapply(record,function(x) x$neg),xlab='Posterior of Ne*g',main='Prior is Exp(1)')
#hist(sapply(record,function(x) x$R),xlab='Posterior of R',main='Prior is Exp(1)')

#Compare inferred and theoretical distributions for the number of cases
ncases=rep(0,length(record)/2)
for (i in (length(record)/2):length(record)) {
  ptree=ttreeFromFullTree(record[[i]]$tree)
  ncases[i-length(record)/2+1]=nrow(ptree)
}
par(mfrow=c(2,1))
ncases2=ncases;ncases2[which(ncases2>100)]=100
hist(ncases2,breaks=0.5:100.5,main='Inferred number of cases',ylab='')
plot(dbinom(length(which(simu[,2]==0&simu[,3]==0)),1:100,pi),type='h',main='Theoretical distribution for number of cases',ylab='')

#Show age of transmission events
# locs=c()
# for (i in (length(record)/2):length(record)) {
#   ftree=record[[i]]$tree  
#   for (j in 1:nrow(ftree)) {
#     if (ftree[j,2]>0 && ftree[j,3]==0) {
#       locs=c(locs,ftree[j,1])
#     }
#   }
# }
# hist(locs)
#   