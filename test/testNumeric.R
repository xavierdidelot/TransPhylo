rm(list=ls())
w.shape=10;w.scale=0.1
pi=0.01#0.01
T=2010
off.r=9.574812;off.p=0.5

dt=0.05;L=1000
omega=rep(NA,L);omega[1]=1;omegabar=rep(NA,L);omegabar[1]=1
dgammastore=dgamma(dt*(1:(L-1)),shape=w.shape,scale=w.scale)
coef=c(0.5,rep(1,L-1))
omegaStar <- uniroot(function(x) {x-(1-pi)*((1-off.p)/(1-off.p*x))^off.r},c(0,1))$root #This is Equation (1)
for (k in 1:(L-1)) {
  omegabar[k+1]=min(1,sum(coef[1:k]*dgammastore[seq(k,1,-1)]*omega[1:k]*dt)+1-pgamma(k*dt,shape=w.shape,scale=w.scale))
  omega[k+1]=(1-pi*pgamma(k*dt,shape=w.shape,scale=w.scale))*((1-off.p)/(1-off.p*omegabar[k+1]))^off.r
  if (abs(omegabar[k+1]-omegaStar)<0.01) {omegabar[k+2:L]=omegaStar;omega[k+2:L]=omegaStar;break} 
  if (k==(L-1)) warning('Convergence not reached')
}
fomega   =function(x) {omega   [max(1,min(L,round((T-x)/dt)))] }
fomegabar=function(x) {omegabar[max(1,min(L,round((T-x)/dt)))] }

xs=seq(T-50,T+10,0.01);plot(xs,lapply(xs,fomega),col='blue',type='l',xlab='',ylab='',ylim=c(0,1))
lines(xs,lapply(xs,fomegabar),col='red',type='l')
points(T-50,omegaStar)
legend('topleft',legend=c('omega','omegabar'),col=c('blue','red'),lty=1)


