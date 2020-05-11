#' Calculates the log-probability of a transmission tree
#' @param ttree Transmission tree
#' @param off.r First parameter of the negative binomial distribution for offspring number
#' @param off.p Second parameter of the negative binomial distribution for offspring number
#' @param pi probability of sampling an infected individual
#' @param w.shape Shape parameter of the Gamma probability density function representing the generation time
#' @param w.scale Scale parameter of the Gamma probability density function representing the generation time 
#' @param ws.shape Shape parameter of the Gamma probability density function representing the sampling time
#' @param ws.scale Scale parameter of the Gamma probability density function representing the sampling time 
#' @param dateT Date when process stops (this can be Inf for fully simulated outbreaks)
#' @return Probability of the transmission tree
probTTreeR = function(ttree,off.r,off.p,pi,w.shape,w.scale,ws.shape,ws.scale,dateT)  {
  if (is.list(ttree)) ttree=ttree$ttree
  prob <- 0 
  n <- nrow(ttree)
  
  if (dateT==Inf) {
    # This is the case of a finished outbreak
    omegaStar <- uniroot(function(x) {x-(1-pi)*((1-off.p)/(1-off.p*x))^off.r},c(0,1))$root #This is Equation (2)
    alphaStar <- rep(NA, n+1)
    for (i in (1:n)) { 
      if (is.na(ttree[i,2])) prob<-prob+log(1-pi) #This is the first term in the product in Equation (7)
      else prob<-prob+log(pi)+dgamma((ttree[i,2]-ttree[i,1]),shape=ws.shape,scale=ws.scale,log=TRUE) #This is the second term in the product in Equation (7)
      offspring <- which( ttree[ ,3] == i ) 
      d <- length(offspring)
      if (is.na(alphaStar[d+1])) {
        notinf=max(d,1000)
        alphaStar[d+1]=sum(choose(d:notinf,d)*dnbinom(d:notinf,off.r,1-off.p)*omegaStar^{0:(notinf-d)})#This is Equation (4)
        alphaStar[d+1]=log(alphaStar[d+1])
      }
      prob <- prob + alphaStar[d+1] #This is the third term in the product in Equation (7)
      for (j in offspring) {
        prob <- prob + dgamma((ttree[j,1]-ttree[i,1]),shape=w.shape,scale=w.scale,log=TRUE) #This is the fourth term in the product in Equation (7)
      } 
    } 
    
  } else {
    #This is the case of an ongoing outbreak
    dt=0.01;L=round((dateT-min(ttree[,1]))/dt)
    omegabar=getOmegabarR(L,dt,off.r,off.p,pi,w.shape,w.scale,ws.shape,ws.scale)
    #pit      =function(t) {pi*pgamma((dateT-t),shape=w.shape,scale=w.scale) }#This is Equation (6), but replaced with pi*trunc
    #fomega   =function(x) {omega   [max(1,min(L,round((dateT-x)/dt)))] }
    fomegabar=function(x) {omegabar[max(1,round((dateT-x)/dt))] }
    for (i in (1:n)) { 
      tinf=ttree[i,1]
      ltruncW =pgamma(dateT-tinf,shape= w.shape,scale= w.scale,log.p=T)
       truncWS=pgamma(dateT-tinf,shape=ws.shape,scale=ws.scale)
      if (is.na(ttree[i,2])) prob<-prob+log(1-pi*truncWS) #This is the first term in the product in Equation (11)
      else prob<-prob+log(pi)+dgamma((ttree[i,2]-tinf),shape=ws.shape,scale=ws.scale,log=TRUE) #This is the second term in the product in Equation (11) Note simplification of truncWS/truncWS
      offspring <- which(ttree[ ,3]==i)
      d <- length(offspring)
      notinf=max(1000,d)
      alpha=sum(choose(d:notinf,d)*dnbinom(d:notinf,off.r,1-off.p)*fomegabar(tinf)^{0:(notinf-d)}) #This is Equation (10)
      prob <- prob + log(alpha) #This is the third term in the product in Equation (11)
      prob <- prob + sum(dgamma((ttree[offspring,1]-tinf),shape=w.shape,scale=w.scale,log=TRUE)-ltruncW)#This is the fourth term in the product in Equation (11)
    } 
  }
  return(prob)
} 

getOmegabarR=function(L,dt,off.r,off.p,pi,w.shape,w.scale,ws.shape,ws.scale) {
  omega=rep(NA,L);omega[1]=1;omegabar=rep(NA,L);omegabar[1]=1
  dgammastore=dgamma(dt*(1:(L-1)),shape=w.shape,scale=w.scale)
  coef=c(0.5,rep(1,L-1))
  #omegaStar <- uniroot(function(x) {x-(1-pi)*((1-off.p)/(1-off.p*x))^off.r},c(0,1))$root #This is Equation (2)
  for (k in 1:(L-1)) {
    omegabar[k+1]=min(1,sum(coef[1:k]*dgammastore[seq(k,1,-1)]*omega[1:k]*dt)+1-pgamma(k*dt,shape=w.shape,scale=w.scale))
    omega[k+1]=(1-pi*pgamma(k*dt,shape=ws.shape,scale=ws.scale))*((1-off.p)/(1-off.p*omegabar[k+1]))^off.r #This is Equation (S3)
    #if (abs(omega[k+1]-omegaStar)<0.0001) {omegabar[(k+2):L]=omegaStar;omega[(k+2):L]=omegaStar;break} 
    #if (k==L-1) {warning('Convergence not reached in getOmegabarR')
    #  warning(sprintf('omegaStar=%f,omegabar[L]=%f,w.shape=%f,w.scale=%f,pi=%f,off.r=%f,off.p=%f',omegaStar,omegabar[L],w.shape,w.scale,pi,off.r,off.p))}
  }
  return(omegabar)
}

wstar_rootFinderR=function(pi,off.p,off.r) {
  omegaStar <- uniroot(function(x) {x-(1-pi)*((1-off.p)/(1-off.p*x))^off.r},c(0,1))$root #This is Equation (2)
}