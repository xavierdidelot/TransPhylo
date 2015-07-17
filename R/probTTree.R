#' Calculates the log-probability of a transmission tree
#' @param ttree Transmission tree
#' @param off.r First parameter of the negative binomial distribution for offspring number
#' @param off.p Second parameter of the negative binomial distribution for offspring number
#' @param pi probability of sampling an infected individual
#' @param w.shape Shape parameter of the Gamma probability density function representing the generation length w
#' @param w.scale Scale parameter of the Gamma probability density function representing the generation length w 
#' @param T Date when process stops (this can be Inf for fully simulated outbreaks)
#' @return Probability of the transmission tree
probTTree = function(ttree,off.r,off.p,pi,w.shape,w.scale,T)  {
  prob <- 0 
  n <- nrow(ttree)
  
  if (T==Inf) {
    # This is the case of a finished outbreak
    omegaStar <- uniroot(function(x) {x-(1-pi)*((1-off.p)/(1-off.p*x))^off.r},c(0,1))$root #This is Equation (1)
    alphaStar <- rep(NA, n+1)
    for (i in (1:n)) { 
      if (is.na(ttree[i,2])) prob<-prob+log(1-pi) #This is the first term in the product in Equation (5)
      else prob<-prob+log(pi)+dgamma((ttree[i,2]-ttree[i,1]),shape=w.shape,scale=w.scale,log=TRUE) #This is the second term in the product in Equation (5)
      offspring <- which( ttree[ ,3] == i ) 
      d <- length(offspring)
      if (is.na(alphaStar[d+1])) {
        alphaStar[d+1]=0
        notinf=d+2*off.r*off.p/(1-off.p)
        for (k in d:notinf) alphaStar[d+1]=alphaStar[d+1]+choose(k,d)*dnbinom(k,off.r,off.p)*omegaStar^{k-d} #This is in Equation (2)
        alphaStar[d+1]=log(alphaStar[d+1])
      }
      prob <- prob + alphaStar[d+1] #This is the third term in the product in Equation (5)
      for (j in offspring) {
        prob <- prob + dgamma((ttree[j,1]-ttree[i,1]),shape=w.shape,scale=w.scale,log=TRUE) #This is the fourth term in the product in Equation (5)
      } 
    } 
    
  } else {
    #This is the case of an ongoing outbreak
    dt=0.5;L=1000
    omega=rep(NA,L);omega[1]=1;omegabar=rep(NA,L);omegabar[1]=1
    dgammastore=dgamma(dt*(1:(L-1)),shape=w.shape,scale=w.scale)
    coef=c(0.5,rep(1,L-1))
    omegaStar <- uniroot(function(x) {x-(1-pi)*((1-off.p)/(1-off.p*x))^off.r},c(0,1))$root #This is Equation (1)
    for (k in 1:(L-1)) {
      omegabar[k+1]=sum(coef[1:k]*dgammastore[seq(k,1,-1)]*omega[1:k]*dt)+1-pgamma(k*dt,shape=w.shape,scale=w.scale)
      omega[k+1]=(1-pi*pgamma(k*dt,shape=w.shape,scale=w.scale))*((1-off.p)/(1-off.p*omegabar[k+1]))^off.r
      if (abs(omegabar[k+1]-omegaStar)<0.01) {omegabar[k+2:L]=omegaStar;omega[k+2:L]=omegaStar;break} 
      if (k==L-1) {warning('Convergence not reached')
        print(sprintf('omegaStar=%f,omegabar[L]=%f,w.shape=%f,w.scale=%f,pi=%f,off.r=%f,off.p=%f',omegaStar,omegabar[L],w.shape,w.scale,pi,off.r,off.p))}
    }
    #pit      =function(t) {pi*pgamma((T-t),shape=w.shape,scale=w.scale) }#This is Equation (6), but replaced with pi*trunc
    fomega   =function(x) {omega   [max(1,min(L,round((T-x)/dt)))] }
    fomegabar=function(x) {omegabar[max(1,min(L,round((T-x)/dt)))] }
    
    for (i in (1:n)) { 
      tinf=ttree[i,1]
      trunc=pgamma(T-tinf,shape=w.shape,scale=w.scale)
      ltrunc=log(trunc)
      if (is.na(ttree[i,2])) prob<-prob+log(1-pi*trunc) #This is the first term in the product in Equation (9)
      else prob<-prob+log(pi*trunc)+dgamma((ttree[i,2]-ttree[i,1]),shape=w.shape,scale=w.scale,log=TRUE)-ltrunc #This is the second term in the product in Equation (9)
      offspring <- which( ttree[ ,3] == i ) 
      d <- length(offspring)
      alpha=0
      notinf=d+2*off.r*off.p/(1-off.p)
      for (k in d:notinf) alpha=alpha+choose(k,d)*dnbinom(k,off.r,off.p)*fomegabar(tinf)^{k-d} #This is in Equation (8)
      alpha=log(alpha)
      prob <- prob + alpha #This is the third term in the product in Equation (9)
      for (j in offspring) {
        prob <- prob + dgamma((ttree[j,1]-ttree[i,1]),shape=w.shape,scale=w.scale,log=TRUE)-ltrunc #This is the fourth term in the product in Equation (9)
      } 
    } 
  }
  
  return(prob)
} 
