#' Plot a transmission tree
#' @param ttree Transmission tree
#' @param w.shape Shape parameter of the Gamma probability density function representing the generation length w
#' @param w.scape Scale parameter of the Gamma probability density function representing the generation length w 
plotTTree = function(ttree,w.shape,w.scale) {
  n=nrow(ttree)
  #Determine ys 
  ys <- rep(0, n)
  scale <- rep(1,n)
  todo=c(which(ttree[,3]==0))
  while (length(todo)>0) {
    f=which(ttree[,3]==todo[1])
    o=rank(-ttree[f,1])
    f[o]=f
    for (i in f) {ys[i]=ys[todo[1]]+scale[todo[1]]*which(f==i)/(length(f)+1);scale[i]=scale[todo[1]]/(length(f)+1);todo=c(todo,i)}
    todo=todo[-1]
  }
  ys=rank(ys)
  par(yaxt='n',bty='n')
  mi=min(ttree[,1])
  ma=max(ttree[which(!is.na(ttree[,2])),2])
  xstep=(ma-mi)/2000
  plot(c(),c(),xlim=c(mi-(ma-mi)*0.05,ma+(ma-mi)*0.05),ylim=c(0,n+1),xlab='',ylab='')
  maxcol=max(dgamma(seq(0,ma-mi,xstep),shape=w.shape,scale=w.scale))
  for (i in 1:n) {
    as=seq(ttree[i,1],ma,xstep)
    bs=rep(ys[i],length(as))
    cs=abs((maxcol-dgamma(as-ttree[i,1],shape=w.shape,scale=w.scale))/maxcol)
    cs=gray(cs)
    segments(as,bs,x1=as+xstep,col=cs)
    points(ttree[i,2],ys[i],col = 'red') 
    text(ma+(ma-mi)*0.05,ys[i],i)
    if (ttree[i,3]==0) {next}
    arrows(ttree[i,1],ys[ttree[i,3]],ttree[i,1],ys[i],length=0.1)
  }
}