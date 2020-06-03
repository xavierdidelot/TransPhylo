#' Plot a transmission tree in a detailed format
#' @param ttree Transmission tree
#' @param w.shape Shape parameter of the Gamma probability density function representing the generation time
#' @param w.scale Scale parameter of the Gamma probability density function representing the generation time 
#' @param showLabels Whether or not to show the labels 
#' @param maxTime Maximum value of time to show on x axis
#' @param cex Expansion factor
#' @return Returns invisibly the first parameter
#' @examples 
#' plotTTree(extractTTree(simulateOutbreak()),2,1)
#' @export
plotTTree = function(ttree,w.shape,w.scale,showLabels=TRUE,maxTime=NA,cex=1) {
  nam=ttree$nam
  ttree=ttree$ttree
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
  oldpar <- par('yaxt','bty')
  on.exit(par(oldpar))
  par(yaxt='n',bty='n')
  mi=min(ttree[,1])
  if (is.na(maxTime)) ma=max(ttree[which(!is.na(ttree[,2])),2]) else ma=maxTime
  xstep=(ma-mi)/2000
  plot(c(),c(),xlim=c(mi-(ma-mi)*0.05,ma+(ma-mi)*0.05),ylim=c(0,n+1),xlab='',ylab='')
  maxcol=max(dgamma(seq(0,ma-mi,xstep),shape=w.shape,scale=w.scale))
  for (i in 1:n) {
    as=seq(ttree[i,1],ma,xstep)
    bs=rep(ys[i],length(as))
    cs=abs((maxcol-dgamma(as-ttree[i,1],shape=w.shape,scale=w.scale))/maxcol)
    cs=grDevices::gray(cs)
    segments(as,bs,x1=as+xstep,col=cs)
    points(ttree[i,2],ys[i],col = 'red') 
    if (showLabels) text(ma+(ma-mi)*0.05,ys[i],nam[i],cex=cex)
    if (ttree[i,3]==0) {next}
    arrows(ttree[i,1],ys[ttree[i,3]],ttree[i,1],ys[i],length=0.1)
  }
  return(invisible(ttree))
}
