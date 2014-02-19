#' Plot a transmission tree
#' @param ttree Transmission tree
plotTTree = function(ttree) {
  n=nrow(ttree)
  #Determine ys 
  ys <- matrix(0, n, 1) 
  for (i in 1:n) {
    if (ttree[i,3]==0) {next}
    ys[which(ys>ys[ttree[i,3]])]=ys[which(ys>ys[ttree[i,3]])]+1
    ys[i]=ys[ttree[i,3]]+1
  }
  par(yaxt='n',bty='n')
  mi=min(ttree[,1])
  ma=max(ttree[,2])
  plot(0,0,xlim=c(mi,ma+(ma-mi)*0.05),ylim=c(-1,n+1),xlab='',ylab='')
  for (i in 1:n) {
    lines(c(ttree[i,1],ma),c(ys[i],ys[i]))
    points(ttree[i,2],ys[i],col = 'red') 
    text(ma+(ma-mi)*0.05,ys[i],i)
    if (ttree[i,3]==0) {next}
    arrows(ttree[i,1],ys[ttree[i,3]],ttree[i,1],ys[i],length=0.1)
  }
}