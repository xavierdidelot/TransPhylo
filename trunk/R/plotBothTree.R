#' Plot both phylogenetic and transmission trees using colors on the phylogeny
#' @param tree Combined phylogenetic/transmission tree
plotBothTree = function(tree)  {
  n <- sum(tree[ ,2]+tree[ ,3] == 0) 
  plot.new()
  method <- 1 
  par(yaxt='n',bty='n')
  plot(0,0,xlim=c(min(tree[,1]),max(tree[,1])),ylim=c(0,n+1),xlab='',ylab='')
  lines(c(10,11),c(10,11))
  lines(c(10,11),c(11,10))
  host <- tree[ ,4] 
  palette(rainbow(n))#Need as many unique colors as there are hosts
  
  #Determine ys 
  ys <- matrix(0, n, 1) 
  todo <- cbind(nrow(tree),0,0.5,1);#Matrix of nodes to do,with associated starting x and y coordinates and scale 
  while (nrow(todo) > 0)  { 
    w <- todo[1,1] 
    x <- todo[1,2] 
    y <- todo[1,3] 
    scale <- todo[1,4] 
    if (tree[w,2] == 0 && tree[w,3] == 0)  { 
      #Leaf node 
      ys[w] <- y 
    } else if (tree[w,3] == 0)  { 
      #Transmission node 
      todo <- rbind(todo,cbind(tree[w,2],tree[w,1],y,scale,deparse.level=0)) 
    } else { 
      #Binary node 
      todo <- rbind(todo,cbind(tree[w,2],tree[w,1],y + scale/2,scale/2,deparse.level=0),cbind(tree[w,3],tree[w,1],y-scale/2,scale/2,deparse.level=0))
    } 
    todo <- rbind(todo[-1,]);
  } 
  MySort <- sort(ys,index.return = TRUE); ys <- MySort$ix 
  ys[ys] <- 1:length(ys) 
  ylims <- cbind(ys,ys)
  ylims <- rbind(ylims,matrix(0,nrow(tree)-n,2))
  for (i in ((n+1):nrow(tree))) { 
    f <- 1 + which( tree[i,2:3] > 0 ) 
    ylims[i,1] <- min(ylims[tree[i,f],1])
    ylims[i,2] <- max(ylims[tree[i,f],2])
    ys[i] <- mean(ys[tree[i,f]])
  } 
  
  todo <- cbind(nrow(tree),tree[nrow(tree),1]);
  while (nrow(todo) > 0)  { 
    w <- todo[1,1] 
    x <- todo[1,2] 
    y <- ys[w] 
    col=host[w];
    if (tree[w,2] == 0 && tree[w,3] == 0)  { 
      #Leaf node 
      lines(c(x,tree[w,1]),c(y,y),col=col,lwd=2) 
      text(tree[w,1] + (max(cbind(tree[ ,1]))-min(cbind(tree[ ,1])))/100,y,w)
    } else if (tree[w,3] == 0)  { 
      #Transmission node 
      lines(c(x,tree[w,1]),c(y,y),col=col,lwd=2) 
      points(tree[w,1],y,col = 'red',pch=8) 
      todo <- rbind(todo,cbind(tree[w,2],tree[w,1])) 
    } else { 
      #Binary node 
      lines(c(x,tree[w,1]),c(y,y),col=col,lwd=2)
      lines(c(tree[w,1],tree[w,1]),cbind(ys[tree[w,2]],ys[tree[w,3]]),col=col,lwd=2)
      todo <- rbind(todo,cbind(tree[w,2],tree[w,1]),cbind(tree[w,3],tree[w,1])) 
    } 
    todo <- rbind(todo[-1,]);
  } 
} 
