#' Plot both phylogenetic and transmission trees using colors on the phylogeny
#' @param tree Combined phylogenetic/transmission tree
#' @examples
#' plotBothTree(simulateOutbreak())
plotBothTree = function(tree)  {
  nsam <- sum(tree[ ,2]+tree[ ,3] == 0) 
  nh <- nrow(tree)-3*nsam+1
  ntot <- nsam+nh
  par(yaxt='n',bty='n')
  plot(0,0,type='l',xlim=c(min(tree[,1]),max(tree[,1])),ylim=c(0,nsam+1),xlab='',ylab='')
  host <- tree[ ,4] 
  palette(rainbow(ntot))#Need as many unique colors as there are hosts
  
  #Determine ys for leaves
  root<-which(host==0)
  ys <- matrix(0, nsam, 1) 
  todo <- cbind(root,0,0.5,1);#Matrix of nodes to do,with associated starting x and y coordinates and scale 
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
    todo <- rbind(todo[-1,])
  } 
  MySort <- sort(ys,index.return = TRUE); ys <- MySort$ix 
  ys[ys] <- 1:length(ys) 
  
  #Determine ys for non-leaves
  for (i in ((nsam+1):nrow(tree))) { 
    children <- c()
    todo <- i
    while (length(todo)>0) {
      children=c(children,todo[1])
      todo=c(todo[-1],setdiff(tree[todo[1],2:3],0))
    }
    ys[i] <- mean(ys[children[which(children<=nsam)]])
  } 
  
  todo <- cbind(root,tree[root,1]);
  while (nrow(todo) > 0)  { 
    w <- todo[1,1] 
    x <- todo[1,2] 
    y <- ys[w] 
    col=host[w]
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
