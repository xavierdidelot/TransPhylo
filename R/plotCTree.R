#' Plot both phylogenetic and transmission trees using colors on the phylogeny
#' @param tree Combined phylogenetic/transmission tree
#' @param showLabels Whether or not to show the labels 
#' @param cols Colors to use for hosts
#' @param maxTime Maximum time to show on the x axis
#' @examples
#' plotCTree(simulateOutbreak())
#' @export
plotCTree = function(tree,showLabels=TRUE,cols=NA,maxTime=NA)  {
  nam=tree$nam
  tree=tree$ctree
  nsam <- sum(tree[ ,2]+tree[ ,3] == 0) 
  nh <- nrow(tree)-3*nsam+1
  ntot <- nsam+nh
  par(yaxt='n',bty='n')
  plot(0,0,type='l',xlim=c(min(tree[,1]),ifelse(is.na(maxTime),max(tree[,1]),maxTime)),ylim=c(0,nsam+1),xlab='',ylab='')
  host <- tree[ ,4] 
  if (ntot>1) {
      if (is.na(cols[1])) palette(rainbow(min(1024,ntot)))#Need as many unique colors as there are hosts. If there are more than 1024 hosts, colors are recycled.
    else palette(cols)
    }
  
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
  ys<-rank(ys)
  
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
      if (showLabels) text(tree[w,1] + (max(cbind(tree[ ,1]))-min(cbind(tree[ ,1])))/100,y,nam[w])
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
