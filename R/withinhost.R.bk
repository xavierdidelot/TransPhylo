# Simulates the within-host coalescent model
# @param times times at which N samples are taken(counted forward in time from infection time) 
# @param neg is the product of the within-host effective population size and the generation duration in days 
# @return array of size(2N)*3 where each row is a node,the first column indicate the date of the node and the last two columns indicate the two children. This array has internal nodes sorted in order of most recent to most ancient node(and remains so during the algorithm). The last node corresponds to infection time and only has one child 
.withinhost = function(times,neg)  {
  prob <- 0 
  MySort <- sort(times,decreasing=TRUE,index.return = TRUE); tim <- MySort$x; ind <- MySort$ix 
  n <- length(tim) 
  nodes <- cbind(0,ind[1],0);#Start with one node at time 0 and with the first isolate connected to it 
  i <- 2 
  while (i <= n) {#Graft branches one by one 
    r <- -log(runif(1)) * neg 
    curt <- tim[i];#Current time:start with date of isolate and go back in time until coalescence happens 
    fi <- which( nodes[ ,1] < curt ) ;fi<-fi[1]
    for (j in (.seqML(fi,nrow(nodes))))  {
      if (r > (curt-nodes[j,1]) * (i-j))  { 
        r <- r-(curt-nodes[j,1]) * (i-j) 
        curt <- nodes[j,1] 
      } else { 
        curt <- curt-r/(i-j);#Found the time for grafting
        prob <- prob + log(dexp(r,neg^(-1))) 
        r <- 0 
        break 
      } 
    } 
    if (r>0) {next} 
    #Create new node 
    a <- nodes[ ,2:3];a[a >= j + n] <- a[a >= j + n] + 1;nodes[ ,2:3] <- a;#Renumbering according to table insertion in next line 
    nodes <- rbind(nodes[.seqML(1,j-1), ],c(curt,ind[i],0),nodes[.seqML(j,nrow(nodes)),]) 
    #Now choose on which branch to regraft amongst the branches alive at time curt 
    no <- j 
    side <- 2 
    prob <- prob + log(1/(nrow(nodes)-j))
    w <- 1 + floor(runif(1) * (nrow(nodes)-j)) 
    while (w > 0)  { 
      no <- no + side-1 
      side <- 3-side 
      if (nodes[no,side + 1] <= n ||(nodes[no,side + 1] > n && nodes[nodes[no,side + 1]-n,1] > curt))  { 
        w <- w-1 
      } 
    } 
    nodes[j,3] <- nodes[no,side + 1] 
    nodes[no,side + 1] <- n + j 
    i <- i + 1 
  } 
  nodes <- rbind(matrix(0, nrow = n, ncol = 3),nodes) 
  nodes[1:n,1] <- times 
  return(list(nodes <- nodes,prob <- prob))
} 

.seqML <- function(from, to, by=1) {if (from > to) integer(0) else seq.int(from, to, by)}

