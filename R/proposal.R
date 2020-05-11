#Proposal distribution used in rjMCMC to update the transmission tree
#Returns the proposed tree as well as qr=proposal ratio part of the Metropolis-Hastings ratio
#And also old=list of hosts in previous tree for which subtree changed
#And also new=list of hosts in new tree for which subtree changed
proposal = function(tree)  {
  which=sample.int(3,1)
  if (which==1) return(move1(tree))#Add
  if (which==2) return(move2(tree))#Remove
  if (which==3) return(move3(tree))#Update
}

#MOVE1: Add a new transmission event
move1 = function(tree) {
  fathers <- rep(NA, nrow(tree));fathers[tree[ ,2:3] + 1] <- 1:nrow(tree);fathers <- fathers[-1] 
  totbralen=sum(utils::head(tree[,1],-2)-utils::head(tree[fathers,1],-1))
  if (totbralen==0) return(list(tree=tree,qr=1))
  loc=runif(1)*totbralen
  bra=1
  while (loc>(tree[bra,1]-tree[fathers[bra],1])) {
    loc=loc- (tree[bra,1]-tree[fathers[bra],1])
    bra=bra+1}
  nsam <- sum(tree[ ,2] == 0&tree[ ,3] == 0) 
  old=tree[bra,4]
  nage=tree[fathers[bra],1]+loc
  tree=treeAdd(tree,nage,bra,fathers[bra])
  tree <- cbind(tree[ ,1:3],.computeHost(tree)) 
  ntraeve=sum( tree[ ,2] > 0&tree[ ,3] == 0)
  qr=totbralen/(ntraeve-nsam)
  a=which(tree[,1]==nage)#find the new transmission event
  b=tree[a,2]
  return(list(tree=tree,qr=qr,old=old,new=c(tree[a,4],tree[b,4])))
}

#MOVE2: remove a transmission event if possible
move2 = function(tree) {
  nsam <- sum(tree[ ,2] == 0&tree[ ,3] == 0) 
  ntraeve=sum(tree[ ,2]  > 0&tree[ ,3] == 0)
  if (nsam==ntraeve) return(list(tree=tree,qr=1,old=as.numeric(c()),new=as.numeric(c())))#Nothing to remove
  
  #Choose a transmission event that can be removed
  host <- tree[ ,4]
  fathers <- rep(NA, nrow(tree));fathers[tree[ ,2:3] + 1] <- 1:nrow(tree);fathers <- fathers[-1] 
  totbralen=sum(utils::head(tree[,1],-2)-utils::head(tree[fathers,1],-1))
  w=0
  while (w==0||w==nrow(tree)||(infector<=nsam && infected<=nsam)){
    w<-which( tree[ ,2] > 0&tree[ ,3] == 0)
    if (length(w)>1) w <- sample(w ,1)
    tmp=tree[w,2]
    infector <- host[w] 
    infected <- host[tree[w,2]]}
  
  #Remove it
  tree=treeRem(tree,w,fathers[w])
  tree <- cbind(tree[ ,1:3],.computeHost(tree)) 
  qr=(ntraeve-nsam)/totbralen
  return(list(tree=tree,qr=qr,old=c(infector,infected),new=tree[tmp,4]))
}

#MOVE3: update a transmission event
move3 = function(tree) {
  nsam <- sum(tree[ ,2] == 0&tree[ ,3] == 0) 
  host <- tree[ ,4]
  fathers <- rep(NA, nrow(tree));fathers[tree[ ,2:3] + 1] <- 1:nrow(tree);fathers <- fathers[-1] 
  
  #Choose a transmission event
  w<-which( tree[ ,2] > 0&tree[ ,3] == 0 ) 
  if (length(w)>1) w <- sample(w,1)
  infector <- host[w] 
  infected <- host[tree[w,2]]   
  
  if (w == nrow(tree))  { 
    #Update age of transmission to the index case 
    r <- (runif(1)-0.5)/1
    mini <- tree[nrow(tree),1]-tree[nrow(tree)-1,1]
    if (r < mini)  { 
      r <- mini + (mini-r) 
    } 
    #tree[1:(nrow(tree)-1),1] <- tree[1:(nrow(tree)-1),1] + r 
    tree[nrow(tree),1]<-tree[nrow(tree),1]-r
    if (tree[nrow(tree),1]>tree[nrow(tree)-1,1]) warning('Warning with root age')
    return(list(tree=tree,qr=1,old=tree[nrow(tree)-1,4],new=tree[nrow(tree)-1,4]))          
  } 
    
  #Remove the transmission node 
  tree=treeRem(tree,w,fathers[w])
  host=host[-w]
  fathers <- rep(NA, nrow(tree));fathers[tree[ ,2:3] + 1] <- 1:nrow(tree);fathers <- fathers[-1] 
  
  #Choose where to add it back
  islocs <- rep(0, nrow(tree)) #Where is the star allowed to be moved to
  
  if (infector<=nsam && infected<=nsam) {
  #If transmission event is from a sampled case to a sampled case, it has to stay on the path from one leaf to another
  path <- infector 
  islocs[path] <- 1 
  while (path < nrow(tree))  { 
    path <- fathers[path] 
    islocs[path] <- 1 
  } 
  path <- infected 
  islocs[path] <- 1 
  while (path < nrow(tree))  { 
    path <- fathers[path] 
    islocs[path] <- 1-islocs[path] 
  } 
  } else {
  #Otherwise it can move anywhere within infector or infected (except on the root)
    islocs[which(host==infected|host==infector)]<-1
    islocs[length(islocs)-1]<-0
  }
  
  #Choose where to add it back in the set of possible locations
  locs <- which(islocs==1) 
  lens <- rep(0,length(locs)) 
  for (i in 1:length(locs)) { 
    lens[i] <- tree[locs[i],1]-tree[fathers[locs[i]],1] 
  } 
  r <- runif(1) * sum(lens) 
  i=1
  while (r>lens[i]) {r<-r-lens[i];i=i+1}
  
  #Add it back
  nage=tree[locs[i],1]-r
  tree=treeAdd(tree,nage,locs[i],fathers[locs[i]])
  tree <- cbind(tree[ ,1:3],.computeHost(tree)) 
  
  a=which(tree[,1]==nage)#find the new transmission event
  b=tree[a,2]
  return(list(tree=tree,qr=1,old=c(infector,infected),new=c(tree[a,4],tree[b,4])))  
}

#.reordernodes = function(tree) {
#  #Reorder nodes chronologically 
#  nsam <- sum(tree[ ,2] == 0&tree[ ,3] == 0) 
#  MySort <- sort(tree[(nsam+1):nrow(tree),1],decreasing=TRUE,index.return = TRUE); ind <- MySort$ix 
#  for (i in (nsam+1):nrow(tree)) for (j in (2:3)) if (tree[i,j] > nsam) tree[i,j] <- nsam + which( ind == tree[i,j]-nsam ) 
#  tree <- tree[c(1:nsam,nsam + ind), ] 
#  tree <- cbind(tree[ ,1:3],.computeHost(tree)) 
#}

#Add a transmission node onto a tree
treeAdd = function(tree,age,child,father) {
  nsam <- sum(tree[ ,2] == 0&tree[ ,3] == 0) 
  w=nsam+min(which(tree[(nsam+1):nrow(tree),1]<age))
  tree[which(tree[,2]>=w),2]=tree[which(tree[,2]>=w),2]+1
  tree[which(tree[,3]>=w),3]=tree[which(tree[,3]>=w),3]+1
  tree[father,1+which(tree[father,2:3]==child)]=w
  tree=rbind(tree[1:(w-1),],c(age,child,0,0),tree[w:nrow(tree),])
}

#Remove a transmission node from a tree
treeRem = function(tree,w,father) {
  tree[father,1+which(tree[father,2:3]==w)]=tree[w,2]
  tree[which(tree[,2]>w),2]=tree[which(tree[,2]>w),2]-1
  tree[which(tree[,3]>w),3]=tree[which(tree[,3]>w),3]-1
  tree=tree[-w,]
}

