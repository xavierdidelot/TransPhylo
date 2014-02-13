makeTTree <-function(R) { 
#Creates a transmission tree, and returns the result as a N*3 matrix in the
#following format:
#One row per infected host
#First column is time of infection
#Second column is time of sampling
#Third column is infector
#R is the basic reproduction number
#Assume w is Gamma(2,1)
ttree<-matrix(0,1,3)
prob<-0
todo<-1
while (length(todo)>0) {
    draw<-rgamma(1,2,1)
    ttree[todo[1],2]<-ttree[todo[1],1]+draw
    prob<-prob+log(dgamma(draw,2,1))
    offspring<-rpois(1,R)
    prob<-prob+log(dpois(offspring,R))
    if (offspring>0) {
    for (i in 1:offspring) {
        draw<-rgamma(1,2,1)
        prob<-prob+log(dgamma(draw,2,1))
        ttree<-rbind(ttree,c(ttree[todo[1],1]+draw,0,todo[1]))
        todo<-c(todo,nrow(ttree))
        if (nrow(ttree)>100) {return(list(NULL,NULL))}
    }
    }
    todo<-todo[-1] 
}
return(list(ttree,prob))
}