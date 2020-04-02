#' Simultaneously infer transmission trees given phylogenetic trees 
#' User can specify any subset of parameters that will be shared by providing a character vector of parameter names to 
#' the argument "share".  
#' @param ptree_lst List of phylogenetic tree
#' @param w.shape Shape parameter of the Gamma probability density function representing the generation time
#' @param w.scale Scale parameter of the Gamma probability density function representing the generation time 
#' @param ws.shape Shape parameter of the Gamma probability density function representing the sampling time
#' @param ws.scale Scale parameter of the Gamma probability density function representing the sampling time 
#' @param mcmcIterations Number of MCMC iterations to run the algorithm for
#' @param thinning MCMC thinning interval between two sampled iterations
#' @param startNeg Starting value of within-host coalescent parameter Ne*g
#' @param startOff.r Starting value of parameter off.r
#' @param startOff.p Starting value of parameter off.p
#' @param startPi Starting value of sampling proportion pi
#' @param prior_pi_a First shape parameter of Beta prior for pi
#' @param prior_pi_b Second shape parameter of Beta prior for pi
#' @param updateNeg Whether of not to update the parameter Ne*g
#' @param updateOff.r Whether or not to update the parameter off.r
#' @param updateOff.p Whether or not to update the parameter off.p
#' @param updatePi Whether or not to update the parameter pi
#' @param share Character vector of parameters to be shared. For example, share = c("off.r", "off.p") would 
#' share the offspring distribution. Allowed parameter names are "neg", "off.r", "off.p" and "pi". 
#' @param startCTree_lst Optional combined list of trees to start from
#' @param updateTTree Whether or not to update the transmission tree
#' @param optiStart Type of optimisation to apply to MCMC start point (0=none, 1=slow, 2=fast)
#' @param dateT Date when process stops (this can be Inf for fully simulated outbreaks)
#' @param verbose Whether or not to use verbose mode (default is false)
#' @return list the same size as input, each element contains posterior transmission trees inferred from
#' corresponding phylogenetic tree
#' @author Yuanwei Xu
#' @export
infer_multittree_share_param = function(ptree_lst,w.shape=2,w.scale=1,ws.shape=w.shape,ws.scale=w.scale,mcmcIterations=1000,
                      thinning=1,startNeg=100/365,startOff.r=1,startOff.p=0.5,startPi=0.5,prior_pi_a=1,prior_pi_b=1,
                      updateNeg=TRUE,updateOff.r=TRUE,updateOff.p=FALSE,updatePi=TRUE,share=NULL,
                      startCTree_lst=rep(NA,length(ptree_lst)),updateTTree=TRUE,optiStart=2,dateT=Inf,verbose=F) {

  if (length(ptree_lst)>=1) for (i in 1:length(ptree_lst)) class(ptree_lst[[i]])<-'list'
  ptree_lst <- purrr::map(ptree_lst, function(x) within(x, ptree[,1] <- ptree[,1]+runif(nrow(ptree))*1e-10))
  #MCMC algorithm
  neg <- startNeg
  off.r <- startOff.r
  off.p <- startOff.p
  pi <- startPi
  
  ctree_lst <- vector("list", length(ptree_lst))
  for(k in seq_along(ptree_lst)){ # starting ctree
    if(is.na(startCTree_lst[[k]]))
      ctree_lst[[k]] <- makeCTreeFromPTree(ptree_lst[[k]],off.r,off.p,neg,pi,w.shape,w.scale,ws.shape,ws.scale,dateT,optiStart)
    else
      ctree_lst[[k]] <- startCTree_lst[[k]]
  }
  ntree <- length(ptree_lst)
  neg_lst <- as.list(rep(neg, ntree))
  off.r_lst <- as.list(rep(off.r, ntree))
  off.p_lst <- as.list(rep(off.p, ntree))
  pi_lst <- as.list(rep(pi, ntree))
  not_share <- setdiff(c("neg", "off.r", "off.p", "pi"), share)

  one_update <- function(ctree, pTTree, pPTree, neg, off.r, off.p, pi, not_share){
    # Get a copy of current ttree
    ttree <- extractTTree(ctree)
    
    if (updateTTree) {
      #Metropolis update for transmission tree 
      if (verbose) message("Proposing ttree update")
      prop <- proposal(ctree$ctree)
      ctree2 <- list(ctree=prop$tree,nam=ctree$nam)
      class(ctree2)<-'ctree'
      ttree2 <- extractTTree(ctree2)
      pTTree2 <- probTTree(ttree2$ttree,off.r,off.p,pi,w.shape,w.scale,ws.shape,ws.scale,dateT) 
      pPTree2 <- probPTreeGivenTTree(ctree2$ctree,neg) 
      if (log(runif(1)) < log(prop$qr)+pTTree2 + pPTree2-pTTree-pPTree)  { 
        ctree <- ctree2 
        ttree <- ttree2
        pTTree <- pTTree2 
        pPTree <- pPTree2 
      } 
    }
    
    if (("neg" %in% not_share) && updateNeg) {
      #Metropolis update for Ne*g, assuming Exp(1) prior 
      neg2 <- abs(neg + (runif(1)-0.5)*0.5)
      if (verbose) message(sprintf("Proposing unshared Ne*g update %f->%f",neg,neg2))
      pPTree2 <- probPTreeGivenTTree(ctree$ctree,neg2) 
      if (log(runif(1)) < pPTree2-pPTree-neg2+neg)  {neg <- neg2;pPTree <- pPTree2} 
    }
    
    if (("off.r" %in% not_share) && updateOff.r) {
      #Metropolis update for off.r, assuming Exp(1) prior 
      off.r2 <- abs(off.r + (runif(1)-0.5)*0.5)
      if (verbose) message(sprintf("Proposing unshared off.r update %f->%f",off.r,off.r2))
      pTTree2 <- probTTree(ttree$ttree,off.r2,off.p,pi,w.shape,w.scale,ws.shape,ws.scale,dateT) 
      if (log(runif(1)) < pTTree2-pTTree-off.r2+off.r)  {off.r <- off.r2;pTTree <- pTTree2}
    }
    
    if (("off.p" %in% not_share) && updateOff.p) {
      #Metropolis update for off.p, assuming Unif(0,1) prior 
      off.p2 <- abs(off.p + (runif(1)-0.5)*0.1)
      if (off.p2>1) off.p2=2-off.p2
      if (verbose) message(sprintf("Proposing unshared off.p update %f->%f",off.p,off.p2))
      pTTree2 <- probTTree(ttree$ttree,off.r,off.p2,pi,w.shape,w.scale,ws.shape,ws.scale,dateT) 
      if (log(runif(1)) < pTTree2-pTTree)  {off.p <- off.p2;pTTree <- pTTree2}
    }
    
    if (("pi" %in% not_share) && updatePi) {
      # use a beta prior
      pi2 <- pi + (runif(1)-0.5)*0.1
      if (pi2<0.01) pi2=0.02-pi2
      if (pi2>1) pi2=2-pi2
      if (verbose) message(sprintf("Proposing unshared pi update %f->%f",pi,pi2))
      pTTree2 <- probTTree(ttree$ttree,off.r,off.p,pi2,w.shape,w.scale,ws.shape,ws.scale,dateT) 
      log_beta_ratio <- (prior_pi_a - 1) * log(pi2) + (prior_pi_b - 1) * log(1 - pi2) -
        (prior_pi_a - 1) * log(pi) - (prior_pi_b - 1) * log(1 - pi)
      if (log(runif(1)) < pTTree2-pTTree + log_beta_ratio)  {pi <- pi2;pTTree <- pTTree2}       
    }
    
    list(ctree=ctree, pTTree=pTTree, pPTree=pPTree, neg=neg, off.r=off.r, off.p=off.p, pi=pi)
  }
  
  one_update_share <- function(ctree_lst, pTTree_lst, pPTree_lst, neg_lst, off.r_lst, off.p_lst, pi_lst, share){
    ttree_lst <- purrr::map(ctree_lst, extractTTree)
    ttree_lst <- purrr::map(ttree_lst, "ttree") # list of ttree matrices

    if(("neg" %in% share) && updateNeg){
      neg <- neg_lst[[1]]
      #Metropolis update for Ne*g, assuming Exp(1) prior 
      neg2 <- abs(neg + (runif(1)-0.5)*0.5)
      if (verbose) message(sprintf("Proposing shared Ne*g update %f->%f",neg,neg2))
      pPTree <- purrr::flatten_dbl(pPTree_lst)
      pPTree2 <- purrr::map_dbl(ctree_lst, function (x,neg) probPTreeGivenTTree(x$ctree,neg), neg = neg2) 
      if (log(runif(1)) < sum(pPTree2)-sum(pPTree)-neg2+neg){
        neg_lst <- as.list(rep(neg2, ntree))
        pPTree_lst <- as.list(pPTree2)
      }
    }
    
    if(("off.r" %in% share) && updateOff.r) {
      off.r <- off.r_lst[[1]]
      #Metropolis update for off.r, assuming Exp(1) prior 
      off.r2 <- abs(off.r + (runif(1)-0.5)*0.5)
      if (verbose) message(sprintf("Proposing shared off.r update %f->%f",off.r,off.r2))
      pTTree <- purrr::flatten_dbl(pTTree_lst)
      pTTree2 <- purrr::pmap_dbl(list(ttree=ttree_lst, pOff=off.p_lst, pi=pi_lst), probTTree, rOff=off.r2,  
                                 shGen=w.shape, scGen=w.scale, shSam=ws.shape, scSam=ws.scale, dateT=dateT)
      if (log(runif(1)) < sum(pTTree2)-sum(pTTree)-off.r2+off.r){
        off.r_lst <- as.list(rep(off.r2, ntree))
        pTTree_lst <- as.list(pTTree2)
      }
    }
    
    if(("off.p" %in% share) && updateOff.p) {
      off.p <- off.p_lst[[1]]
      #Metropolis update for off.p, assuming Unif(0,1) prior 
      off.p2 <- abs(off.p + (runif(1)-0.5)*0.1)
      if (off.p2>1) off.p2=2-off.p2
      if (verbose) message(sprintf("Proposing shared off.p update %f->%f",off.p,off.p2))
      pTTree <- purrr::flatten_dbl(pTTree_lst)
      pTTree2 <- purrr::pmap_dbl(list(ttree=ttree_lst, rOff=off.r_lst, pi=pi_lst), probTTree, pOff=off.p2, 
                                 shGen=w.shape, scGen=w.scale, shSam=ws.shape, scSam=ws.scale, dateT=dateT)
      if (log(runif(1)) < sum(pTTree2)-sum(pTTree)){
        off.p_lst <- as.list(rep(off.p2, ntree))
        pTTree_lst <- as.list(pTTree2)
      }
    }
    
    if(("pi" %in% share) && updatePi){
      pi <- pi_lst[[1]]
      # use a beta prior
      pi2 <- pi + (runif(1)-0.5)*0.1
      if (pi2<0.01) pi2=0.02-pi2
      if (pi2>1) pi2=2-pi2
      if (verbose) message(sprintf("Proposing shared pi update %f->%f",pi,pi2))
      pTTree <- purrr::flatten_dbl(pTTree_lst)
      pTTree2 <- purrr::pmap_dbl(list(ttree=ttree_lst, rOff=off.r_lst, pOff=off.p_lst), probTTree, pi=pi2, 
                                 shGen=w.shape, scGen=w.scale, shSam=ws.shape, scSam=ws.scale, dateT=dateT)
      log_beta_ratio <- (prior_pi_a - 1) * log(pi2) + (prior_pi_b - 1) * log(1 - pi2) -
        (prior_pi_a - 1) * log(pi) - (prior_pi_b - 1) * log(1 - pi)
      if (log(runif(1)) < sum(pTTree2) - sum(pTTree) + ntree * log_beta_ratio){
        pi_lst <- as.list(rep(pi2, ntree))
        pTTree_lst <- as.list(pTTree2)
      } 
    }
    
    list(ctree=ctree_lst, pTTree=pTTree_lst, pPTree=pPTree_lst, neg=neg_lst, off.r=off.r_lst, off.p=off.p_lst, pi=pi_lst)
  }
  
  # Create outpout data structure: nested list
  record <- vector('list', ntree)
  for(i in seq_along(record)){
    record[[i]] <- vector("list", mcmcIterations/thinning)
  }
  
  # Initialize MCMC state in multi-tree space
  pTTree_lst <- vector("list", ntree)
  pPTree_lst <- vector("list", ntree)
  for(k in 1:ntree){
    ttree <- extractTTree(ctree_lst[[k]])
    pTTree_lst[[k]] <- probTTree(ttree$ttree,off.r,off.p,pi,w.shape,w.scale,ws.shape,ws.scale,dateT) 
    pPTree_lst[[k]] <- probPTreeGivenTTree(ctree_lst[[k]]$ctree,neg)  
  }
  mcmc_state <- list(ctree=ctree_lst, pTTree=pTTree_lst, pPTree=pPTree_lst, 
                     neg=neg_lst, off.r=off.r_lst, off.p=off.p_lst, pi=pi_lst)
  
  #Main MCMC loop
  if (verbose==F) pb <- utils::txtProgressBar(min=0,max=mcmcIterations,style = 3)
  for (i in 1:mcmcIterations) {
    
    # Update shared parameters
    out <- with(mcmc_state, one_update_share(ctree, pTTree, pPTree, neg, off.r, off.p, pi, share))
    mcmc_state[["ctree"]] <- out[["ctree"]]
    mcmc_state[["pTTree"]] <- out[["pTTree"]]
    mcmc_state[["pPTree"]] <- out[["pPTree"]]
    mcmc_state[["neg"]] <- out[["neg"]]
    mcmc_state[["off.r"]] <- out[["off.r"]]
    mcmc_state[["off.p"]] <- out[["off.p"]]
    mcmc_state[["pi"]] <- out[["pi"]]

    # Update unshared parameters
    state_new <- purrr::pmap(mcmc_state, one_update, not_share=not_share)
    
    if (i%%thinning == 0) {
      #Record things 
      if (verbose==F) utils::setTxtProgressBar(pb, i)
      if (verbose==T) message(sprintf('it=%d',i))
      for(k in seq_along(ctree_lst)){
        record[[k]][[i/thinning]]$ctree <- state_new[[k]]$ctree
        record[[k]][[i/thinning]]$pTTree <- state_new[[k]]$pTTree 
        record[[k]][[i/thinning]]$pPTree <- state_new[[k]]$pPTree 
        record[[k]][[i/thinning]]$neg <- state_new[[k]]$neg 
        record[[k]][[i/thinning]]$off.r <- state_new[[k]]$off.r
        record[[k]][[i/thinning]]$off.p <- state_new[[k]]$off.p
        record[[k]][[i/thinning]]$pi <- state_new[[k]]$pi
        record[[k]][[i/thinning]]$w.shape <- w.shape
        record[[k]][[i/thinning]]$w.scale <- w.scale
        record[[k]][[i/thinning]]$ws.shape <- ws.shape
        record[[k]][[i/thinning]]$ws.scale <- ws.scale
        record[[k]][[i/thinning]]$source <- with(state_new[[k]]$ctree, ctree[ctree[which(ctree[,4]==0),2],4])
        if (record[[k]][[i/thinning]]$source<=length(state_new[[k]]$ctree$nam)) 
          record[[k]][[i/thinning]]$source=state_new[[k]]$ctree$nam[record[[k]][[i/thinning]]$source] 
        else record[[k]][[i/thinning]]$source='Unsampled'
      }
    }
    # Assign updated state to current state
    mcmc_state <- purrr::transpose(state_new)

  }#End of main MCMC loop
  class(record)<-'resTransPhylo'
  return(record)
}
