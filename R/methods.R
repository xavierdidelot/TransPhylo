#' Plotting for resTransPhylo
#' @param x Output from inferTTree
#' @param ... Additional parameters are passed on
#' @return Plot of TransPhylo results
#' @export
plot.resTransPhylo = function(x,...) {
  stopifnot(inherits(x, "resTransPhylo"))
  plotTraces(x)
}

#' Plotting for ctree
#' @param x Object of class ctree, ie a colored phylogenetic tree
#' @param ... Additional parameters are passed on
#' @return Plot of ctree
#' @examples plot(simulateOutbreak())
#' @export
plot.ctree = function(x,...) {
  stopifnot(inherits(x, "ctree"))
  plotCTree(x,...)
}

#' Plotting for ttree
#' @param x Object of class ttree, ie  a transmission tree
#' @param type Type of plot to display, can be 'detailed' or 'summarised' (default)
#' @param w.shape Shape parameter of the generation time, needed for detailed plot only
#' @param w.scale Scale parameter of the generation time, needed for detailed plot only
#' @param ... Additional parameters are passed on
#' @return Plot of ttree
#' @examples 
#' plot(extractTTree(simulateOutbreak()))
#' @export
plot.ttree = function(x,type='summarised',w.shape=NA,w.scale=NA,...) {
  stopifnot(inherits(x, "ttree"))
  if (type=='summarised') plotTTree2(x,...) 
  if (type=='detailed') {
    if (is.na(w.shape)||is.na(w.scale)) stop('You need to specify w.shape and w.scale to display this plot.')
    else plotTTree(x,w.shape=w.shape,w.scale=w.scale,...)
  }
}

#' Plotting for ptree
#' @param x Object of class ptree, ie  a phylogenetic tree
#' @param ... Additional parameters are passed on to ape::plot.phylo
#' @return Plot of ptree
#' @examples 
#' plot(ptreeFromPhylo(ape::rtree(5),2020))
#' @export
plot.ptree = function(x,...) {
  stopifnot(inherits(x, "ptree"))
  phy=phyloFromPTree(x)
  ape::plot.phylo(phy,...)
  ape::axisPhylo(backward = F)
}

#' Print function for resTransPhylo objects
#' @param x output from inferTTree
#' @param ... Additional parameters are passed on
#' @return Print out details of TransPhylo results
#' @export
print.resTransPhylo <- function(x, ...)
{
  stopifnot(inherits(x, "resTransPhylo"))
  cat( 'Result from TransPhylo analysis\n')
  coda=convertToCoda(x,0.5)
  for (nam in colnames(coda)) {
    v=coda[,nam]
    v=sort(v)
    vals=c(mean(v),v[pmax(1,floor(length(v)*c(0.025,0.975)))])
    cat(sprintf('%s=%.2e [%.2e;%.2e]\n',nam,vals[1],vals[2],vals[3]))
  }
  invisible(x)
}

#' Print function for ctree objects
#' @param x Object of class ctree, ie a colored phylogenetic tree
#' @param ... Additional parameters are passed on
#' @return Print out details of the ctree
#' @examples 
#' print(simulateOutbreak())
#' @export
print.ctree <- function(x, ...)
{
  stopifnot(inherits(x, "ctree"))
  cat( 'Colored phylogenetic tree\n')
  cat(sprintf('Number of sampled individuals=%d\n',length(x$nam)))
  cat(sprintf('Total number of hosts=%d\n',max(x$ctree[,4])))
  invisible(x)
}

#' Print function for ttree objects
#' @param x Object of class ttree, ie a transmission tree
#' @param ... Additional parameters are passed on
#' @return Print out details of the ttree
#' @examples 
#' print(extractTTree(simulateOutbreak()))
#' @export
print.ttree <- function(x, ...)
{
  stopifnot(inherits(x, "ttree"))
  cat( 'Transmission tree\n')
  cat(sprintf('Number of sampled individuals=%d\n',length(x$nam)))
  cat(sprintf('Total number of hosts=%d\n',nrow(x$ttree)))
  invisible(x)
}

#' Print function for ptree objects
#' @param x Object of class ptree, ie a phylogenetic tree
#' @param ... Additional parameters are passed on
#' @return Print out details of the ptree
#' @examples 
#' print(extractPTree(simulateOutbreak()))
#' @export
print.ptree <- function(x, ...)
{
  stopifnot(inherits(x, "ptree"))
  cat( 'Phylogenetic tree\n')
  cat(sprintf('Number of sampled individuals=%d\n',length(x$nam)))
  invisible(x)
}

#' Summary function for resTransPhylo objects
#' @param object output from inferTTree
#' @param ... Passed on to print.phylo
#' @return Print out details of TransPhylo results
#' @export
summary.resTransPhylo <- function(object, ...){
  stopifnot(inherits(object, "resTransPhylo"))
  cat( 'Result from TransPhylo analysis\n')
  invisible(object)
}

#' Convert to coda mcmc format
#' @param x Output from inferTTree
#' @param burnin Proportion of the MCMC output to be discarded as burnin
#' @return mcmc object from coda package
#' @export
as.mcmc.resTransPhylo <- function(x,burnin=0.5) {
  return(convertToCoda(x,burnin))
}

#' Return the date of last sample from a ttree or ctree or ptree
#' @param x A transmission tree or colored tree or phylogenetic tree
#' @return date of the last sample
#' @export
dateLastSample <- function(x) {
  if (inherits(x,'ctree')) return(max(x$ctree[,1]))
  if (inherits(x,'ttree')) return(max(x$ttree[,2],na.rm = T))
  if (inherits(x,'ptree')) return(max(x$ptree[,1]))
}

#' Return the combined tree corresponding to a given iteration of the TransPhylo results
#' @param res Output from inferTTree command
#' @param iteration Number of the iteration to be extracted
#' @return The colored tree at the specified iteeatino
#' @export
extractCTree <- function(res,iteration) {
  stopifnot(inherits(res, "resTransPhylo"))
  return(res[[iteration]]$ctree)
}
