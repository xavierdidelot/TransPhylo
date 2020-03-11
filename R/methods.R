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
#' @export
plot.ctree = function(x,...) {
  stopifnot(inherits(x, "ctree"))
  plotCTree(x,...)
}

#' Plotting for ttree
#' @param x Object of class ttree, ie  a transmission tree
#' @param ... Additional parameters are passed on
#' @return Plot of ttree
#' @export
plot.ttree = function(x,...) {
  stopifnot(inherits(x, "ttree"))
  plotTTree(x,...)
}

#' Plotting for ptree
#' @param x Object of class ptree, ie  a phylogenetic tree
#' @param ... Additional parameters are passed on
#' @return Plot of ptree
#' @export
plot.ptree = function(x,...) {
  stopifnot(inherits(x, "ptree"))
  phy=phyloFromPTree(x)
  ape::plot.phylo(phy,...)
  ape::axisPhylo()
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

