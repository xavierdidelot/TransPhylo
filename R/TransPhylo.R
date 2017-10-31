#' @name TransPhylo
#' @title Inference of Transmission Tree under a Within-Host Evolution Model
#'
#' @description Inference of transmission tree under a within-host evolution model.
#'
#' The main functions of the package are:
#' \itemize{
#'
#' \item inferTTree
#' \item simulateOutbreak
#' }
#'
#' @author Xavier Didelot \email{xavier.didelot@gmail.com}
#'
#' @references Didelot et al (2014,2017) Molecular Biology an Evolution
#' @seealso https://github.com/xavierdidelot/TransPhylo

#' @importFrom Rcpp evalCpp
#' @importFrom utils setTxtProgressBar txtProgressBar head
#' @importFrom grDevices gray gray.colors palette rainbow
#' @useDynLib TransPhylo
#' @import stats
#' @import graphics
NULL
