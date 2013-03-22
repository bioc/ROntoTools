#' Pathway-Express result class
#' 
#' This class is used to encode the results of the pathway analysis performed by the function \code{\link{pe}}.
#' 
#' @details
#' 
#' The slots \code{input} and \code{ref} record global information related to the whole analysis, 
#' while the \code{pathways} slot records the specific results as \code{\link{pePathway}} for each one of the pathways used in the analysis.
#' 
#' 
#' @section Slots: 
#' 
#' \describe{
#'     \item{\code{pathways}:}{A list of \code{\link{pePathway}} objects.}
#'     \item{\code{input}:}{named vector of fold changes used for the analysis. The names of the vector are the IDs originaly used.}
#'     \item{\code{ref}:}{character vector containing the IDs used as reference in the analysis.}
#'     \item{\code{cutOffFree}:}{boolean value indicating if a cut-of-free analysis has been performed.}
#'   }
#' 
#' @author
#' 
#' Calin Voichita and Sorin Draghici
#'
#' @seealso \code{\link{pe}} \code{\link{pePathway}}
#'
#' @aliases peRes-class
#' @exportClass peRes
setClass("peRes",
         representation(pathways = "list",
                        input = "numeric",
                        ref = "character",
                        cutOffFree = "logical"),
         prototype(pathways = list(), cutOffFree = FALSE)
)

#' Class that encodes the result of Pathway-Express for a single pathway
#' 
#' 
#' @section Slots:
#' 
#' \describe{
#'    \item{map}{an object of type graph (e.g., \code{\link{graphNEL}}).}
#'    \item{input}{named vector of fold changes for genes on this pathway. The names of the genes are the orignal IDS used in the analysis}
#'    \item{ref}{vector of reference IDs on this pathway}
#'    \item{boot}{an object of class \code{boot} encoding the bootstrap information.}
#'    \item{PF}{the gene perturbation factors for all genes on the pathway, as computed by Pathway-Express.}
#'    \item{Acc}{the gene accumulations for all genes on the pathway, as computed by Pathway-Express.}
#' }
#' 
#' @author
#' 
#' Calin Voichita and Sorin Draghici
#'
#'             
#' @seealso \code{\link{pe}} \code{\link{peRes}}
#' 
#' @aliases pePathway-class
#' @import graph
#' @exportClass pePathway
setClass("pePathway", 
         representation(map = "graph",
                        input = "numeric",
                        ref = "character",
                        boot = "ANY",
                        PF = "numeric",
                        Acc = "numeric"
         ), 
         prototype(map = new("graphNEL")
         )
)
