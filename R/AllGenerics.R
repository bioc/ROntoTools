#' Retrieve the node weights of a graph
#' 
#' @description
#' 
#' A generic function that returns the node weights of a graph. If \code{index} is specified, only the weights of the specified nodes are returned. The user can control which node attribute is interpreted as the weight.
#' 
#' @param object A graph, any object that inherits the \code{graph} class.
#' @param index If supplied, a character or numeric vector of node names or indices.
#' @param ... Unused.
#' @param attr The name of the node attribute to use as a weight. You can view the list of defined node attributes and their default values using nodeDataDefaults.
#' @param default The value to use if \code{object} has no node attribute named by the value of \code{attr}. The default is the value 1.
#' 
#' @details
#' 
#' The weights of all nodes identified by the \code{index} are returned. If \code{index} is not supplied, the weights of all nodes are returned.
#' 
#' By default, \code{nodeWeights} looks for an node attribute with name "weight" and, if found, uses these values to construct the node weight vector.
#' You can make use of attributes stored under a different name by providing a value for the \code{attr} argument.
#' For example, if \code{object} is a graph instance with an node attribute named "WTS", then the call \code{nodeWeights(object, attr="WTS")} will attempt to use those values.
#' 
#' If the graph instance does not have an node attribute with name given by the value of the \code{attr} argument, \code{default} will be used as the weight for all nodes.
#' Note that if there is an attribute named by \code{attr}, then its default value will be used for nodes not specifically customized.
#' See nodeData and nodeDataDefaults for more information.
#' 
#' @return
#' 
#' A named vector with the node weights. The names of the vector are the names of the specified \code{index}, or all nodes if \code{index} was not provided.
#' 
#' @author
#' 
#' Calin Voichita and Sorin Draghici
#' 
#' @seealso
#' 
#' \link{nodes}, \link{nodeData}
#' 
#' @examples
#' 
#' library(graph)
#' V <- LETTERS[1:4]
#' g <- graphNEL(nodes = V, edgemode = "directed")
#' nodeWeights(g)
#' nodeWeights(g, "B")
#' nodeWeights(g, attr = "WT", default = 3)
#' 
#' @rdname nodeWeights
#' 
#' @export
setGeneric("nodeWeights",
           function(object, index, ..., attr="weight", default=1)
             standardGeneric("nodeWeights")
)
