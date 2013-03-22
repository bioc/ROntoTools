
#' Set node weights
#' 
#' 
#' @param graphList a list of \code{graph} (e.g., \code{\link{graphNEL}}) objects
#' @param weights named vector or matrix; if vector, the node is going to have the same weight in all graphs it appears;
#'                         if matrix, the rows represent nodes and columns represent graphs and the node will have different weights in each pathway  
#' @param defaultWeight the default weight for all nodes not set by the parameter \code{weights}
#' 
#' @return
#' 
#' The \code{graphList} with the node weights set.
#' 
#' @author
#' 
#' Calin Voichita and Sorin Draghici
#' 
#' @examples
#' 
#' # load the set of pathways
#' kpg <- keggPathwayGraphs("hsa")
#' 
#' kpg <- setNodeWeights(kpg)
#' 
#' nodeWeights(kpg[["path:hsa04110"]])
#' 
#' 
#' @export
setNodeWeights <- function(graphList, weights = NULL, defaultWeight = 1)
{
  if(is.null(weights))
  {
    graphList <- lapply(graphList,
                        function(g)
                        {                      
                          nodeDataDefaults(g, "weight") <- defaultWeight
                          return(g)
                        })
    return(graphList)
  }  
  if(is.vector(weights))
  {
    graphList <- lapply(graphList,
                        function(g)
                        {
                          i <- nodes(g)[nodes(g) %in% names(weights)]
                          
                          nodeDataDefaults(g, "weight") <- defaultWeight
                          
                          if(length(i) != 0)
                            nodeData(g, i, "weight") <- weights[i]
                          
                          return(g)
                        })
    return(graphList)
  }
  if(is.matrix(weights))
  {
    graphList <- Recall(graphList, numeric(0), defaultWeight)
    
    
    
    gi <- names(graphList)[names(graphList) %in% colnames(weights)]
    
    if(length(gi) == 0)
      return(graphList)
    
    for(i in 1:length(gi))
    {
      
      w <- weights[,gi[i]]
      names(w) <- rownames(weights)
      
      
      g <- graphList[[gi[i]]]
      
      j <- nodes(g)[nodes(g) %in% names(w)]
      
      nodeDataDefaults(g, "weight") <- defaultWeight
      
      if(length(j) != 0)
        nodeData(g, j, "weight") <- w[j]
      
      graphList[gi[i]] <- g
    }  
    
    return(graphList)
  }
  return(NULL)
}


#' Set gene weights based on edge type
#' 
#' \code{setEdgeWeights}
#' 
#' @param graphList a list of \code{\link{graphNEL}} objects
#' @param edgeTypeAttr edge attribute to be considered as the edge type. If the edge has multiple types, 
#'                     the edge type attribute is considered as a comma separeted list of types
#' @param edgeWeightByType named list of weigths, where the names of the list are the 
#'                         edge type (values of the attribute defined by \code{edgeTypeAttr})
#' @param defaultWeight default value for an edge with a type not defined in \code{edgeWeightByType}
#' @param combineWeights for the edges with multiple types, the function to be applied on the vector of weights
#' @param nodeOnlyGraphs boolean value marking if graphs with no edges should be returned or not; note that graphs with all edge weights equal to 0 are considered node only graphs
#' 
#' @return
#' 
#' The \code{graphList} with the edge weights set.
#' 
#' @author
#' 
#' Calin Voichita and Sorin Draghici
#' 
#' @examples
#' 
#' # load the set of pathways
#' kpg <- keggPathwayGraphs("hsa")
#' 
#' kpg <- setEdgeWeights(kpg)
#' 
#' edgeWeights(kpg[["path:hsa04110"]])
#' 
#' @export
setEdgeWeights <- function(graphList,
                           edgeTypeAttr = "subtype",
                           edgeWeightByType = list(activation = 1, inhibition = -1, expression = 1, repression = -1),
                           defaultWeight = 0,
                           combineWeights = sum,
                           nodeOnlyGraphs = FALSE)
{
  graphList <- lapply(graphList,
                      function(g)
                      {
                        gftM <- graph2ftM(g)
                        
                        if(nrow(gftM) == 0)
                          return(g)
                        
                        weigths <- sapply(strsplit(unlist(edgeData(g, gftM[,1], gftM[,2], edgeTypeAttr)), ","),
                                          function(x){
                                            weightNames <- names(edgeWeightByType)[match(x, names(edgeWeightByType))]
                                            weightNames <- weightNames[!is.na(weightNames)]
                                            
                                            if(length(weightNames) == 0)
                                              return(defaultWeight)
                                            
                                            return(combineWeights(unlist(edgeWeightByType[weightNames])))
                                          })
                        
                        suppressWarnings(g <- addEdge(gftM[,1], gftM[,2], g, weigths))
                        
                        return(g)
                      })
  
  if (!nodeOnlyGraphs)
    return(graphList[sapply(graphList, function(x) length(unlist(edges(x)))) != 0])
  
  return(graphList)
}
