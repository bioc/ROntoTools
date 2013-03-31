#' Extract edge render information from a \code{pePathway} object
#' 
#' @param x an object of class \code{\link{pePathway}}
#' @param pos.col color of the edges with possitive weight
#' @param pos.lty line type of the edges with possitive weight
#' @param pos.ah arrow head of the edges with possitive weight 
#' @param neg.col color of the edges with negative weight
#' @param neg.lty line type of the edges with negative weight 
#' @param neg.ah arrow head of the edges with negative weight 
#' @param zero.col color of the edges with zero weight 
#' @param zero.lty color of the edges with zero weight
#' @param zero.ah color of the edges with zero weight
#' 
#' @value a named list as expected by \code{\link{edgeRenderInfo}}
#' 
#' @author Calin Voichita and Sorin Draghici
#' 
#' @seealso \code{\link{edgeRenderInfo}},\code{\link{par}}
#' 
#' @examples
#' 
#' # load experiment
#' load(system.file("extdata/E-GEOD-21942.topTable.RData", package = "ROntoTools"))
#' fc <- top$logFC[top$adj.P.Val <= .01]
#' names(fc) <- top$entrez[top$adj.P.Val <= .01]
#' ref <- top$entrez
#' 
#' # load the set of pathways
#' kpg <- keggPathwayGraphs("hsa")
#' kpg <- setEdgeWeights(kpg)
#' kpg <- setNodeWeights(kpg, defaultWeight = 1)
#' 
#' # perform the pathway analysis
#' peRes <- pe(fc, graphs = kpg, ref = ref, nboot = 100, verbose = TRUE)
#' 
#' p <- peRes@@pathways[[50]]
#' g <- layoutGraph(p@@map, layoutType = "dot")
#' graphRenderInfo(g) <- list(fixedsize = FALSE)
#' edgeRenderInfo(g) <- peEdgeRenderInfo(p)
#' nodeRenderInfo(g) <- peNodeRenderInfo(p)
#' # notice the different type of edges in the graph (solid/dashed/dotted) 
#' renderGraph(g)
#' 
#' @export
peEdgeRenderInfo <- function(x,
                             pos.col = "black", pos.lty = "solid", pos.ah = "vee",
                             neg.col = "black", neg.lty = "dashed", neg.ah = "tee",
                             zero.col = "lightgray", zero.lty = "dotted", zero.ah = "none")
{
  stopifnot(class(x) == "pePathway")
  
  ew <- unlist(edgeWeights(x@map))
  
  aHead <- rep(zero.ah, length(edgeNames(x@map)))
  names(aHead) <- edgeNames(x@map)
  aHead[ew > 0] <- pos.ah
  aHead[ew < 0] <- neg.ah
  aHead <- aHead[setdiff(seq(along=ew), removedEdges(x@map))]
  
  eCol <- rep(zero.col, length(edgeNames(x@map)))
  names(eCol) <- edgeNames(x@map)
  eCol[ew > 0] <- pos.col
  eCol[ew < 0] <- neg.col
  eCol <- eCol[setdiff(seq(along=ew), removedEdges(x@map))]
  
  eStyle <- rep(zero.lty, length(edgeNames(x@map)))
  names(eStyle) <- edgeNames(x@map)  
  eStyle[ew > 0] <- pos.lty
  eStyle[ew < 0] <- neg.lty
  eStyle <- eStyle[setdiff(seq(along=ew), removedEdges(x@map))]
  
  return(list(
    arrowhead = aHead,
    col = eCol,
    lty = eStyle
  ))
}

#' Extract node render information from a \code{pePathway} object
#' 
#' @param x an object of class \code{\link{pePathway}}
#' @param y a string representing the factor to be represented (\code{Pert, Acc} or \code{input}; see \code{\link{pePathway}})
#' @param input.shape shape of nodes that have measured expression change
#' @param default.shape shape of all other nodes
#' @param pos.col color of nodes with a positive \code{y} factor
#' @param neg.col color of nodes with a negative \code{y} factor
#' @param zero.col color of nodes with the \code{y} factor equal to zero
#' 
#' @value a named list as expected by \code{\link{nodeRenderInfo}}
#' 
#' @author Calin Voichita and Sorin Draghici
#' 
#' @seealso \code{\link{nodeRenderInfo}},\code{\link{par}}
#' 
#' @examples
#' 
#' # load experiment
#' load(system.file("extdata/E-GEOD-21942.topTable.RData", package = "ROntoTools"))
#' fc <- top$logFC[top$adj.P.Val <= .01]
#' names(fc) <- top$entrez[top$adj.P.Val <= .01]
#' ref <- top$entrez
#' 
#' # load the set of pathways
#' kpg <- keggPathwayGraphs("hsa")
#' kpg <- setEdgeWeights(kpg)
#' kpg <- setNodeWeights(kpg, defaultWeight = 1)
#' 
#' # perform the pathway analysis
#' peRes <- pe(fc, graphs = kpg, ref = ref, nboot = 100, verbose = TRUE)
#' 
#' p <- peRes@@pathways[[50]]
#' g <- layoutGraph(p@@map, layoutType = "dot")
#' graphRenderInfo(g) <- list(fixedsize = FALSE)
#' edgeRenderInfo(g) <- peEdgeRenderInfo(p)
#' nodeRenderInfo(g) <- peNodeRenderInfo(p)
#' # notice the different type of nodes in the graph (box/circle)
#' # the color of each node represents the perturbation (red = positive, blue = negative)
#' # the shade represents the stregth of the perturbation 
#' renderGraph(g)
#' 
#' nodeRenderInfo(g) <- peNodeRenderInfo(p, "Acc")
#' # now, the color of each node represents the accumulation (red = positive, blue = negative)
#' # notice that square nodes with no parents have no accumulation
#' renderGraph(g)
#' 
#' @export
peNodeRenderInfo <- function(x, y = "Pert",
                             input.shape = "box",
                             default.shape = "ellipse",
                             pos.col = "red",
                             neg.col = "blue",
                             zero.col = "white")
{
  stopifnot(class(x) == "pePathway")
  stopifnot(y %in%  c("input", "Pert", "Acc"))
  
  nShape <- rep(default.shape, length(nodes(x@map)))
  names(nShape) <- nodes(x@map)
  nShape[names(x@input)] <- input.shape
  
  nFillColor <- rep(zero.col, length(nodes(x@map)))
  names(nFillColor) <- nodes(x@map)
  pf <- slot(x, y)
  nFillColor[pf <= 0] <- colorRampPalette(c(zero.col,neg.col))(256)[as.numeric(cut(abs(pf[pf<=0]), 256))]                                  
  nFillColor[pf >= 0] <- colorRampPalette(c(zero.col,pos.col))(256)[as.numeric(cut(abs(pf[pf>=0]), 256))]
  
  return(list(
    shape = nShape,
    fill = nFillColor
  ))
}