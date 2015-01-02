#' @import KEGGREST
#' @import KEGGgraph
loadKEGGpathwayDataREST <- function(organism = "hsa",
                                    updateCache = FALSE,
                                    verbose = TRUE)
{ 
  dfUnparsed <- paste(system.file("extdata",package="ROntoTools"),
                      "/KEGGRESTunparsed_",organism,".RData", sep = "") 
  
  if (!file.exists(dfUnparsed) || updateCache)
  {
    pathList <- keggList("pathway", organism)
   
    if (verbose)
    {
      message("Downloading pathway data:")
      pb <- txtProgressBar(min = 0, max = length(pathList)-1, style = 3)
    }
    
    tmpDir <- tempfile("ROntoTools", tempdir())
    dir.create(tmpDir)
    
    
    allPathwayInfo <- lapply(names(pathList), 
                             function(pathID) {
                               p <- keggGet(pathID, "kgml")
                               pKgml <- file.path(tmpDir, paste(strsplit(pathID, ":")[[1]][2], ".kgml", sep = ""))
                               write(p, pKgml)
                               pathData <- parseKGML(pKgml)
                               file.remove(pKgml)
                               if(verbose)
                                 setTxtProgressBar(pb, getTxtProgressBar(pb) + 1)
                               return(pathData)
                             })
    dbInfo <- keggInfo("pathway")
    
    file.remove(tmpDir)
    
    save(allPathwayInfo, dbInfo, file = dfUnparsed)
  }else{    
    load(dfUnparsed)
    message(paste("Using cached pathway data. Database info:\n", dbInfo, sep =""))
  }
  
  return(allPathwayInfo)
}


#' Download and parse KEGG pathway data
#' 
#' @param organism organism code as defined by KEGG
#' @param targRelTypes target relation types
#' @param relPercThresh percentage of the number of relation types over all possible realtions in the pathway
#' @param nodeOnlyGraphs allow graphs with no edges
#' @param updateCache re-download KEGG data
#' @param verbose show progress of downloading and parsing
#' 
#' @return
#' 
#' A list of \code{\link{graphNEL}} objects encoding the pathway information.
#' 
#' @author
#' 
#' Calin Voichita and Sorin Draghici
#' 
#' @seealso \code{\link{keggPathwayNames}}
#' 
#' @examples
#'
#' # The pathway cache provided as part of the pathway contains only the 
#' # pathways that passed the default filtering. We recommend, re-downloading
#' # the pathways using the updateCache parameter
#' kpg <- keggPathwayGraphs("hsa")
#' 
#' # to update the pathway cache for human run:
#' # kpg <- keggPathwayGraphs("hsa", updateCache = TRUE)
#' # this is time consuming and depends on the available bandwith.
#' 
#' head(names(kpg))
#' 
#' kpg[["path:hsa04110"]]
#' head(nodes(kpg[["path:hsa04110"]]))
#' head(edges(kpg[["path:hsa04110"]]))
#' 
#' @importMethodsFrom KEGGgraph
#' @export
keggPathwayGraphs <- function(organism = "hsa", 
                              targRelTypes = c("GErel","PCrel","PPrel"),
                              relPercThresh = 0.9,
                              nodeOnlyGraphs = FALSE,
                              updateCache = FALSE,
                              verbose = TRUE)
{
  defaultParameters <- FALSE
  if ((organism == "hsa") & all(targRelTypes %in% c("GErel","PCrel","PPrel")) &
        (relPercThresh == 0.9) & (nodeOnlyGraphs == FALSE))
    defaultParameters <- TRUE  
  
  allPathwayInfo <- loadKEGGpathwayDataREST(organism, updateCache, verbose)
  
  if (defaultParameters & !updateCache)
  {
    message("Default parameters detected. Using pre-parsed data.")
    load(paste(system.file("extdata",package="ROntoTools"), "/kpgDefault.RData", sep = ""))
    return(pathwayGraphs)
  }
  
  l <- lapply(allPathwayInfo, function(path) {
    l <- sapply(edges(path), getType)
    if(length(l) == 0)
      return(0)
    t <- table(l)
    return(t)
  })
  
  allRelTypes <- unique(unlist(lapply(l, names)))
  
  counts <- do.call(rbind,lapply(l, function(x) as.vector(x[allRelTypes])))
  colnames(counts) <- allRelTypes
  
  counts[is.na(counts)] <- 0
  
  accIndex <- rowSums(counts[,targRelTypes], na.rm=T) / rowSums(counts, na.rm=T) >= relPercThresh
  accIndex[is.na(accIndex)] <- nodeOnlyGraphs
  allPathwayInfo <- allPathwayInfo[accIndex]
  
  names(allPathwayInfo) <- sapply(allPathwayInfo, getName)
  
  if (verbose)
  {
    message("Parsing pathway data:")
    pb <- txtProgressBar(min = 0, max = length(allPathwayInfo)-1, style = 3)
  }
  
  pathwayGraphs <- lapply(allPathwayInfo, function(g) 
  {
    
#     tryCatch(
      g <- KEGGpathway2Graph(g, expandGenes=TRUE)
#       , error = function(e) {})    
    if (class(g) != "graphNEL") {
      return(NULL)
    }
    
    kg <- new("graphNEL", nodes(g), edges(g), edgemode = "directed")

    if (length(getKEGGedgeData(g)) == 0)
    {
      if (verbose)
        setTxtProgressBar(pb, getTxtProgressBar(pb) + 1)
      if (nodeOnlyGraphs)
        return(kg)
      else
        return(NULL)
    }
    edgeDataDefaults(kg, "subtype") <- NA
    

    cf <- cbind(
      do.call(rbind, strsplit(names(getKEGGedgeData(g)), '~')),
      sapply(getKEGGedgeData(g), function(e) 
        paste(lapply(getSubtype(e), getName), collapse=","))
    );

    if (nrow(cf) < 2) {
      ucf <- cf
    } else {
      ucf <- cf[unique(rownames(cf)),]
      ucf[,3] <- tapply(cf[,3], rownames(cf), function(ll) return(paste( unique(ll), collapse = ',')))[rownames(ucf)]
    }
    

    relGeneTable <- data.frame(ucf, stringsAsFactors = FALSE)
    names(relGeneTable) <- c("from","to","subtype")

    edgeData(kg, relGeneTable$from, relGeneTable$to, "subtype") <- relGeneTable$subtype

    if (verbose)
      setTxtProgressBar(pb, getTxtProgressBar(pb) + 1)
    
    return(kg)
  })

  pathwayGraphs <- pathwayGraphs[!sapply(pathwayGraphs, is.null)]
  
  if (defaultParameters)
    save(pathwayGraphs, file = paste(system.file("extdata",package="ROntoTools"), "/kpgDefault.RData", sep = ""))
  
  return(pathwayGraphs)
}

#' Obtain KEGG pathway titles
#' 
#' @param organism organism code as defined by KEGG
#' @param updateCache re-download KEGG data
#' @param verbose show progress of downloading and parsing
#' 
#' @return
#' 
#' A named vector of pathway titles. The names of the vector are the pathway KEGG IDs.
#' 
#' @author
#' 
#' Calin Voichita and Sorin Draghici
#' 
#' @seealso \code{\link{keggPathwayGraphs}}
#' 
#' @examples
#'
#' kpn <- keggPathwayNames("hsa")
#' 
#' # to update the pathway cache for human run:
#' # kpn <- keggPathwayNames("hsa", updateCache = TRUE)
#' # this is time consuming and depends on the available bandwidth.
#' 
#' head(kpn)
#' 
#' @import KEGGgraph
#' @export
keggPathwayNames <- function(organism = "hsa", 
                              updateCache = FALSE,
                              verbose = TRUE)
{
  allPathwayInfo <- loadKEGGpathwayDataREST(organism, updateCache, verbose)
  
  allNames <- sapply(allPathwayInfo, getTitle)
  
  names(allNames) <- sapply(allPathwayInfo, getName)
  
  return(allNames)
}

#' Modified version of the same function from KEGGgraph
#' 
#' @keywords internal
KEGGpathway2Graph <- function (pathway, genesOnly = TRUE, expandGenes = TRUE) 
{
  stopifnot(is(pathway, "KEGGPathway"))
  pathway <- splitKEGGgroup(pathway)
  if (expandGenes) {
    pathway <- expandKEGGPathway(pathway)
  }
  knodes <- nodes(pathway)
  kedges <- unique(edges(pathway))
  node.entryIDs <- getEntryID(knodes)
  edge.entryIDs <- getEntryID(kedges)
  V <- node.entryIDs
  edL <- vector("list", length = length(V))
  names(edL) <- V
  if (is.null(nrow(edge.entryIDs))) {
    for (i in seq(along = edL)) {
      edL[[i]] <- list()
    }
  }
  else {
    for (i in 1:length(V)) {
      id <- node.entryIDs[i]
      hasRelation <- id == edge.entryIDs[, "Entry1ID"]
      if (!any(hasRelation)) {
        edL[[i]] <- list(edges = NULL)
      }
      else {
        entry2 <- unname(edge.entryIDs[hasRelation, "Entry2ID"])
        edL[[i]] <- list(edges = entry2)
      }
    }
  }
  gR <- new("graphNEL", nodes = V, edgeL = edL, edgemode = "directed")
  names(kedges) <- sapply(kedges, function(x) paste(getEntryID(x), 
                                                    collapse = "~"))
  env.node <- new.env()
  env.edge <- new.env()
  assign("nodes", knodes, envir = env.node)
  assign("edges", kedges, envir = env.edge)
  nodeDataDefaults(gR, "KEGGNode") <- env.node
  edgeDataDefaults(gR, "KEGGEdge") <- env.edge
  if (genesOnly) {    
    gR <- subGraphByNodeType(gR, "gene")
  }
  return(gR)
}

#' Modified version of the same function from KEGGgraph
#' 
#' @keywords internal
subGraphByNodeType <- function (graph, type = "gene") 
{
  kegg.node.data <- getKEGGnodeData(graph)
  types <- sapply(kegg.node.data, getType)
  isType <- grep(type, types)
  if (!any(isType)) {
    stop("No '", type, "' type found in the file, maybe it is a map file. Please try parsing the file with 'genesOnly=FALSE'\n")
  }
  new.nodes <- names(types[isType])
  new.edges <- edges(graph)
  new.edges <- new.edges[names(new.edges) %in% new.nodes]
  new.edges <- lapply(new.edges, function(eL) return(eL[eL %in% new.nodes]))
  
  subgraph <- new("graphNEL", new.nodes, new.edges, edgemode = "directed")
    
  subnodes <- new.nodes
  subedges <- getRgraphvizEdgeNames(subgraph)
  keggnodes <- get("nodes", nodeDataDefaults(graph, "KEGGNode"))
  keggedges <- get("edges", edgeDataDefaults(graph, "KEGGEdge"))
  subkeggnodes <- keggnodes[subnodes]
  subkeggedges <- keggedges[subedges]
  env.node <- new.env()
  env.edge <- new.env()
  assign("nodes", subkeggnodes, envir = env.node)
  assign("edges", subkeggedges, envir = env.edge)
  nodeDataDefaults(subgraph, "KEGGNode") <- env.node
  edgeDataDefaults(subgraph, "KEGGEdge") <- env.edge
  
  return(subgraph)
}
