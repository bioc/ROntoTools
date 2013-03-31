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
  if ((organism == "hsa") & all(targRelTypes == c("GErel","PCrel","PPrel")) &
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
  
  accIndex <- rowSums(counts[,targRelTypes], na.rm=T) / rowSums(counts, na.rm=T) >= relPercThresh
  accIndex[is.na(accIndex)] <- FALSE  
  allPathwayInfo <- allPathwayInfo[accIndex]
  
  names(allPathwayInfo) <- sapply(allPathwayInfo, getName)
  
  if (verbose)
  {
    message("Parsing pathway data:")
    pb <- txtProgressBar(min = 0, max = length(allPathwayInfo)-1, style = 3)
  }
  
  pathwayGraphs <- lapply(allPathwayInfo, function(g) 
  {
    g <- KEGGgraph::KEGGpathway2Graph(g)
    kg <- new("graphNEL", nodes(g), edges(g), edgemode = "directed")
    
    if (length(getKEGGedgeData(g)) == 0)
    {
      if (verbose)
        setTxtProgressBar(pb, getTxtProgressBar(pb) + 1)
      return(NULL)
    }
    edgeDataDefaults(kg, "subtype") <- NA  
        
    relGeneTable <- data.frame(cbind(
      do.call(rbind, strsplit(names(getKEGGedgeData(g)), '~')),
      sapply(getKEGGedgeData(g), function(e) 
        paste(lapply(getSubtype(e), getName), collapse=","))
    ), stringsAsFactors = FALSE)
    names(relGeneTable) <- c("from","to","subtype")
        
    edgeData(kg, relGeneTable$from, relGeneTable$to, "subtype") <- relGeneTable$subtype

    if (verbose)
      setTxtProgressBar(pb, getTxtProgressBar(pb) + 1)
    
    return(kg)
  })

  pathwayGraphs <- pathwayGraphs[!sapply(pathwayGraphs, is.null)]
  
  if (defaultParameters)
    save(pathwayGraphs, paste(system.file("extdata",package="ROntoTools"), "/kpgDefault.RData", sep = ""))
  
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

