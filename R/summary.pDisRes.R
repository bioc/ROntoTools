#' Summarize the results of a primary dis-regulation (pDis) analysis
#' 
#' @usage summary.pDisRes(object, ..., pathNames = NULL, totalpDis = TRUE, normalize = TRUE, 
#'  ppDis = TRUE, pORA = TRUE, 
#'  comb.pv = c("ppDis", "pORA"), comb.pv.func = compute.fisher,
#'  order.by = "pComb", adjust.method = "fdr")
#' 
#' @param object pDis analysis result object obtained using \code{\link{pDis}}
#' @param ... ignored
#' @param pathNames named vector of pathway names; the names of the vector are the IDs of the pathways
#' @param totalpDis boolean value indicating if the total primary dis-regulation should be computed
#' @param normalize boolean value indicating if normalization with regards to the boostrap simulations should be performed on totalpDis
#' @param ppDis boolean value indicating if the significance of the total primary dis-regulation in regards to the bootstrap permutations should be computed
#' @param pORA boolean value indicating if the over-represtation p-value should be computed
#' @param comb.pv vector of the p-value names to be combine (any of the above p-values)
#' @param comb.pv.func the function to combine the p-values; takes as input a vector of p-values and returns the combined p-value
#' @param order.by the name of the p-value that is used to order the results
#' @param adjust.method the name of the method to adjust the p-value (see \link{p.adjust})
#' 
#' @seealso \code{\link{pDis}}
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
#' pDisRes <- pDis(fc, graphs = kpg, ref = ref, nboot = 100, verbose = TRUE)
#' 
#' # obtain summary of results
#' head(summary(pDisRes))
#' 
#' kpn <- keggPathwayNames("hsa")
#' 
#' head(summary(pDisRes))
#' 
#' head(summary(pDisRes, pathNames = kpn, totalpDis = FALSE, 
#'             pORA = FALSE, comb.pv = NULL, order.by = "pDis"))
#' 
#' @export
summary.pDisRes <- function(object, ..., pathNames = NULL, totalpDis = TRUE, normalize = TRUE, 
                   ppDis = TRUE, pORA = TRUE, 
                   comb.pv = c("ppDis", "pORA"), comb.pv.func = compute.fisher,
                   order.by = "pComb", adjust.method = "fdr")
{  
  ifelse <- function(test, trueCase, falseCase){
    if(test) return(trueCase)
    else return(falseCase)
  }
  
  pathStats <- function(pDisPath)
  {
    pStats <- NULL
    
    
    pStats$totalpDis <- ifelse(totalpDis, ifelse(!pDisPath@asGS, get.totalpDis(pDisPath), NA), NULL)
    
    pStats$totalpDisNorm <- ifelse(totalpDis & normalize, ifelse(!pDisPath@asGS, get.totalpDisNorm(pDisPath), NA), NULL)
    
    pStats$ppDis <- ifelse(ppDis, ifelse(!pDisPath@asGS, compute.ppDis(pDisPath), NA), NULL)
    
    pStats$pORA <- ifelse(pORA & !object@cutOffFree, compute.pORA(pDisPath, length(object@input), length(object@ref)), NULL)
    
    pStats$pComb <- ifelse(!is.null(comb.pv) & !any(is.null(pStats[comb.pv])), 
                           ifelse(!any(is.na(pStats[comb.pv])), as.numeric(comb.pv.func(unlist(pStats[comb.pv]))), NA), NULL)
    
    return(unlist(pStats))
  }
  
  if (pORA & object@cutOffFree)
  {
    pORA <- FALSE
    if ("pORA" %in% comb.pv)
    {
      order.by <- setdiff(comb.pv, "pORA")[1]
      comb.pv <- NULL
    }
    message("The over-representaion p-value is not defined for cut-off free analysis and will not be computed!")  
  }
  
  if(!is.null(comb.pv))
  {
    if(!all(comb.pv %in% c("ppDis","pORA")))
    {
      warning("The p-value to be combined are not specified correctly. No combination p-value will be calculated!")
      comb.pv <- NULL
      if(order.by == "pComb")
        order.by <- NULL
    }else{
      for(i in 1:length(comb.pv))
        assign(comb.pv[i], TRUE)
    }
  }
  
  topStats <- data.frame(do.call(rbind, lapply(object@pathways, pathStats)))
  
  if(!is.null(pathNames))
  {
    pathNames <- pathNames[rownames(topStats)]
    topStats <- cbind(pathNames, topStats)
  }
  
  if(order.by %in% colnames(topStats))
  {
    topStats <- topStats[order(topStats[,order.by]),]
  }
  
  allPVs <- c("ppDis","pORA", "pComb")
  
  lapply(allPVs[allPVs %in% colnames(topStats)],
         function(pv)
           topStats[[paste(pv, "." , adjust.method, sep = "")]] <<- p.adjust(topStats[[pv]], adjust.method)
  )
  return(topStats)
}