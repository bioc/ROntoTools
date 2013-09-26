#' Summarize the results of a Pathway-Express analysis
#' 
#' 
#' @usage Summary(x, pathNames = NULL, totalAcc = TRUE, totalPert = TRUE, normalize = TRUE, 
#'  pPert = TRUE, pAcc = TRUE, pORA = TRUE, 
#'  comb.pv = c("pPert", "pORA"), comb.pv.func = compute.fisher,
#'  order.by = "pComb", adjust.method = "fdr")
#' 
#' @param x Pathways-Express result object obtained using \code{\link{pe}}
#' @param pathNames named vector of pathway names; the names of the vector are the IDs of the pathways
#' @param totalAcc boolean value indicating if the total accumulation should be computed
#' @param totalPert boolean value indicating if the total perturbation should be computed
#' @param normalize boolean value indicating if normalization with regards to the boostrap simulations should be performed on totalAcc and totalPert
#' @param pPert boolean value indicating if the significance of the total perturbation in regards to the bootstrap permutations should be computed
#' @param pAcc boolean value indicating if the significance of the total accumulation in regards to the bootstrap permutations should be computed
#' @param pORA boolean value indicating if the over-represtation p-value should be computed
#' @param comb.pv vector of the p-value names to be combine (any of the above p-values)
#' @param comb.pv.func the function to combine the p-values; takes as input a vector of p-values and returns the combined p-value
#' @param order.by the name of the p-value that is used to order the results
#' @param adjust.method the name of the method to adjust the p-value (see \link{p.adjust})
#' 
#' @seealso \code{\link{pe}}
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
#' # obtain summary of results
#' head(Summary(peRes))
#' 
#' kpn <- keggPathwayNames("hsa")
#' 
#' head(Summary(peRes))
#' 
#' head(Summary(peRes, pathNames = kpn, totalAcc = FALSE, totalPert = FALSE, 
#'              pAcc = FALSE, pORA = FALSE, comb.pv = NULL, order.by = "pPert"))
#' 
#' @rdname Summary-methods
#' 
#' @aliases Summary.peRes
#' @aliases Summary,peRes-method
#' @export
setMethod("Summary", c("x" = "peRes"),
          function(x, pathNames = NULL, totalAcc = TRUE, totalPert = TRUE, normalize = TRUE, 
                   pPert = TRUE, pAcc = TRUE, pORA = TRUE, 
                   comb.pv = c("pPert", "pORA"), comb.pv.func = compute.fisher,
                   order.by = "pComb", adjust.method = "fdr")
          {  
            ifelse <- function(test, trueCase, falseCase){
              if(test) return(trueCase)
              else return(falseCase)
            }
            
            pathStats <- function(pePath)
            {
              pStats <- NULL
              
              pStats$totalAcc <- ifelse(totalAcc, ifelse(!pePath@asGS, get.totalAcc(pePath), NA), NULL)
              pStats$totalPert <- ifelse(totalPert, ifelse(!pePath@asGS, get.totalPert(pePath), NA), NULL)
              
              pStats$totalAccNorm <- ifelse(totalAcc & normalize, ifelse(!pePath@asGS, get.totalAccNorm(pePath), NA), NULL)
              pStats$totalPertNorm <- ifelse(totalPert & normalize, ifelse(!pePath@asGS, get.totalPertNorm(pePath), NA), NULL)
              
              pStats$pPert <- ifelse(pPert, ifelse(!pePath@asGS, compute.pPert(pePath), NA), NULL)
              pStats$pAcc <- ifelse(pAcc, ifelse(!pePath@asGS, compute.pAcc(pePath), NA), NULL)
              
              pStats$pORA <- ifelse(pORA & !x@cutOffFree, compute.pORA(pePath, length(x@input), length(x@ref)), NULL)
              
              pStats$pComb <- ifelse(!is.null(comb.pv) & !any(is.null(pStats[comb.pv])), 
                                     ifelse(!any(is.na(pStats[comb.pv])), as.numeric(comb.pv.func(unlist(pStats[comb.pv]))), NA), NULL)
              
              return(unlist(pStats))
            }
            
            if (pORA & x@cutOffFree)
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
              if(!all(comb.pv %in% c("pPert","pAcc","pORA")))
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
            
            topStats <- data.frame(do.call(rbind, lapply(x@pathways, pathStats)))
            
            if(!is.null(pathNames))
            {
              pathNames <- pathNames[rownames(topStats)]
              topStats <- cbind(pathNames, topStats)
            }
            
            if(order.by %in% colnames(topStats))
            {
              topStats <- topStats[order(topStats[,order.by]),]
            }
            
            allPVs <- c("pPert","pAcc","pORA", "pComb")
            
            lapply(allPVs[allPVs %in% colnames(topStats)],
                   function(pv)
                     topStats[[paste(pv, "." , adjust.method, sep = "")]] <<- p.adjust(topStats[[pv]], adjust.method)
            )
            return(topStats)
          }
)