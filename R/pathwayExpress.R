#' Pathway-Express: Pathway analysis of signaling pathways
#' 
#' 
#' @param x named vector of log fold changes for the differentially expressed genes; \code{names(x)} must use the same id's as \code{ref} and the nodes of the \code{graphs}
#' @param graphs list of pathway graphs as objects of type \code{graph} (e.g., \code{\link{graphNEL}}); the graphs must be weighted graphs (i.e., have an atrribute \code{weight} for both nodes and edges)
#' @param ref the reference vector for all genes in the analysis; if the reference is not provided or it is identical to \code{names(x)} a cut-off free analysis is performed
#' @param nboot number of bootstrap iterations
#' @param verbose print progress output
#' @param cluster a cluster object created by makeCluster for parallel computations
#' @param seed an integer value passed to set.seed() during the boostrap permutations
#' 
#' @details
#' 
#' See details in the cited articles.
#' 
#' @return
#' 
#' An object of class \code{\link{peRes}}.
#' 
#' @author
#' 
#' Calin Voichita and Sorin Draghici
#' 
#' @references
#' 
#' Voichita C., Donato M., Draghici S.: "Incorporating gene significance in the impact analysis of signaling pathways", IEEE Machine Learning and Applications (ICMLA), 2012 11th International Conference on, Vol. 1, p.126-131, 2012
#' 
#' Tarca AL., Draghici S., Khatri P., Hassan SS., Kim J., Kim CJ., Kusanovic JP., Romero R.: "A Signaling Pathway Impact Analysis for Microarray Experiments", 2008, Bioinformatics, 2009, 25(1):75-82.
#' 
#' Khatri P., Draghici S., Tarca AL., Hassan SS., Romero R.: "A system biology approach for the steady-state analysis of gene signaling networks". Progress in Pattern Recognition, Image Analysis and Applications, Lecture Notes in Computer Science. 4756:32-41, November 2007. 
#' 
#' Draghici S., Khatri P., Tarca A.L., Amin K., Done A., Voichita C., Georgescu C., Romero R.: "A systems biology approach for pathway level analysis". Genome Research, 17, 2007. 
#' 
#' @seealso \code{\link{Summary}} \code{\link{keggPathwayGraphs}} \code{\link{setNodeWeights}} \code{\link{setEdgeWeights}}
#' 
#' @examples
#' 
#' # load a multiple sclerosis study (public data available in Array Express 
#' # ID: E-GEOD-21942)
#' # This file contains the top table, produced by the limma package with 
#' # added gene information. All the probe sets with no gene associate to them,
#' # have been removed. Only the most significant probe set for each gene has been
#' # kept (the table is already ordered by p-value)
#' # The table contains the expression fold change and signficance of each  
#' # probe set in peripheral blood mononuclear cells (PBMC) from 12 MS patients
#' # and 15 controls.
#' load(system.file("extdata/E-GEOD-21942.topTable.RData", package = "ROntoTools"))
#' head(top)
#' 
#' # select differentially expressed genes at 1% and save their fold change in a 
#' # vector fc and their p-values in a vector pv
#' fc <- top$logFC[top$adj.P.Val <= .01]
#' names(fc) <- top$entrez[top$adj.P.Val <= .01]
#' 
#' pv <- top$P.Value[top$adj.P.Val <= .01]
#' names(pv) <- top$entrez[top$adj.P.Val <= .01]
#' 
#' # alternativly use all the genes for the analysis
#' # NOT RUN: 
#' # fc <- top$logFC
#' # names(fc) <- top$entrez
#' 
#' # pv <- top$P.Value
#' # names(pv) <- top$entrez
#' 
#' # get the reference
#' ref <- top$entrez
#' 
#' # load the set of pathways
#' kpg <- keggPathwayGraphs("hsa")
#' 
#' # set the beta information (see the citated documents for meaning of beta)
#' kpg <- setEdgeWeights(kpg)
#' 
#' # inlcude the significance information in the analysis (see Voichita:2012 
#' # for more information)
#' # set the alpha information based on the pv with one of the predefined methods
#' kpg <- setNodeWeights(kpg, weights = alphaMLG(pv), defaultWeight = 1)
#' 
#' # perform the pathway analysis
#' # in order to obtain accurate results the number of boostraps, nboot, should 
#' # be increase to a number like 2000
#' peRes <- pe(fc, graphs = kpg, ref = ref, nboot = 100, verbose = TRUE)
#' 
#' # obtain summary of results
#' head(Summary(peRes))
#' 
#' @export
pe <- function(x, graphs, ref = NULL, nboot = 2000, verbose = TRUE, cluster = NULL, seed = NULL)
{
  cutOffFree <- FALSE
  if (is.null(ref))
  {
    ref <- names(x)
    cutOffFree <- TRUE
  }
  else
  {
    if (!any(is.na(match(ref,names(x)))))
    {
      warning("All the reference IDs are part of the input. Cut-off free analysis is performed.") 
      ref <- names(x)
      cutOffFree <- TRUE
    }
  }
  
  preservedSeed <- NULL
  if (exists(".Random.seed"))
    preservedSeed <- .Random.seed
  
  if(!is.null(seed))
    set.seed(seed)
  
  if (any(is.na(match(names(x), ref))))
  {
    warning("There are input IDs not available in the reference. These will be excluded from analysis.")
    x <- x[!is.na(match(names(x), ref))]
  }
  
  peRes <- pf.helper(x = x, ref = ref, graphs = graphs, nboot = nboot, verbose = verbose, cluster = cluster, seed = seed)
  peRes@cutOffFree <- cutOffFree
  
  if(is.null(preservedSeed))
  {
    if(!is.null(seed))
      rm(.Random.seed)
  }
  else
    .Random.seed <- preservedSeed
    
  return(peRes)
}

#' @import parallel
pf.helper <- function(x, graphs, ref = NULL, nboot = 2000, verbose = TRUE, cluster = NULL, seed = NULL)
{
  if(verbose)
  {
    message("Performing pathway analysis...")
    if(is.null(cluster))
    {
      pb <- txtProgressBar(min = 0, max = length(graphs), style = 3)  
    }
  }
  
  t1 <- Sys.time()
  if (is.null(cluster))
  {
    allBoot <- lapply(graphs, function(g, seed) {
      if(!is.null(seed))
        set.seed(seed)
      ret <- pe.boot(g, x = x, ref = ref, nboot = nboot)
      if(verbose)
        setTxtProgressBar(pb, getTxtProgressBar(pb) + 1)
      return(ret)
    }, seed = seed)
  }
  else
  {
    clusterExport(cluster, c("pe.boot", "compute.inverse", "compute.B"))
    clusterEvalQ(cluster, library(ROntoTools))
    
    allBoot <- parLapply(cluster, graphs, function(g, seed) {
      if(!is.null(seed))
        set.seed(seed)
      ret <- pe.boot(g, x = x, ref = ref, nboot = nboot)
      return(ret)
    }, seed = seed)
    
  }
  t2 <- Sys.time()
  
  if(verbose)
  {
    message("Analysis completed in ", format(t2-t1), ".")
  }
  
  allBoot <- allBoot[!sapply(allBoot, is.null)]
  
  peRes <- new("peRes", pathways = allBoot, input = x, ref = ref)
  
  return(peRes)
}

#' @import boot
#' @keywords internal
pe.boot <- function(g, x, ref, nboot, all.genes = F)
{
  inv <- compute.inverse(compute.B(g))
    
  pePath <- new("pePathway", 
                      map = g, 
                      input = x[names(x) %in% nodes(g)],
                      ref = ref[ref %in% nodes(g)])

  if (is.null(inv) | (length(pePath@input) == 0))
    return(NULL)
  
  # same number of DE genes at any position in the pathway 
  # (given by the gene from the pathway in the reference)
  ran.gen.de <- function(x, l)  {
    y <- sample(l$fc, length(x))
    names(y) <- sample(l$ref, length(x))
    return(y)
  }
    
  pePath@boot <- boot(pePath@input, 
    function(x, inv) {
      xx <- rep(0, nrow(inv));  names(xx) <- rownames(inv);
      xx[names(x)] <- x
      xx <- xx * nodeWeights(pePath@map, names(xx))
      tt = inv %*% xx;
      ret <- c(sum(abs(tt-xx)), sum(abs(tt)))
      names(ret) <- c("tAcc", "tPert")
      return(ret)
    }, 
            nboot,
            "parametric", ran.gen = ran.gen.de, mle = list(ref = pePath@ref, fc = as.numeric(x)),
            inv = inv
    )
  colnames(pePath@boot$t) <- names(pePath@boot$t0)
  
  xx <- rep(0, nrow(inv))
  names(xx) <- rownames(inv)
  xx[names(pePath@input)] <- pePath@input  
  pePath@PF = (inv %*% xx)[,1];
  pePath@Acc = pePath@PF - xx
    
  return(pePath) 
}

compute.inverse <- function(M, eps = 1e-5)
{
  if ( abs(det(M)) >= eps )
  {
    s = svd(M);
    inv = s$v %*% diag(1/s$d) %*% t(s$u)
    rownames(inv) <- colnames(inv) <- rownames(M)
    return(inv)
  }else{
    return(NULL)
  }
}

compute.B <- function(g, non.zero = TRUE)
{
  if (non.zero)
    # number of downstream genes (like in SPIA)
    nds <- sapply(edgeWeights(g), function(x) sum(x != 0 ))
  else
    nds <- sapply(edges(g), length)
  
  # add 1 for all genes that do not have downstream genes to avoid division by 0
  # this does not affect the computation
  nds[nds == 0] <- 1
  
  # compute B = (I - beta/nds)
  #B <- t(diag(length(nodes(g))) - as(g, "matrix") / nds)
  B <- diag(length(nodes(g))) - t(as(g, "matrix")) / matrix(nds, byrow=TRUE, nrow = length(nds), ncol = length(nds))
  
  return(B)
}

  