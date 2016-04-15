#' Plot pathway level statistics
#' 
#' @description Display graphical representation of pathway level statistic like:
#' i) two way comparison between the measured expression change and one of the 
#' factors computed by Pathway-Express (\code{\link{pe}}) or ii) the boostrap 
#' statistics of the same factors.
#' 
#' 
#' @param x an object of type \code{\link{pePathway-class}}
#' @param y if provided, the factor to be ploted (either \code{Acc} (default) or \code{Pert}; see \code{\link{pePathway-class}})
#' @param main title
#' @param ... Arguments to be passed to methods, such as \code{\link{par}}
#' @param type type of plot (either \code{two.way} (default) or \code{boot})
#' @param eps any value smaller than this will be ploted as 0
#' 
#' @author Calin Voichita and Sorin Draghici
#' 
#' @seealso \code{\link{pe}}, \code{\link{plot,peRes,missing-method}}, \code{\link{peNodeRenderInfo}}, \code{\link{peEdgeRenderInfo}}
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
#' # perform the pathway analysis (for more accurate results use nboot = 2000)
#' peRes <- pe(fc, graphs = kpg, ref = ref, nboot = 100, verbose = TRUE)
#' 
#' plot(peRes@@pathways[[50]])
#' 
#' plot(peRes@@pathways[[50]], "Pert", main = "Perturbation factor")
#' 
#' plot(peRes@@pathways[[50]], type = "boot")
#' 
#' plot(peRes@@pathways[[50]], "Pert", type = "boot", main = "Perturbation factor")
#' 
#' @rdname plot.pePathway-methods
#' 
#' @export
setMethod("plot", signature(x="pePathway", y="missing"),
          function(x, y, ..., type = "two.way", eps = 1e-6)
          {
            plot(x, y = "Acc", ..., type = type, eps = eps)           
          }
)

#' @rdname plot.pePathway-methods
#' 
#' @export
setMethod("plot", signature(x="pePathway", y="character"),
          function(x, y, main = "", ... , type = "two.way", eps = 1e-6)
          {
            if (!(y %in% c("Acc", "Pert")))
              stop("Undefined slot selected: ", y,".")
                        
            switch(type,
                   two.way={        
                     iy <- y
                     
                     extInput <- rep(0, length(slot(x, iy)))
                     names(extInput) <- names(slot(x, iy))
                     extInput[names(x@input)] <- x@input
                     
                     cl <- rep("black", length(slot(x, iy)))
                     cl[abs(slot(x, iy)) >= eps] <- "green"
                     cl[abs(extInput) >= eps] <- "blue"
                     cl[abs(slot(x, iy)) >= eps & abs(extInput) >= eps] <- "red"
                     
                     plot(slot(x, iy), extInput, pch = 16, xlab = y, ylab = "Log2 FC", main = main, ...)
                     abline(v=0,h=0, lwd = .5)
                     points(slot(x, iy), extInput, pch = 16, col = cl)
                     
                     return(invisible())
                   },
                   boot={
                     iy <- paste("t", y, sep = "")
                     
                     tB <- x@boot$t0[iy]
                     allB <- x@boot$t[,iy]
                     
                     tB <- (tB - mean(allB)) / sd(allB)
                     allB <- (allB - mean(allB)) / sd(allB)
                     
                     xlim <- c(min(allB, tB), max(allB, tB)) * 1.10
                     
                     plot(density(allB, from = xlim[1], to = xlim[2]), xlab = y, main = main, ...)
                     abline(v=0, lwd = .5)
                     abline(v=tB, lwd = 1, col = "red")
                     
                     return(invisible())
                   }
            )
            stop(type, " is not a valid plot type.")
            
          }
)

#' Plot Pathway-Express result
#' 
#' @description Display a two-way plot using two of the p-values from the Pathway-Express analysis.
#' 
#' @param x an object of type \code{\link{peRes-class}}
#' @param y vector of two p-values names to be combined using \code{comb.pv.func} (default: \code{c("pAcc", "pORA")}).
#' @param ... Arguments to be passed to methods, such as \code{\link{par}}.
#' @param comb.pv.func the function to combine the p-values - takes as input a vector of p-values 
#' and returns the combined p-value (default: \code{\link{compute.fisher}}).
#' @param adjust.method the name of the method to adjust the p-value (see \code{\link{p.adjust}})
#' @param threshold corrected p-value threshold
#' @param eps any value smaller than this will be considered as \code{eps} (default: \code{1e-6}).
#' 
#' @author Calin Voichita and Sorin Draghici
#' 
#' @seealso \code{\link{pe}}, \code{\link{summary.peRes}}, \code{\link{plot,pePathway,missing-method}}
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
#' # perform the pathway analysis (for more accurate results use nboot = 2000)
#' peRes <- pe(fc, graphs = kpg, ref = ref, nboot = 100, verbose = TRUE)
#' 
#' plot(peRes)
#' 
#' plot(peRes, c("pPert","pORA"), comb.pv.func = compute.normalInv, threshold = .01)
#' 
#' @rdname plot.peRes-methods
#' 
#' @export
setMethod("plot", signature(x="peRes", y="missing"),
          function(x, y, ... , comb.pv.func = compute.fisher, adjust.method = "fdr", threshold = .05, eps = 1e-6)
          {
            plot(x, y = c("pAcc", "pORA"), ... , comb.pv.func = comb.pv.func, adjust.method = adjust.method,
                 threshold = threshold, eps = eps)
          }
)

#' @rdname plot.peRes-methods
#'
#' @export
setMethod("plot", signature(x="peRes", y="character"),
          function(x, y, ... , comb.pv.func = compute.fisher, adjust.method = "fdr", threshold = .05, eps = 1e-6)
          {
            
            st <- Summary(x, comb.pv = y, comb.pv.func = comb.pv.func, adjust.method = adjust.method)
            st <- st[!is.na(st[, paste("pComb", adjust.method, sep = ".")]),]
            
            st[,y[1]][st[,y[1]] <= eps] <- eps
            st[,y[2]][st[,y[2]] <= eps] <- eps
            
            i <- st[, paste("pComb", adjust.method, sep = ".")] <= threshold
            thr.comb <- mean(min(st[!i,"pComb"]),max(st[i,"pComb"]))
            
            xrange <- c(min(-log(st[,y[1]])), max(-log(st[,y[1]])))
            xrange[2] <- xrange[2] + (xrange[2]-xrange[1]) * .1
            yrange <- c(min(-log(st[,y[2]])), max(-log(st[,y[2]]))) 
            yrange[2] <- yrange[2] + (yrange[2]-yrange[1]) * .1
            
            i <- seq(xrange[1], xrange[2], length.out=200)
            j <- seq(yrange[1], yrange[2], length.out=200)
            expGrid <- expand.grid(i,j)
            z <- apply(1/exp(expGrid), 1, comb.pv.func) <= thr.comb
            
            plot(c(min(-log(st[,y[1]])),min(-log(st[,y[2]]))), 
                 xlab = y[1], ylab = y[2], 
                 xlim = xrange,
                 ylim = yrange,
                 col = "white", ...)
            nonSig <- expGrid[!z,][chull(expGrid[!z,]),]
            sig <- expGrid[z,][chull(expGrid[z,]),]
            
            polygon(nonSig, col="gray90")
            polygon(sig, col="lightcyan")
            
            i <- st[, paste("pComb", adjust.method, sep = ".")] <= threshold
            
            points(-log(st[,y[1]]), -log(st[,y[2]]),  xlab = y[1], ylab = y[2], pch = 19)
            if (any(i))
            {
              points(-log(st[,y[1]])[i], -log(st[,y[2]])[i], pch = 21, bg = "red")
              text(-log(st[,y[1]])[i], -log(st[,y[2]])[i] - .5, labels = rownames(st)[i], 
                   cex = .75)
            }
            
            i <- st[, paste(y[1], adjust.method, sep = ".")] <= threshold
            if (any(i))
            {              
              thr1 <- mean(min(st[!i,y[1]]),max(st[i,y[1]]))            
              abline(v = -log(thr1), col = "red", lwd = 2, lty = 2)
            }
            
            i <- st[, paste(y[2], adjust.method, sep = ".")] <= threshold
            if (any(i))
            {
              thr2 <- mean(min(st[!i,y[2]]),max(st[i,y[2]]))                        
              abline(h = -log(thr2), col = "red", lwd = 2, lty = 2)
            }
          }
)
