#' @rdname nodeWeights
#' @aliases nodeWeights,graph,character-method
#' @import graph
#' @export
setMethod("nodeWeights", signature(object="graph", index="character"),
          function (object, index, attr, default) 
          {
            if (!is.character(attr) || length(attr) != 1) 
              stop("'attr' must be character(1)")
            if (!attr %in% names(nodeDataDefaults(object))) {
              nodeDataDefaults(object, attr) <- default
            }
            nw <- nodeData(object, index, attr = attr)
            if (length(nw) != 0) 
              return(unlist(nw))
            else
              return(numeric(0))
          }
)

#' @rdname nodeWeights
#' @aliases nodeWeights,graph,numeric-method
#' @import graph
#' @export
setMethod("nodeWeights", signature(object="graph", index="numeric"),
          function (object, index, attr, default) 
          {
            index <- nodes(object)[index]
            nodeWeights(object, index, attr = attr, default = default)
          })

#' @rdname nodeWeights
#' @aliases nodeWeights,graph,missing-method
#' @import graph
#' @export
setMethod("nodeWeights", signature(object="graph", index="missing"),
          function (object, index, attr, default) 
          {
            index <- nodes(object)
            nodeWeights(object, index, attr = attr, default = default)
          })