#' Autodetect region neighbours
#'
#' This function allows you to automatically identify neighbours of a selected region.
#'
#' One way of using method this is to find spots surrounding a certain cluster. First, you need to make sure
#' the identity of the Seurat object is set to the meta.data column that you want to use, so for example
#' `se <- SetIdent(se, value = "seurat_clusters")` if you want to use the default seurat clusters.
#' Then you select the label that defined the region of interest using the `id` parameter, so for example
#' `ìd = "1"` will use cluster 1 as the region. If you set the `keep.idents` parameter to TRUE, the cluster ids
#' of the neighbouring spots will be kept in the result, otherwise they will be returned as one single goup.
#' You can also activate the `keep.within.id` parameter to include all spots of the selected region in the output,
#' otherwise only the spots along the region border will be kept.
#'
#' @param object Seurat object.
#' @param id Group label used to define region, e.g. a cluster id. The region will be selected using the
#' active identity (see `Idents`).
#' @param column.key Sets the name of the column in the meta.data slot where the output is stored. [default: 'nbs_']
#' @param keep.idents If set to TRUE, the identities of the neighbours are kept in the output, otherwise all
#' neighbours will be named 'nbs_id' where id is the group label defined by the aprameter `ìd` [default: FALSE]
#' @param keep.within.id If set to TRUE, all id spots are kept, otherwise only the spots with outside neighbours are kept
#' @param verbose Print messages
#' 
#' @importFrom Seurat Idents
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Cluster data, find neighbours to cluster 10 and then plot the results
#' se <- FindNeighbours(se) %>% FindClusters()
#' se <- RegionNeighbours(se, id = 10)
#' ST.FeaturePlot(se, features = "nbs_10")
#'
#' # Find neighbours to cluster 10 and keep all id spots
#' se <- RegionNeighbours(se, id = 10)
#' ST.FeaturePlot(se, features = "nbs_10", keep.within.id = TRUE)
#' }
#'
RegionNeighbours <- function (
  object,
  id,
  column.key = "nbs_",
  keep.idents = FALSE,
  keep.within.id = FALSE,
  verbose = FALSE
) {

  if (verbose) cat("Creating Connected Netork using KNN ... \n")
  spatnet <- GetSpatNet(object)
  spatnet <- do.call(rbind, lapply(seq_along(spatnet), function(i) {
    spnet <- spatnet[[i]]
    spnet$sample <- paste0(i)
    return(spnet)
  }))

  spatnet$label_from <- Idents(object)[spatnet$from] %>% as.character()
  spatnet$label_to <- Idents(object)[spatnet$to] %>% as.character()
  spatnet <- subset(spatnet, label_from == paste0(id))
  if (verbose) cat(paste0("Found ", nrow(spatnet), " neighbours for id ", id, " ... \n"))

  if (!keep.within.id) {
    if (verbose) cat("Excluding neighbours from the same group ... \n")
    spatnet <- subset(spatnet, label_from != label_to)
    if (verbose) cat(paste0(nrow(spatnet), " neighbours left ... \n"))
  }

  column.key <- paste0(column.key, id)
  column.value <- setNames(rep(NA, ncol(object)), nm = colnames(object))
  column.value[unique(union(spatnet$from, spatnet$to))] <- Idents(object)[unique(union(spatnet$from, spatnet$to))] %>% as.character()
  if (keep.idents) {
    if (verbose) cat("Naming neighbours to id_nb_to* ... \n")
    column.value[column.value != paste0(id) & !is.na(column.value)] <- paste0(column.value[column.value != paste0(id) & !is.na(column.value)], "_nb_to_", id)
  } else {
    if (verbose) cat(paste0("Naming all neighbours ", paste0("nbs_", id), " ... \n"))
    column.value[column.value != paste0(id) & !is.na(column.value)] <- paste0("nbs_", id)
  }

  if (verbose) cat(paste0("Saving neighbour ids to column '", column.key, "' ... \n"))
  object@meta.data[, column.key] <- column.value
  if (verbose) cat(paste0("Finished. \n\n"))
  return(object)
}
