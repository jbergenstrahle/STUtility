#' Create Spatial Networks
#'
#' Create spatial networks based on spot center to center distances from a Seurat object.
#'
#' @param object Seurat object
#' @param nNeighbours Number of nearest neighbours to calculate for each spot. The default
#' number of neighbours is 6 for the 'Visium' platform and 4 for the '1k' and '2k' platforms.
#' @param maxdist Distance cut-off for nearest neighbours to consider. The default is 1.5 for the
#' 'Visium' and '2k' platforms and 2 for the '1k' platform.
#' @param minK Minimum nearest neigbhours if maxdist is not provided [default: 0]
#'
#' @importFrom dplyr group_by mutate ungroup
#' @importFrom dbscan kNN
#' @importFrom igraph graph_from_data_frame
#'
#' @export
#' @examples
#' spatial.networks <- GetSpatNet(se)

GetSpatNet <- function (
  object,
  nNeighbours = NULL,
  maxdist = NULL,
  minK = 0
) {

  # Check to see if Staffli object is present
  if (!"Staffli" %in% names(object@tools)) stop("Staffli object is missing from Seurat object... \n", call. = FALSE)
  st.object <- object@tools$Staffli

  # Check if images are loaded
  if (length(x = st.object@rasterlists) == 0) stop("He images need to be loaded before running this function ... \n")

  # spatial information
  xys = setNames(st.object@meta.data[, c("pixel_x", "pixel_y", "sample")], c("x", "y", "sample"))

  # Split x, y, s
  xys.list <- split(xys, xys$sample)

  # Obtain platforms
  platforms <- st.object@platforms

  # Compute network
  knn_spatial.norm.list <- lapply(seq_along(xys.list), function(i) {
    xys <- xys.list[[i]]

    # vector matching spot_ID and order
    spotnames <- rownames(xys)
    names(spotnames) <- c(1:nrow(xys)) %>% paste0()

    # Get spot distances
    sdist <- st.object@pixels.per.um[i]

    nNeighbours <- nNeighbours %||% ifelse(platforms[i] == "Visium", 6, 4)
    maxdist <- maxdist %||% ifelse(platforms[i] == "1k", 270*sdist, 150*sdist)

    knn_spatial <- kNN(x = xys[, c("x", "y")] %>% as.matrix(), k = nNeighbours)
    knn_spatial.norm <- data.frame(from = rep(1:nrow(knn_spatial$id), nNeighbours),
                                 to = as.vector(knn_spatial$id),
                                 weight = 1/(1 + as.vector(knn_spatial$dist)),
                                 distance = as.vector(knn_spatial$dist))
    #nw_spatial.norm = igraph::graph_from_data_frame(knn_spatial.norm, directed = FALSE)
    #CN <- igraph::as_adjacency_matrix(nw_spatial.norm)

    # create network for coordinates
    spatnet <- knn_spatial.norm
    spatnet$from <- spotnames[spatnet$from]
    spatnet$to <- spotnames[spatnet$to]
    spatnet <- spatnet %>% group_by(from) %>% mutate(rnk = rank(distance)) %>% ungroup()
    spatnet =  subset(spatnet, distance <= maxdist | rnk <= minK)

    # Add coordinates
    spatnet <- cbind(spatnet, setNames(xys[spatnet$from, 1:2], paste0("start_", c("x", "y"))))
    spatnet <- cbind(spatnet, setNames(xys[spatnet$to, 1:2], paste0("end_", c("x", "y"))))

    return(spatnet)
  })

  return(knn_spatial.norm.list)
}


# TODO: fix imported libs (spdep not loaded), fix factor in output

#' Find genes with high spatial autocorrelation
#'
#' This function can be used to find genes with spatial structure in ST datasets.
#' A more detailed decription of the algorithm is outlined in the Details section below.
#'
#' overview of method:
#' \itemize{
#'    \item{Build a connection network from the array x,y coordinates for each sample. For a 'Visium' array, this would typically be 6 neighbours
#'    because of the hexagonal structure of spots.}
#'    \item{Combine connection networks from multiple samples}
#'    \item{Compute the lag vector for each feature}
#'    \item{Compute the correlation between the lag vector and the original vector}
#' }
#' The connection network is build by defining edges between each spot and its `nNeighborurs` closest
#' neighbours that are within a maximum distance defined by `maxdist`. This is to make sure that spots
#' along the tissue edges or holes have the correct number of neighbours. A connection network is built for
#' each section separately but they are then combined into one large connection network so that the
#' autocorrelation can be computed for the whole dataset.
#'
#' Now that we have a neighbour group defined for
#' each spot, we can calculate the lag vector for each feature. The lag vector of a features is essentially
#' the summed expression of that feature in the neighbour groups, computed for all spots and can be thought
#' of as a "smoothing" estimate.
#'
#' If we consider a spot A and its neighbours nbA, a feature with high spatial
#' corelation should have similar expression levels in both groups. We can therefore compute the a
#' correlation score between the lag vector and the "normal" expression vector to get an estimate of
#' the spatial autocorrelation.
#'
#' @param object Seurat object
#' @param assay Name of assay the function is being run on
#' @param slot Slot to use as input [default: 'scale.data']
#' @param features Features to rank by spatial autocorrelation. If no features are provided, the
#' features are selected using the `VariableFeatures` function in Seurat, meaning that the top variable genes
#' will be used.
#' @param nNeighbours Number of neighbours to find for each spot, For Visium data, this parameter is set to
#' 6 because of the spots are arranged in a hexagonal pattern and should have maximum 6 neighbors.
#' @param maxdist Maximum allowed distance to define neighbouring spots [default: 1.5]. If not provided, a
#' maximum distance is automatically selected depending on the platform. For Visium data, this maximum distance
#' is set to 150 microns.
#'
#' @return data.frame with gene names and correlation scores
#'
#' @importFrom Matrix bdiag
#' @importFrom spdep mat2listw
#'
#' @export
#'
CorSpatialGenes <- function (
  object,
  assay = NULL,
  slot = 'scale.data',
  features = NULL,
  nNeighbours = NULL,
  maxdist = NULL
) {

  # Check if adespatial is installed
  #if (!requireNamespace("adespatial")) stop("R package adespatial is required to run this function ... \n")
  #if (!requireNamespace("spdep")) stop("R package spdep is required to run this function ... \n")

  # Collect Staffli object
  if (!"Staffli" %in% names(object@tools)) stop("There is no 'Staffli' object present in this 'Seurat' object ...", call. = FALSE)
  st.object <- GetStaffli(object)

  # Obtain data
  features <- features %||% VariableFeatures(object)
  assay <- assay %||% DefaultAssay(object)
  data.use <- GetAssayData(object, slot = slot, assay = assay)
  if (length(data.use) == 0) stop(paste0(slot, " is empty in assay ", assay, ". Please normalise the data first or chose another slot  ... \n"), call. = FALSE)
  data.use <- data.use[features, ]

  # Create a combined network for the samples
  CN <- do.call(rbind, GetSpatNet(object = object, nNeighbours = nNeighbours, maxdist = maxdist))
  resCN <- as.matrix(data.frame(reshape2::dcast(CN, formula = from ~ to, value.var = "distance", fill = 0), row.names = 1))
  resCN[resCN > 0] <- 1
  empty.CN <- matrix(0, nrow = ncol(data.use), ncol = ncol(data.use), dimnames = list(colnames(data.use), colnames(data.use)))
  colnames(resCN) <- gsub(pattern = "\\.", replacement = "-", x = colnames(resCN))
  colnames(resCN) <- gsub(pattern = "^X", replacement = "", x = colnames(resCN))
  empty.CN[rownames(resCN), colnames(resCN)] <- resCN
  listw <- mat2listw(empty.CN)
  fun <- function (x) lag.listw(listw, x, TRUE)

  # Calculate the lag matrix from the network
  tablag <- apply(t(data.use), 2, fun)

  sp.cor <- unlist(lapply(1:nrow(data.use), function(i) {
    cor(data.use[i, ], tablag[, i])
  }))

  res <- data.frame(gene = rownames(data.use), cor = sp.cor, stringsAsFactors = F)
  res <- res[order(sp.cor, decreasing = T), ]
  rownames(res) <- res$gene
  return(res)
}

