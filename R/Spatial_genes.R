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

  # spatial information
  xys = st.object@meta.data[, c("x", "y", "sample")]

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

    nNeighbours <- nNeighbours %||% ifelse(platforms[i] == "Visium", 6, 4)
    maxdist <- maxdist %||% ifelse(platforms[i] == "1k", 2, 1.5)

    knn_spatial <- dbscan::kNN(x = xys[, c("x", "y")] %>% as.matrix(), k = nNeighbours)
    knn_spatial.norm <- data.frame(from = rep(1:nrow(knn_spatial$id), nNeighbours),
                                 to = as.vector(knn_spatial$id),
                                 weight = 1/(1 + as.vector(knn_spatial$dist)),
                                 distance = as.vector(knn_spatial$dist))
    nw_spatial.norm = igraph::graph_from_data_frame(knn_spatial.norm, directed = FALSE)

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


#' Find genes with spatial structure in multiple datasets
#'
#' Rapid computation of genes that are spatially clustered in multiple samples.
#'
#' @section Acknowledgement:
#' The function is based on the `binGetSpatialGenes` function from the `Giotto` R package,
#' but have been modified to work with `Seurat` objects and with multiple samples.
#'
#' @param object Seurat object
#' @param spots Character vector with spot IDs to plot [default: all spots]
#' @param slot Which slot to pull the data from? [default: 'data']
#' @param features.include Features to include in the detection computation. By default,
#' most variable features are selected.
#' @param nstart kmeans: nstart parameter
#' @param iter_max kmeans: iter.max parameter
#' @param rank_percentage Percentage of top spots for binarization
#' @param do_fisher_test Perform Fisher test
#' @param community_expectation Degree expectation in spatial communities
#' @param verbose Print messages
#'
#' @return List of data.tables with results (see details)
#'
#' @details Two ways are provided to identify spatial genes based on gene expression binarization.
#' Both methods are identicial except for how binarization is performed.
#' \itemize{
#'   \item{1. binarize: }{Each gene is binarized (0 or 1) in each cspot with \bold{kmeans} (k = 2) or based on \bold{rank} percentile}
#'   \item{2. network: }{Alll spots are connected through a k-nearest neighbor network}
#'   \item{3. contingency table: }{A contingency table is calculated based on all pairwise spot-spot interactions (0-0, 0-1, 1-0 or 1-1)}
#'   \item{4. For each gene an odds-ratio (OR) and fisher.test (optional) is calculated}
#' }
#' Additionally 2 other statistics are provided:
#' \itemize{
#'   \item{Number of spots with high expression (binary = 1)}
#'   \item{total and ratio of highly connected spots: }{ Spots with a connectivity higher than community_expectation}
#' }
#' By selecting a subset of likely spatial genes (e.g. highly variable genes) the function will be much faster.
#'
#' @importFrom magrittr %>%
#' @importFrom zeallot %<-%
#' @importFrom data.table as.data.table dcast.data.table
#' @importFrom reshape2 melt
#' @importFrom dplyr summarize group_by
#'
#' @export
#'
#' @examples
#' spatgenes <- BinSpatialGenes(se)
#'
#' @return list of data.frames with spatial scores for each gene

BinSpatialGenes = function (
  object,
  spots = NULL,
  bin_method = c('kmeans', 'rank'),
  slot = 'data',
  features.include = NULL,
  nstart = 3,
  iter_max = 10,
  percentage_rank = 10,
  do_fisher_test = F,
  community_expectation = 5,
  verbose = FALSE
) {

  # Check to see if Staffli object is present
  if (!"Staffli" %in% names(object@tools)) stop("Staffli object is missing from Seurat object... \n", call. = FALSE)
  st.object <- object@tools$Staffli

  # Collect data
  spots <- spots %||% colnames(x = object)
  features.include <- features.include %||% object@assays[[DefaultAssay(object)]]@var.features
  if (length( x = features.include) == 0) stop("No variable features found in the Seurat object ... \n")
  data <- GetAssayData(object = object, slot = slot)
  data <- data[features.include, spots]

  # set binarization method
  bin_method = match.arg(bin_method, choices = c('kmeans', 'rank'))

  # spatial network
  spatial_networks = GetSpatNet(object)

  # binarize matrix
  if(bin_method == 'kmeans') {
    bin_matrix = t(apply(X = data, MARGIN = 1, FUN = kmeans_binarize, nstart = nstart, iter.max = iter_max))
  } else if(bin_method == 'rank') {
    max_rank = (ncol(data)/100)*percentage_rank
    bin_matrix = t(apply(X = data, MARGIN = 1, FUN = function(x) rank_binarize(x = x, max_rank = max_rank)))
  }

  if (verbose) cat("Finished binarization of expression data ... \n")

  # extra info: average expression of high expression group
  sel_data = data * bin_matrix
  av_expr = apply(sel_data, MARGIN = 1, FUN = function(x) {
    mean(x[x > 0])
  })
  av_expr_df = data.frame(genes = names(av_expr), av_expr = av_expr, stringsAsFactors = FALSE) %>% arrange(genes)

  # dcast
  bin_matrix_df = as.data.table(melt(bin_matrix, varnames = c('genes', 'spots'), value.name = 'value'))
  nr_high_spots <- bin_matrix_df %>% group_by(genes) %>% subset(value == 1) %>% summarize(count = n()) %>% as.data.table()

  # extra info: nr of spots with high expression
  #nr_high_spots = bin_matrix_DT[, .N, by = .(genes, value)][value == 1]
  #nr_high_spots = nr_high_spots[,.(genes,N)]

  # Split binarized matrix
  bin_matrix_df.list <- list()
  for (s in unique(st.object@samplenames)) {
    bin_matrix_df.list[[s]] <- subset(bin_matrix_df, spots %in% rownames(subset(st.object@meta.data, sample == s)))
  }

  # combine binarized matrix with spatial network
  combined.list <- lapply(seq_along(bin_matrix_df.list), function(i) {
    spatial_network <- spatial_networks[[i]] %>% as.data.table(); bin_matrix_df <- bin_matrix_df.list[[i]] %>% as.data.table()
    spatial_network_min = spatial_network[, .(to, from)]
    spatial_network_min = data.table:::merge.data.table(x = spatial_network_min, by.x = 'to', y = bin_matrix_df, by.y = 'spots', allow.cartesian = T)
    setnames(spatial_network_min, 'value', 'to_value')
    spatial_network_min[bin_matrix_df, from_value := value, on = c(genes = 'genes', from = 'spots')]
    spatial_network_min[, comb := paste0(to_value,'-',from_value)]
    tablecomb = spatial_network_min[, .N, by = .(genes, comb)]
    setorder(tablecomb, genes, comb)
    dtablecomb = dcast.data.table(tablecomb, formula = genes ~ comb, value.var = 'N')

    ## fisher test or odds-ratio only ##

    if(do_fisher_test == TRUE) {
      dtablecomb = dtablecomb[, fish_function2(A = `0-0`, B = `0-1`, C = `1-0`, D = `1-1`), by = genes]
    } else {
      # OR only
      dtablecomb = dtablecomb[, OR_function2(A = `0-0`, B = `0-1`, C = `1-0`, D = `1-1`), by = genes]
      #dtablecomb[, OR := ((`0-0`* `1-1`)/(`0-1`*`1-0`)), by = genes]
    }
    return(list(dtablecomb,  spatial_network_min))
  })

  merged.DTs <- lapply(seq_along(combined.list), function(i) {
    c(dtablecomb, spatial_network_min) %<-% combined.list[[i]]
    ## estimate for community ##
    # create count table for individual spots for all conditions
    tospots = spatial_network_min[, .(to, genes, comb)]
    setnames(tospots, 'to', 'spots')
    fromspots = spatial_network_min[, .(from, genes, comb)]
    setnames(fromspots, 'from', 'spots')
    allspots = rbind(tospots, fromspots)
    counttable_spots = allspots[, .N, by = .(genes, comb, spots)]

    # uniq spots per combination (0-0, 1-1, ...)
    count_uniq_spots = counttable_spots[, length(unique(spots)), by = .(genes, comb)]
    setorder(count_uniq_spots, genes, comb)

    # spots with higher connectivity per combination
    count_comm_spots = counttable_spots[, sum(N >= community_expectation), by = .(genes, comb)]
    setorder(count_comm_spots, genes, comb)
    setnames(count_comm_spots, 'V1', 'comm')

    count_comm_spots[, total := count_uniq_spots$V1]
    count_comm_spots[, ratio := round(comm/total, 2)]
    count_comm_spots = count_comm_spots[comb == '1-1']

    # merge different information
    mergeDT = merge(av_expr_df, nr_high_spots, by = 'genes') %>% as.data.table()
    mergeDT = merge(mergeDT, dtablecomb, by = 'genes')
    mergeDT = merge(mergeDT, count_comm_spots[,.(genes, ratio)], by = 'genes')

    mergeDT[, total_score := av_expr*ratio*log2(OR+1)]
    setorder(mergeDT, -total_score)

    return(mergeDT)
  })

  return(merged.DTs)
}



#' Find genes with spatial structure
#'
#' Runs a distance based method to detect genes with spatial structure.
#'
#' @section Acknowledgement:
#' The function is based on the `calculate_spatial_genes_python()` available from
#' the `Giotto` R package, but have been modified to work with Seurat objects
#' and on multiple samples.
#'
#' @param object Seurat object
#' @param spots Character vector with spot IDs to plot [default: all spots]
#' @param slot Which slot to pull the data from? [default: 'data']
#' @param metric The distance method used for `scipy.spatial.distance.pdist`
#' The distance function can be ‘braycurtis’, ‘canberra’, ‘chebyshev’, ‘cityblock’,
#' ‘correlation’, ‘cosine’, ‘dice’, ‘euclidean’, ‘hamming’, ‘jaccard’,
#' ‘jensenshannon’, ‘kulsinski’, ‘mahalanobis’, ‘matching’, ‘minkowski’,
#' ‘rogerstanimoto’, ‘russellrao’, ‘seuclidean’, ‘sokalmichener’,
#' ‘sokalsneath’, ‘sqeuclidean’ or ‘yule’.
#' @param features.include Features to include in the detection computation
#' @param rbp_p Fractional binarization method
#' @param examine_top Top fraction to evaluate with silhouette
#' @param python_path Specify specific path to python if required
#'
#' @export
#'
#' @examples
#' spatgenes <- SpatialGenes(se)
#'
#' @return list of data.frames with spatial scores for each gene

SpatialGenes <- function (
  object,
  spots = NULL,
  slot = "data",
  metric = "euclidean",
  features.include = NULL,
  rbp_p = 0.95,
  examine_top = 0.3,
  python_path = NULL
) {

  # Check that Giotto is installed
  if (!all(c("reticulate", "Giotto") %in% rownames(installed.packages()))) stop("R package Giotto or reticulate is missing ... \n")

  # Check to see if Staffli object is present
  if (!"Staffli" %in% names(object@tools)) stop("Staffli object is missing from Seurat object... \n", call. = FALSE)
  st.object <- object@tools$Staffli

  # Collect data
  spots <- spots %||% colnames(x = object)
  features.include <- features.include %||% object@assays[[DefaultAssay(object)]]@var.features
  data <- GetAssayData(object = object, slot = slot)
  data <- data[features.include, spots]

  # Spatial locations
  spatlocs <- st.object@meta.data[spots, c("x", "y", "sample")]

  # Split data and spatlocs
  data.list <- list()
  spatlocs.list <- list()
  for (s in unique(spatlocs$sample)) {
    data.list[[s]] <- data[, spatlocs$sample == s]
    spatlocs.list[[s]] <- spatlocs[spatlocs$sample == s, ]
  }

  # Python path
  python_path <- python_path %||% system("which python3", intern = TRUE)

  ## prepare python path and louvain script
  #reticulate::use_python(required = TRUE, python = python_path)
  python_leiden_function = system.file("python", "python_spatial_genes.py", package = 'Giotto')
  reticulate::source_python(file = python_leiden_function)

  cap <- reticulate::py_capture_output({
    output_python.list <- lapply(seq_along(data.list), function(i) {
      output_python <- python_spatial_genes(spatial_locations = spatlocs,
                          expression_matrix = as.data.frame(as.matrix(data)),
                          metric = metric,
                          rbp_p = rbp_p,
                          examine_top = examine_top)
    })
  })

  # unlist output
  spatial_python_df.list <- lapply(seq_along(output_python.list), function(i) {
    output_python <- output_python.list[[i]]
    genes = unlist(lapply(output_python, FUN = function(x) {
      y = x[1][[1]]
    }))
    scores = unlist(lapply(output_python, FUN = function(x) {
      y = x[2][[1]]
    }))
    spatial_python_df = data.frame(genes = genes, scores = scores, stringsAsFactors = FALSE)
  })
}


#' kmeans_binarize
#'
#' Create binarized scores using kmeans
#'
#' @param x Expression vector
#' @param nstart Number of random sets
#' @param iter.max Maximum number of allowed iterations
kmeans_binarize = function (
  x,
  nstart = 3,
  iter.max = 10
) {

  sel_gene_km = stats::kmeans(x, centers = 2, nstart = nstart, iter.max = iter.max)$cluster
  mean_1 = mean(x[sel_gene_km == 1])
  mean_2 = mean(x[sel_gene_km == 2])

  if(mean_1 > mean_2) {
    mean_1_value = 1
    mean_2_value = 0
  } else {
    mean_1_value = 0
    mean_2_value = 1
  }

  sel_gene_bin = x
  sel_gene_bin[sel_gene_km == 1] = mean_1_value
  sel_gene_bin[sel_gene_km == 2] = mean_2_value

  return(sel_gene_bin)

}


#' rank_binarize
#'
#' Create binarized scores using arbitrary rank of top genes
#'
#' @param x Expression vector
#' @param max_rank Maximum gene rank to consider
rank_binarize = function(x, max_rank = 200) {

  sel_gene_rank = rank(-x, ties.method = 'average')

  sel_gene_bin = x
  sel_gene_bin[sel_gene_rank <= max_rank] = 1
  sel_gene_bin[sel_gene_rank > max_rank] = 0

  return(sel_gene_bin)

}

#' fish_function
#'
#' Perform fisher exact test
fish_function = function(x_to, x_from) {

  fish_table = table(x_to == '1',
                     x_from == '1')

  fish_res = stats::fisher.test(fish_table)

  return(list(pval = fish_res$p.value, OR = fish_res$estimate))
}


#' fish_function2
#'
#' Perform fisher exact test
fish_function2 = function(A, B, C, D) {

  # set NA's to 0
  A = ifelse(is.na(A), 0, A)
  B = ifelse(is.na(B), 0, B)
  C = ifelse(is.na(C), 0, C)
  D = ifelse(is.na(D), 0, D)

  fish_matrix = matrix(c(A, B, C, D), nrow = 2)

  fish_res = stats::fisher.test(fish_matrix)

  return(list(pval = fish_res$p.value, OR = fish_res$estimate))
}


#' OR_function2
#'
#' Calculate odds-ratio
OR_function2 = function(A, B, C, D) {

  fish_matrix = matrix(c(A, B, C, D), nrow = 2)
  fish_matrix = fish_matrix/1000 # to prevent overflow

  OR = ((fish_matrix[1]*fish_matrix[4]) / (fish_matrix[2]*fish_matrix[3]))
  return(list(OR = OR))
}

