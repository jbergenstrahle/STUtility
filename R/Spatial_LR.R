#' Detect ligand-receptor interactions
#'
#' This function can be used to find ligand-receptor pairs that are expressed
#' in the same area of the tissue.
#'
#' @param object Seurat object
#' @param LR A data.frame with ligand-receptor pairs. Columns should be named
#' "ligand" and "receptor".
#' @param method "binarize" or "lee". See details for more information.
#' @param smooth.data Logical indicating whether or not the data should be smoothened
#' or not. Smoothing of expression data is used by averaging the expression across
#' neighboring spots.
#' @param nNeighbors Number of neighbors to use for smoothing. More neighbors will
#' lead to more smoothing of the data.
#' @param return.interaction.assay Returns an assay with the interactions.
#' @param do.parallel Run the function on multiple cores. The number of cores can be
#' set by the MC_CORES envornment variable, e.g:
#' `Sys.setenv("MC_CORES" = 7L)`
#' @param verbose Print messages
#'
#' @return data.frame with ligand-receptor interactions
#'
#' @importFrom spdep lee
#' @importFrom pbmcapply pbmclapply
#' @importFrom parallel detectCores
#' @importFrom stats fisher.test p.adjust
#' @importFrom zeallot %<-%
#'
SpatialLR <- function (
  object,
  LR,
  method = "binarize",
  slot = "data",
  smooth.data = TRUE,
  group.by = NULL,
  nNeighbors = 20,
  return.interaction.assay = FALSE,
  return.smooth.assay = FALSE,
  do.parallel = TRUE,
  verbose = TRUE
) {

  # Check column names
  if (!all(c("ligand", "receptor") %in% colnames(LR))) stop("Invalid LRdb ... \n")

  # Get expression matrix
  data.use <- as.matrix(GetAssayData(object, slot = slot))

  # Specify lapply function
  if (verbose & do.parallel) {
    n.cores <- ifelse(Sys.getenv("MC_CORES") == "", 2, Sys.getenv("MC_CORES"))
    cat("Running parallel using ", n.cores, " cores ... \n", sep = "", file = stderr())
  }
  fl <- ifelse(do.parallel, get("pbmclapply"), get("lapply"))

  # Create a combined network for the samples
  if (smooth.data | (method == "lee")) {
    if (verbose) cat("Constructing spatial network ... \n", file = stderr())
    CN <- do.call(rbind, GetSpatNet(object = object, nNeighbours = 20, maxdist = Inf))
    resCN <- as.matrix(data.frame(reshape2::dcast(CN, formula = from ~ to, value.var = "distance", fill = 0), row.names = 1))
    diag(resCN) <- 1
    resCN[resCN > 0] <- 1
    empty.CN <- matrix(0, nrow = ncol(data.use), ncol = ncol(data.use), dimnames = list(colnames(data.use), colnames(data.use)))
    colnames(resCN) <- gsub(pattern = "\\.", replacement = "-", x = colnames(resCN))
    colnames(resCN) <- gsub(pattern = "^X", replacement = "", x = colnames(resCN))
    empty.CN[rownames(resCN), colnames(resCN)] <- resCN
    listw <- mat2listw(empty.CN)
  }

  # Subset data
  data <- data.use[rowSums(data.use) > 0, ]
  shared.genes <- intersect(union(LR$ligand, LR$receptor), rownames(data))
  if (length(shared.genes) == 0) stop("No L/R pairs found in data ... \n")
  data <- data[shared.genes, ]
  spots <- colnames(data)
  LRdb.subset <- LR[(LR$ligand %in% shared.genes) & (LR$receptor %in% shared.genes), ]

  # Calculate the lag matrix from the network
  if (smooth.data) {
    if (verbose) cat("Spatial smoothing ... \n", file = stderr())
    fun <- function (x) lag.listw(listw, x, TRUE)
    data <- do.call(rbind, fl(seq_along(shared.genes), function(i) {
      fun(x = data[shared.genes[i], ])
    }))
    rownames(data) <- shared.genes
    colnames(data) <- spots
  }

  # Run test
  if (method == "binarize") {

    if (verbose) cat("Checking for ligand-receptor interactions using the 'binarize' method ... \n", file = stderr())
    dataLR <- apply(data, 1, binarizer)
    tests.LR <- do.call(rbind, fl(1:nrow(LRdb.subset), function(i) {
      L <- LRdb.subset[i, "ligand"]; R <- LRdb.subset[i, "receptor"]
      x <- dataLR[, L]
      y <- dataLR[, R]
      test = paste0(x, '-', y)

      if(length(unique(test)) < 4) {
        possibs = c("1-1", "0-1", "1-0", "0-0")
        missings_possibs = possibs[!possibs %in% unique(test)]
        test = c(test, missings_possibs)
      }

      fish_res = fisher.test(matrix(table(test), byrow = T, nrow = 2))[c('p.value', 'estimate')]
      return(data.frame(R, L,
                        p_value = fish_res$p.value,
                        odds_ratio = fish_res$estimate))
    }))
    tests.LR <- tests.LR[!duplicated(tests.LR[, 1:2]), ]
    if (verbose) cat(nrow(tests.LR), "interactions returned ... \n", file = stderr())
    tests.LR$adj_p_val <- p.adjust(tests.LR$p_value)
    tests.LR <- tests.LR %>% arrange(adj_p_val)

  } else if (method == "lee") {

    if (verbose) cat("Checking for ligand-receptor interactions using the 'lee' method ... \n", file = stderr())
    if (verbose) cat("Running test for all L/R pairs ... \n", file = stderr())
    tests <- fl(1:nrow(LRdb.subset), function(i) {
      sp <- lee(x = data[LRdb.subset$ligand[i], ],
                y = data[LRdb.subset$receptor[i], ],
                listw = listw, n = ncol(object),
                zero.policy = TRUE)
      v <- sp$localL
      d <- data.frame(R = LRdb.subset$receptor[i], L = LRdb.subset$ligand[i], lee = sp$L)
      return(list(d, v))
    })
    tests.LR <- do.call(rbind, lapply(tests, function(ls) {
      c(d, v) %<-% ls
      return(d)
    }))
    tests.LR <- tests.LR[!duplicated(tests.LR[, 1:2]), ]
    if (verbose) cat(nrow(tests.LR), "interactions returned ... \n", file = stderr())
    tests.LR <- tests.LR %>% arrange(-lee)

  }

  # Calculate enrichment score
  if (!is.null(group.by)) {
    if (!group.by %in% colnames(object[[]])) stop(paste0("'", group.by, "' not present in meta.data slot ... \n"))
    groups <- object[[]][, group.by]
    if (!class(groups) %in% c("character", "factor")) stop(paste0("Invalid class for column '", group.by, "' ... \n"))

    if (class(groups) == "factor") {
      lvls <- levels(groups)
    } else {
      lvls <- unique(groups)
    }

    # Get expression matrix
    data.use <- as.matrix(GetAssayData(object, slot = slot))
    #data.use <- data
    pairs <- tests.LR[, c("L", "R")]
    pairs <- pairs[!duplicated(pairs), ]
    data.use <- GetAssayData(object, slot = "data")[union(pairs$L, pairs$R), ]
    all.scores <- do.call(rbind, lapply(lvls, function(cl) {
      d1 <- data.use[, groups %in% cl]
      d2 <- data.use[, !groups %in% cl]
      scores <- ((rowMeans(d1) + 1)/(rowMeans(d2) + 1))*((rowMeans(d1 > 0) + 1)/(rowMeans(d2 > 0) + 1))
      cl.scores <- data.frame(LR.pair = paste0(pairs$L, "-", pairs$R),
                              score = scores[pairs$L]*scores[pairs$R],
                              cluster = cl) %>%
        arrange(-score)
    }))
    all.scores <- all.scores %>%
      group_by(LR.pair) %>%
      top_n(wt = score, n = 1)
    all.scores <- setNames(all.scores, nm = c("LR_pair", "enrichment_score", group.by))
    tests.LR$LR_pair <- paste0(tests.LR$L, "-", tests.LR$R)
    tests.LR <- merge(tests.LR, all.scores, by = "LR_pair")
  }

  if (return.interaction.assay | return.smooth.assay) {
    if (return.interaction.assay) {
      if (verbose) cat("Computing interaction matrix ... \n", file = stderr())
      interMat <- do.call(rbind, lapply(1:nrow(tests.LR), function(i) {
        data[tests.LR$R[i], ]*data[tests.LR$L[i], ]
      }))
      rownames(interMat) <- paste0(tests.LR$L, "-", tests.LR$R)
      int.assay <- CreateAssayObject(data = interMat)
      object@assays[["int"]] <- int.assay
    }
    if (return.smooth.assay) {
      smooth.assay <- CreateAssayObject(data = data)
      object@assays[["smooth"]] <- smooth.assay
    }
    return(list(tests.LR, object))
  } else {
    return(tests.LR)
  }
}



binarizer <- function (
  x,
  seed = 123
) {
  set.seed(seed)
  a <- as.logical(kmeans(x, 2)$cluster - 1)
  ma <- mean(x[a])
  mb <- mean(x[!a])
  res <- ifelse(rep(ma > mb, length(a)), a, !a)
  return(as.numeric(res))
}
