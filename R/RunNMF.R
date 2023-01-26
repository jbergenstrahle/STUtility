#' @include generics.R Staffli.R
NULL


#' Plot a bar chart with gene loadings
#'
#' This function takes a `Seurat`object as input with a 'NMF' reduction slot
#' and plots a bar chart of gene weights (loadings) for the top most driving genes.
#'
#' @param object Seurat object
#' @param factor Factor to plot gene loadings for [default: 1]
#' @param topn Top genes to show [default: 20]
#' @param dark.theme Use a dark theme for the plot
#' 
#' @importFrom ggplot2 geom_bar coord_flip 
#' @importFrom Seurat DarkTheme
#' @importFrom stats reorder
#'
#' @export
#'
FactorGeneLoadingPlot <- function (
  object,
  factor = 1,
  topn = 20,
  dark.theme = FALSE
) {

  # require NNLM
  if (!requireNamespace("NNLM")) {
    devtools::install_github('linxihui/NNLM')
  }
  
  # Check if NMF has been computed
  if (!"NMF" %in% names(object@reductions)) stop("NMF has not been computed ... \n", call. = FALSE)

  ftr <- paste0("factor_", factor)
  nmf <- object@reductions$NMF@feature.loadings[, ftr]
  gene <- names(nmf)
  df <- data.frame(gene, val = nmf, stringsAsFactors = F)
  df <- df[order(df$val, decreasing = T), ]
  df <- df[1:topn, ]
  df$gene <- factor(df$gene, levels = df$gene)
  p <- ggplot(df[1:topn, ], aes(reorder(gene, val), val)) +
    geom_bar(stat = "identity", fill = ifelse(dark.theme, "dark gray", "lightgray"), color = ifelse(dark.theme, "lightgray", "black"), width = 0.7) +
    coord_flip() +
    labs(x = "gene", y = "value")

  if (dark.theme) p <- p + DarkTheme() else p <- p + theme_minimal()

  return(p)
}


#' Run Non-negative Matrix Factorization
#'
#' Decompose an expression matrix A with non-negative elements into matrices WxH, also with
#' non-negative elements. W is the feature loading matrix (features x factors) and H is the
#' low dimensional embedding of the spots (factors x spots).
#'
#' @param object Seurat object
#' @param assay Assay Name of Assay NMF is being run on
#' @param slot Slot to pull data from.
#' @param features Features to compute the NMF for. Note that these features must be present in the
#' slot used to compute the NMF. By default, the `features` is set to `VariableFeatures(object)`
#' to include the most variable features selected in the normalization step.
#' @param nfactors Total Number of factors to compute and store (20 by default)
#' @param rescale Rescale data to make sure that values of the input matrix are non-n
#' @param reduction.name Dimensional reduction name, "NMF" by default
#' @param reduction.key Dimensional reduction key, specifies the prefix of the factor ids, e.g.
#' "factor_1", "factor_2", etc.
#' @param n.cores Number of threads to use in computation
#' @param order.by.spcor Order factors by spatial correlation
#' @param sort.spcor.by.var Sort factors by decreasing variance
#' @param ... Additional parameters
#'
#' @importFrom parallel detectCores
#' @importFrom Seurat CreateDimReducObject DefaultAssay VariableFeatures GetAssayData
#'
#' @export
#'
RunNMF <- function (
  object,
  assay = NULL,
  slot = "scale.data",
  features = NULL,
  nfactors = 20,
  rescale = TRUE,
  reduction.name = "NMF",
  reduction.key = "factor_",
  n.cores = NULL,
  order.by.spcor = FALSE,
  sort.spcor.by.var = FALSE,
  ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  var.genes <- features %||% VariableFeatures(object)
  norm.counts <- GetAssayData(object, slot = slot, assay = assay)
  if (rescale) {
    norm.counts <- t(apply(norm.counts, 1, function(x) (x - min(x))/(max(x) - min(x))))
  }
  if (min(norm.counts) < 0) stop("Negative values are not allowed")
  nmf.results <- rnmf(A = norm.counts[var.genes, ], k = nfactors)
  #nmf.results$W <- swne::ProjectFeatures(norm.counts, nmf.results$H, n.cores = n.cores)
  feature.loadings <- nmf.results$W
  cell.embeddings <- t(nmf.results$H)

  # Set cores
  n.cores <- n.cores %||% {detectCores() - 1}

  # Order factors based on spatial correlation
  if (order.by.spcor) {
    CN <- do.call(rbind, GetSpatNet(object = object, nNeighbours = NULL, maxdist = NULL))
    resCN <- as.matrix(data.frame(reshape2::dcast(CN, formula = from ~ to, value.var = "distance", fill = 0), row.names = 1))
    resCN[resCN > 0] <- 1
    empty.CN <- matrix(0, nrow = nrow(cell.embeddings), ncol = nrow(cell.embeddings), dimnames = list(rownames(cell.embeddings), rownames(cell.embeddings)))
    colnames(resCN) <- gsub(pattern = "\\.", replacement = "-", x = colnames(resCN))
    colnames(resCN) <- gsub(pattern = "^X", replacement = "", x = colnames(resCN))
    empty.CN[rownames(resCN), colnames(resCN)] <- resCN
    listw <- spdep::mat2listw(empty.CN)
    fun <- function (x) spdep::lag.listw(listw, x, TRUE)

    # Calculate the lag matrix from the network
    tablag <- apply(cell.embeddings, 2, fun)

    # Split sp.cor by sample
    if (sort.spcor.by.var) {
      sp.cor.split <- do.call(rbind, lapply(unique(GetStaffli(object)@meta.data$sample), function(s) {
        tablag.split <- tablag[GetStaffli(object)@meta.data$sample == s, ]
        cell.embeddings.split <- cell.embeddings[GetStaffli(object)@meta.data$sample == s, ]
        unlist(lapply(1:ncol(cell.embeddings.split), function(i) {
          cor(tablag.split[, i], cell.embeddings.split[, i])
        }))
      }))
      order.vec <- order(apply(sp.cor.split, 2, var))
    } else {
      sp.cor <- unlist(lapply(1:ncol(cell.embeddings), function(i) {
        cor(cell.embeddings[, i], tablag[, i])
      }))
      order.vec <- order(sp.cor, decreasing = TRUE)
    }

    cell.embeddings <- cell.embeddings[, order.vec]
    colnames(cell.embeddings) <- paste0(reduction.key, 1:ncol(cell.embeddings))
  }

  rownames(x = feature.loadings) <- var.genes
  colnames(x = feature.loadings) <- paste0(reduction.key, 1:nfactors)
  rownames(x = cell.embeddings) <- colnames(x = object)
  colnames(x = cell.embeddings) <- colnames(x = feature.loadings)
  reduction.data <- CreateDimReducObject (
    embeddings = cell.embeddings,
    loadings = feature.loadings,
    assay = assay,
    key = reduction.key
  )
  object[[reduction.name]] <- reduction.data
  return(object)
}


#' Run NMF with ica init
#' 
#' @param A Expression matrix
#' @param k = Number of factors
#' @param alpha Alpha parameter value
#' @param init Method to initiate NMF by
#' @param n.cores Number of threads to run
#' @param loss Loss method
#' @param max.iter Maximum number of iterations for NMF
#' @param ica.fast Should fast ica be run?
#' 
#' @importFrom stats runif
#' 
rnmf <- function (
  A,
  k,
  alpha = 0,
  init = "ica",
  n.cores = 1,
  loss = "mse",
  max.iter = 500,
  ica.fast = F
) {
  if (any(A < 0))
    stop("The input matrix contains negative elements !")
  if (k < 3)
    stop("k must be greater than or equal to 3 to create a viable SWNE plot")
  if (!init %in% c("ica", "nnsvd", "random")) {
    stop("Invalid initialization method")
  }
  A <- as.matrix(A)
  if (any(A < 0)) {
    stop("Input matrix has negative values")
  }
  if (init == "ica") {
    nmf.init <- ica_init(A, k, ica.fast = ica.fast)
  }
  else if (init == "nnsvd") {
    nmf.init <- nnsvd_init(A, k, LINPACK = T)
  }
  else {
    nmf.init <- NULL
  }
  if (is.null(nmf.init)) {
    nmf.res <- NNLM::nnmf(A, k = k, alpha = alpha, n.threads = n.cores,
                          loss = loss, max.iter = max.iter)
  }
  else {
    A.mean <- mean(A)
    zero.eps <- 1e-06
    nmf.init$W[nmf.init$W < zero.eps] <- 0
    nmf.init$H[nmf.init$H < zero.eps] <- 0
    zero.idx.w <- which(nmf.init$W == 0)
    zero.idx.h <- which(nmf.init$H == 0)
    nmf.init$W[zero.idx.w] <- runif(length(zero.idx.w), 0,
                                    A.mean/100)
    nmf.init$H[zero.idx.h] <- runif(length(zero.idx.h), 0,
                                    A.mean/100)
    nmf.res <- NNLM::nnmf(A, k = k, alpha = alpha, init = nmf.init,
                          n.threads = n.cores, loss = loss, max.iter = max.iter)
  }
  colnames(nmf.res$W) <- rownames(nmf.res$H) <- sapply(1:ncol(nmf.res$W),
                                                       function(i) paste("factor", i, sep = "_"))
  return(nmf.res)
}



#' Summarize features associated with cselected factors
#'
#' Extracts the top driving features per factor and returns
#'
#' @param object Seurat object
#' @param dims Factors to use
#' @param features.return Number of features to return per factor
#' @param features.use Select features (genes) to subset the data on
#'
#' @export
#'
SummarizeAssocFeatures <- function (
  object,
  dims = NULL,
  features.return = 10,
  features.use = NULL
) {

  if (!"NMF" %in% names(object@reductions)) stop(paste0("No factors available in Seurat object. Run RunNMF() first "), call. = FALSE)

  feature.factor.assoc <- object@reductions[["NMF"]]@feature.loadings
  if (!is.null(features.use)) {
    feature.factor.assoc <- feature.factor.assoc[features.use, ]
  }
  if (!is.null(dims)) {
    feature.factor.assoc <- feature.factor.assoc[, dims]
  }
  factor.features.df <- do.call("rbind", lapply(1:ncol(feature.factor.assoc),
                                                function(i) {
                                                  features.df <- data.frame(assoc_score = feature.factor.assoc[, i])
                                                  features.df$feature <- rownames(feature.factor.assoc)
                                                  features.df$factor <- colnames(feature.factor.assoc)[[i]]
                                                  features.df <- features.df[order(features.df$assoc_score, decreasing = T), ]
                                                  head(features.df, n = features.return)
                                                }))
  rownames(factor.features.df) <- NULL
  gene.loadings.selected <- feature.factor.assoc[unique(factor.features.df$feature), ]
  return(list(factor.features.df, gene.loadings.selected))
}


# Set a default value if an object is null
#
# @param lhs An object to set if it's null
# @param rhs The value to provide if x is null
#
# @return rhs if lhs is null, else lhs
#
# @author Hadley Wickham
# @references https://adv-r.hadley.nz/functions.html#missing-arguments
#
`%||%` <- function(lhs, rhs) {
  if (!is.null(x = lhs)) {
    return(lhs)
  } else {
    return(rhs)
  }
}

#' Initiate NMF using ICA
#'
#' @param A input matrix
#' @param k number of components to compute
#' @param ica.fast Should a fast implementation of ICA be used?
#'

ica_init <- function (
    A, 
    k, 
    ica.fast = F
) {
  
  if (!requireNamespace("irlba")) install.packages("irlba")
  if (!requireNamespace("ica")) install.packages("ica")
  
  if (ica.fast) {
    pc.res.h <- irlba::irlba(t(A), nv = 50, maxit = 100,
                             center = rowMeans(A))
    ica.res.h <- ica::icafast(pc.res.h$u, nc = k, maxit = 25,
                              tol = 1e-04)
    return(list(W = (A - Matrix::rowMeans(A)) %*% ica.res.h$S,
                H = t(ica.res.h$S)))
  }
  else {
    ica.res <- ica::icafast(t(A), nc = k, maxit = 25, tol = 1e-04)
    return(list(W = ica.res$M, H = t(ica.res$S)))
  }
}

nnsvd_init <- function (A, k, LINPACK)
{
  size <- dim(A)
  m <- size[1]
  n <- size[2]
  W <- matrix(0, m, k)
  H <- matrix(0, k, n)
  s = svd(A, k, k, LINPACK = LINPACK)
  U <- s$u
  S <- s$d
  V <- s$v
  W[, 1] = sqrt(S[1]) * abs(U[, 1])
  H[1, ] = sqrt(S[1]) * abs(t(V[, 1]))
  for (i in seq(2, k)) {
    uu = U[, i]
    vv = V[, i]
    uup = .pos(uu)
    uun = .neg(uu)
    vvp = .pos(vv)
    vvn = .neg(vv)
    n_uup = .norm(uup)
    n_vvp = .norm(vvp)
    n_uun = .norm(uun)
    n_vvn = .norm(vvn)
    termp = n_uup %*% n_vvp
    termn = n_uun %*% n_vvn
    if (termp >= termn) {
      W[, i] = sqrt(S[i] * termp) * uup/n_uup
      H[i, ] = sqrt(S[i] * termp) * vvp/n_vvp
    }
    else {
      W[, i] = sqrt(S[i] * termn) * uun/n_uun
      H[i, ] = sqrt(S[i] * termn) * vvn/n_vvn
    }
  }
  return(list(W = W, H = H))
}

