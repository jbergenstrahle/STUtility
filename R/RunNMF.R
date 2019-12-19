#' @include generics.R Staffli.R
NULL

#' Plot spatial heatmap of NMF factors
#'
#' @param object Seurat object
#' @param factors Specify factors to plot (integer vector)
#'
NMFPlot <- function (
  object,
  factors = 1:5
) {
  p.list <- lapply(paste0("factor_", factors), function(f) {
    data <- setNames(data.frame(object[[]], object@reductions[["NMF"]]@cell.embeddings[, f]), nm = c(colnames(object[[]]), f))
    ggplot(data, aes_string("ads_x", "36 - ads_y", color = f)) +
      geom_point(size = 0.3) +
      scale_x_continuous(limits = c(0, 33)) +
      scale_y_continuous(limits = c(0, 35)) +
      scale_color_gradientn(colours = c("black", "dark blue", "cyan", "yellow", "red", "dark red")) +
      facet_wrap(~ids, ncol = 3) +
      theme_void() +
      labs(title = f) +
      theme(plot.background = element_rect(fill = "black"),
            legend.text = element_text(colour = "white"),
            legend.title = element_text(colour = "white"),
            strip.text = element_text(colour = "white"))
  })
  cowplot::plot_grid(plotlist = p.list, ncol = 1)
}

#' Plot density distribution of factor values
#'
#' @param object Seurat object
#' @param factors Specify factors to display (default 1-5)
#' @param top.n.features Number of features to include in bar plot
#'
NMFRidgePlot <- function (
  object,
  factors = 1:5,
  top.n.features = 10
) {
  factors <- paste0("factor_", factors)
  data <- setNames(data.frame(object[[]], object@reductions[["NMF"]]@cell.embeddings[, factors]), nm = c(colnames(object[[]]), factors))
  datam <- reshape2::melt(data, measure.vars = factors)

  bw_theme <- theme(plot.background = element_rect(fill = "black", color = "black"),
                    panel.grid = element_blank(),
                    axis.text = element_text(colour = "white"),
                    strip.background = element_blank(),
                    rect = element_rect(fill = "black"),
                    panel.background = element_rect(fill = "black", color = "black"),
                    legend.text = element_text(color = "white"),
                    legend.title = element_text(color = "white"))

  p1 <- ggplot(datam, aes(value, ids, fill = factor(..quantile..))) +
    stat_density_ridges(
      geom = "density_ridges_gradient", calc_ecdf = TRUE,
      quantiles = 4, quantile_lines = TRUE
    ) +
    facet_wrap(~variable, scales = "free_x", ncol = 1) +
    viridis::scale_fill_viridis(option = "E", discrete = TRUE, name = "Quartiles") +
    labs(fill = "value") +
    bw_theme

  feature.loadings <- test@reductions[["NMF"]]@feature.loadings[, factors]
  genes <- rownames(feature.loadings)
  gg <- do.call(rbind, lapply(factors, function(f) {
    x <- feature.loadings[, f]
    genes.rank <- order(x, decreasing = T)[1:top.n.features]
    data.frame(gene = genes[genes.rank], value = x[genes.rank], variable = f)
  })) %>% arrange(variable, value) %>% mutate(order = row_number())

  p2 <- ggplot(gg, aes(x = order, y = value)) +
    geom_bar(stat = "identity", color = "black", alpha = .7) +
    coord_flip() +
    labs(x = "") +
    facet_wrap(~variable, scales = "free", ncol = 1) +
    theme_minimal() +
    scale_x_continuous(
      breaks = gg$order,
      labels = gg$gene,
      expand = c(0,0)
    ) +
    bw_theme

  cowplot::plot_grid(p1, p2)
}

#' Run Non-negative Matrix Factorization
#'
#' @param object Seurat object
#' @param assay assay Name of Assay PCA is being run on
#' @param features features to compute the NMF for
#' @param nfactors Total Number of factors to compute and store (20 by default)
#' @param nfactors.print factorss to print genes for
#' @param nfeatures.print Number of genes to print for each factor
#' @param reduction.name dimensional reduction name,  pca by default
#' @param reduction.key dimensional reduction key, specifies the string before
#' @param n.cores Number of threds to use in computation
#' @param order.by.spatcor Order factors by spatial correlation
#'
#' @export
#'
RunNMF <- function (
  object,
  assay = NULL,
  features = NULL,
  nfactors = 20,
  ndims.print = 1:5,
  nfeatures.print = 30,
  reduction.name = "NMF",
  reduction.key = "factor_",
  n.cores = 7,
  order.by.spcor = FALSE,
  ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  var.genes <- features %||% VariableFeatures(object)
  norm.counts <- GetAssayData(object, slot = "scale.data", assay = assay)
  norm.counts <- t(apply(norm.counts, 1, function(x) (x - min(x))/(max(x) - min(x))))
  nmf.results <- rnmf(A = norm.counts[var.genes, ], k = nfactors)
  #nmf.results$W <- swne::ProjectFeatures(norm.counts, nmf.results$H, n.cores = n.cores)
  feature.loadings <- nmf.results$W
  cell.embeddings <- t(nmf.results$H)

  # Order factors based on spatial correlation
  if (order.by.spcor) {
    st.object <- GetStaffli(object)
    m.list <- lapply(unique(st.object[[, "sample", drop = T]]), function(s) {
      resCN <- adespatial::chooseCN(xy = xy, ask = FALSE, type = 6, plot.nb = FALSE, edit.nb = FALSE, result.type = "listw", k = 6)
      m <- listw2mat(resCN)
    })
    nb <- Matrix::bdiag(m.list) %>% as.matrix()
    listw <- mat2listw(nb)
    fun <- function (x) lag.listw(listw, x, TRUE)
    tablag <- apply(cell.embeddings, 2, fun)
    sp.cor <- unlist(lapply(1:ncol(cell.embeddings), function(i) {
      cor(cell.embeddings[, i], tablag[, i])
    }))
    cell.embeddings <- cell.embeddings[, order(sp.cor, decreasing = TRUE)]
    feature.loadings <- feature.loadings[, order(sp.cor, decreasing = TRUE)]
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



#' Summarize features associated with certain factors
#'
#' @param object Seurat object
#' @param features.return Number of features to return per factor
#' @param features.use Select features to use for summary
#'
SummarizeAssocFeatures <- function (
  object,
  features.return = 10,
  features.use = NULL
) {

  if (!"NMF" %in% names(object@reductions)) stop(paste0("No factors available in Seurat object. Run RunNMF() first "), call. = FALSE)

  feature.factor.assoc <- object@reductions[["NMF"]]@feature.loadings
  if (!is.null(features.use)) {
    feature.factor.assoc <- feature.factor.assoc[features.use, ]
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


ica_init <- function (A, k, ica.fast = F)
{
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
