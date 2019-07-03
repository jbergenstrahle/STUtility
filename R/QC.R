#' @import dplyr magrittr
#' @import Matrix
"_PACKAGE"

#' Quality metrics for spots
#'
#' @param exprMat expression matrix with spots in columns and genes in rows
#' @param delim delimiter for expression matrix headers
#' @return data.frame with x, y coordinates and quality metrics

qcSpots <- function(
  exprMat,
  delim = "x"
) {
  qc.table <- GetCoords(colnames(exprMat), delim)
  qc.table <- cbind(qc.table, data.frame(nCount_RNA = colSums(exprMat),
                                         nFeature_RNA = apply(exprMat, 2, function(x) sum(x > 0))))
  qc.table$nCount_per_nFeature_RNA <- qc.table$nCount_RNA/qc.table$nFeature_RNA
  return(qc.table)
}


#' Quality metrics for genes
#'
#' @param object expression matrix with spots in columns and genes in rows or Seurat object
#' @param delim delimiter
#' @param labels character vector with labels corresponding to dataset IDs in merged expression matrix [optional]
#' @return data.frame with gene names and quality metrics

qcGenes <- function(
  object,
  delim = "x",
  labels = NULL
) {
  stopifnot(class(object) %in% c("matrix", "data.frame", "Seurat"))

  if (class(object) == "Seurat") {
    coords <- GetCoords(colnames(object), delim)
    exprMat <- GetAssayData(object = object, slot = 'counts')
    exprMat <- exprMat[rownames(object), colnames(object)]
  } else {
    coords <- GetCoords(colnames(object), delim)
    exprMat <- object
  }

  qc.table <- data.frame(gene = rownames(object))
  if ("sample" %in% colnames(coords)) {
    samples <- unique(coords$sample)
    for (i in seq_along(samples)) {
      qc.table[, ifelse(!is.null(labels), labels[i], samples[i])] <- rowSums(exprMat[, coords$sample == samples[i]])
    }
  } else {
    qc.table[, "count"] <- rowSums(exprMat)
  }
  return(qc.table)
}


#' Quality metrics for samples
#'
#' @param object expression matrix with spots in columns and genes in rows or Seurat object
#' @param delim delimiter
#' @param labels character vector with labels corresponding to dataset IDs in merged expression matrix [optional]
#' @return data.frame with quality metrics per sample

qcSamples <- function(
  object,
  delim = "x",
  labels = NULL
) {
  stopifnot(class(object) %in% c("matrix", "data.frame", "Seurat"))

  if (class(object) == "Seurat") {
    exprMat <- GetAssayData(object = object, slot = 'counts')
    exprMat <- exprMat[rownames(object), colnames(object)]
  } else {
    exprMat <- object
  }

  qc.table <- qcSpots(exprMat, delim)
  summarize_fkn <- function(grouped_df) {
    grouped_df %>% summarize(avg.nCount_RNA = round(mean(nCount_RNA), digits = 2),
                             max.nCount_RNA = max(nCount_RNA),
                             min.nCount_RNA = min(nCount_RNA),
                             avg.nFeature_RNA = round(mean(nFeature_RNA), digits = 2),
                             max.nFeature_RNA = max(nFeature_RNA),
                             min.nFeature_RNA = min(nFeature_RNA)) %>% as.data.frame()
  }
  if ("sample" %in% colnames(qc.table)) {
    qc.table <- qc.table %>%
      group_by(sample) %>%
      summarize_fkn
  } else {
    qc.table <- qc.table %>%
      summarize_fkn
  }
  return(qc.table)
}

#' function(x, y, ...) used to plot the contents of each panel of the display in pairs function
#'
#' @param x
#' @param y
#' @return panel for pairs function
#' @export internal

qc.scatter <- function(x,y){
    dns <- densCols(x,y);
    points(x,y, col=dns, pch=".", panel.first=grid());
    abline(a=0, b=1, col="brown")
}


