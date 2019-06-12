#' @import dplyr magrittr
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
  qc.table <- cbind(qc.table, data.frame(nTranscripts = colSums(exprMat),
                                         nGenes = apply(exprMat, 2, function(x) sum(x > 0))))
  qc.table$nTranscripts_per_gene <- qc.table$nTranscripts/qc.table$nGenes
  return(qc.table)
}


#' Quality metrics for genes
#'
#' @param exprMat expression matrix with spots in columns and genes in rows
#' @param delim delimiter
#' @param labels character vector with labels corresponding to dataset IDs in merged expression matrix [optional]
#' @return data.frame with gene names and quality metrics

qcGenes <- function(
  exprMat,
  delim = "x",
  labels = NULL
) {
  coords <- GetCoords(colnames(exprMat), delim)
  qc.table <- data.frame(gene = rownames(exprMat))
  if ("sampleID" %in% colnames(coords)) {
    samples <- unique(coords$sampleID)
    for (i in seq_along(samples)) {
      qc.table[, ifelse(!is.null(labels), labels[i], samples[i])] <- rowSums(exprMat[, coords$sampleID == samples[i]])
    }
  } else {
    qc.table[, "count"] <- rowSums(exprMat)
  }
  return(qc.table)
}


#' Quality metrics for samples
#'
#' @param exprMat expression matrix with spots in columns and genes in rows
#' @param delim delimiter
#' @param labels character vector with labels corresponding to dataset IDs in merged expression matrix [optional]
#' @return data.frame with quality metrics per sample

qcSamples <- function(
  exprMat,
  delim = "x",
  labels = NULL
) {
  qc.table <- qcSpots(exprMat, delim)
  summarize_fkn <- function(grouped_df) {
    grouped_df %>% summarize(avg.nTranscripts = round(mean(nTranscripts), digits = 2),
                             max.nTranscripts = max(nTranscripts),
                             min.nTranscripts = min(nTranscripts),
                             avg.nGenes = round(mean(nGenes), digits = 2),
                             max.nGenes = max(nGenes),
                             min.nGenes = min(nGenes)) %>% as.data.frame()
  }
  if ("sampleID" %in% colnames(qc.table)) {
    qc.table <- qc.table %>%
      group_by(sampleID) %>%
      summarize_fkn
  } else {
    qc.table <- qc.table %>%
      summarize_fkn
  }
  return(qc.table)
}
