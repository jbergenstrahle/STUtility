#' @importFrom data.table fread
"_PACKAGE"

#' load count files
#'
#' @param path Path to count files
#' @param delim delimiter
#' @param row.names set column as row.names
#' @keywords internal

st.load.matrix.2 = function(path, delim="\t", row.names=1, ...) {
  x = c()
  tmp = suppressWarnings({try({x = data.frame(fread(input = path, integer64 = "character",
                                   sep = delim),
                             row.names = row.names,
                             check.names = FALSE)})})
  if(inherits(tmp, 'try-error')) {
    return(as.matrix(c()))
  } else {
    return(as.matrix(x))
  }
}

#' Merge expression matrices
#'
#' @description This function takes a list of expression matrices as input with ST spots in columns
#' and genes in rows. All expression matrices will be merged into one matrix where all genes that are
#' present in any of the expression matrices will be included. If a there are no counts available for
#' a gene in one or more of the expression matrices, the counts for that gene will be set to 0 in those
#' expression matrices. By default, a unique number will be prefixed to the headers of each expression
#' matrix.
#' @param exp.list list of expression matrices
#' @param delim delimiter used to separate coordinates in expression matrix headers
#' @param labels labels to use as suffix for the headers of each expression matrix

expr.merger <- function(exp.list, delim = "_", labels = NULL) {
  # Check labels
  if (!is.null(labels)) {
    stopifnot(length(labels) == length(ls) & class(labels) == "character")
  }

  # collect all unique genes
  genes <- unique(unlist(lapply(exp.list, rownames)))

  # Obtain
  cols <- unlist(lapply(exp.list, colnames))

  merged_exprMat <- do.call(cbind, lapply(seq_along(exp.list), function(i) {
    exprMat <- exp.list[[i]]
    A <- as.data.frame(exprMat)
    A <- A[genes, ]
    rownames(A) <- genes
    A[is.na(A)] <- 0
    colnames(A) <- paste(ifelse(is.null(labels), i, labels[i]), colnames(A), sep = delim)
    as(as.matrix(A), "dgCMatrix")
  }))
  return(merged_exprMat)
}


