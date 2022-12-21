#' @import Matrix
"_PACKAGE"

#' load count files
#'
#' @param path Path to count files
#' @param delim Delimiter used in gene count matrix file
#' @param row.names set column as row.names
#' @param visium Load 10x visium output files
#'
#' @importFrom data.table fread
#'
#' @keywords internal

st.load.matrix = function (
  path,
  delim = "\t",
  row.names = 1,
  visium = F
) {
  stopifnot(file.exists(path))
  x = c()
  if(visium == F){ #Original ST loading
    tmp = suppressWarnings({try({x = data.frame(fread(input = path, integer64 = "character", sep = delim), row.names = row.names, check.names = FALSE)})})
  } else { #10X Visium loading
    suff <- getExtension(path)
    if (dir.exists(path) & suff == basename(path)) {
      tmp = try({x = Seurat::Read10X(data.dir = path)})
    } else {
      tmp = try({x = Seurat::Read10X_h5(filename = path)})
    }
  }
  if(inherits(tmp, 'try-error')) {
    stop("Failed to read data.")
  } else {
    return(as(as.matrix(x), "dgCMatrix"))
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
#' @param sparse.matrix.fmt return mergedexpression matrix in dgCMatrix format to save memory
#'
#' @return merged expression matrix

Merger <- function(
  exp.list,
  delim = "x",
  labels = NULL,
  sparse.matrix.fmt = F
) {

  # Check labels
  if (!is.null(labels)) {
    stopifnot(length(labels) == length(ls) & class(labels) == "character")
  }

  # List all genes
  genes <- Reduce(union, lapply(exp.list, rownames))

  matrix.list <- lapply(seq_along(exp.list), function(i) {
    count <- exp.list[[i]]
    nspots <- ncol(count)
    curgenes <- rownames(count)
    m <- matrix(0, nrow=length(genes), ncol=nspots)
    rownames(m) <- genes
    m[curgenes,] <- count
    colnames(m) <- paste(ifelse(!is.null(labels), labels[i], i), colnames(count), sep = delim)

    # Save in sparse format?
    if (sparse.matrix.fmt) {
      m <- as(m, "dgCMatrix")
    }
    return(m)
  })
  return(do.call(cbind, matrix.list))
}

#' Obtain file extension
#'
#' @param file Path to file
#'
getExtension <- function (
  file
){
  ex <- strsplit(basename(file), split = "\\.")[[1]]
  if (ex[length(ex)] == "gz") suff <- paste(ex[length(ex) - 1], ex[length(ex)], sep = ".") else suff <- ex[length(ex)]
  return(suff)
}
