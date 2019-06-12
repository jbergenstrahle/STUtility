#' Extract table of x, y coordinates from ST count matrix headers
#'
#' @param spotnames x, y coordinates separated by delimiter (default "x")
#' @param delim delimiter
#' @return data.frame with x, y coordinates and optionally a sampleID column

GetCoords <- function(spotnames, delim = "x") {
  stopifnot(class(spotnames) == "character")
  if (sum(duplicated(spotnames)) > 0) {
    stop("Duplicate names are not allowed ...")
  }
  coords <- do.call(rbind, strsplit(spotnames, split = delim))
  if (ncol(coords) == 3) {
    coords <- setNames(data.frame(coords, stringsAsFactors = F), nm = c("x", "y", "sampleID"))
  } else if (ncol(coords) == 2) {
    coords <- setNames(data.frame(coords, stringsAsFactors = F), nm = c("x", "y"))
  } else {
    stop("Invalid spotnames ...")
  }
  coords$x <- as.numeric(coords$x)
  coords$y <- as.numeric(coords$y)
  rownames(coords) <- spotnames
  return(coords)
}
