#' Extract table of x, y coordinates from ST count matrix headers
#'
#' @param spotnames x, y coordinates separated by delimiter (default "x")
#' @param delim delimiter
#' @return data.frame with x, y coordinates and optionally a sampleID column
#' @export

GetCoords <- function(
  spotnames,
  delim = "x"
) {
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


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


#' Get the names of objects within a Seurat object that are of a certain class
#'
#' @param object A Seurat object
#' @param classes.keep A vector of names of classes to get
#'
#' @return A vector with the names of objects within the Seurat object that are of class \code{classes.keep}

FilterObjects <- function(object, classes.keep = c('Assay', 'DimReduc')) {
  slots <- na.omit(object = Filter(
    f = function(x) {
      return(class(x = slot(object = object, name = x)) == 'list')
    },
    x = slotNames(x = object)
  ))
  slots <- grep(pattern = 'tools', x = slots, value = TRUE, invert = TRUE)
  slots <- grep(pattern = 'misc', x = slots, value = TRUE, invert = TRUE)
  slots.objects <- unlist(
    x = lapply(
      X = slots,
      FUN = function(x) {
        return(names(x = slot(object = object, name = x)))
      }
    ),
    use.names = FALSE
  )
  object.classes <- sapply(
    X = slots.objects,
    FUN = function(i) {
      return(class(x = object[[i]]))
    }
  )
  object.classes <- object.classes[object.classes %in% classes.keep]
  return(names(x = object.classes))
}
