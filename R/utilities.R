#' Extract table of x, y coordinates from ST count matrix headers
#'
#' @param spotnames x, y coordinates separated by delimiter (default "x")
#' @param delim delimiter
#' @return data.frame with x, y coordinates and optionally a sample column
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
    coords <- setNames(data.frame(coords, stringsAsFactors = F), nm = c("x", "y", "sample"))
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

#' Palette selection
#'
#' @param palette Palette choice for plotting spatial expression histology heatmap
#' @importFrom viridis viridis cividis magma plasma
#' @importFrom RColorBrewer brewer.pal
#' @keywords internal

palette.select <- function(palette, info = F) {
  palettes <- list (
    Reds = colorRampPalette(brewer.pal(n = 9, name = "Reds")),
    LgBu = colorRampPalette(c("grey", "lightblue", "blue")),
    GnRd = colorRampPalette(c("green", "black", "red")),
    BuRd = colorRampPalette(c("blue", "black", "red")),
    GrRd = colorRampPalette(c("lightgray", "mistyrose", "red")),
    magma = colorRampPalette(magma(9)),
    GnBu = colorRampPalette(rev(brewer.pal(9,"GnBu"))),
    the.cols = colorRampPalette(c(rgb(255,255,217, maxColorValue=255),
                                  rgb(65,182,196, maxColorValue=255),
                                  rgb(8, 29, 88, maxColorValue=255)),
                                space="Lab"),
    offwhite.to.black = colorRampPalette(c(rgb(220,220,220, maxColorValue=255),
                                           rgb(0, 0, 0, maxColorValue=255)),
                                         space="Lab"),
    viridis = colorRampPalette(viridis(9)),
    cividis = colorRampPalette(cividis(9)),
    plasma = colorRampPalette(plasma(9)),
    heat = colorRampPalette(c("dark blue", "cyan", "yellow", "red")),
    spectral = colorRampPalette(rev(brewer.pal(9,"Spectral"))),
    RdBu = colorRampPalette(rev(brewer.pal(9,"RdBu"))),
    MaYl = colorRampPalette(c("#FF00FF", "black", "#FFFF00")),
    RdYlBu = colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))
  )
  if (info) {
    return(data.frame(palette = names(palettes), category = c("seq", "seq", "seq", "seq", "seq", "seq", "seq", "seq", "seq", "seq", "seq", "seq", "seq", "div", "div", "div", "div"), stringsAsFactors = F))
  }
  if (!palette %in% names(palettes)) {
    stop("Invalid palette name: ", palette, call. = FALSE)
  }
  return(palettes[[palette]])
}

#' Convert gene names
#'
#' @param exprMat expression matrix
#' @param annotation data.frame with one id.column containing the current gene ids and a matching replace
#' column containing the new gene ids
#' @param id.column character vector of unique gene ids of the same type as in the Seurat object
#' @param replace.column character vector of unique gene ids to replace the old ids
#'
#' @return A Seurat object with covnerted gene ids
#' @export
ConvertGeneNames <- function(
  exprMat,
  annotation,
  id.column,
  replace.column
) {

  if (!id.column %in% colnames(annotation) | !replace.column %in% colnames(annotation)) stop(paste0(id.column, " and ", replace.column, " must be present in the annotation table"), call. = F)
  if (sum(duplicated(annotation[, id.column])) > 0 | sum(duplicated(annotation[, id.column])) > 0) stop("Duplicated ids are not allowed", call. = F)

  rownames(annotation) <- annotation[, id.column]
  rownames(exprMat) <- annotation[rownames(exprMat), ]$gene_name
  return(exprMat)
}
