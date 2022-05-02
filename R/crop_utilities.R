#' Get crop geometries
#' 
#' This function can be used to define "crop windows" based on sets 
#' of target spots. 
#' 
#' @param object A Seurat object created with STutility
#' @param group.by A meta.data column to group the spots by. Only character 
#' vector or factors are allowed. Selects "labels" by default.
#' @param groups.to.keep Provide a caharcter vector of groups to keep. 
#' @param keep.all.spots Keep all spots within selected crop windows and not just the spots 
#' belonging to the predefined groups. All spots that are present in multiple crop windows will be removed.
#' @param xy_padding Increase the crop area in all directions. Given in pixels.
#' @return A list of "crop geometries" to be used with \code{\link{CropImages}})
#' 
#' @export

GetCropWindows <- function (
  object,
  group.by = NULL,
  groups.to.keep = NULL,
  symmetric = TRUE,
  keep.all.spots = FALSE,
  xy_padding = 50
) {
  
  # Check to see if Staffli object is present
  if (!"Staffli" %in% names(object@tools)) stop("Staffli object is missing from Seurat object. ", call. = FALSE)
  st.object <- object@tools$Staffli
  
  # Check group.by
  group.by <- group.by %||% "labels"
  group.by.keep <- group.by
  if (!group.by %in% colnames(object@meta.data)) stop("group.by is not available in the Seurat object meta.data slot. ", call. = FALSE)
  group.by <- object@meta.data[, group.by]
  if (!class(group.by) %in% c("character", "factor")) stop(paste0("Invalid group.by column of class ", class(group.by)), call. = FALSE)
  group.by <- as.character(group.by)
  #if (length(unique(group.by)) == 1) stop("Only one group provided.")
  
  coords <- data.frame(st.object@meta.data[, c("pixel_x", "pixel_y")], Var = group.by, section = st.object@meta.data$sample)
  
  # Reduce groups
  if (!is.null(groups.to.keep)) {
    coords <- subset(coords, Var %in% groups.to.keep) 
  }
  
  # Get groups
  if (!is.null(groups.to.keep)) {
    grps <- groups.to.keep
  } else {
    if (is.null(group.by)) {
      grps <- setdiff(unique(coords$Var), "Default")
    } else {
      grps <- unique(coords$Var)
    }
  }
  
  crop.geoms.all.sections <- Reduce(c, lapply(unique(coords$section), function(sid) {
    coords.subset <- subset(coords, section %in% sid)
    crop.geoms <- lapply(grps, function(grp) {
      coords.subset.grp <- subset(coords.subset, Var %in% grp)
      if (nrow(coords.subset.grp) == 0) return(NULL)
      minxy <- apply(coords.subset.grp[, c("pixel_x", "pixel_y")], 2, range)
      minxy[1, ] <- minxy[1, ] - xy_padding
      minxy[2, ] <- minxy[2, ] + xy_padding
      minxy <- round(minxy)
      wh <- apply(minxy, 2, diff)
      if (symmetric) {
        wh <- rep(max(wh), 2)
        mid.points <- apply(minxy, 2, mean)
        minxy <- matrix(c(mid.points[1] - max(wh)/2, 
                          mid.points[1] + max(wh)/2,
                          mid.points[2] - max(wh)/2,
                          mid.points[2] + max(wh)/2), ncol = 2)
        geom <- setNames(data.frame(paste0(wh[1], "x", wh[2], "+", minxy[1, 1], "+", minxy[1, 2]), grp, group.by.keep, keep.all.spots), nm = c("geom", "group", "group.by", "all.spots"))
      } else {
        geom <- setNames(data.frame(paste0(wh[1], "x", wh[2], "+", minxy[1, 1], "+", minxy[1, 2]), grp, group.by.keep, keep.all.spots), nm = c("geom", "group", "group.by", "all.spots"))
      }
    })
    crop.geoms <- crop.geoms[!sapply(crop.geoms, is.null)]
    names(crop.geoms) <- rep(sid, length(crop.geoms))
    return(crop.geoms)
  }))

  return(crop.geoms.all.sections)
  
}
