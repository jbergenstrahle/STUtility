#' @importFrom Seurat CreateSeuratObject
"_PACKAGE"

#' load count files
#'
#' @param path Path to count files
#' @param suffix add suffix to path
#' @param row.names set column as row.names
#' @keywords internal

#' Parse the output from the ST spot detector tool
#'
#' @param path Path to spot files
#' @param delim delimiter
#' @keywords internal

parse.spot.file = function(path, delim = "\t", ...) {
  x = c()
  tmp = suppressWarnings({try({x = data.frame(fread(input = path, integer64 = "character",
                                                    sep = delim),
                                              check.names = FALSE)})})
  if(inherits(tmp, 'try-error')) {
    return(as.matrix(c()))
  } else {
    return(as.matrix(x))
  }
}


#' Create S4 object from info table.
#' This function is a wrapper to create a complete S4 object with all the samples and metadata
#'
#' @param infotable table with paths to count files and metadata
#' @param min.gene.count filter away genes that has a total read count below this threshold
#' @param min.gene.spots filter away genes that is not expressed below this number of capture spots
#' @param min.spot.count filter away capture spots that contains a total read count below this threshold
#' @param topN OPTIONAL: Filter out the top most expressed genes
#' @param scaleVisium 10X visium scale factor for pixel coordinates, default set to 0.1039393
#' @param pattern.remove Regex pattern used to filter out certain genes
#' @param verbose Print messages
#' @param ... parameters passed to \code{\link{CreateSeuratObejct}}
#'
#' @inheritParams ConvertGeneNames
#'
#' @export

InputFromTable <- function (
  infotable,
  topN = 0,
  min.gene.count = 0,
  min.gene.spots = 0,
  min.spot.count = 0,
  annotation = NULL,
  id.column = NULL,
  replace.column = NULL,
  platform = "Visium",
  scaleVisium = 0.1039393,
  pattern.remove = NULL,
  verbose = TRUE,
  ...
){

  # Generate empty lists
  counts <- list()
  spotFileData <- list()
  countPaths <- infotable[, "samples"]
  rownames(infotable) <- paste0(1:nrow(infotable))

  # Check if spotfiles are present
  if ("spotfiles" %in% colnames(infotable)){
    print("Removing all spots outside of tissue")
  }

  # Specify platform
  if (!"platform" %in% colnames(infotable)) {
    platforms <- rep(platform, nrow(infotable))
  } else {
    platforms <- infotable$platform
  }

  # Parse data files and store in counts and spotFileData
  for (i in seq_along(countPaths)) {
    path <- countPaths[i]
    if (verbose) cat(paste0("Loading ", path, " count matrix from a '", platforms[i], "' experiment\n"))

    if (platforms[i] == "Visium") {
      # Load data
      if (length(grep(path, pattern = ".h5")) == 1) {
        counts[[path]] <- st.load.matrix(path, visium = T)
      } else if (length(grep(path, pattern = ".tsv")) == 1) {
        counts[[path]] <- t(st.load.matrix(path))
      } else {
        stop("Currently only .h5 and .tsv formats are supported for Visium samples")
      }


      # Load spotdata
      # ------------------------------------------------
      # Check that spotfiles are provided
      #if (!"spotfiles" %in% colnames(infotable)) stop("Spotfiles are required for 10X Visium samples", call. = FALSE)
      if ("spotfiles" %in% colnames(infotable)) {
        spotsData <- data.frame(parse.spot.file(infotable[which(infotable$samples == path), "spotfiles"], delim = ","), stringsAsFactors = F)
        if (ncol(spotsData) == 1) {
          spotsData <- setNames(data.frame(parse.spot.file(infotable[which(infotable$samples == path), "spotfiles"], delim = "\t"), stringsAsFactors = F),
                                nm = c("x", "y", "adj_x", "adj_y", "pixel_x", "pixel_y"))
          rownames(spotsData) <- paste(spotsData[, "x"], spotsData[, "y"], sep = "x")
          spotsData <- spotsData[intersect(rownames(spotsData), colnames(counts[[path]])), ]
          counts[[path]] <- counts[[path]][, intersect(rownames(spotsData), colnames(counts[[path]]))]
          spotFileData[[i]] <- spotsData
        } else {
          rownames(spotsData) <- as.character(spotsData[, 1])
          colnames(spotsData) <- c("barcode", "visium", "adj_y", "adj_x", "pixel_y", "pixel_x") #OBS, what is column nr2,3,4?
          spotsData[, c("adj_y", "adj_x", "pixel_y", "pixel_x")] <- apply(spotsData[, c("adj_y", "adj_x", "pixel_y", "pixel_x")], 2, as.numeric)
          spotsData[, c("x", "y")] <- spotsData[, c("adj_x", "adj_y")]
          spotsData <- spotsData[,  c("x", "y", "adj_x", "adj_y", "pixel_x", "pixel_y", "barcode", "visium")]
          spotsData <- spotsData[intersect(rownames(spotsData), colnames(counts[[path]])), ]
          counts[[path]] <- counts[[path]][, intersect(rownames(spotsData), colnames(counts[[path]]))]
          spotsData$pixel_x <- as.numeric(spotsData$pixel_x) * scaleVisium
          spotsData$pixel_y <- as.numeric(spotsData$pixel_y) * scaleVisium
          spotFileData[[i]] <- spotsData
        }
      } else {
        warning(paste0("Extracting spot coordinates from gene count matrix headers. It is highly recommended to use spotfiles."), call. = FALSE)

        # Check if headers can be extracted
        spotsData <- GetCoords(colnames(counts[[path]]))
        if (ncol(spotsData) != 2) stop("Headers are not valid. You have to provide spotfiles or make sure that headers contains (x, y) coordinates", call. = FALSE)
        spotFileData[[i]] <- spotsData
      }
    } else {
      counts[[path]] <- t(st.load.matrix(path))

      # Load spotdata
      # ------------------------------------------------
      if ("spotfiles" %in% colnames(infotable)){
        spotsData <- as.data.frame(parse.spot.file(infotable[which(infotable$samples == path), "spotfiles"]))
        if ("selected" %in% colnames(spotsData)) {
          spotsData <- subset(spotsData, selected == 1)
          spotsData$selected <- NULL
        }
        spotsData <- setNames(spotsData, nm = c("x", "y", "adj_x", "adj_y", "pixel_x", "pixel_y"))
        rownames(spotsData) <- paste(spotsData$x, spotsData$y, sep = "x")
        intersecting.spots <- intersect(rownames(spotsData), colnames(counts[[path]]))
        spotsData <- spotsData[intersecting.spots, ]
        counts[[path]] <- counts[[path]][, intersecting.spots]
        spotFileData[[i]] <- spotsData #Save pixel coords etc
      } else {
        # Obtain x/y coordinates from headers
        spotsData <- GetCoords(colnames(counts[[path]]))
        if (ncol(spotsData) != 2) stop("No spotfiles provided and the headers are invalid. Please make sure that the count matrices are correct.")
        spotFileData[[i]] <- spotsData
      }
    }
  }

  # ---- Merge counts
  genes <- c()
  for (count in counts)
    genes <- union(genes, rownames(count))

  # Check that some genes were found
  if (length(x = Reduce(intersect, lapply(counts, rownames))) == 0) stop("No intersecting genes found across experiment. Please make sure that the count matrices have the same symbols")

  # Generate empty merged count matrix
  samples <- c()
  cnt <- matrix(0, nrow = length(genes), ncol = 0)
  rownames(cnt) <- genes

  # Merge counts and add unique id
  for (i in seq_along(counts)) {
    count <- counts[[i]]
    nspots <- ncol(count)
    curgenes <- rownames(count)
    m <- matrix(0, nrow = length(genes), ncol = nspots)
    rownames(m) <- genes
    m[curgenes,] <- count
    colnames(m) <- paste(colnames(count), "_", i, sep = "")
    cnt <- cbind(cnt, m)
    samples <- c(samples, rep(paste0(i), nspots))
  }

  # Collect meta data if available
  intersecting_columns <- Reduce(intersect, lapply(spotFileData, colnames))
  if (length(x = intersecting_columns) < 2) stop("No spot coordinates found. Aborting ... \n", call. = FALSE)
  meta_data_staffli <- do.call(rbind, lapply(seq_along(spotFileData), function(i) {
    x <- spotFileData[[i]][, intersecting_columns]
    rownames(x) <- paste(rownames(x), "_", i, sep = "")
    return(x)
  }))

  # Add column samples to meta_data used for Staffli object
  meta_data_staffli[, "sample"] <- samples

  # Create an empty meta data table for Seurat
  meta_data_seurat <- data.frame(row.names = rownames(meta_data_staffli))

  #----- Add metadata from infotable to Seurat meta data table
  metaData <- infotable[, -which(colnames(infotable) %in% c("samples", "spotfiles", "imgs")), drop = F]
  if (ncol(metaData) >= 1){
    for (column in colnames(metaData)){
      meta_data_seurat[, column] <- metaData[samples, column]
    }
  }

  # Convert gene symbols if an annotation file is provided
  if (!is.null(annotation)) {

    # Check that annotation is a data.frame
    if (class(annotation) != "data.frame") stop("annotation table must be a data.frame", call. = FALSE)
    id.column <- id.column %||% "gene_id"
    replace.column <- id.column %||% "gene_name"

    if (verbose) cat(paste0("Using provided annotation table with id.column '", id.column, "' and replace column '", replace.column, "' to convert gene names \n"))

    # Check that all genes are present in annotation table
    genes.in.ann <- which(rownames(cnt) %in% annotation[, id.column])

    if (length(x = genes.in.ann) == 0) stop(paste0("None of the genes present in the merged data could be found in the '", id.column, "' column of the annotation table"), call. = FALSE)

    if (sum(rownames(cnt) %in% annotation[, id.column]) != nrow(cnt)) {
      warning(paste0(nrow(cnt) - length(x = genes.in.ann), " are missing from the annotation table and will be removed from the count matrices"), call. = FALSE)
      cnt <- cnt[genes.in.ann, ]
    }

    cnt <- ConvertGeneNames(cnt, annotation, id.column, replace.column)
  }

  # ---- pre filtering
  keep.genes <- rowSums(cnt) >= min.gene.count & rowSums(cnt > 0) > min.gene.spots
  keep.spots <- colSums(cnt) >= min.spot.count

  before <- dim(cnt)
  cnt <- cnt[keep.genes, keep.spots]
  meta_data_staffli <- meta_data_staffli[keep.spots, ]
  meta_data_seurat <- meta_data_seurat[keep.spots, , drop = FALSE]
  after <- dim(cnt)

  if (verbose) {
    cat("\n")
    cat("------------- Filtering (not including images based filtering) -------------- \n")
    cat(paste("  Spots removed: ", before[2] - after[2], " \n"))
    cat(paste("  Genes removed: ", before[1] - after[1], " \n"))
  }

  # Skip meta data if it is empty
  if (ncol(meta_data_seurat) == 0) {
    m <- CreateSeuratObject(counts = cnt)
  } else {
    m <- CreateSeuratObject(counts = cnt, meta.data = meta_data_seurat)
  }

  # Filter top genes
  if (topN > 0){
    m <- m[order(rowSums(as.matrix(m[["RNA"]]@counts)), decreasing = TRUE), ]
    m <- m[1:topN, ]
  }

  # Filter using regex
  if (!is.null(pattern.remove)) {
    remove.genes <- grep(pattern = pattern.remove, x = rownames(cnt))
    if (length(remove.genes) > 0) {
      if (verbose) cat(paste0("Removing ", length(x = remove.genes), " genes matching '", pattern.remove, "' regular expression \n"))
      cnt <- cnt[-remove.genes, ]
    } else {
      if (verbose) cat(paste0("Provided regular expression '", pattern.remove, "' to filter genes but no matches were found \n"))
    }
  }


  # ---- Add image paths
  m@tools$Staffli <- CreateStaffliObject(imgs = ifelse(rep("imgs" %in% colnames(infotable), nrow(infotable)), infotable[, "imgs"], NULL),
                                           meta.data = meta_data_staffli,
                                           platforms = platforms)

  #Add info for manual annotation tool:
  m@meta.data$id <- seq(1:dim(m)[2])
  m@meta.data$labels <- "Default"

  if (verbose) {
    cat(paste0("After filtering the dimensions of the experiment is: [", paste(dim(m), collapse = " genes, "), " spots] \n"))
  }

  return(m)

}
