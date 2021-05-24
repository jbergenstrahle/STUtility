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

parse.spot.file = function(path, delim = "\t") {
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

# TODO: load data without spotfiles Visium?

#' Create Seurat object from Spatial Transcriptomics data
#'
#' This function is a wrapper to create a complete S4 Seurat object with all the samples and metadata.
#' The input is a data.frame containing paths to all relevant files, s.a. gene count matrices, HE images and
#' spot selection files. The function can also perform some basic filtering and gene conversion.
#'
#' @details
#' This wrapper function has been written to make it easier to various types of Spatial Transcriptomics data, both for the
#' 10X Visium platform and the older array types, here called '1k' and '2k' (with 1000 and 2000 spots respectively). The infotable should
#' at minimum contain paths to the gene count matrices in a column named 'samples', which can be provided as:
#' \itemize{\item{'.tsv' or '.tsv.gz' if the platform is set to either '1k' or '2k'}}
#' \itemize{\item{'.h5' or path to the output folder containing 'barcodes.tsv', 'features.tsv' and 'matrix.mtx' files if the platform is set to 'either 'Visium'}}
#'
#' It is however recommended to include a column named 'imgs' in the `infotable` with paths to the HE images, otherwise none of of the image related
#' processing and visualization methods can be utilized. In addition to the images, a third column named 'spotfiles' should be provided with paths
#' to the tables specifying what spatial coordinates are located under the tissue and which also specifies the pixel coorindates of the spots in the HE images.
#' It is important that the 'spotfiles' pixel coordinates match the HE images provided in the 'imgs' column. If you have run SpaceRanger, the HE images would
#' would correspond to the 'tissue_hires_image.png' files and the spotfiles to the 'tissue_positions_list.csv' files in the output folder. For '1k' and '2k'
#' data, you have to make sure that whatever HE images ('imgs') you are using have a matched selection table in tsv format ('spotfiles'). See
#' [spotdetector tool github](https://github.com/SpatialTranscriptomicsResearch/st_spot_detector/wiki/ST-Spot-Detector-Usage-Guide) for more information.
#'
#' NOTE that for Visium samples, you also need to provide scaling factors which is described in the section below
#'
#' @section samples:
#' If you are analyzing 'Visium' samples, you should provide the '.h5' files in the 'samples' column of infotable.
#' If you are analyzing '1k' or '2k' samples processed with the ST_Pipeline, you can provide the '.tsv' files as 'samples'.
#' We recommend you to use the output '.tsv' file from the ST Pipeline without any modifications and make sure that the
#' matrix is transposed correctly. The matrix will be transposed by defualt, but this can be deactivated by setting `transpose = FALSE`.
#'
#' @section imgs:
#' If you are analyzing 'Visium' samples, you should provide the 'tissue_hires_image.png' files in the 'imgs' column of infotable.
#' If you are analyzing '1k' or '2k' samples you should provide the HE images that were processed with the ST Spot Detector. It is
#' important that these images match the spotfiles!
#'
#' @section spotfiles:
#' If you are analyzing 'Visium' samples, you should provide the 'tissue_positions_list.csv' files in the 'spotfiles' column of infotable.
#' If you are analyzing '1k' or '2k' samples you should provide the HE images that were processed with the ST Spot Detector, typically called
#' 'spot_data-selection...' or 'spot_data-all...'. It is important that these spotfiles are matched with the HE images!
#'
#' @section scaleVisium:
#' If you are analyzing 'Visium' samples, you need to provide a scaling factor to align the spots to the correct positions on the array.
#' These scaling factors are provided in the 'scalefactors_json.json' file in the SpaceRanger output folder. The scaling factor needed for STutility
#' is called "tissue_hires_scalef" and can be provided with the `scaleVisium` argument in `InputFromTable` if the scaling factor is the same
#' for all of your samples. If they are not the same, you can either specify the scaling factors manually into a column named 'scaleVisium'
#' or you can provide the paths to the json files into a column named 'json' in of the infotable data.frame.
#'
#' @section gene id conversion:
#' If you need to convert the gene ids of your expression matrices, you can provide a data.frame with an `ìd.column`
#' with gene symbols matching the symbols of your input matrices and a `replace.column` with the gene symbols that you
#' want to convert to. NOTE that any genes not found in the annotation data.frame will be discarded.
#'
#' @section platform:
#' If you are not analyzing 'Visium' data (default platform), you need to specify what other platform you have used, i.e. '1k' or '2k'.
#' You can also provide a column named 'platform' in the infotable data.frame with a charcter vector specifying the platforms
#' used for each sample.
#'
#' @section meta data:
#' You can also add additional meta data into the infotable data.frame which will be included in the `@meta.data` slot
#' of the returned Seurat object.
#'
#' @section notes:
#' Make sure to check that the paths are correct and preferably absolute paths. If you change the working directory
#' and want to reload the HE images into your Seurat object, you need to make sure that these files can be found on your
#' system.
#'
#' When creating the infotable data.frame, set the parameter `stringsAsFactors = FALSE`.
#'
#' @param infotable table with paths to count files and metadata. See details below for more information.
#' @param transpose transposes the gene expression matrices [default: TRUE]. Only active when reading data
#' from platforms '1k' or '2k'ld
#' @param minUMICountsPerGene filter away genes that has a total UMI count below this threshold (previously called "min.gene.count")
#' @param minSpotsPerGene filter away genes that is not expressed below this number of capture spots (previously called "min.gene.spots)
#' @param minGenesPerSpot ilter away capture spots that contains a total gene count below this threshold (previously called "min.spot.feature.count")
#' @param minUMICountsPerSpot filter away capture spots that contains a total UMI count below this threshold (previously called "min.spot.count")
#' @param topN OPTIONAL: Filter out the top most expressed genes
#' @param annotation data.frame containing columns needed for gene id conversion. See the gene id conversion section
#' for more information.
#' @param id.column column name in annotation data.frame providing the gene ids of the input matrices
#' @param replace.column column name in annotation data.frame providing the gene ids for the conversion
#' @param platform name of the platform used to generate the data [options: 'Visium', '1k', '2k']
#' @param scaleVisium 10X visium scale factor for pixel coordinates matching the "tissue_hires_image.png" files [required].
#' If a numeric value is given, it is assumed that all samples have the same scaling factor. Alternatively, an additional
#' column named "scaleVisium" can be provided with a scaling factor for each sample, or a column named "json" with paths to
#' the "scalefactors_json.json" files.
#' @param verbose Print messages
#' @param ... additional parameters
#'
#' @inheritParams ConvertGeneNames
#'
#' @importFrom jsonlite read_json
#'
#' @export

InputFromTable <- function (
  infotable,
  transpose = TRUE,
  minUMICountsPerGene = 0,
  minSpotsPerGene = 0,
  minGenesPerSpot = 0,
  minUMICountsPerSpot = 0,
  topN = 0,
  annotation = NULL,
  id.column = NULL,
  replace.column = NULL,
  platform = "Visium",
  scaleVisium = NULL,
  disable.subset = FALSE,
  verbose = TRUE,
  ...
){
  
  # Check deprecated arguments
  additional_args <- list(...)
  if ("min.gene.count" %in% names(additional_args)) {
    warning("min.gene.counts has been deprecated, please use minUMICountsPerGene instead")
    minUMICountsPerGene <- additional_args[["min.gene.count"]]
  }
  if ("min.gene.spots" %in% names(additional_args)) {
    warning("min.gene.spots has been deprecated, please use minUMICountsPerGene instead")
    minSpotsPerGene  <- additional_args[["min.gene.spots"]]
  }
  if ("min.spot.feature.count" %in% names(additional_args)) {
    warning("min.spot.feature.count has been deprecated, please use minGenesPerSpot instead")
    minGenesPerSpot <- additional_args[["min.spot.feature.count"]]
  }
  if ("min.spot.count" %in% names(additional_args)) {
    warning("min.spot.count has been deprecated, please use minUMICountsPerSpot instead")
    minUMICountsPerSpot <- additional_args[["min.spot.count"]]
  }

  if (!"samples" %in% colnames(infotable)) stop("No samples have been provided ... \n", call. = FALSE)

  # Generate empty lists
  counts <- list()
  spotFileData <- list()
  ranges.list <- NULL
  spot.diamater.list <- NULL
  countPaths <- infotable[, "samples"]
  rownames(infotable) <- paste0(1:nrow(infotable))

  #
  if (!is.null(scaleVisium) & length(scaleVisium) == 1) scaleVisium <- rep(scaleVisium, nrow(infotable))

  # Check if spotfiles are present
  if ("spotfiles" %in% colnames(infotable) & !disable.subset){
    cat("Using spotfiles to remove spots outside of tissue\n")
  }

  # Specify platform
  if (!"platform" %in% colnames(infotable)) {
    platforms <- rep(platform, nrow(infotable))
  } else {
    platforms <- infotable$platform
  }

  # Check if Visium platform i Visium, in which case the scalefactors are required
  if ("imgs" %in% colnames(infotable)) {
    if (any(platforms == "Visium")) {
      if ("json" %in% colnames(infotable)) {
        suffs <- sapply(infotable[, "json"], getExtension)
        if (!all(suffs == "json")) stop("Incorrect format of json files in infotable ...", call. = FALSE)
        scaleVisium <- scaleVisium %||% sapply(infotable[, "json"], function(f) {read_json(f)$tissue_hires_scalef})
        spot.diamater.list <- sapply(infotable[, "json"], function(f) {read_json(f)$spot_diameter_fullres})*scaleVisium
        infotable[, "json"] <- NULL
      } else if ("scaleVisium" %in% colnames(infotable)) {
        if (!class(infotable[, "scaleVisium"]) == "numeric") stop("Column scaleVisium is not numeric ... \n", call. = FALSE)
        scaleVisium <- infotable[, "scaleVisium"]
        infotable[, "scaleVisium"] <- NULL
      } else if (class(scaleVisium) == "numeric" & length(scaleVisium) == 1) {
        scaleVisium <- rep(scaleVisium, nrow(infotable))
        warning("Only 1 scaleVisium provided. Using this scalefactor for all samples ... \n", call. = FALSE)
      } else {
        warning("No scalefactors provided for Visium samples ... \n", call. = FALSE)
      }
    }
  } else {
    warning("No images provided. You will not be able to use image related functions... \n", call. = FALSE)
  }

  # Parse data files and store in counts and spotFileData
  for (i in seq_along(countPaths)) {
    path <- countPaths[i]
    if (verbose) cat(paste0("Loading ", path, " count matrix from a '", platforms[i], "' experiment\n"))

    if (platforms[i] == "Visium") {
      # Load data
      if (getExtension(path) %in% c("h5", "mtx") | dir.exists(path)) {
        counts[[i]] <- st.load.matrix(path, visium = TRUE)
      } else if (getExtension(path) %in% c("tsv", "tsv.gz")) {
        counts[[i]] <- t(st.load.matrix(path))
      } else {
        stop("Currently only .h5, .mtx and .tsv formats are supported for Visium samples")
      }


      # Load spotdata
      # ------------------------------------------------
      # Check that spotfiles are provided
      #if (!"spotfiles" %in% colnames(infotable)) stop("Spotfiles are required for 10X Visium samples", call. = FALSE)
      if ("spotfiles" %in% colnames(infotable)) {
        spotsData <- data.frame(parse.spot.file(infotable[i, "spotfiles"], delim = ","), stringsAsFactors = F)
        if (ncol(spotsData) == 1) {
          spotsData <- data.frame(parse.spot.file(infotable[i, "spotfiles"], delim = "\t"), stringsAsFactors = F)
          if (ncol(spotsData) == 6) nms <-  c("x", "y", "adj_x", "adj_y", "pixel_x", "pixel_y") else nms <- c("x", "y", "adj_x", "adj_y", "pixel_x", "pixel_y", "selected")
          spotsData <- setNames(spotsData, nm = nms)
          if (ncol(spotsData) == 7 & !disable.subset) {
            spotsData <- subset(spotsData, selected == 1)
          }
          rownames(spotsData) <- paste(spotsData[, "x"], spotsData[, "y"], sep = "x")
          spotsData <- spotsData[intersect(rownames(spotsData), colnames(counts[[i]])), ]
          counts[[i]] <- counts[[i]][, intersect(rownames(spotsData), colnames(counts[[i]]))]
          spotFileData[[i]] <- spotsData
        } else {
          rownames(spotsData) <- as.character(spotsData[, 1])
          if (getExtension(path) %in% c("mtx") | dir.exists(path)){
            if (length(grep(pattern = "\\-1$", x = rownames(spotsData))) > 0) {
                rownames(spotsData) <- gsub(pattern = "\\-1$", replacement = "", x = rownames(spotsData))
            }
          }
          if (ncol(spotsData) == 6) {
            colnames(spotsData) <- c("barcode", "selection", "adj_y", "adj_x", "pixel_y", "pixel_x")
            #if (!disable.subset) {
            #  spotsData <- subset(spotsData, selection == 1)
            #}
          } else if (ncol(spotsData) == 7 & !getExtension(path) %in% c("tsv", "tsv.gz")) {
            colnames(spotsData) <- c("barcode", "visium", "adj_y", "adj_x", "pixel_y", "pixel_x")
          } else if (ncol(spotsData) %in% c(6, 7) & getExtension(path) %in% c("tsv", "tsv.gz")) {
            colnames(spotsData) <- c("x", "y", "adj_x", "adj_y", "pixel_x", "pixel_y", "selection")
            #if (!disable.subset) {
            #  spotsData <- subset(spotsData, selected == 1)
            #}
          } else {
            stop("Spotfiles format not recognized ... \n", call. = FALSE)
          }

          # Subset data
          spotsData[, c("adj_y", "adj_x", "pixel_y", "pixel_x")] <- apply(spotsData[, c("adj_y", "adj_x", "pixel_y", "pixel_x")], 2, as.numeric)
          spotsData$pixel_x <- spotsData$pixel_x * scaleVisium[i]
          spotsData$pixel_y <- spotsData$pixel_y * scaleVisium[i]
          # Save ranges
          ranges.list[[i]] <- sapply(spotsData[, c("pixel_x", "pixel_y")], range)
          if (!disable.subset & "selection" %in% colnames(spotsData)) {
            spotsData <- subset(spotsData, selection == 1)
          }

          spotsData[, c("x", "y")] <- spotsData[, c("adj_x", "adj_y")]
          spotsData <- spotsData[,  c("x", "y", "adj_x", "adj_y", "pixel_x", "pixel_y")]
          spotsData <- spotsData[intersect(rownames(spotsData), colnames(counts[[i]])), ]
          counts[[i]] <- counts[[i]][, intersect(rownames(spotsData), colnames(counts[[i]]))]
          spotFileData[[i]] <- spotsData
        }
      } else {
        if (getExtension(path) %in% c("h5", "mtx") | dir.exists(path)) {
          if (platforms[i] == "Visium") stop("Spotfiles are required for Visium data provided in .h5 or .mtx format.", call. = FALSE)
        }
        warning(paste0("Extracting spot coordinates from gene count matrix headers. It is highly recommended to use spotfiles."), call. = FALSE)

        # Check if headers can be extracted
        spotsData <- GetCoords(colnames(counts[[i]]), delim = "x|_")
        if (ncol(spotsData) != 2) stop("Headers are not valid. You have to provide spotfiles or make sure that headers contains (x, y) coordinates", call. = FALSE)
        spotFileData[[i]] <- spotsData
      }
    } else {
      if (transpose) {
        counts[[i]] <- t(st.load.matrix(path))
      } else{
        counts[[i]] <- st.load.matrix(path)
      }

      # Load spotdata
      # ------------------------------------------------
      if ("spotfiles" %in% colnames(infotable)){
        spotsData <- as.data.frame(parse.spot.file(infotable[i, "spotfiles"]))
        if ("selected" %in% colnames(spotsData)) {
          spotsData <- subset(spotsData, selected == 1)
          spotsData$selected <- NULL
        }
        spotsData <- setNames(spotsData, nm = c("x", "y", "adj_x", "adj_y", "pixel_x", "pixel_y"))
        rownames(spotsData) <- paste(spotsData$x, spotsData$y, sep = "x")
        intersecting.spots <- intersect(rownames(spotsData), colnames(counts[[i]]))
        spotsData <- spotsData[intersecting.spots, ]
        counts[[i]] <- counts[[i]][, intersecting.spots]
        spotFileData[[i]] <- spotsData #Save pixel coords etc
      } else {
        # Obtain x/y coordinates from headers
        spotsData <- GetCoords(colnames(counts[[i]]), delim = "x|_")
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
  cnt <- as(matrix(0, nrow = length(genes), ncol = 0), "dgCMatrix")
  rownames(cnt) <- genes

  # Merge counts and add unique id
  for (i in seq_along(counts)) {
    count <- counts[[i]]
    nspots <- ncol(count)
    curgenes <- rownames(count)
    m <- matrix(0, nrow = length(genes), ncol = nspots)
    rownames(m) <- genes
    m[curgenes, ] <- as.matrix(count)
    colnames(m) <- paste(colnames(count), "_", i, sep = "")
    cnt <- cbind(cnt, as(m, "dgCMatrix"))
    samples <- c(samples, rep(paste0(i), nspots))
  }

  # Collect meta data if available
  intersecting_columns <- Reduce(intersect, lapply(spotFileData, colnames))
  if (length(x = intersecting_columns) < 2) stop("No spot coordinates found. Data cannot be read ... \n", call. = FALSE)
  meta_data_staffli <- do.call(rbind, lapply(seq_along(spotFileData), function(i) {
    x <- spotFileData[[i]][, intersecting_columns]
    rownames(x) <- paste(rownames(x), "_", i, sep = "")
    return(x)
  }))

  # Add column samples to meta_data used for Staffli object
  meta_data_staffli[, "sample"] <- samples

  # Add uncropped coordinates
  meta_data_staffli$original_x <- meta_data_staffli$pixel_x
  meta_data_staffli$original_y <- meta_data_staffli$pixel_y

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
  keep.genes <- rowSums(cnt) >= minUMICountsPerGene & rowSums(cnt > 0) > minSpotsPerGene
  keep.spots <- colSums(cnt) >= minUMICountsPerSpot & colSums(cnt > 0) > minGenesPerSpot

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


  # ---- Add image paths
  if ("imgs" %in% colnames(infotable)) {
    imgs <- infotable[, "imgs"]
  } else {
    imgs <- NULL
  }
  
  m@tools$Staffli <- CreateStaffliObject(imgs = imgs,
                                           meta.data = meta_data_staffli,
                                           platforms = platforms)

  # Add ranges.list if available
  if (!is.null(ranges.list) & !is.null(spot.diamater.list)) {
    if (verbose) cat("Saving capture area ranges to Staffli object \n")
    dims.list <- lapply(seq_along(ranges.list), function(i) {
      rs <- ranges.list[[i]]
      data.frame(min_x = rs[1, 1], max_x = rs[2, 1], min_y = rs[1, 2], max_y = rs[2, 2], spot_diameter = spot.diamater.list[i])
    })
    m@tools$Staffli@dims <- dims.list
  }

  # Add info for manual annotation tool:
  m@meta.data$id <- seq(1:dim(m)[2])
  m@meta.data$labels <- "Default"

  if (verbose) {
    cat(paste0("After filtering the dimensions of the experiment is: [", paste(dim(m), collapse = " genes, "), " spots] \n"))
  }

  return(m)

}
