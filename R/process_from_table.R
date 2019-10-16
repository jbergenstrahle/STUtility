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

parse.spot.file = function(path, delim="\t", ...) {
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
#' @param transpose set TRUE if count files are in the form of genes as columns and spatial
#' coordiantes as rows.
#' @param min.gene.count filter away genes that has a total read count below this threshold
#' @param min.gene.spots filter away genes that is not expressed below this number of capture spots
#' @param min.spot.count filter away capture spots that contains a total read count below this threshold
#' @param topN OPTIONAL: Filter out the top most expressed genes
#' @param spot.file OPTIONAL: Turn this off (FALSE) if images or/and spotfiles are not used as input
#' @param visium For 10X Visium input. Does not need to be specified if using .h5 files as input
#' @param scaleVisium 10X visium scale factor for pixel coordinates, default set to 0.1039393
#' @param ... parameters passed to \code{\link{CreateSeuratObejct}}
#'
#' @inheritParams ConvertGeneNames
#'
#' @export

InputFromTable <- function(
  infotable,
  transpose=TRUE,
  topN=0,
  min.gene.count = 0,
  min.gene.spots = 0,
  min.spot.count = 0,
  spot.file = TRUE,
  annotation = NULL,
  id.column = NULL,
  replace.column = NULL,
  visium=F,
  scaleVisium=0.1039393,
  ...
){
  counts <- list()
  spotFileData <- list()
  countPaths <- infotable[,"samples"]
  i=1
  rownames(infotable) <- paste0(1:nrow(infotable))


  if(spot.file != FALSE){
    print("Removing all spots outside of tissue - turn this off by parameter spot.file=FALSE")
  }

  for(path in countPaths) {
    cat(paste0("Loading ", path, "\n"))
    if(transpose==FALSE){

      if(length(grep(path, pattern=".h5")) == 1){ #If --> 10x file path
        counts[[path]] <- st.load.matrix(path, visium=T)
        visium=T
        }else{
          counts[[path]] <- st.load.matrix(path)
            }
    }else{

      if(length(grep(path, pattern=".h5")) == 1){
        counts[[path]] <- t(st.load.matrix(path, visium=T))
        visium=T
      }else{
        counts[[path]] <- t(st.load.matrix(path))
          }
      }

    if(spot.file != FALSE){ #Remove spots outside of tissue

      if(visium==T){
        spotsData <- data.frame(parse.spot.file(infotable[which(infotable$samples==path), "spotfiles"], delim=","), stringsAsFactors = F)
        rownames(spotsData) <- as.character(spotsData[,1])
        colnames(spotsData) <- c("barcode", "visium", "new_y", "new_x", "pixel_y", "pixel_x") #OBS, what is column nr2,3,4?
        spotsData$new_x <- as.numeric(spotsData$new_x)
        spotsData$new_y <- as.numeric(spotsData$new_y)
        spotsData <- spotsData[intersect(rownames(spotsData), colnames(counts[[path]])), ]
        counts[[path]] <- counts[[path]][, intersect(rownames(spotsData), colnames(counts[[path]]))]
        spotsData$pixel_x <- as.numeric(spotsData$pixel_x) * scaleVisium
        spotsData$pixel_y <- as.numeric(spotsData$pixel_y) * scaleVisium
        spotFileData[[i]] <- spotsData #Save pixel coords etc

      }else{
        spotsData <- as.data.frame(parse.spot.file(infotable[which(infotable$samples==path), "spotfiles"]))
        if ("selected" %in% colnames(spotsData)) {
          spotsData <- subset(spotsData, selected == 1)
        }
        rownames(spotsData) <- paste(spotsData$x, spotsData$y, sep="x")
        spotsData <- spotsData[intersect(rownames(spotsData), colnames(counts[[path]])), ]
        counts[[path]] <- counts[[path]][, intersect(rownames(spotsData), colnames(counts[[path]]))]
        spotFileData[[i]] <- spotsData #Save pixel coords etc
      }
    }

    i=i+1
    visium=F #If mutiple files from both orginial ST and 10x Visium is used as input
  }

  # ---- Merge counts

  genes <- c()
  for(count in counts)
    genes <- union(genes, rownames(count))

  samples <- c()
  cnt <- matrix(0, nrow=length(genes), ncol=0)
  rownames(cnt) <- genes
  pixelx <- c()
  pixely <- c()
  adj_x <- c()
  adj_y <- c()

  idx <- 1
  for(count in counts) {
    nspots <- ncol(count)
    curgenes <- rownames(count)
    m <- matrix(0, nrow=length(genes), ncol=nspots)
    rownames(m) <- genes
    m[curgenes,] <- count
    colnames(m) <- paste(colnames(count), "_", idx, sep="")
    cnt <- cbind(cnt, m)
    samples <- c(samples, rep(paste0(idx), nspots))
    if(spot.file != FALSE){
      pixelx <- c(pixelx, spotFileData[[idx]]$pixel_x)
      pixely <- c(pixely, spotFileData[[idx]]$pixel_y)
      adj_x <- c(adj_x, spotFileData[[idx]]$new_x)
      adj_y <- c(adj_y, spotFileData[[idx]]$new_y)
    }
    idx <- idx + 1
  }

  if (spot.file) {
    meta_data <- data.frame(sample = samples, ads_x = adj_x, ads_y = adj_y, pixel_x = pixelx, pixel_y = pixely, row.names = colnames(cnt), stringsAsFactors = F)
  } else {
    meta_data <- data.frame(sample = samples, row.names = colnames(cnt), stringsAsFactors = F)
  }

  #----- Add metadata
  metaData <- infotable[, -which(colnames(infotable) %in% c("samples", "spotfiles", "imgs")), drop = F]

  if(ncol(metaData) >= 1){
    for(column in colnames(metaData)){
      meta_data[, column] <- metaData[samples, column]
    }
  }

  if (!is.null(annotation)) {
    if (class(annotation) != "data.frame") stop("annotation table must be a data.frame", call. = FALSE)
    id.column <- id.column %||% "gene_id"
    replace.column <- id.column %||% "gene_name"
    cnt <- ConvertGeneNames(cnt, annotation, id.column, replace.column)
  }

  # ---- pre filtering

  keep.genes <- rowSums(cnt)>=min.gene.count &
    apply(cnt, 1, function(i) sum(i > 0)) > min.gene.spots
  keep.spots <- colSums(cnt)>=min.spot.count

  before <- dim(cnt)
  cnt <- cnt[keep.genes, keep.spots]
  meta_data <- meta_data[keep.spots, ]
  after <- dim(cnt)

  print("------------- Filtering (not including images based filtering) --------------")
  print(paste("Spots removed: ", before[2] - after[2] ))
  print(paste("Genes removed: ", before[1] - after[1]))

  if(is.vector(meta_data)){
    meta_data <- data.frame(sample=meta_data, row.names = colnames(cnt), stringsAsFactors = F)
  }
  m <- CreateSeuratObject(counts = cnt, meta.data = meta_data, ...)

  #Filter top genes
  if (topN > 0){
    m <- m[order(rowSums(as.matrix(m[["RNA"]]@counts)), decreasing=TRUE), ]
    m <- m[1:topN, ]
  }

  # ---- Add image paths
  if(spot.file == FALSE){
    m@misc$spotfile = FALSE
  }else{m@misc$spotfile = TRUE
  m@tools <- list(imgs = infotable$imgs)}

  #Add info for manual annotation tool:
  m@meta.data$id <- seq(1:dim(m)[2])
  m@meta.data$labels <- "Default"

  cat(paste("After filtering the dimensions of the experiment is: "))
  print(dim(m))

  return(m)

}
