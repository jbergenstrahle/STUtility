#' @importFrom SingleCellExperiment SingleCellExperiment
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
#' @param sampleTable table with paths to count files and metadata
#' @param transpose set TRUE if count files are in the form of genes as columns and spatial
#' coordiantes as rows.
#' @param topN OPTIONAL: Filter out the top most expressed genes
#' @param th.gene OPTIONAL: Filter genes with total counts across all samples below threshold
#' @param th.spot OPTIONAL: Filter `"capture-spots"` with total count values below threshold
#' @export

prep.from.table <- function(sampleTable, transpose=TRUE,
                           topN=0, th.gene=0, th.spot=0,
                           object.type = "Seurat",
                           spot.file = TRUE,  ...){
  outObj = NULL
  counts <- list()
  spotFileData <- list()
  countPaths <- sampleTable[,"samples"]
  i=1

  if(spot.file != FALSE){
    print("Removing all spots outside of tissue - turn this off by parameter spot.file=FALSE")
  }

  for(path in countPaths) {
    cat(paste0("Loading ", path, "\n"))
    if(transpose==FALSE){
      counts[[path]] <- st.load.matrix(path)
    }else{counts[[path]] <- t(st.load.matrix(path))}

    if(spot.file != FALSE){ #Remove spots outside of tissue
      spotsData <- as.data.frame(parse.spot.file(sampleTable[which(sampleTable$samples==path), "spotfiles"]))
      spotsDataKeep <- paste(spotsData[which(spotsData$selected==1),]$x, spotsData[which(spotsData$selected==1),]$y, sep="x")
      counts[[path]] <- counts[[path]][, c(spotsDataKeep)]
      spotFileData[[i]] <- spotsData[which(spotsData$selected==1), ] #Save pixel coords etc
    }

    #pre-filter
    counts[[path]] <- counts[[path]][rowSums(counts[[path]])>th.gene, colSums(counts[[path]])>th.spot]
    i=i+1
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

  if(object.type=="Seurat"){
    m <- CreateSeuratObject(counts = cnt)
  }else{
    m <- SingleCellExperiment(assays=list(counts=cnt))
  }

  #add meta
  m$sample <- samples
  if(spot.file != FALSE){
    m$pixel_x <- pixelx
    m$pixel_y <- pixely
    m$ads_x <- adj_x
    m$ads_y <- adj_y
  }

  #Filter top genes
  if (topN > 0){
    m <- m[order(rowSums(as.matrix(m[["RNA"]]@counts)), decreasing=TRUE),]
    m <- m[1:topN, ]
  }

  #----- Add metadata
  metaData <- sampleTable[, -which(colnames(sampleTable) %in% c("samples", "spotfiles", "imgs"))]
  if(!is.null(ncol(metaData)) & ncol(metaData)>1){
    for(column in colnames(metaData)){
      #print(column)
      m[[column]] <- metaData[samples, column]
    }
  }
  # ---- Add image paths
  m@tools <- list(imgs = sampleTable$imgs)

  cat(paste("After filtering the dimensions of the experiment is: "))
  print(dim(m))

  return(m)

}
