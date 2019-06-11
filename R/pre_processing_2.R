#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom Seurat CreateSeuratObject
"_PACKAGE"

#' load count files
#'
#' @param path Path to count files
#' @param suffix add suffix to path
#' @param row.names set column as row.names
#' @keywords internal

st.load.matrix = function(path, suffix="", row.names=1, ...) {
  x = c()
  tmp = try({ x = read.delim(paste(path, suffix, sep=""),
                             header=T,
                             row.names=row.names,
                             sep=",",
                             check.names=F,
                             ...)})
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

prepare.counts <- function(sampleTable, transpose=TRUE,
                           topN=0, th.gene=0, th.spot=0,
                           spot.file,
                           object.type = "Seurat", ...){
  outObj = NULL
  counts <- list()
  paths <- sampleTable[,1]
  for(path in paths) {
    cat(paste0("Loading ", path, "\n"))
    if(transpose==FALSE){
      counts[[path]] <- st.load.matrix(path)
    }else{counts[[path]] <- t(st.load.matrix(path))}

    #pre-filter
    counts[[path]] <- counts[[path]][rowSums(counts[[path]])>th.gene, colSums(counts[[path]])>th.spot]

  }

  # ---- Merge counts

  genes <- c()
  for(count in counts)
    genes <- union(genes, rownames(count))

  samples <- c()
  cnt <- matrix(0, nrow=length(genes), ncol=0)
  rownames(cnt) <- genes


  idx <- 1
  for(count in counts) {
    nspots <- ncol(count)
    curgenes <- rownames(count)
    m <- matrix(0, nrow=length(genes), ncol=nspots)
    rownames(m) <- genes
    m[curgenes,] <- count
    colnames(m) <- paste(colnames(count), "_", idx, sep="")
    cnt <- cbind(cnt, m)
    samples <- c(samples, rep(idx, nspots))
    idx <- idx + 1
  }

  if(is.na(spot.file)){
    print("No spot file provided")
  }else{
    #Parse spot file
  }

  if(object.type=="Seurat"){
    m <- CreateSeuratObject(counts = cnt)
    m$sample <- samples
    m <- m[order(rowSums(as.matrix(m[["RNA"]]@counts)), decreasing=TRUE),]
  }else{
    m <- SingleCellExperiment(assays=list(counts=cnt))
    m$sample <- samples
    m <- m[order(rowSums(assays(m)$counts), decreasing=TRUE),]
  }
  if (topN > 0){
    m <- m[1:topN, ]
  }


  #----- Add metadata
  if(!is.null(ncol(sampleTable)) & ncol(sampleTable)>1){
    for(column in colnames(sampleTable[2:ncol(sampleTable)])){
      #print(column)
      m[[column]] <- sampleTable[samples, column]
    }
  }


  cat(paste("After filtering the dimensions of the experiment is: "))
  print(dim(m))

  return(m)

}
