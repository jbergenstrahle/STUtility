#' @importFrom alphahull ahull
#' @importFrom dplyr filter %>%
#' @importFrom deldir deldir tile.list mid.in
#' @importFrom viridis viridis cividis plasma magma
#' @importFrom RColorBrewer brewer.pal
#' @importFrom Rtsne Rtsne
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom Seurat CreateSeuratObject CreateDimReducObject
"_PACKAGE"

#' Palette selection
#'
#' @param palette Palette choice for plotting spatial expression histology heatmap
#' @keywords internal

palette.select <- function(palette) {
  palettes <- list(
    GnBu = colorRampPalette(RColorBrewer::brewer.pal(9,"GnBu")),
    the.cols = colorRampPalette(c(rgb(255,255,217, maxColorValue=255),
                                  rgb(65,182,196, maxColorValue=255),
                                  rgb(8, 29, 88, maxColorValue=255)),
                                space="Lab"),
    spectral = colorRampPalette(RColorBrewer::brewer.pal(9,"Spectral")),
    offwhite.to.black = colorRampPalette(c(rgb(220,220,220, maxColorValue=255),
                                           rgb(0, 0, 0, maxColorValue=255)),
                                         space="Lab"),
    viridis = colorRampPalette(viridis::viridis(9)),
    cividis = colorRampPalette(viridis::cividis(9)),
    magma = colorRampPalette(viridis::magma(9)),
    plasma = colorRampPalette(viridis::plasma(9)),
    heat = colorRampPalette(c("dark blue", "cyan", "yellow", "red")),
    RdBu = colorRampPalette(RColorBrewer::brewer.pal(9,"RdBu")),
    MaYl = colorRampPalette(c("#FF00FF", "black", "#FFFF00")),
    RdYlBu = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 9, name = "RdYlBu")))
  )
  if (!palette %in% names(palettes)) {
    stop("Invalid palette name: ", palette, call. = FALSE)
  }
  return(palettes[[palette]])
}

#'  ==============  Spatial Expression Histology ====================================
#'
#' Covariance based spatial expression Histology
#'
#' @param se s4 object
#' @param n.clust number of clusters
#' @param sample.wise.z sample wize z-score calculation
#' @param global.z global z-score calculation
#' @param m.var variance?
#' @param log.freq use log of freq
#' @param use.cov use covarance
#' @param n.start number of starts kMeans
#' @param m.iter number of iterations kMeans
#' @export

ST.seh <- function(se, n.clust=10,
                   n.top = 0,
                   log.freq=TRUE, n.start=10, n.iter=10, use.moduleScore=FALSE, use.scale=TRUE
                   ...) {

  # OBS only seurat work atm - add to include singlecellExp also

  # ================================= Clustering of relative counts =======================================

  rnaAssay <- se[["RNA"]]@counts
  rnaAssay <- rnaAssay[order(Matrix::rowSums(rnaAssay), decreasing = T), ]

  if(n.top>0){
    se[["RNA"]]@counts <- se[["RNA"]]@counts[order(Matrix::rowSums(se[["RNA"]]@counts), decreasing = T), ]
    se <- se[1:n.top, ]
    rnaAssay <- se[["RNA"]]@counts
  }

  #  ----- logarithms of relative frequencies in the spots
  if(log.freq==TRUE){
    print("Calculating relative frequences of log(1+counts)")
    rnaAssay<- prop.table(log(1 + rnaAssay, 2))
  } else {
    print("Calculating lo(1+counts)")
    rnaAssay <- log(1 + rnaAssay)
  }


  #Covariance
  print("Calculating Covariance for clustering")
  rnaCov <- cov(x=t(Matrix::as.matrix(rnaAssay)))


  #kMeans clustering
  print("kMeans clustering ...")
  km <- kmeans(rnaCov, centers = n.clust, iter.max=n.iter, nstart=n.start)
  clusters <- km$cluster
  #Save this in s4 obj
  se@assays$cluster <- clusters

  # =============================== Using clusters to sum relative counts ================================

 # for(cluster in unique(sort(clusters))){ #can this be removed?
  #  name = paste("c_", cluster, sep="")
  #  se@meta.data[[name]] <- 0
  #}

 # se@assays$geneContribution <- rep(0, length(rownames(se[["RNA"]]))) #Storage for indivudal gene contribution to cluster
 #instead of subsettting - doing this now:


  plotValues <- matrix(0, ncol=length(unique(clusters)), nrow=dim(se)[[2]])
  colnames(plotValues) <- unique(sort(clusters))
  rownames(plotValues) <- sapply(strsplit(colnames(se), "_"), "[[", 2)

<<<<<<< HEAD
  print("Obtaining SEH values...")

  if(use.moduleScore!=TRUE){
    for(samp in unique(sort(se$sample))) { # Per sample
        x <- se[, se$sample==samp]
        x[["RNA"]]@counts <- prop.table(Matrix::as.matrix(x[["RNA"]]@counts), 2)

      for(cluster in unique(sort(x@assays$cluster))) { #Per gene cluster
        clusterName = paste("c_", cluster, sep="")
        z <- x[x@assays$cluster == cluster, ]
        if(use.scale==TRUE){
          plotValues[which(rownames(plotValues)==samp),
                     cluster] <- scale(Matrix::colSums(z[["RNA"]]@counts), scale=FALSE)
        }else{
          plotValues[which(rownames(plotValues)==samp),
                     cluster] <- Matrix::colSums(z[["RNA"]]@counts)
        }

      }
=======
  for(samp in unique(sort(cm$sample))) { # Per sample
    x <- se[, cm$sample==samp]
    x[["RNA"]]@counts <- prop.table(Matrix::as.matrix(x[["RNA"]]@counts), 2)
    for(cluster in unique(sort(x@assays$cluster))) { #Per gene cluster
      clusterName = paste("c-", cluster, sep="")
      z <- x[x@assays$cluster == cluster, ]
     # se[se@assays$cluster == cluster,
      #   se$sample == samp]@assays$geneContribution <- Matrix::rowSums(z[["RNA"]]@counts)
      #se[se@assays$cluster == cluster,
       #  se$sample == samp]@meta.data[[clusterName]] <- Matrix::colSums(z[["RNA"]]@counts) #values for plots
        # -------------
      plotValues[which(rownames(plotValues)==samp)  , cluster] <- Matrix::colSums(z[["RNA"]]@counts)
>>>>>>> origin/master
    }
  }else{ # ---------- Using Seurats AddModuleScore
    features <- split(names(clusters), clusters)
    clusterScore <- AddModuleScore(se, features=features, assay="RNA", name="CovClusters")

    }

  rownames(plotValues) <- colnames(se)
  se[["SEH"]] <- CreateDimReducObject(embeddings  = plotValues, key = "c_", assay="RNA")
  print("Completed")
  return(se)
}


#' Plot spatial expression histology heatmap
#'
#' @param se Summarized experiment.
#' @param filename Name of output file.
#' @param bg Backsampound color.
#' @param col.scale.cluster Colors will be scaled per gene cluster across all samples.
#' @param disable.voronoi Plot points instead of voronoi
#' @param marBot Bottom plot margindfsdfdsf
#' @param marLeft Left plot margin
#' @param marRight Right plot margin
#' @param marTop Top plot margin
#' @export

seh.plot <- function(se, filename, bg="black", col.scale.cluster=TRUE,
                     disable.voronoi=FALSE, marBot=0, marLeft=0, marRight=0, marTop=0){

  print(paste("Creating factor images based on:"))
  print(paste(dim(se)[[1]], " genes and", dim(se)[[2]], " spots"))

  ns = length(unique(se$sample))
  ncol = length(grep(colnames(se@meta.data), pattern = "c_")) # number of clusters saved in SE object

  plot.asp = (ns + 2) / (ncol + 2)
  filename <- paste(filename, "_expressionHeatmap.pdf", sep="")
  pdf(file = filename, width = 6*ncol, height = 6*ns*plot.asp)

  par(mfrow=c(ns, ncol), mar=c(marBot,marLeft,marTop,marRight), bg=bg)

  #pal <- palette.select("spectral")(9)
  #cols <- rev(rgb(pal(seq(0, 1, length.out = 9)), maxColorValue = 255))

  #cols <- RColorBrewer::brewer.pal(9,"Spectral")
  cols <- palette.select("spectral")(9)

  #Sort by total factor precence
  sums <- NULL
  for(cluster in unique(sort(se@assays$cluster))){
    clusterName = paste("c_", cluster, sep="")
    cSum <- sum(as.data.frame(se[["SEH"]]@cell.embeddings[clusterName]))
    cSum <- as.data.frame(cbind(sum=cSum, cluster = cluster))
    sums <- rbind(sums, cSum)
  }
  clusterOrder <- order(sums[, 1], decreasing=TRUE)
  print("The clusters are printed in order of overall contribution, in order: ")
  print(clusterOrder)
  for(samp in unique(unique(se$sample))){

    for(cluster in clusterOrder) {
      spotNames <- colnames(se[,se$sample==samp])
      spotNames <- sapply(strsplit(spotNames, "_"), "[[", 1)
      xcoord = as.numeric(sapply(strsplit(spotNames, "x"), "[[", 1))
      ycoord = as.numeric(sapply(strsplit(spotNames, "x"), "[[", 2))
      coord_df = as.data.frame(cbind(x=xcoord, y=ycoord))

      clusterName = paste("c_", cluster, sep="")
      if(col.scale.cluster==TRUE){
        factorS <- as.data.frame(se[["SEH"]]@cell.embeddings[,clusterName])
        factorS$sample <- se$sample
        factorS$colorRamp <- colorRampPalette(cols)(nrow(factorS))[rank(factorS[, 1])]
        factorS <- factorS %>% filter(sample==samp)

      }else{
        factorS <- as.data.frame(se[se[["SEH"]]@cell.embeddings==cluster, se$sample==samp]@meta.data[clusterName])
        factorS$colorRamp <- colorRampPalette(cols)(nrow(factorS))[rank(factorS[, 1])]
      }

      if(disable.voronoi==FALSE){
        vtess <- deldir(xcoord, ycoord)
        pts <- cbind(x=xcoord, y=ycoord)
        hull <- ahull(pts, alpha=5)
        indx=hull$arcs[,"end1"]
        the.hull <- list(x = xcoord[indx], y = ycoord[indx])
        tiles <- tile.list(vtess)

        plotJ(tiles, verbose = FALSE, close = FALSE, pch = 1,
              fillcol = factorS$colorRamp, col.pts=NULL,
              col.num=NULL,border=NULL, lwd=5,
              add = FALSE, asp = 1, xlab = "",
              ylab = "", main = paste(cluster), warn=FALSE,
              number=FALSE,adj=NULL,
              showpoints=FALSE,
              clipp=the.hull)
      }else{
        factorS <- as.data.frame(se[se[["SEH"]]@cell.embeddings, se$sample==samp]@meta.data[clusterName])
        factorS$colorRamp <- colorRampPalette(cols)(nrow(factorS))[rank(factorS[, 1])]
        plot(x=xcoord, y=ycoord, col=factorS$colorRamp, lwd=5, asp=1,
             ylab="", xlab="", pch=7,
             xaxt="n", yaxt="n", bty="n")
      }
    }
  }
  dev.off()

}


#' Create tSNE based on spatial expression histology
#'
#' @param se s4 object
#' @param perpelxity perplexity paramter tSNE
#' @param priorPCA use PCA as input to tSNE
#' @param bg plot background color
#' @export

seh.tsne <- function(se, filename, perplexity=30, priorPCA=FALSE, bg="black"){

  ns = length(unique(colData(se)$X))
  ncol = 5

  plot.asp = (ns + 2) / (ncol + 2)
  pdf(file = paste(filename, "_tSNE.pdf",sep=""), width = 6*ncol, height = 6*ns*plot.asp)

  par(mfrow=c(ns, ncol), mar=c(0,0,0,0), bg=bg)

  tsneDf = as.data.frame(colData(se)["c-1"])

  for(cluster in 2:length(unique(sort(rowData(se)$cluster)))){
    clusterName = paste("c-", cluster, sep="")
    tsneDf = cbind(tsneDf, colData(se)[clusterName])
  }

  tsne <- Rtsne(as.matrix(tsneDf), dims=3, pca=priorPCA, perplexity=perplexity)
  range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  newVals <- round(range01(tsne$Y) * 255,0)
  rgbs <- as.factor(rgb(newVals, maxColorValue = 255))
  colData(se)$rgbs <- rgbs

  for(samp in unique(sort(colData(se)$X))){ #plot per sample
    spotNames <- colnames(se[,colData(se)$X==samp])
    spotNames <- sapply(strsplit(spotNames, "_"), "[[", 1)
    xcoord = as.numeric(sapply(strsplit(spotNames, "x"), "[[", 1))
    ycoord = as.numeric(sapply(strsplit(spotNames, "x"), "[[", 2))
    coord_df = as.data.frame(cbind(x=xcoord, y=ycoord))
    vtess <- deldir(xcoord, ycoord)
    pts <- cbind(x=xcoord, y=ycoord)
    hull <- ahull(pts, alpha=5)
    indx=hull$arcs[,"end1"]
    the.hull <- list(x = xcoord[indx], y=ycoord[indx])

    clusterName = paste("c-", cluster, sep="")

    colors <- colData(se[, colData(se)$X==samp])$rgbs
    tiles <- tile.list(vtess)
    plotJ(tiles, verbose = FALSE, close = FALSE, pch = 1,
          fillcol = colors, col.pts=NULL,
          col.num=NULL,border=NULL, lwd=5,
          add = FALSE, asp = 1, xlab = "",
          ylab = "", main = "", warn=FALSE,
          number=FALSE,adj=NULL,
          showpoints=FALSE,
          clipp=the.hull)

  }

  dev.off()

}

#' Gene cluster tables
#'
#' Creates tables of the gene clusters, both with raw count data,
#' and the relative freq. data used in the heatmap plots.
#' @param se s4 object
#' @param filename filename
#' @export

save.table <- function(se, filename){


  df <- data.frame(matrix(ncol = 4, nrow = 0))
  x <- c("cluster", "sample", "spot", "contribution")
  colnames(df) <- x

  df2 <- data.frame(matrix(ncol = 3, nrow = 0))
  x <- c("cluster", "gene", "contribution")   #, "qval")
  colnames(df2) <- x

  #Sort by total factor precence
  sums <- NULL
  for(cluster in unique(sort(rowData(se)$cluster))){
    clusterName = paste("c-", cluster, sep="")
    cSum <- sum(as.data.frame(colData(se)[clusterName]))
    cSum <- as.data.frame(cbind(sum=cSum, cluster = cluster))
    sums <- rbind(sums, cSum)
  }
  clusterOrder <- order(sums[, 1], decreasing=TRUE)

  for(samp in unique(sort(colData(se)$X))){

    for(cluster in clusterOrder) {

      clusterName = paste("c-", cluster, sep="")
      factorS <- as.data.frame(colData(se[rowData(se)$cluster==cluster, colData(se)$X==samp])[clusterName])
      names(factorS) = NULL

      x <- se[rowData(se)$cluster==cluster, colData(se)$X==samp]
      dfNew = as.data.frame(cbind(cluster=cluster, sample=samp, spot=colnames(x), contribution=as.vector(factorS)))

      df = rbind(df, dfNew)

    }

  }

  for(cluster in clusterOrder) {
    genes <- rownames(se[rowData(se)$cluster == cluster, ])
    contri <- rowData(se[rowData(se)$cluster == cluster, ])$geneContribution
    dfNew = as.data.frame(cbind(cbind(cluster=cluster, gene=genes), contribution=contri))
    df2 = rbind(df2, dfNew)
  }

  write.table(df, file=paste(filename, "_t1.csv", sep=""), row.names=FALSE, append=FALSE)
  write.table(df2, file=paste(filename, "_t2.csv", sep=""), row.names=FALSE, append=FALSE)

}

#' ======================================= INTERNAL ====================================
#'  plot tile list for voronoi
#'
#' @keywords internal

plotJ <- function (x, verbose = FALSE, close = FALSE, pch = 1, fillcol = getCol(x,
                                                                                warn = warn), col.pts = NULL, col.num = NULL, border = NULL,
                   showpoints = !number, add = FALSE, asp = 1, clipp = NULL,
                   xlab = "x", ylab = "y", main = "", warn = FALSE, number = FALSE,
                   adj = NULL, axes=FALSE, ...)
{
  object <- x
  if (!inherits(object, "tile.list"))
    stop("Argument \"object\" is not of class tile.list.\n")
  clip <- !is.null(clipp)
  n <- length(object)
  rw <- attr(object, "rw")
  rx <- rw[1:2]
  ry <- rw[3:4]
  x.pts <- unlist(lapply(object, function(w) {
    w$pt[1]
  }))
  y.pts <- unlist(lapply(object, function(w) {
    w$pt[2]
  }))
  if (!add)
    plot(0, 0, type = "n", asp = asp, xlim = rx, ylim = ry,
         xlab = xlab, ylab = ylab, main = main, axes=axes)
  fillcol <- apply(col2rgb(fillcol, TRUE), 2, function(x) {
    do.call(rgb, as.list(x/255))
  })
  fillcol <- rep(fillcol, length = length(object))
  hexbla <- do.call(rgb, as.list(col2rgb("black", TRUE)/255))
  hexwhi <- do.call(rgb, as.list(col2rgb("white", TRUE)/255))
  if (is.null(col.pts)) {
    col.pts <- ifelse(fillcol == hexbla, hexwhi, hexbla)
  }
  else {
    col.pts <- apply(col2rgb(col.pts, TRUE), 2, function(x) {
      do.call(rgb, as.list(x/255))
    })
    col.pts <- rep(col.pts, length = length(object))
  }
  if (is.null(col.num)) {
    col.num <- ifelse(fillcol == hexbla, hexwhi, hexbla)
  }
  else {
    col.num <- apply(col2rgb(col.num, TRUE), 2, function(x) {
      do.call(rgb, as.list(x/255))
    })
    col.num <- rep(col.num, length = length(object))
  }
  if (is.null(border))
    border <- if (all(fillcol == hexbla))
      hexwhi
  else hexbla
  else if (length(border) > 1)
    stop("Argument \"border\" must be a scalar or NULL.\n")
  lnwid <- if (all(fillcol == hexbla))
    2
  else 1
  ptNums <- sapply(x, function(u) {
    u$ptNum
  })
  Adj <- adj
  if (is.null(Adj))
    Adj <- if (showpoints)
      -1
  else 0
  pch <- rep(pch, n)
  okn <- logical(n)
  for (i in 1:n) {
    if (clip) {
      if (requireNamespace("polyclip", quietly = TRUE)) {
        pgon <- polyclip::polyclip(object[[i]], clipp)
        ok <- length(pgon) > 0
      }
      else {
        stop("Cannot clip the tiles; package \"polyclip\" not available.\n")
      }
    }
    else {
      pgon <- list(object[[i]])
      ok <- TRUE
    }
    okn[i] <- ok
    inner <- !any(object[[i]]$bp)
    for (ii in seq(along = pgon)) {
      ptmp <- pgon[[ii]]
      polygon(ptmp, col = fillcol[i], border = NA)
      if (close | inner) {
        polygon(ptmp, col = NA, border = border, lwd = lnwid)
      }
      else {
        x <- ptmp$x
        y <- ptmp$y
        ni <- length(x)
        for (j in 1:ni) {
          jnext <- if (j < ni)
            j + 1
          else 1
          do.it <- mid.in(x[c(j, jnext)], y[c(j, jnext)],
                          rx, ry)
          if (do.it)
            segments(x[j], y[j], x[jnext], y[jnext],
                     col = border, lwd = lnwid)
        }
      }
    }
    if (ok & verbose) {
      if (showpoints)
        points(object[[i]]$pt[1], object[[i]]$pt[2],
               pch = pch[i], col = col.pts[i], ...)
      if (number)
        text(object[[i]]$pt[1], object[[i]]$pt[2], labels = ptNums[i],
             col = col.num[i], adj = Adj, ...)
      if (i < n)
        readline(paste("i = ", i, "; Go? ", sep = ""))
      if (i == n)
        cat("i = ", i, "\n", sep = "")
    }
  }
  if (showpoints & !verbose)
    points(x.pts[okn], y.pts[okn], pch = pch[okn], col = col.pts[okn],
           ...)
  if (number & !verbose)
    text(x.pts[okn], y.pts[okn], labels = ptNums[okn], col = col.num[okn],
         adj = Adj, ...)
  invisible()
}




