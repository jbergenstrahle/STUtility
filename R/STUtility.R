#' @importFrom alphahull ahull
#' @importFrom dplyr filter %>%
#' @importFrom deldir deldir tile.list mid.in
#' @importFrom viridis cividis viridis
#' @importFrom RColorBrewer brewer.pal
#' @importFrom Rtsne Rtsne
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom Seurat CreateSeuratObject
"_PACKAGE"

#' Palette selection
#'
#' @param palette Palette choice for plotting spatial expression histology heatmap
#' @keywords internal

palette.select <- function(palette) {
  palettes <- list(
    GnBu = colorRampPalette(brewer.pal(9,"GnBu")),
    the.cols = colorRampPalette(c(rgb(255,255,217, maxColorValue=255),
                           rgb(65,182,196, maxColorValue=255),
                           rgb(8, 29, 88, maxColorValue=255)),
                         space="Lab"),
    spectral = colorRampPalette(brewer.pal(9,"Spectral")),
    offwhite.to.black = colorRampPalette(c(rgb(220,220,220, maxColorValue=255),
                                    rgb(0, 0, 0, maxColorValue=255)),
                                  space="Lab"),
    viridis = colorRampPalette(viridis(9)),
    cividis = colorRampPalette(cividis(9)),
    magma = colorRampPalette(magma(9)),
    plasma = colorRampPalette(plasma(9)),
    heat = colorRampPalette(c("dark blue", "cyan", "yellow", "red")),
    RdBu = colorRampPalette(brewer.pal(9,"RdBu")),
    MaYl = colorRampPalette(c("#FF00FF", "black", "#FFFF00")),
    RdYlBu = colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))
  )
  if (!palette %in% names(palettes)) {
    stop("Invalid palette name: ", palette, call. = FALSE)
  }
  return(palettes[[palette]])
}

#' Plot spatial expression histology heatmap
#'
#' @param se Summarized experiment.
#' @param filename Name of output file.
#' @param bg Backsampound color.
#' @param col.scale.cluster Colors will be scaled per gene cluster across all samples.
#' @param disable.voronoi Plot points instead of voronoi
#' @param marBot Bottom plot margin
#' @param marLeft Left plot margin
#' @param marRight Right plot margin
#' @param marTop Top plot margin
#' @export

seh.plot <- function(se, filename, bg="black", col.scale.cluster=TRUE, disable.voronoi=FALSE, marBot=0, marLeft=0, marRight=0, marTop=0){

  print(paste("Creating factor images based on:"))
  print(paste(dim(se)[[1]], " genes and", dim(se)[[2]], " spots"))

  ns = length(unique(colData(se)$X))
  ncol = length(unique(colData(se)))-1 # number of clusters saved in SE object

  plot.asp = (ns + 2) / (ncol + 2)
  filename <- paste(filename, "_expressionHeatmap.pdf", sep="")
  pdf(file = filename, width = 6*ncol, height = 6*ns*plot.asp)

  par(mfrow=c(ns, ncol), mar=c(marBot,marLeft,marTop,marRight), bg=bg)

  pal <- palette.select("spectral")
  cols <- rev(rgb(pal(seq(0, 1, length.out = 9)), maxColorValue = 255))

  #Sort by total factor precence
  sums <- NULL
  for(cluster in unique(sort(rowData(se)$cluster))){
    clusterName = paste("c-", cluster, sep="")
    cSum <- sum(as.data.frame(colData(se)[clusterName]))
    cSum <- as.data.frame(cbind(sum=cSum, cluster = cluster))
    sums <- rbind(sums, cSum)
  }
  clusterOrder <- order(sums[, 1], decreasing=TRUE)
  print("The clusters are printed in order of overall contribution, in order: ")
  print(clusterOrder)
  for(samp in unique(sort(colData(se)$X))){

    for(cluster in clusterOrder) {
      spotNames <- colnames(se[,colData(se)$X==samp])
      spotNames <- sapply(strsplit(spotNames, "_"), "[[", 1)
      xcoord = as.numeric(sapply(strsplit(spotNames, "x"), "[[", 1))
      ycoord = as.numeric(sapply(strsplit(spotNames, "x"), "[[", 2))
      coord_df = as.data.frame(cbind(x=xcoord, y=ycoord))

      clusterName = paste("c-", cluster, sep="")
      if(col.scale.cluster==TRUE){
        factorS <- as.data.frame(colData(se)[clusterName])
        factorS$sample <- colData(se)$X
        factorS$colorRamp <- colorRampPalette(cols)(nrow(factorS))[rank(factorS[, 1])]
        factorS <- factorS %>% filter(sample==samp)

      }else{
        factorS <- as.data.frame(colData(se[rowData(se)$cluster==cluster, colData(se)$X==samp])[clusterName])
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
        factorS <- as.data.frame(colData(se[rowData(se)$cluster==cluster, colData(se)$X==samp])[clusterName])
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

#' Spatial Expression Histology
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

runSEH <- function(se, n.clust=10, sample.wise.z = FALSE, global.z = FALSE, m.var=0,
                   log.freq=TRUE, use.cov = TRUE, n.start=10, m.iter=10,
                   db.scan = FALSE, db.sca.eps = 0.15, db.scan.MinPts = 5) {

  # --------------- OBS OBS ----- TAnk pa TOP genes EJ SORTED BY DEFAULT HAR

  if(sample.wise.z){
    useCov=FALSE
  }
  if(global.z){
    useCov=FALSE
  }
  # ---- filter out most variable
  #mostVar: e.g. 0.25 would give top 25% variable
  if(m.var>0){
    #Convert to log counts
    logX <- log(assays(se)$counts+1)
    #Calculate coefficient of variation
    compute_cv <- function(x) sd(x) / mean(x)
    cv <- apply(logX, 1, compute_cv)
    #Selection of genes with highest CV.
    cutoff <- m.var
    variable <- logX[rank(cv) / length(cv) > 1 - cutoff, ]
    se <- se[rownames(variable), ]
  }

  #  ----- logarithms of relative frequencies in the spots
  if(log.freq==TRUE){
    print("Relative frequences of log")
    assays(se)$y <- prop.table(log(1 + assays(se)$counts, 2))
  }else{
    print("Log of counts")
    assays(se)$y <- log(1 + assays(se)$counts)
  }

  myscale = function(x) {
    cm <- colMeans(x)
    y <- t(t(x) - cm)
    sds <- apply(y, 2, sd)
    sds[sds==0] = 1
    t(t(y)/sds)
  }

  if(use.cov==TRUE){
    #Covariance
    print("Using Covariance matrix for clustering")
    y <- cov(x=t(assays(se)$y))
    metadata(se)$cov <- y

  }else{

    if(sample.wise.z){
      print("Performing sample wise zscore calculations")
      # z-transform the genes across the spots within the samples
      for(samp in unique(sort(colData(se)$X))) { # Per sample
        assays(se[, colData(se)$X==samp])$y <- myscale(assays(se[, colData(se)$X==samp])$y)
      }
      y <- assays(se)$y
      metadata(se)$y <- y
    }else {
      print("Performing global zscore calculation")
      assays(se)$y <- myscale(assays(se)$y)
      y <- assays(se)$y
      metadata(se)$y <- y
    }
  }


  #if(db.scan == FALSE){
  print("kMeans clustering ...")
  km <- kmeans(y, centers = n.clust, iter.max=m.iter, nstart=n.start)
  rowData(se)$cluster <- km$cluster
  #}else{
  #  print("DBSCAN clustering ...")
  #  db <- fpc::dbscan(y, eps=db.scan.eps, MinPts = db.scan.MinPts)
  #}

  for(cluster in unique(sort(rowData(se)$cluster))){
    name = paste("c-", cluster, sep="")
    colData(se)[name] = 0
  }

  rowData(se)$geneContribution <- 0 #Storage for indivudal gene contribution to cluster

  print("Preparing data for vizualisation")
  for(samp in unique(sort(colData(se)$X))) { # Per sample
    x <- se[, colData(se)$X==samp]
    assays(x)$counts <- prop.table(assays(x)$counts, 2)
    for(cluster in unique(sort(rowData(se)$cluster))) { #Per gene cluster
      clusterName = paste("c-", cluster, sep="")
      z <- x[rowData(x)$cluster==cluster, ]
      rowData(se[rowData(se)$cluster==cluster, colData(se)$X==samp])$geneContribution <- rowSums(assays(z)$counts)
      colData(se[rowData(se)$cluster==cluster, colData(se)$X==samp])[clusterName] <- colSums(assays(z)$counts) #values for plots

    }
  }
  print("Completed")
  return(se)
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

#' plot tile list for voronoi
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



