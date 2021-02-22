#' Squeeze 2 or 3 column feature data into the unit cube and converts into RGB space
#'
#' @param data data.frame containing feature values and coordinates
#' @param channels.use Select channels to use for blending. Default is red, green and blue but the order can be shuffled.
#' For 2 features, the default is red and green. (options: "red", "green" and "blue")

ColorBlender <- function (
  data,
  channels.use = NULL
) {
  rgb.order <- setNames(1:3, c("red", "green", "blue"))
  if (!length(channels.use) == ncol(data)) {
    stop(paste0("channels.use must be same length as number of features or dimensions"))
  } else if (!all(channels.use %in% names(rgb.order))) {
    stop("Invalid color names in channels.use. Valid options are: 'red', 'green' and 'blue'")
  } else if (sum(duplicated(channels.use))){
    stop("Duplicate color names are not allowed in channels.use")
  }
  col.order <- rgb.order[channels.use]

  if (ncol(data) == 2) {
    first_vec <- data[, 1]
    second_vec <- data[, 2]
    data <- matrix(data = 0, nrow = nrow(data), ncol = 3)
    data[, col.order[1]] <- first_vec; data[, col.order[2]] <- second_vec
  } else if (ncol(data) == 3) {
    data <- data[, col.order]
  }
  color.codes <- rgb(data)
}


#' Find the quantile of a data
#'
#' Converts a quantile in character form to a number regarding some data
#' String form for a quantile is represented as a number prefixed with 'q'
#' For example, 10th quantile is 'q10' while 2nd quantile is 'q2'
#'
#' Will only take a quantile of non-zero data values
#'
#' @param cutoff The cutoff to turn into a quantile
#' @param data The data to find the quantile of
#'
#' @references \url{https://github.com/satijalab/seurat/blob/master/R/visualization.R}
#'
#' @return The numerical representation of the quantile
#'
#' @importFrom stats quantile

SetQuantile <- function (
  cutoff,
  data
) {
  if (grepl(pattern = '^q[0-9]{1,2}$', x = as.character(x = cutoff), perl = TRUE)) {
    this.quantile <- as.numeric(x = sub(
      pattern = 'q',
      replacement = '',
      x = as.character(x = cutoff)
    )) / 100
    data <- unlist(x = data)
    data <- data[data > 0]
    cutoff <- quantile(x = data, probs = this.quantile)
  }
  return(as.numeric(x = cutoff))
}

#' Function used to scale numerical features
#'
#' @param data data.frame containing x, y coordinates and columns with numerical features
#' @param features feature names
#' @param min.cutoff,max.cutoff Vector of minimum and maximum cutoff values for each feature,
#'  may specify quantile in the form of 'q##' where '##' is the quantile (eg, 'q1', 'q10')

feature.scaler <- function (
  data,
  features,
  min.cutoff,
  max.cutoff
) {
  min.cutoff <- setNames(mapply(
    FUN = function(cutoff, feature) {
      return(ifelse(
        test = is.na(x = cutoff),
        yes = min(data[, feature]),
        no = cutoff
      ))
    },
    cutoff = min.cutoff,
    feature = features
  ), nm = features)
  max.cutoff <- setNames(mapply(
    FUN = function(cutoff, feature) {
      return(ifelse(
        test = is.na(x = cutoff),
        yes = max(data[, feature]),
        no = cutoff
      ))
    },
    cutoff = max.cutoff,
    feature = features
  ), nm = features)
  check.lengths <- unique(x = vapply(
    X = list(features, min.cutoff, max.cutoff),
    FUN = length,
    FUN.VALUE = numeric(length = 1)
  ))
  if (length(x = check.lengths) != 1) {
    stop("There must be the same number of minimum and maximum cuttoffs as there are features")
  }

  data[, features] <- sapply(X = features,
                             FUN = function(index) {
                               data.feature <- as.vector(x = data[, index])
                               min.use <- SetQuantile(cutoff = min.cutoff[index], data.feature)
                               max.use <- SetQuantile(cutoff = max.cutoff[index], data.feature)
                               data.feature[data.feature < min.use] <- min.use
                               data.feature[data.feature > max.use] <- max.use
                               return(data.feature)
                             }) %>% as.data.frame()
  return(data)
}

#' Generates a dark theme for STPlot
#' @importFrom ggplot2 element_rect element_text theme

dark_theme <- function() {
  theme(plot.background = element_rect(fill = "black", color = "black"),
        panel.background = element_rect(fill = "black", color = "black"),
        plot.title = element_text(colour = "white"),
        strip.text = element_text(colour = 'white'),
        legend.title = element_text(colour = "white"),
        legend.text = element_text(colour = "white"))
}

#' Store rgb values after ST.Dimplot
#' OBS! could be merged with ST.Dimplot?
#'
#' @param plotObj A ST.DimPlot output object
#'
#' @importFrom plotly ggplotly
#' @importFrom ggplot2 ggplot_build
#'
#' @export

ST.rgbDimPlot <- function(
  plotObj){
  ggBuild <- ggplot_build(gg)
  rgbTransform <- as.data.frame(t(col2rgb(ggBuild$data[[1]]$colour)))
  plotObj$layers[[1]]$mapping$rgbs <- paste("R:", round(rgbTransform$red/2.55),
                                            "% G:", round(rgbTransform$green/2.55),
                                            "% B:", round(rgbTransform$blue/2.55), "%",
                                            sep="")

  plot <- ggplotly(gg , tooltip="rgbs")
  return(plot)

}


#' Select spots by rgb values after ST.DimPlot
#'
#' @param object A Seurat object
#' @param dimPlot A 'ggplot' object created with \code{ST.DimPlot}
#' @param label1,label2 Optional names of the 1st and 2nd groups of capture-spots
#'
#' @export

labelRGB <- function (
  object,
  dimPlot,
  label1 = "label1",
  label2 = NULL,
  red1 = c(0, 255),
  green1 = c(0, 255),
  blue1 = c(0, 255),
  red2 = c(0, 255),
  green2 = c(0, 255),
  blue2 = c(0, 255)
){
  # --------------------- Add this part to parameter in ST.Dimplot
  ggBuild <- ggplot_build(dimPlot)
  rgbTransform <- as.data.frame(t(col2rgb(ggBuild$data[[1]]$colour)))
  ggBuild$layers[[1]]$mapping$rgbs <- paste("R:", round(rgbTransform$red/2.55),
                                            "% G:", round(rgbTransform$green/2.55),
                                            "% B:", round(rgbTransform$blue/2.55), "%",
                                            sep="")

  object$rgb.red <- round(rgbTransform$red/2.55)
  object$rgb.green <- round(rgbTransform$green/2.55)
  object$rgb.blue <- round(rgbTransform$blue/2.55)
  # -------------------------------------------------------------
  object$rgbInfo <- 0
  object$rgbInfo[which(object$rgb.red > red1[1] & object$rgb.red < red1[2]
                       & object$rgb.green > green1[1] & object$rgb.green < green1[2]
                       & object$rgb.blue > blue1[1] & object$rgb.blue < blue1[2])] <- label1

  if(!is.null(label2)){
    object$rgbInfo[which(object$rgb.red > red2[1] & object$rgb.red < red2[2]
                         & object$rgb.green > green2[1] & object$rgb.green < green2[2]
                         & object$rgb.blue > blue2[1] & object$rgb.blue < blue2[2])] <- label2
  }

  object$rgbInfo  <- as.factor(object$rgbInfo)

  #seuratobject <- AddMetaData(object, metadata=object$rgbInfo)
  Idents(object = object) <- object$rgbInfo #skriver over?

  return(object)
}


#' Create legend
#'
#' Creates a colorscale legend for a dataset in grob format
#'
#' @param data data.frame with values that should be mapped onto colorscale
#' @param data.type input data type, e.g. "factor", "numeric", "character", ...
#' @param variable character string specifying tha name of the column (variable) containing values
#' @param center.zero Specifies whether or not the colorscale should be centered at zero
#' @param cols Character vector with color ids to be used in colorscale
#' @param val.limits Specifies the limits for values in colorscale
#' @param dark.theme Should a dark theme be used?
#'
#' @importFrom ggplot2 ggplot geom_point aes_string scale_fill_gradient2 labs scale_fill_gradientn element_rect element_text ggplot_gtable ggplot_build theme

g_legend <- function (
  data,
  data.type = "numeric",
  variable,
  center.zero,
  cols,
  val.limits,
  dark.theme = TRUE
) {
  lg <- ggplot() + geom_point(data = data, aes_string("x", "y", fill = paste0("`", variable, "`")), shape = 21, alpha = 1, size = 8)
  if (center.zero & !any(data.type %in% c("character", "factor"))) {
    lg <- lg +
      scale_fill_gradient2(low = cols[1], mid = cols[median(seq_along(cols))], high = cols[length(cols)], midpoint = 0, na.value = ifelse(dark.theme, "#000000", "#FFFFFF"), limits = val.limits)
  } else if (any(data.type %in% c("character", "factor"))) {
    lg <- lg +
      labs(fill = variable) +
      scale_fill_manual(values = cols)
  } else {
    lg <- lg +
      scale_fill_gradientn(colours = cols, na.value = ifelse(dark.theme, "#000000", "#FFFFFF"), limits = val.limits)
  }
  lg <- lg + theme(legend.background = element_rect(fill = ifelse(dark.theme, "#000000", "#FFFFFF")),
                   legend.key = element_rect(fill = ifelse(dark.theme, "#000000", "#FFFFFF"), colour = ifelse(dark.theme, "#000000", "#FFFFFF")),
                   legend.text = element_text(color = ifelse(dark.theme, "#FFFFFF", "#000000")),
                   legend.title = element_text(colour = ifelse(dark.theme, "#FFFFFF", "#000000")))
  tmp <- ggplot_gtable(ggplot_build(lg))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


#' Obtain array coordinates
#'
#' @param object A Seurat object
#' @param data Data.frame with array coordinates
#' @param image.type One of "raw", "masked" or "processed"
#' @param spots Spots to subset data on

obtain.array.coords <- function (
  st.object,
  data,
  image.type,
  spots
) {

  if (all(c("warped_x", "warped_y") %in% colnames(st.object[[]]))) {
    data <- cbind(data, setNames(st.object[[spots, c("warped_x", "warped_y")]], nm = c("x", "y")))
    image.type <- "processed"
  } else if ("raw" %in% names(st.object@rasterlists)) {
    data <- cbind(data, setNames(st.object[[spots, c("pixel_x", "pixel_y")]], nm = c("x", "y")))
    image.type <- "raw"
  } else if (all(c("adj_x", "adj_y") %in% colnames(st.object[[]]))) {
    data <- cbind(data, setNames(st.object[[spots, c("adj_x", "adj_y")]], nm = c("x", "y")))
  } else if (all(c("x", "y") %in% colnames(st.object[[]]))) {
    data <- cbind(data, st.object[[spots, c("x", "y")]])
  }
  return(list(data, image.type))
}

#' Draws a scalebar for ST plots
#'
#' Takes a ggplot object as input together with some additional parameters
#' specifying the position and widht of a scalebar representing the actual
#' width of corresponding 500um in the plot. The actual width of the scale
#' bar is determined by the x, xend parameters and has to be calulated
#' separately. To do this, you can use the pixels.per.um slot of the
#' `Staffli` object.
#'
#' @param p An object of class 'ggplot'
#' @param x,xend The start and positions of the scalebar along the x axis.
#' @param y The position of the scalebar along the y axis.
#' @param pxum A data.frame object with columns for 'x', 'xend', 'y' and 'sample'
#' used for facetted plots.
#' @param sb.size Size of scalebar
#' @param dark.theme Switches color of scalebar to 'white'
#'
#' @return An object of class 'ggplot' with a scalebar drawn on top of it
#'
#' @importFrom ggplot2 geom_segment annotate aes
#'

draw_scalebar <- function (
  p,
  x = NULL,
  xend = NULL,
  y = NULL,
  pxum = NULL,
  sb.size = 2.5,
  dark.theme = FALSE
) {

  if (sb.size == 0) {
    return(p)
  }
  end.width <- sb.size/25
  mid.width <- sb.size/50
  sb.color = ifelse(dark.theme, "white", "black")

  if (!is.null(pxum)) {
    p <- p + geom_segment(data = pxum, aes(x = x, xend = xend, y = y, yend = y, group = sample), color = sb.color) +
      geom_segment(data = pxum, aes(x = x, xend = x, y = y - (xend - x)*end.width, yend = y + (xend - x)*end.width, group = sample), color = sb.color) +
      geom_segment(data = pxum, aes(x = x + (xend - x)/5, xend = x + (xend - x)/5, y = y - (xend - x)*mid.width, yend = y + (xend - x)*mid.width, group = sample), color = sb.color) +
      geom_segment(data = pxum, aes(x = x + (2*(xend - x))/5, xend = x + (2*(xend - x))/5, y = y - (xend - x)*mid.width, yend = y + (xend - x)*mid.width, group = sample), color = sb.color) +
      geom_segment(data = pxum, aes(x = x + (3*(xend - x))/5, xend = x + (3*(xend - x))/5, y = y - (xend - x)*mid.width, yend = y + (xend - x)*mid.width, group = sample), color = sb.color) +
      geom_segment(data = pxum, aes(x = x + (4*(xend - x))/5, xend = x + (4*(xend - x))/5, y = y - (xend - x)*mid.width, yend = y + (xend - x)*mid.width, group = sample), color = sb.color) +
      geom_segment(data = pxum, aes(x = xend, xend = xend, y = y - (xend - x)*end.width, yend = y + (xend - x)*end.width, group = sample), color = sb.color) +
      geom_text(data = pxum, aes(x = x + (xend - x)/2, y = y + (xend - x)/3, label = paste0("500 ", "\u03bcm"), group = sample), size = sb.size, color = sb.color)
  } else {
    p <- p + geom_segment(aes(x = x, xend = xend, y = y, yend = y), color = sb.color) +
      geom_segment(aes(x = x, xend = x, y = y - (xend - x)*end.width, yend = y + (xend - x)*end.width), color = sb.color) +
      geom_segment(aes(x = x + (xend - x)/5, xend = x + (xend - x)/5, y = y - (xend - x)*mid.width, yend = y + (xend - x)*mid.width), color = sb.color) +
      geom_segment(aes(x = x + 2*(xend - x)/5, xend = x + 2*(xend - x)/5, y = y - (xend - x)*mid.width, yend = y + (xend - x)*mid.width), color = sb.color) +
      geom_segment(aes(x = x + 3*(xend - x)/5, xend = x + 3*(xend - x)/5, y = y - (xend - x)*mid.width, yend = y + (xend - x)*mid.width), color = sb.color) +
      geom_segment(aes(x = x + 4*(xend - x)/5, xend = x + 4*(xend - x)/5, y = y - (xend - x)*mid.width, yend = y + (xend - x)*mid.width), color = sb.color) +
      geom_segment(aes(x = xend, xend = xend, y = y - (xend - x)*end.width, yend = y + (xend - x)*end.width), color = sb.color) +
      annotate(geom = "text", x = x + (xend - x)/2, y = y + (xend - x)/3, label = paste0("500 ", "\u03bcm"), size = sb.size, color = sb.color)
  }
  return(p)
}


#' Prep scalebar
#'
#' @param st.object An object of class 'Staffli'
#' @param data data.frame for plotting
#' @param data.type String specifying the class of the variable used for splitting
#' @param indices Integer vector specifying sample indices to include in the plot
#' @param split.labels Split plots
#' @param features Features included in data
#' @param dims List of image dimensions
#' @param show.sb Return NULL if FALSE
#'
#' @return A data.frame containing positions of a scale bar used for plotting
#'

prep.sb <- function (
  st.object,
  data,
  data.type,
  indices,
  split.labels,
  features,
  dims,
  show.sb
) {
  if (show.sb) {
    if (length(x = st.object@pixels.per.um) > 0) {
      pixels.per.um <- st.object@pixels.per.um[indices]
      if (split.labels) {
        if (length(x = features) > 1) stop(paste0("Only one feature allowed when splitting labels ..."), call. = FALSE)
        if (!data.type %in% c("character", "factor")) stop(paste0("Only categorical variables can be used for splitting ..."), call = FALSE)
        pxum <- data.frame(pixels.per.um = rep(pixels.per.um, length(unique(data[, features]))), sample = unique(data[, features]))
        hewidths <- rep(dims[[1]][1], length(unique(data[, features])))
        heheights <- rep(dims[[1]][2], length(unique(data[, features])))
      } else {
        pxum <- data.frame(pixels.per.um, sample = factor(paste0(indices), levels = unique(paste0(indices))))
        hewidths <- lapply(dims, function(x) x[1]) %>% unlist()
        heheights <- lapply(dims, function(x) x[2]) %>% unlist()
      }

      pxum$sb500 <- pxum$pixels.per.um*500
      pxum$hewidth <-  hewidths
      pxum$x <- 7*pxum$hewidth/9
      pxum$xend <- 7*pxum$hewidth/9 + pxum$sb500
      pxum$y <- heheights - heheights/8

      # Check that the scale bar is inside plot
      dx <- pxum$hewidth - pxum$xend
      pxum$x <- ifelse(dx < 0, pxum$x + dx*2, pxum$x)
      pxum$xend <- ifelse(dx < 0, pxum$xend + dx*2, pxum$xend)

    } else {
      pxum <- NULL
    }
  } else {
    pxum <- NULL
  }
  return(pxum)
}
