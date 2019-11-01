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
    data <- cbind(data, rep(0, nrow(data)))
    col.order <- c(col.order, setdiff(1:3, col.order))
    data <- data[, col.order]
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
#' @param min.cutoff,max.cutoff Vector of minimum and maximum cutoff values for each feature,
#'  may specify quantile in the form of 'q##' where '##' is the quantile (eg, 'q1', 'q10')

feature.scaler <- function (
  data,
  features,
  min.cutoff,
  max.cutoff,
  spots
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
                             })
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
#' @param object S4 object
#' @param dimPlot A ST.DimPlot output object
#' @param label1 Optional name of the 1st group of capture-spots
#' @param label2 Optional name of the 2nd group of capture-spots
#'
#' @export

labelRGB <- function(object,
                     dimPlot,
                     label1 = "label1",
                     label2 = NULL,
                     red1=c(0,255),
                     green1=c(0,255),
                     blue1=c(0,255),
                     red2=c(0,255),
                     green2=c(0,255),
                     blue2=c(0,255)){
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
#' @param variable character string specifying tha name of the column (variable) containing values
#' @param center.zero Specifies whether or not the colorscale should be centered at zero
#' @param cols Character vector with color ids to be used in colorscale
#' @param val.limits Specifies the limits for values in colorscale
#'
#' @importFrom ggplot2 ggplot geom_point aes_string scale_fill_gradient2 labs scale_fill_gradientn element_rect element_text ggplot_gtable ggplot_build theme

g_legend <- function (
  data,
  variable,
  center.zero,
  cols,
  val.limits
) {
  lg <- ggplot() + geom_point(data = data, aes_string("x", "y", fill = paste0("`", variable, "`")), shape = 22, alpha = 0)
  if (center.zero) {
    lg <- lg +
      scale_fill_gradient2(low = cols[1], mid = cols[median(seq_along(cols))], high = cols[length(cols)], midpoint = 0, na.value = "#000000", limits = val.limits)
  } else if (any(data.type %in% c("character", "factor"))) {
    lg <- lg +
      labs(fill = variable)
  } else {
    lg <- lg +
      scale_fill_gradientn(colours = cols, na.value = "#000000", limits = val.limits)
  }
  lg <- lg + theme(legend.background = element_rect(fill = "#000000"), legend.key = element_rect(fill = "#FFFFFF"), legend.text = element_text(color = "#FFFFFF"), legend.title = element_text(colour = "#FFFFFF"))
  tmp <- ggplot_gtable(ggplot_build(lg))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


#' Obtain array coordinates
#'
#' @param object Seurat object
#'
obtain.array.coords <- function (
  st.object,
  data,
  image.type
) {

  if (all(c("warped_x", "warped_y") %in% colnames(st.object[[]]))) {
    data <- cbind(data, setNames(st.object[[, c("warped_x", "warped_y")]], nm = c("x", "y")))
    image.type <- "processed"
  } else if ("raw" %in% names(st.object@rasterlists)) {
    data <- cbind(data, setNames(st.object[[, c("pixel_x", "pixel_y")]], nm = c("x", "y")))
    image.type <- "processed"
  } else if (all(c("adj_x", "adj_y") %in% colnames(st.object[[]]))) {
    data <- cbind(data, setNames(st.object[[, c("adj_x", "adj_y")]], nm = c("x", "y")))
  } else if (all(c("x", "y") %in% colnames(st.object[[]]))) {
    data <- cbind(data, st.object[[, c("x", "y")]])
  }
  return(list(data, image.type))
}

