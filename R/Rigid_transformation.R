#' Obtain edges of binaryh mask stores in masked.masks list
#'
#' @param object Seurat object
#' @param index Sample index
#' @param verbose Print messages
#'
#' @importFrom imager imgradient add map_il

get.edges <- function (
  object,
  index,
  verbose = FALSE,
  type = "masked.masks"
) {
  if (verbose) cat(paste0(" Detecting edges of sample ", index, "\n"))
  im <- object@tools[[type]][[index]]

  grad <- imgradient(as.cimg(im))
  grad.sq <- grad %>% map_il(~ .^2)

  grad.sq <- imager::add(grad.sq)
  grad.sq <- apply(grad.sq, c(1, 2), max)
  return(grad.sq/max(grad.sq))
}



#' Match two point sets using iterative closest point search
#'
#' @param set1 Point set from image to be aligned with reference
#' @param set2 Point set from reference image
#' @param iterations Number of iterations
#' @param mindist Minimum distance that definces valid points
#' @param type Select type of transform to be applied
#' @param threads number of trheads to use
#' @return list with transformed x, y coordinates and list of transformation matrices
#'
#' @importFrom Rvcg vcgCreateKDtree vcgSearchKDtree
#' @importFrom Morpho computeTransform applyTransform

icpmat <- function (
  set1,
  set2,
  iterations,
  mindist = 1e+15,
  type = c("rigid", "similarity", "affine"),
  threads = 1
) {
  set1 <- cbind(set1, 0)
  set2 <- cbind(set2, 0)
  xtmp <- set1
  yKD <- vcgCreateKDtree(set2)
  transformations <- list()

  for (i in 1:iterations) {
    clost <- vcgSearchKDtree(yKD, xtmp, 1, threads = threads)
    good <- which(clost$distance < mindist)
    trafo <- computeTransform(set2[clost$index[good], ], xtmp[good, ], type = type, weights = NULL, centerweight = FALSE)
    xtmp <- applyTransform(xtmp[, ], trafo)
    transformations[[i]] <- trafo
  }
  xtmp <- xtmp[, 1:2]
  return(list(xy = xtmp, map =  Reduce(`%*%`, transformations)))
}



#' Finds optimal transform based on RMSE
#'
#' Tests different types of reflection settings and return the optimal solution
#' based on RMSE between the transformed points and thre reference set
#'
#' @param set1 Point set from image to be aligned with reference
#' @param set2 Point set from reference image
#' @param xdim Width of image
#' @param ydim Height of image
#' @return list with the list of tranformation matrices, reflection coordinates and rmse score
#' for the optimal transformation
#'
#' @importFrom Rvcg vcgKDtree

find.optimal.transform <- function (
  set1,
  set2,
  xdim,
  ydim
) {
  os <- matrix(c(c(0, 1, 0, 1)*xdim, c(0, 0, 1, 1)*ydim), ncol = 2)
  trf <- lapply(1:nrow(os), function(i) {
    p <- set1
    px_dims <- os[i, ]
    p <- t(abs(t(p) - px_dims))
    icpr <- icpmat(p, set2, iterations = 10)
    RMSE <- sqrt(sum(vcgKDtree(set2, icpr$xy, k = 1)$distance^2))
    return(list(icp = icpr, os = px_dims, rmse = RMSE))
  })

  rmses <- unlist(lapply(trf, function(x) x$rmse))
  return(trf[[which.min(rmses)]])
}

#' Apply rigid transformation to a set of points
#'
#' Takes a list of obtained with \code{\link{FindOptimalTransform}} and
#' a matrix of x, y coordinates and returns the transformed x, y coordinates
#'
#' @param icp List containing transformation matrices
#' @param set Matrix of x, y coordinates to be transformed
#' @return Matrix of transformed x, y coordinates

apply.transform <- function (
  map,
  set
) {
  set.new <- cbind(as.matrix(set), 0)
  set.new <- applyTransform(x = set.new, map)
  return(set.new[, 1:2])
}


#' Generate a map function
#'
#' Given two sets of points (2D) this function will create a rigid transformation
#' function that minimized the rmse of the two points sets. Because we use the backward
#' transformation function for image warping, we need to compute the forward transformation
#' function for pixel coordinate warping. If set1 is (x, y) and set2 is the target (x', y')
#' point set, we want to find the backward transformation function M(x', y') -> (x, y).
#' The forward transformation will be the inverse of M, i.e. M^-1(x, y) -> (x', y').
#' The iterative closest point algorithm used to find the best rigid transformation function
#' may not find the best solution if the image is reflected. For this reason we calculate
#' the transofmration function for all different combinations of reflections and select the
#' function with the lowest rmse between the aligned set and reference set.
#'
#' @param icps results obtained with \code{\link{find.optimal.transform}}
#' @return A transformation function that takes x and y coordinates as input and outputs a
#' list of warped x, y coordinates

generate.map.affine <- function (
  icps,
  forward = FALSE
) {
  #icps <- find.optimal.transform(set2, set1, xdim, ydim)
  if (forward) {
    map.affine <- function (x, y) {
      p <- cbind(x, y)
      os <- icps$os
      xy <- apply.transform(map = solve(icps$icp$map), p)
      xy <- t(abs(t(xy) - os))
      list(x = xy[, 1], y = xy[, 2])
    }
  } else {
    map.affine <- function (x, y) {
      p <- cbind(x, y)
      p <- t(abs(t(p) - icps$os))
      xy <- apply.transform(map = icps$icp$map, p)
      list(x = xy[, 1], y = xy[, 2])
    }
  }
  return(map.affine)
}


#' Automatic alignment of HE stained tissue images
#'
#' Image alignment or image registration consists in finding a rigid tranformation
#' function that remaps pixels between two images so that two images are aligned.
#' AlignImages allows you to align all images to a reference (any image present in the Seurat object)
#' which simplifies the interpretation of spatial heatmaps and can be useful for creating 3D
#' models. The transformation function is learned using the ICP (Iterative Closest Point) on two point
#' sets which defines the edges of the tissues in the HE images. Note that this alignment works best
#' for tissue sections that are intact, i.e. not cropped or folded. Also, because the method only used
#' the edges of the tissue for alignment, you might end up with strange results if the tissue
#' shape is symmetrical.
#'
#' @param object A Seurat object
#' @param indices Integer: sample indices of images to align with reference
#' @param reference.index Integer: sample index of referece image
#' @param verbose Print messages
#'
#' @importFrom imager as.cimg imwarp
#' @importFrom grDevices as.raster
#'
#' @export

AlignImages <- function (
  object,
  indices = NULL,
  reference.index = NULL,
  verbose = FALSE
) {

  if (!"masked" %in% names(object@tools)) stop(paste0("Masked images are not present in Seurat object"), call. = FALSE)
  if (!any(c("pixel_x", "pixel_y") %in% colnames(object[[]]))) stop(paste0("Pixel coordinates are missing in Seurat object"), call. = FALSE)

  reference.index <- reference.index %||% 1
  if (verbose) cat(paste0("Selecting image ", reference.index, " as reference for alignment. \n"))

  reference.edge <- get.edges(object, index = reference.index)
  indices <- indices %||% (1:length(object@tools$imgs))[-reference.index]
  edge.list <- lapply(indices, function(i) {
    get.edges(object, index = i, verbose = verbose)
  })

  xyset.ref <- which(reference.edge > 0, arr.ind = T)
  colnames(xyset.ref) <- c("x", "y")
  xyset <- setNames(lapply(edge.list, function(edge) {
    xy <- which(edge > 0, arr.ind = T)
    colnames(xy) <- c("x", "y")
    return(xy)
  }), nm = indices)

  # Obtain reference image
  im.ref <- as.cimg(object@tools$raw[[reference.index]])

  # Create empty lists
  transformations <- setNames(ifelse(rep("transformations" %in% names(object@tools), length(object@tools$imgs)), object@tools$transformations, lapply(1:length(object@tools$imgs), function(i) {diag(c(1, 1, 1))})), nm = names(object@tools$masked))
  processed.images <- setNames(ifelse(rep("processed" %in% names(object@tools), length(object@tools$imgs)), object@tools$processed, object@tools$masked), nm = names(object@tools$masked))
  processed.masks <- object@tools$masked.masks
  warped_coords <- object[[c("pixel_x", "pixel_y")]]

  for (i in indices) {
    ima <- as.cimg(object@tools$raw[[i]])
    ima.msk <- as.cimg(object@tools$masked.masks[[i]])


    if (verbose) cat(paste0("Processing image ", i, " \n Estimating transformation function ... \n"))
    xdim <- object@tools$xdim
    width <- as.numeric(object@tools$dims[[i]][2]); height <- as.numeric(object@tools$dims[[i]][3])
    ydim <- round(height/(width/xdim))

    # Obtain optimal transform and create map functions
    icps <- find.optimal.transform(xyset.ref, xyset[[paste0(i)]], xdim, ydim)
    tr <- icps$icp$map
    tr <- tr[-3, -3]
    transformations[[i]] <- tr%*%transformations[[i]]
    map.affine.backward <- generate.map.affine(icps)
    map.affine.forward <- generate.map.affine(icps, forward = T)

    # Warp images
    if (verbose) cat(paste0(" Applying rigid transformation ... \n"))
    imat = imwarp(ima, map = map.affine.backward, dir = "backward", interpolation = "cubic")
    imat.msk = imwarp(ima.msk, map = map.affine.backward, dir = "backward", interpolation = "linear")
    inds <- which(imat.msk != 255)

    # Obtain scale factors
    dims.raw <- as.numeric(object@tools$dims[[i]][2:3])
    dims.scaled <- dim(object@tools$raw[[i]])
    sf.xy <- dims.raw/dims.scaled
    pixel_xy <- subset(object[[]], sample == paste0(i))[, c("pixel_x", "pixel_y")]/sf.xy

    # Warp coordinates
    warped_xy <- map.affine.forward(pixel_xy[, 1], pixel_xy[, 2])
    warped_coords[rownames(pixel_xy), ] <- sapply(setNames(as.data.frame(t(t(do.call(cbind, warped_xy))*sf.xy)), nm = c("x", "y")), round, 2)

    if (verbose) cat(paste0(" Cleaning up background ... \n"))
    imrst <- as.raster(imat)
    imat[inds] <- 255
    imrst <- as.raster(imat)
    tab.im <- table(imrst)
    if (length(tab.im) > 2) {
      imrst[imrst == names(which.max(tab.im))] <- "#FFFFFF"
    }

    if (verbose) cat(paste0(" Image ", i, " alignment complete. \n\n"))
    processed.images[[i]] <- imrst
    processed.masks[[i]] <- as.raster(imat.msk)
  }

  object@tools$processed <- processed.images
  object@tools$processed.masks <- processed.masks
  object@tools$transformations <- transformations
  object[[c("warped_x", "warped_y")]] <- warped_coords

  return(object)
}


#' Manual alignment of images
#'
#' Creates an interactive shiny application to align images manually
#'
#' @param object Seurat object
#' @param Image type used for alignment
#' @param reference.index Specifies reference sample image for alignment(default: 1)
#' @param edges Uses the tissue edges as points set for alignment
#' @param verbose Print messages
#'
#' @importFrom shiny runApp fluidPage fluidRow column sliderInput checkboxInput selectInput actionButton plotOutput reactive renderPlot eventReactive observe stopApp
#' @importFrom shinyjs useShinyjs reset
#'
#' @export

ManualAlignImages <- function (
  object,
  type = NULL,
  reference.index = 1,
  verbose = FALSE,
  edges = TRUE
) {

  # use processed images as input if available
  type <- type %||% {
    if ("processed" %in% names(object@tools)) {
      "processed"
    } else if ("masked" %in% names(object@tools)) {
      "masked"
    } else {
      stop(paste0("No masked images are available in the Seurat object"), call. = FALSE)
    }
  }
  if (verbose) cat(paste0("Using ", type, " images as input for alignment ... \n"))

  # Obtain point sets from each image
  scatters <- grid.from.seu(object, type = type, edges = TRUE)
  fixed.scatter <- scatters[[reference.index]]$scatter
  counter <- NULL
  coords.ls <- NULL
  tr.matrices <- lapply(seq_along(object@tools[[type]]), function(i) diag(c(1, 1, 1)))


  ui <- fluidPage(
    useShinyjs(),
    fluidRow(
      column(3,
                    sliderInput(
                      inputId = "angle",
                      label = "Rotation angle",
                      value = 0, min = -120, max = 120, step = 0.1
                    ),
                    sliderInput(
                      inputId = "shift_x",
                      label = "Move along x axis",
                      value = 0, min = -200, max = 200, step = 1
                    ),
                    sliderInput(
                      inputId = "shift_y",
                      label = "Move along y axis",
                      value = 0, min = -200, max = 200, step = 1
                    ),
                    sliderInput(
                      inputId = "size",
                      label = "Change point size",
                      value = 0.5, min = 0.1, max = 6, step = 0.1
                    ),
                    checkboxInput(inputId = "flip_x",
                                         label = "Mirror along x axis",
                                         value = FALSE),
                    checkboxInput(inputId = "flip_y",
                                         label = "Mirror along y axis",
                                         value = FALSE),
                    selectInput(inputId = "sample", choices = 2:6, label = "Select sample", selected = 2),
                    actionButton("myBtn", "Return aligned data")
      ),

      column(8, plotOutput("scatter")
      )
    )
  )

  server <- function(input, output) {

    rotation_angle <- reactive({
      input$angle
    })

    translation_xy <- reactive({
      trxy <- c(input$shift_x, input$shift_y)
      return(trxy)
    })

    mirror_xy <- reactive({
      mirrxy <- c(input$flip_x, input$flip_y)
      return(mirrxy)
    })

    pt_size <- reactive({
      input$size
    })

    coords_list <- reactive({

      # Obtain point set and spot pixel coordinates
      ls <- scatter.coords()
      scatter.t <- ls[[1]]; coords.t <- ls[[2]]

      # Set transformation parameters
      xt.yt <- translation_xy()
      xy.alpha <- rotation_angle()
      mirrxy <-  mirror_xy()

      # Apply reflections
      center <- apply(scatter.t, 2, mean)
      tr.mirror <- mirror(mirror.x = mirrxy[1], mirror.y = mirrxy[2], center.cur = center)

      # Apply rotation
      tr.rotate <- rotate(angle = -xy.alpha, center.cur = center)

      # Apply translation
      tr.translate <- translate(translate.x = xt.yt[1], translate.y = -xt.yt[2])

      # Combine transformations
      tr <- tr.translate%*%tr.rotate%*%tr.mirror


      # Apply transformations
      scatter.t <- t(tr%*%rbind(t(scatter.t), 1))[, 1:2]
      coords.t <- t(tr%*%rbind(t(coords.t), 1))[, 1:2]

      return(list(scatter = scatter.t, coords = coords.t, tr = tr))
    })

    output$scatter <- renderPlot({

      coords.ls <<- coords_list()
      scatter.t <- coords.ls[[1]]; coords.t <- coords.ls[[2]]

      d <- round((sqrt(400^2 + 400^2) - 400)/2)

      plot(fixed.scatter[, 1], 400 - fixed.scatter[, 2], xlim = c(-d, 400 + d), ylim = c(-d, 400 + d))
      points(scatter.t[, 1], 400 - scatter.t[, 2], col = "gray")
      points(coords.t[, 1], 400 - coords.t[, 2], col = "red", cex = pt_size())

    }, height = 800, width = 800)

    scatter.coords <- eventReactive(input$sample, {
      reset("angle"); reset("shift_x"); reset("shift_y"); reset("flip_x"); reset("flip_y")
      if (!is.null(counter)) {
        scatters[[counter]] <<- coords.ls[c(1, 2)]
        if (!is.null(tr.matrices[[counter]])) {
          tr.matrices[[counter]] <<- coords.ls[[3]]%*%tr.matrices[[counter]]
          #cat("Sample:", counter, "\n",  tr.matrices[[counter]][1, ], "\n", tr.matrices[[counter]][2, ], "\n", tr.matrices[[counter]][3, ], "\n\n")
        } else {
          tr.matrices[[counter]] <<- coords.ls[[3]]
        }
      }
      scatter <- scatters[[as.numeric(input$sample)]]$scatter
      coords <- scatters[[as.numeric(input$sample)]]$coords
      counter <<- as.numeric(input$sample)
      return(list(scatter, coords))
    })

    observe({
      if(input$myBtn > 0){
        if (!is.null(counter)) {
          scatters[[counter]] <<- coords.ls[c(1, 2)]
          if (!is.null(tr.matrices[[counter]])) {
            tr.matrices[[counter]] <<- coords.ls[[3]]%*%tr.matrices[[counter]]
            cat("Sample:", counter, "\n",  tr.matrices[[counter]][1, ], "\n", tr.matrices[[counter]][2, ], "\n", tr.matrices[[counter]][3, ], "\n\n")
          } else {
            tr.matrices[[counter]] <<- coords.ls[[3]]
          }
        }
        stopApp(tr.matrices)
      }
    })

  }

  # Returned transformation matrices
  alignment.matrices <- runApp(list(ui = ui, server = server))
  if (verbose) cat(paste("Finished image alignment. \n\n"))
  processed.ids <- which(unlist(lapply(alignment.matrices, function(tr) {!all(tr == diag(c(1, 1, 1)))})))
  print(alignment.matrices)

  # Create lists for transformation
  transformations <- setNames(ifelse(rep("transformations" %in% names(object@tools), length(object@tools$imgs)), object@tools$transformations, lapply(1:length(object@tools$imgs), function(i) {diag(c(1, 1, 1))})), nm = names(object@tools$masked))
  processed.images <- setNames(ifelse(rep("processed" %in% names(object@tools), length(object@tools$imgs)), object@tools$processed, object@tools$masked), nm = names(object@tools$masked))
  masks <- object@tools[[paste0(type, ".masks")]]
  processed.masks <- setNames(ifelse(rep("processed.masks" %in% names(object@tools), length(object@tools$imgs)), object@tools$processed.masks, object@tools$masked.masks), nm = names(object@tools$masked))
  warped_coords <- object[[c("pixel_x", "pixel_y")]]

  for (i in processed.ids) {

    if (verbose) cat(paste0("Loading masked image for sample ", i, " ... \n"))
    m <- object@tools$masked[[i]]

    # Obtain alignment matrix
    tr <- alignment.matrices[[i]]
    transformations[[i]] <- tr%*%transformations[[i]]

    map.rot.backward <- generate.map.rot(tr)
    map.rot.forward <- generate.map.rot(tr, forward = TRUE)

    # Obtain scale factors
    dims.raw <- as.numeric(object@tools$dims[[i]][2:3])
    dims.scaled <- dim(object@tools$raw[[i]])
    sf.xy <- dims.raw/rev(dims.scaled)
    pixel_xy <- subset(object[[]], sample == paste0(i))[, c("pixel_x", "pixel_y")]/sf.xy

    # Warp pixels
    if (verbose) cat(paste0("Warping pixel coordinates for ", i, " ... \n"))
    warped_xy <- sapply(setNames(as.data.frame(do.call(cbind, map.rot.forward(pixel_xy$pixel_x, pixel_xy$pixel_y))), nm = c("warped_x", "warped_y"))*sf.xy, round, digits = 1)
    warped_coords[rownames(pixel_xy), 1:2] <- warped_xy

    if (verbose) cat(paste0("Warping image for ", i, " ... \n"))
    processed.images[[i]] <- Warp(m, map.rot.backward)
    msk <- masks[[i]]
    if (verbose) cat(paste0("Warping image mask for ", i, " ... \n"))
    processed.masks[[i]] <- Warp(msk, map.rot.backward, mask = T)
    if (verbose) cat(paste0("Finished alignment for sample ", i, " \n\n"))
  }

  object@tools$transformations <- transformations
  object@tools$processed <- processed.images
  object@tools$processed.masks <- processed.masks
  object[[c("warped_x", "warped_y")]] <- warped_coords
  return(object)

}


#' Function used to generate a 2D set of points from an HE image
#'
#' The defualt segmentation uses the HE image as input and defines any pixel with an intensity value
#' below a threshold to be a point. The number of points can be downsampeld to limit the maximum number.
#'
#' @param object Seurat object
#' @param type Sets the image type to run segmentation on
#' @param sammple.index Integer value specifying the index of the sample to be analyzed
#' @param limit Sets the intensity threshold in the interval [0, 1]
#' @param maxnum Integer value specifying the maximum number of points
#' @param edges Extracts the coordinates of the edges instead

scatter_HE <- function (
  object,
  type = "masked",
  sample.index = NULL,
  limit = 0.5,
  maxnum = 5e4,
  edges = FALSE
){
  sample.index <- sample.index %||% 1

  if (!type %in% names(object@tools)) stop(paste0(type, " images not fount in Seurat obejct"), call. = FALSE)

  if (edges) {
    bw.image <- get.edges(object, sample.index, type = type)
    xyset = which(bw.image > limit, arr.ind = TRUE)
  } else {
    bw.image = grayscale(as.cimg(object@tools[[type]][[sample.index]]))
    xyset = which(bw.image < limit*255, arr.ind = TRUE)
  }
  set.seed(1)
  if (maxnum < nrow(xyset)) {
    xyset <- xyset[sample(1:nrow(xyset), size = maxnum, replace = FALSE), ]
  }
  xyset <- xyset[, 1:2] %>% as.data.frame() %>% setNames(nm = c("x", "y"))
  return(xyset)
}


#' Creates 2D point patterns for a set of images
#'
#' @param obejct Seurat object
#'
#' @inheritParams scatter_HE

grid.from.seu <- function (
  object,
  type,
  limit = 0.3,
  maxnum = 5e4,
  edges = FALSE
) {

  if (!type %in% names(object@tools)) stop(paste0(type, " images not fount in Seurat obejct"), call. = FALSE)

  setNames(lapply(1:length(object@tools[[type]]), function(sample.index) {
    scatter <- scatter_HE(object = se, sample.index = sample.index, maxnum = 1e3, limit = limit, type = type, edges = edges)
    if (type == "processed") {
      xy.names <- c("warped_x", "warped_y")
    } else {
      xy.names <- c("pixel_x", "pixel_y")
    }
    coords <- subset(se[[]], sample == sample.index)[, xy.names]
    dims.raw <- as.numeric(se@tools$dims[[sample.index]][2:3])
    dims.scaled <- dim(se@tools$raw[[sample.index]])
    sf.xy <- dims.raw/rev(dims.scaled)
    coords <- coords/sf.xy
    return(list(scatter = scatter, coords = coords))
  }), nm = names(object@tools$masked))
}
