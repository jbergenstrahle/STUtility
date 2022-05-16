#' @include generics.R Staffli.R
NULL

#' @importFrom shiny debounce observeEvent reactive
#' @importFrom data.table as.data.table
#' @importFrom ggplot2
#' aes_ aes_string coord_fixed element_blank geom_segment ggplot ggtitle guides
#' guide_legend
#' labs
#' annotation_custom
#' theme theme_bw theme_minimal
#' scale_color_manual scale_fill_manual scale_size
#' scale_x_continuous scale_y_continuous
#' @importFrom shiny pageWithSidebar headerPanel sidebarPanel mainPanel textInput strong actionButton radioButtons sliderInput
#' reactiveValues observeEvent observe hr submitButton runApp HTML modalDialog withProgress showModal
#' @importFrom purrr map
#' @importFrom grid rasterGrob
#' @importFrom magick geometry_size_pixels image_read
#' @importFrom zeallot %<-%
#'
"_PACKAGE"

#' Manual annotation tool via shiny
#' This function takes a seurat object with stored image locations and opens up the manual selection tool in the default browser
#'
#' @param object Seurat object
#' @param type Input image type, e.g. "masked" or "raw"
#' @param res Tissue HE image width in pixels
#' @param verbose Print messages
#'
#' @export
#'


ManualAnnotation <- function (
  object,
  type = NULL,
  res = 1000,
  verbose = FALSE
) {
  
  if (!requireNamespace("ggiraph")) install.packages("ggiraph")
  if (!requireNamespace("tibble")) install.packages("tibble")

  if (!"Staffli" %in% names(object@tools)) stop("Staffli object not present in Seurat object", call. = FALSE)
  st.object <- GetStaffli(object)
  #img <- st.object@imgs[1]
  type <- type %||% {
    if (is.null(rasterlists(st.object))) stop("There are no images present in the Seurat object. Run LoadImages() first.", call. = FALSE)
    choices <- c("processed", "masked", "raw", "processed.masks", "masked.masks")
    match.arg(choices, rasterlists(st.object), several.ok = T)[1]
  }
  img <- st.object@rasterlists[[type]][[1]]
  sampleChoice <- unique(st.object[[, "sample", drop = T]])

  # ===================== UI =======================
  ui <-  pageWithSidebar(
    headerPanel("Manual selection"),



    sidebarPanel(width = 3,
      actionButton(inputId="info", label="Instructions"),  #icon = shiny::icon("info", lib="glyphicon")),
      shiny::hr(),
      selectInput(inputId = "sampleInput", label = "Select sample", choices = sampleChoice, selected = "1"),
      shiny::hr(),
      textInput(inputId = "labelInput", label = "Choose label name", value="", placeholder = "Default"),
      shiny::hr(),
      sliderInput(inputId="alphaValue", label="Opacity [0-1]", min=0, max=1, value=0.7, step=0.2),
      shiny::hr(),
      sliderInput(inputId="spotSize", label="Capture-Spot Size [1-5]", min=0, max=5, value=2.5, step=0.5),
      shiny::hr(),
      actionButton(inputId = "confirm", label="Confirm selection"),
      shiny::hr(),
      actionButton(inputId = "stopApp", label="Quit annotation tool")
    ),
    mainPanel(
      #tags$img(src = jsonlite::base64_enc(img)),
      ggiraph::ggiraphOutput("Plot1", width = "100%", height = paste0(res, "px"))
    )
  )

  # ===================== Server ================================
  server <- function(input, output, session){

    df <- reactiveValues(label = object@meta.data$labels,
                         id = object@meta.data$id,
                         sample = st.object[[, "sample", drop = T]])

    rv <- reactiveValues(sNr = "1", ann = NULL)
    observeEvent(input$sampleInput, {
      rv$sNr = input$sampleInput
      rv$ann <- Create_annotation(object, type, input$sampleInput, res)
    })

    output$Plot1 <- ggiraph::renderGirafe({

     # print(table(df$label[which(df$sample == input$sampleInput)]))
      withProgress(message = "Updating plot", value = 0, {
      x <- ggiraph::girafe(ggobj = make.plot(object,
                                             sampleNr = rv$sNr,
                                             spotAlpha = input$alphaValue,
                                             Labels = df$label[which(df$sample == rv$sNr)],
                                             res = res,
                                             SpotSize = input$spotSize,
                                             ann = rv$ann), width_svg = 12, height_svg = 10)
      x <- ggiraph::girafe_options(x,
                                   ggiraph::opts_zoom(max = 5),
                                   ggiraph::opts_selection(type = "multiple",
                                                           css = "fill:cyan;stroke:black;opacity:0.7;"))
      x
      })
    })

    observeEvent(input$confirm, {
      ids.selected <- as.numeric(input$Plot1_selected)
      df$label[which(df$id %in% ids.selected)] <- input$labelInput
      session$sendCustomMessage(type = 'Plot1_set', message = character(0))
    })

    observe({
      if(input$stopApp > 0){
        print("Stopped")
        object@meta.data$labels <- df$label
        stopApp(returnValue = object)
      }
    })

    observeEvent(input$info, {
      showModal(modalDialog(
        title = "Instructions",
        HTML("1. Select sample (might take a while to load)<br>",
             "2. Specifiy label you want to use<br>",
             "3. Use the blue(select) lasso to label the caputure-spots<br>",
             "4. Press Confirm to set the labels<br>",
             "5. Repeat 1-4 until all labels are set<br>",
             "6. Close the shiny tool to return"),
        easyClose = TRUE,
        footer = NULL
      ))
    })
  }

  runApp(list(ui = ui, server = server), launch.browser = T)
}

#' Used for the manual annotation tool, returns plot and coordinate IDs for selected sample
#'
#' @param object Seurat object
#' @param sampleNr Sample to be plotted
#' @param spotAlpha geom_point opacity
#' @param Labels labels included in the input seurat object
#' @param SpotSize geom_point size
#' @param res resolution of the image
#'
#' @importFrom zeallot %<-%
#' @importFrom magick image_read image_scale image_info geometry_size_pixels
#' @importFrom ggplot2 ggplot aes scale_x_continuous scale_y_continuous theme element_rect element_blank theme_minimal
#'
#' @keywords internal

make.plot <- function (
  object,
  sampleNr,
  spotAlpha,
  Labels,
  SpotSize,
  res,
  ann
) {

  c(ann, coordinates, xmin, xmax, ymin, ymax) %<-% ann

  coordinates$Labels <- Labels

  gg <- ggplot(coordinates, aes(x = x, y = y, data_id = id)) +
    ann +
    ggiraph::geom_point_interactive(size = SpotSize, alpha = spotAlpha, aes(col = Labels)) +
    coord_fixed() +
    scale_x_continuous(expand = c(0, 0), limits = c(xmin, xmax)) +
    scale_y_continuous(expand = c(0, 0), limits = c(ymin, ymax)) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "transparent", colour = NA),
      legend.position = "bottom",
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank()
    )
  return(gg)
}

#' Create an annotation
#' 
#' @param object Seurat object
#' @param type Image type, e.g. "raw" or "masked"
#' @param sampleNr Selected section/sample number
#' @param res Image resolution
#' 
#' @importFrom stats setNames dist
#' 
Create_annotation <- function (
  object,
  type,
  sampleNr,
  res
) {

  st.object <- GetStaffli(object)
  px.ids <- ifelse(rep(type %in% c("raw", "masked", "masked.masks"), 2), c("pixel_x", "pixel_y"), c("warped_x", "warped_y"))
  st.meta <- subset(cbind(st.object@meta.data[, c(px.ids, "sample")], id = object[[]][, "id", drop = T]), sample %in% paste0(sampleNr))
  coordinates <- setNames(st.meta[, c(px.ids, "id")], nm = c("x", "y", "id"))

  image <- image_read(st.object@rasterlists[[type]][[as.integer(sampleNr)]])
  old_width <- st.object@dims[[as.integer(sampleNr)]]$width
  cur_width <- image_info(image)$width
  #old_width <- image_info(image)$width
  #image <- image_scale(image, geometry_size_pixels(width = res, preserve_aspect = T))
  coordinates[, c("x","y")] <- coordinates[, c("x","y")]*(cur_width/old_width)

  r <- min(dist(coordinates[, c("x","y")])) / 2

  c(ymin, ymax) %<-% range(coordinates$y)
  c(xmin, xmax) %<-% range(coordinates$x)

  xmin <- ymin <- min(c(xmin, ymin))
  xmax <- ymax <- max(c(xmax, ymax))

  c(ymin, xmin) %<-% { c(ymin, xmin) %>% map(~. - 3 * r) }
  c(ymax, xmax) %<-% { c(ymax, xmax) %>% map(~. + 3 * r) }

  image <- as.raster(image)

  image <- grid::rasterGrob(
    image,
    width = unit(1, "npc"),
    height = unit(1, "npc"),
    interpolate = TRUE
  )

  if (!is.null(image)) {
    ymin <- max(ymin, 1)
    ymax <- min(ymax, nrow(image$raster))
    xmin <- max(xmin, 1)
    xmax <- min(xmax, ncol(image$raster))
    image$raster <- image$raster[ymin:ymax, xmin:xmax]
    annotation <- ggplot2::annotation_custom(image, -Inf, Inf, -Inf, Inf)
  } else {
    annotation <- NULL
  }

  coordinates$y <- ymax - coordinates$y + ymin
  return(list(annotation, coordinates, xmin, xmax, ymin, ymax))
}
