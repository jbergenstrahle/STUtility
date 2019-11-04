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
#' @importFrom ggiraph girafe renderGirafe
#' geom_point_interactive ggiraphOutput
#' girafe_options
#' @importFrom shiny pageWithSidebar headerPanel sidebarPanel mainPanel textInput strong actionButton radioButtons sliderInput
#' reactiveValues observeEvent observe hr submitButton
#' runApp
#' @importFrom tibble column_to_rownames rownames_to_column
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
#'
#' @export
#'


ManualAnnotation <- function (
  object,
  res = 1500,
  verbose = FALSE
) {

  if (!"Staffli" %in% names(object@tools)) stop("Staffli object not present in Seurat object", call. = FALSE)
  st.object <- GetStaffli(object)
  sampleChoice <- unique(st.object[[, "sample", drop = T]])

  # ===================== UI =======================
  ui <-  pageWithSidebar(
    headerPanel("Manual selection"),



    sidebarPanel(
    actionButton(inputId="info", label="Instructions"),  #icon = shiny::icon("info", lib="glyphicon")),
    shiny::hr(),
    selectInput(inputId = "sampleInput", label = "Select sample", choices = sampleChoice, selected = "1"),
    shiny::hr(),
    textInput(inputId = "labelInput", label = "Choose label name", value="", placeholder = "Default"),
    shiny::hr(),
    sliderInput(inputId="alphaValue", label="Opacity [0-1]", min=0, max=1, value=0.2, step=0.2),
    shiny::hr(),
    sliderInput(inputId="spotSize", label="Capture-Spot Size [1-5]", min=0, max=5, value=2, step=1),
    shiny::hr(),
    actionButton(inputId = "confirm", label="Confirm selection"),
    shiny::hr(),
    actionButton(inputId = "stopApp", label="Quit annotation tool")
    ),
    mainPanel(
      ggiraphOutput("Plot1", height = res)
    )
  )

  # ===================== Server ================================
    server <- function(input, output, session){

      df <- reactiveValues(label = object@meta.data$labels,
                           id = object@meta.data$id,
                           sample = st.object[[, "sample", drop = T]])

      output$Plot1 <- ggiraph::renderGirafe({

        x <- ggiraph::girafe(ggobj = make.plot(object,
                                               sampleNr = input$sampleInput,
                                               spotAlpha = input$alphaValue,
                                               Labels = df$label[which(df$sample == input$sampleInput)],
                                               res = res,
                                               SpotSize = input$spotSize))
        x <- ggiraph::girafe_options(x,
                                     ggiraph::opts_zoom(max=5),
                                     ggiraph::opts_selection(type = "multiple",
                                                             css = "fill:red;stroke:black;r:2pt;" ))
        x
      })

      observeEvent(input$confirm, {
        ids.selected <- as.numeric(input$Plot1_selected)
        df$label[which(df$id %in% ids.selected)] <- input$labelInput
        session$sendCustomMessage(type = 'Plot1_set', message = character(0))
      })

      observe({
        if(input$stopApp > 0){
          print("Stopped")
          object@meta.data$labels <-df$label
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
#' @keywords internal

make.plot <- function (
  object,
  sampleNr,
  spotAlpha,
  Labels,
  SpotSize,
  res
) {
  object.use <- colnames(object[, which(st.object[[, "sample", drop = T]] == sampleNr)])
  object <- SubsetSTData(object, spots = object.use)
  st.object <- GetStaffli(object)
  coordinates <- data.frame(x = st.object[[, "pixel_x", drop = T]],
                            y = st.object[[, "pixel_y", drop = T]],
                            id = object@meta.data$id,
                            Labels = object@meta.data$labels)
  image <- image_read(st.object@imgs)
  old_width <- image_info(image)$width
  image <- image_scale(image, geometry_size_pixels(width = res, preserve_aspect = T))
  coordinates[, c("x","y")] <- coordinates[, c("x","y")]*(res/old_width)

  r <- min(dist(coordinates[, c("x","y")])) / 2

  c(ymin, ymax) %<-% range(coordinates$y)
  c(xmin, xmax) %<-% range(coordinates$x)
  c(ymin, xmin) %<-% { c(ymin, xmin) %>% map(~. - 3 * r) }
  c(ymax, xmax) %<-% { c(ymax, xmax) %>% map(~. + 3 * r) }

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

  coordinates$y <- ymax - coordinates$y+ ymin
  gg <- ggplot(coordinates, aes(x = x, y = y, data_id = id)) +
    annotation +
    ggiraph::geom_point_interactive(size = SpotSize, alpha = spotAlpha, aes(col = Labels)) +
    #coord_fixed() +
    scale_x_continuous(expand = c(0, 0), limits = c(xmin,xmax)) +
    scale_y_continuous(expand = c(0, 0), limits = c(ymin,ymax)) +
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





