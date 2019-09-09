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
#' @importFrom shiny pageWithSidebar headerPanel sidebarPanel mainPanel textInput strong actionButton radioButtons
#' reactiveValues observeEvent observe hr
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


ST.annotation <- function (
  object,
  res = 1500,
  verbose = FALSE
) {

  sampleChoice <- unique(object@meta.data$sample)

  # ===================== UI =======================
  ui <-  pageWithSidebar(
    headerPanel("Manual selection"),

    sidebarPanel(
    selectInput(inputId = "sampleInput", label = "1. Select sample", choices = sampleChoice, selected = "1"),
    shiny::hr(),
    textInput(inputId = "labelInput", label = "2. Choose label name", value="", placeholder = "Default"),
    strong("3. Use lasso tool to select regions"),
    shiny::hr(),
    #shiny::radioButtons(inputId = "selectionColor", label = "Color", choices = c("Red",
     #                                                                    "Green",
      #                                                                   "Blue",
       #                                                                  "Yellow",
        #                                                                 "Black",
         #                                                                "White"), width = NULL),
    numericInput(inputId="alphaValue", label="opacity[0-1]", min=0, max=1, value=0.2, step=0.2),
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

      df <- reactiveValues(label=object@meta.data$labels, id=object@meta.data$id,
                           sample=object@meta.data$sample)

      output$Plot1 <- ggiraph::renderGirafe({

        x <- ggiraph::girafe(ggobj = make.plot(object,
                                               sampleNr = as.numeric(input$sampleInput),
                                               spotAlpha = input$alphaValue,
                                               Labels = df$label[which(df$sample == input$sampleInput)]))
        x <- ggiraph::girafe_options(x,
                                     ggiraph::opts_zoom(max=5),
                                     ggiraph::opts_selection(type = "multiple",
                                                             css = "fill:red;stroke:black;r:2pt;" ))
        x
      })

      observeEvent(input$confirm, {
        ids.selected <- as.numeric(input$Plot1_selected)
        df$label[which(df$id %in% ids.selected)] <- input$labelInput

        #print("we have the following labels: ")
        #print(unique(df$label))
        session$sendCustomMessage(type = 'Plot1_set', message = character(0))
      })

      observe({
        if(input$stopApp > 0){
          print("Stopped")
          object@meta.data$labels <-df$label
          stopApp(returnValue = object)
        }
      })
    }

  runApp(list(ui = ui, server = server), launch.browser = T)
}

#' Used for the manual annotation tool, returns plot and coordinate IDs for selected sample
#'
#' @param object Seurat object
#' @param sampleNr Sample to be plotted
#'
#' @keywords internal

make.plot <- function(
  object,
  sampleNr,
  spotAlpha,
  Labels,
  res=1500
) {
  object.use <- colnames(object[, which(object$"sample" == sampleNr)])
  object <- subset(object, cells = object.use)
  coordinates <- data.frame(x=object@meta.data$pixel_x,
                            y=object@meta.data$pixel_y,
                            id=object@meta.data$id,
                            label=object@meta.data$labels)
  image <- image_read(object@tools$imgs[sampleNr])
  old_width <- image_info(image)$width
  image <- image_scale(image, geometry_size_pixels(width=res, preserve_aspect = T))
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
  gg <- ggplot(coordinates, aes(x=x, y=y, data_id=id)) +
    annotation+
    ggiraph::geom_point_interactive(size = 3, alpha=spotAlpha, aes(col=Labels)) +
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





