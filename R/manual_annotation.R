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
#' @importFrom ggiraph girafe renderGirafe geom_point_interactive ggiraphOutput
#' @importFrom shiny pageWithSidebar headerPanel sidebarPanel mainPanel
#' textInput strong actionButton
#' reactiveValues observeEvent observe
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
  verbose = FALSE #ADD ability to change resolution!
) {

  sampleChoice <- unique(object@meta.data$sample)

  # ===================== UI =======================
  ui <-  pageWithSidebar(
    # Set theme
    #theme = shinytheme("spacelab"),

    # Title text
    headerPanel("Manual selection"),

    sidebarPanel(
    selectInput(inputId = "sampleInput", label = "1. Select sample", choices = sampleChoice, selected = "1"),
    actionButton(inputId = "loadIMG", label="Load Image"),
    hr(),
    textInput(inputId = "labelInput", label = "2. Choose label name", value="", placeholder = "Default"),
    strong("3. Use lasso tool to select regions"),
    hr(),
    actionButton(inputId = "confirm", label="Confirm selection"),
    hr(),
    actionButton(inputId = "stopApp", label="Quit annotation tool")
    ),
    mainPanel(
      ggiraphOutput("Plot1", height = "1000px")
    )
  )

  # ===================== Server ================================
    server <- function(input, output, session){

      df <- reactiveValues()
      df$annotation <- data.frame(label=object@meta.data$labels, id=object@meta.data$id, stringsAsFactors=F)
      # Observes the second feature input for a change
      observeEvent(input$loadIMG,{
        gg <- get.anno.plot(object, sampleNr = as.numeric(input$sampleInput))
        #df$cords <- anno[[1]]
        output$Plot1 <- ggiraph::renderGirafe({

          #gg <- anno[[2]]
          x <- girafe(ggobj = gg)
          x <- girafe_options(x,
                              opts_zoom(max=5),
                              opts_selection(type = "multiple") )
          x

        })
      })


      observeEvent(input$Plot1_selected, {
        print(paste("Looking at sample: ", input$sampleInput))
        print(paste("Setting label: ", input$labelInput))
        print(paste("Nr of spots selected: ", length(input$Plot1_selected)))
        ids.selected <- as.numeric(input$Plot1_selected)
        df$annotation[which(df$annotation$id %in% ids.selected), ]$label <- input$labelInput
        })

      observe({
        if(input$stopApp > 0){
          print("Stopped")
          object@meta.data$labels <-df$annotation$label
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

get.anno.plot <- function(
  object,
  sampleNr
) {
  outObs <- list()
  object.use <- colnames(object[, which(object$"sample" == sampleNr)])
  object <- subset(object, cells = object.use)
  coordinates <- data.frame(x=object@meta.data$pixel_x,
                            y=object@meta.data$pixel_y,
                            id=object@meta.data$id)
  image <- image_read(object@tools$imgs[sampleNr])
  old_width <- image_info(image)$width
  image <- image_scale(image, geometry_size_pixels(width=1000, preserve_aspect = T))
  coordinates[, c("x","y")] <- coordinates[, c("x","y")]*(1000/old_width) #kan andra ordning

  r <- min(dist(coordinates)) / 2

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
  image$raster <- image$raster[ymin:ymax, xmin:xmax]

  if (!is.null(image)) {
    ymin <- max(ymin, 1)
    ymax <- min(ymax, nrow(image$raster))
    xmin <- max(xmin, 1)
    xmax <- min(xmax, ncol(image$raster))

    annotation <- ggplot2::annotation_custom(image, -Inf, Inf, -Inf, Inf)
  } else {
    annotation <- NULL
  }

  coordinates$y <- ymax - coordinates$y + ymin

  gg <- ggplot(coordinates, aes(x=x, y=y, data_id=id)) +
    annotation+
    ggiraph::geom_point_interactive(size = 2, alpha=0) +
   # coord_fixed() +
    scale_x_continuous(expand = c(0, 0), limits = c(xmin, xmax)) +
    scale_y_continuous(expand = c(0, 0), limits = c(ymin, ymax)) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "transparent", colour = NA),
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank()
    )

  coordinates$sample <- sampleNr
  outObs[[1]] <- coordinates #REMOVE ALL THESE AND ONLY RETURN GG
  outObs[[2]] <- gg

  return(gg)
}





