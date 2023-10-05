library(shiny)
library(shinythemes)
library(ggplot2)
library(viridis)
library(dplyr)


huge_merge=readRDS("./huge_merge.RDS")
huge_merge_pro=readRDS("./huge_merge_pro.RDS")



##this is the current final version

unique_celltypes <- unique(huge_merge$celltype.x)
unique_TMACores <- c("All",unique(huge_merge$TMACore.x))
unique_Images <- c("All",unique(huge_merge$Image.x))
unique_markers <-unique(huge_merge$marker)
unique_patient = c("All",unique(huge_merge$Parent))


visualize_cell_expression <- function(data, celltype, source_TMA = NULL, source_Image = NULL, marker_to_plot = "CD8", point_size=2, parent_filter = NULL) {
  
  # Filter data based on celltype and source
  filtered_data <- data
  
  if(!is.null(parent_filter) && parent_filter != "All") {
    filtered_data <- filtered_data %>% filter(Parent == parent_filter)
  } else {
    if(!is.null(source_TMA)) {
      filtered_data <- filtered_data %>% filter(TMACore.x == source_TMA)
    }
    
    if(!is.null(source_Image)) {
      filtered_data <- filtered_data %>% filter(Image.x == source_Image)
    }
  }
  
  if(!is.null(marker_to_plot)) {
    filtered_data <- filtered_data %>% filter(marker == marker_to_plot)
  }
  
  # Define the shape based on the target celltype
  filtered_data$shape <- ifelse(filtered_data$celltype.x == celltype, "Target", "Other")
  
  # Create the plot
  g <- ggplot(filtered_data, aes(x = X, y = Y)) +
    geom_point(aes(color = expression), size = point_size) +
    scale_color_viridis_c(limits = range(filtered_data$expression, na.rm = TRUE), 
                          option = "viridis", end = 0.9) +
    labs(title = paste("Expression of", marker_to_plot, "in", celltype),
         x = "X Coordinate",
         y = "Y Coordinate",
         color = "Expression Level") +
    facet_wrap(~ shape, ncol = 1) +
    theme_light() +
    theme(legend.position = "bottom", legend.title = element_text(face = "bold"), legend.text = element_text(face = "italic"), aspect.ratio = 1)
  
  return(g)
}

visualize_cell_expression_pro <- function(data, celltype, source_TMA = NULL, source_Image = NULL, marker_to_plot = "CD8", point_size = 2, parent_filter = NULL) {
  
  # Filter data based on celltype and source
  filtered_data <- data
  
  if(!is.null(parent_filter) && parent_filter != "All") {
    filtered_data <- filtered_data %>% filter(Parent == parent_filter)
  } else {
    if(!is.null(source_TMA)) {
      filtered_data <- filtered_data %>% filter(TMACore.x == source_TMA)
    }
    
    if(!is.null(source_Image)) {
      filtered_data <- filtered_data %>% filter(Image.x == source_Image)
    }
  }
  
  if(!is.null(marker_to_plot)) {
    filtered_data <- filtered_data %>% filter(marker == marker_to_plot)
  }
  
  # Assign colors based on probability values
  filtered_data$color_class <- cut(filtered_data$expression, 
                                   breaks = c(-Inf, 0.5, 0.6, 0.7, 0.8, Inf), 
                                   labels = c("<=0.5", ">0.5", ">0.6", ">0.7", ">0.8"))
  
  # Define the shape based on the target celltype
  filtered_data$shape <- ifelse(filtered_data$celltype.x == celltype, "Target", "Other")
  
  
  # Create the plot
  g <- ggplot(filtered_data, aes(x = X, y = Y)) +
    geom_point(aes(color = color_class), size = point_size) +
    scale_color_manual(values = c("<=0.5" = "lightblue", ">0.5" = "deepskyblue", ">0.6" = "dodgerblue", ">0.7" = "blue", ">0.8" = "darkblue")) +
    labs(title = paste("Expression of", marker_to_plot, "in", celltype),
         x = "X Coordinate",
         y = "Y Coordinate") +
    facet_wrap(~ shape, ncol = 1) +
    theme_light() +
    theme(legend.position = "bottom", legend.title = element_text(face = "bold"), legend.text = element_text(face = "italic")) +
    coord_fixed(ratio = 1)
  
  return(g)
}

ui <- fluidPage(
  theme = shinytheme("flatly"),  
  titlePanel("Cell Expression Visualization"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput("celltype", "Select Cell Type:", unique_celltypes),
      selectInput("source_TMA", "Select TMA Core:", unique_TMACores),
      selectInput("source_Image", "Select Image:", unique_Images),
      selectInput("parent_filter", "Select Parent:", unique_patient),
      selectInput("marker_to_plot", "Select Marker to Plot:", unique_markers),
      sliderInput("imgScale", "Image Scale:", min = 0.5, max = 2, value = 1, step = 0.1),
      sliderInput("pointSize", "Point Size:", min = 1, max = 5, value = 2, step = 0.1),
      actionButton("plotButton", "Generate Plot")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Expression Visualization", plotOutput("expressionPlot", height = "auto", width = "auto")),
        tabPanel("Expression Probability Visualization", plotOutput("expressionProbPlot", height = "auto", width = "auto"))
      )
    )
  )
)

server <- function(input, output) {
  
  observeEvent(input$plotButton, {
    output$expressionPlot <- renderPlot({
      g <- visualize_cell_expression(huge_merge, input$celltype, 
                                     source_TMA = input$source_TMA, 
                                     source_Image = input$source_Image, 
                                     marker_to_plot = input$marker_to_plot,
                                     point_size = input$pointSize,
                                     parent_filter = input$parent_filter)
      print(g)
    }, width = 600 * input$imgScale, height = 600 * input$imgScale)
  })
  
  observeEvent(input$plotButton, {
    output$expressionProbPlot <- renderPlot({
      g <- visualize_cell_expression_pro(huge_merge_pro, input$celltype, 
                                         source_TMA = input$source_TMA, 
                                         source_Image = input$source_Image, 
                                         marker_to_plot = input$marker_to_plot, 
                                         point_size = input$pointSize,
                                         parent_filter = input$parent_filter)
      print(g)
    }, width = 600 * input$imgScale, height = 600 * input$imgScale)
  })
}

shinyApp(ui = ui, server = server)
