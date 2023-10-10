library(shiny)
library(shinythemes)
library(dplyr)
library(ggplot2)
library(viridis)


##this is the final version

huge_merge=readRDS("huge_merge_final.RDS")
huge_merge_pro=readRDS("huge_merge_pro_final.RDS")

unique_celltypes <- unique(huge_merge$celltype)
unique_TMACores <- c("All",unique(huge_merge$TMACore))
unique_Images <- c("All",unique(huge_merge$Image))
unique_markers <-unique(huge_merge$marker)
unique_patient = c("All",unique(huge_merge$Patient_ID))

visualize_cell_expression <- function(data, celltype, source_TMA = NULL, source_Image = NULL, marker_to_plot = "CD8", point_size=2, parent_filter = NULL, upper=0.95) {
  
  # Filter data based on celltype and source
  filtered_data <- data
  
  if(!is.null(parent_filter) && parent_filter != "All") {
    filtered_data <- filtered_data %>% filter(Patient_ID == parent_filter)
  } 
  
  if(!is.null(source_TMA) && source_TMA != "All") {
    filtered_data <- filtered_data %>% filter(TMACore == source_TMA)
  }
  
  if(!is.null(source_Image) && source_Image != "All") {
    filtered_data <- filtered_data %>% filter(Image == source_Image)
  }
  
  
  if(!is.null(marker_to_plot)) {
    filtered_data <- filtered_data %>% filter(marker == marker_to_plot)
  }
  
  # Define the shape based on the target celltype
  filtered_data$shape <- ifelse(filtered_data$celltype == celltype, "Target", "Other")
  
  # Determine color scale limits
  color_lower_limit <- min(filtered_data$expression, na.rm = TRUE)
  
  # Set the upper limit to the specified quantile (default is 95%)
  color_upper_limit <- quantile(filtered_data$expression, upper, na.rm = TRUE)
  
  # Create the plot
  g <- ggplot(filtered_data, aes(x = X, y = Y)) +
    geom_point(aes(color = expression), size = point_size) +
    scale_color_viridis_c(limits = c(color_lower_limit, color_upper_limit), 
                          option = "viridis", end = 0.9) +
    labs(title = paste("Expression of", marker_to_plot, "in", celltype),
         x = "X Coordinate",
         y = "Y Coordinate",
         color = "Expression Level") +
    facet_wrap(~ shape, ncol = 1) +
    theme_light() +
    theme(legend.position = "bottom", 
          legend.title = element_text(face = "bold"), 
          legend.text = element_text(face = "italic")) +
    coord_fixed(ratio = 1)
  
  return(g)
}


visualize_cell_expression_pro <- function(data, celltype, source_TMA = NULL, source_Image = NULL, marker_to_plot = "CD8", point_size = 2, parent_filter = NULL) {
  
  # Filter data based on celltype and source
  filtered_data <- data
  
  if(!is.null(parent_filter) && parent_filter != "All") {
    filtered_data <- filtered_data %>% filter(Patient_ID == parent_filter)
  } 
  
  
  if(!is.null(source_TMA) && source_TMA != "All") {
    filtered_data <- filtered_data %>% filter(TMACore == source_TMA)
  }
  
  if(!is.null(source_Image) && source_Image != "All") {
    filtered_data <- filtered_data %>% filter(Image == source_Image)
  }
  
  if(!is.null(marker_to_plot)) {
    filtered_data <- filtered_data %>% filter(marker == marker_to_plot)
  }
  
  # Assign colors based on probability values
  filtered_data$color_class <- cut(filtered_data$expression, 
                                   breaks = c(-Inf, 0.5, 0.6, 0.7, 0.8, Inf), 
                                   labels = c("<=0.5", ">0.5", ">0.6", ">0.7", ">0.8"))
  
  # Define the shape based on the target celltype
  filtered_data$shape <- ifelse(filtered_data$celltype == celltype, "Target", "Other")
  
  
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
# UI Definition
ui <- fluidPage(
  theme = shinytheme("flatly"),
  titlePanel("Cell Expression Visualization"),
  
  sidebarLayout(
    sidebarPanel(
      radioButtons("selectMode", "Choose Selection Mode:",
                   choices = c("Patient & TMA Core", "Image"), 
                   selected = "Patient & TMA Core", inline = TRUE),
      conditionalPanel(
        condition = "input.selectMode == 'Patient & TMA Core'",
        selectInput("parent_filter", "Select Patient:", unique_patient),
        selectInput("source_TMA", "Select TMA Core:", unique_TMACores)
      ),
      conditionalPanel(
        condition = "input.selectMode == 'Image'",
        selectInput("source_Image", "Select Image:", unique_Images)
      ),
      selectInput("celltype", "Select Cell Type:", unique_celltypes),
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

# Server Logic
server <- function(input, output, session) {
  
  # Observer to update TMA core choices based on selected patient
  observe({
    patient_selected <- input$parent_filter
    if (patient_selected != "All") {
      tm_cores <- c("All", unique(huge_merge$TMACore[huge_merge$Patient_ID == patient_selected]))
      updateSelectInput(session, "source_TMA", choices = tm_cores)
    } else {
      updateSelectInput(session, "source_TMA", choices = unique_TMACores)
    }
  })
  
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

# Run the Shiny App
shinyApp(ui = ui, server = server)
