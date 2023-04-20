library(shiny)
library(shinydashboard)
library(plotly)
library(ggplot2)
library(dplyr)
library(mclust)
library(Rtsne)

testdata <- read.csv("cluster_genes.csv")
distance_matrix <- read.csv("distance_matrix.csv", row.names = "X")
#replace all nans with 0
distance_matrix[is.na(distance_matrix)] <- 0
distance_matrix <- as.matrix(distance_matrix)
normalised_distance = normalize_input(distance_matrix)
output_tsne <- Rtsne(normalised_distance, theta = 0.0)
needed_output <- as.data.frame(output_tsne$Y)
clustering <- Mclust(needed_output)
tsne1 = output_tsne$Y[,1]
tsne2 = output_tsne$Y[,2]
testing = data.frame(tsne1, tsne2)
testing$mclust = factor(clustering$classification)


ui <- dashboardPage(
  dashboardHeader(title = "Basic dashboard"),
  dashboardSidebar(),
  dashboardBody(
    box("Top 20 enriched genes in cluster",uiOutput("testing"), width = 4, height = 8),
    box(plotlyOutput("mclustplot"), width = 8)
  )
  )

server <- function(input, output) {
  val = reactiveVal(1)
  observe({
    testingclick <- event_data(event = "plotly_click", source = "testing")
    clickval <- testingclick$customdata
    val(clickval)
    
  })
  output$mclustplot <- renderPlotly({
    p <- plot_ly(testing, x = ~tsne1, y = ~tsne2, source = "testing", color = ~mclust,colors = "Set3",
            customdata = ~mclust,
            text = ~paste("Cluster ID:", mclust))
  })
  
 
  output$testing <- renderUI({
    testdata[val(),]
  })
} 

shinyApp(ui, server)
