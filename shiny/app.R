library(shiny)
library(shinydashboard)
library(plotly)
library(ggplot2)
library(dplyr)
library(mclust)
library(Rtsne)

testdata <- read.csv("cluster_genes.csv")
metadata <- read.csv("SraRunTableMod.csv", row.names = 1)
features <-subset(metadata, select = c(Cell_type))

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
combin = data.frame(tsne1, tsne2)
testing <- cbind(combin, features)
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
    plot_ly(testing, x = ~tsne1, y = ~tsne2, 
            color = ~Cell_type,
            colors = "Set3",
            type = "scatter",
            mode = "markers",
            source = "testing",
            customdata = ~mclust,
            legend = FALSE,
            text = ~paste("Cluster ID:", mclust, 
                          "<br>Cell Type:", Cell_type),
            hoverinfo = "text")
    
  })
  
 
  output$testing <- renderUI({
    testdata[val(),]
  })
} 

shinyApp(ui, server)
