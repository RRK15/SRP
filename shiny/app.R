library(shiny)
library(shinydashboard)
library(plotly)
library(ggplot2)

pls <-ggplot(data=iris)+
  geom_point(mapping=aes(x=Petal.Length,y=Petal.Width))
pls2 <- ggplotly(pls)
testdata <- read.csv("cluster_genes.csv")
ui <- dashboardPage(
  dashboardHeader(title = "Basic dashboard"),
  dashboardSidebar(),
  dashboardBody(
    box(uiOutput("testing"), width = 4, height = 8),
    box(
      title = "Select a cluster",
      sliderInput("slider", "Clusters:", 1, 10, 5)
    ),
    box(plotlyOutput("iris"), width = 4)
    
  )
  )

server <- function(input, output) {
  output$testing <- renderUI({
    testdata[input$slider,]
  })
  output$iris <- renderPlotly({
    pls2
  })
}

shinyApp(ui, server)
