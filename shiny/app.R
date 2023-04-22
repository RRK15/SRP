library(shiny)
library(shinydashboard)
library(plotly)
library(ggplot2)
library(dplyr)
library(mclust)
library(Rtsne)



altclust <- read.csv("data/cluster_genes.csv")

metadata <- read.csv("data/SraRunTableMod.csv", row.names = 1)
features <-subset(metadata, select = c(Cell_type))

distance_matrix <- read.csv("data/distance_matrix.csv", row.names = "X")

altdata <- read.csv("data/alternate.csv")

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
  sidebar <- dashboardSidebar(
    sidebarMenu(
      menuItem("Welcome!", tabName = "welcome", icon = icon("dashboard")),
      menuItem("Original pipeline", icon = icon("th"), tabName = "originalpipeline"),
      menuItem("Alternate pipeline", icon = icon("th"), tabName = "alternatepipeline")
    )
  ),
  
  dashboardBody(
    tabItems(
      tabItem(tabName = "originalpipeline",
    fluidRow(
      column(width = 8,
             box(plotlyOutput("mclustplot"), width = NULL)
      ),
      column(width = 4,
             box("Top 10 enriched genes in cluster:",uiOutput("PLS"), width = NULL),
             
      )
    )
  ),
      tabItem(tabName = "alternatepipeline",
        fluidRow(
          column(width = 8,
              box(plotlyOutput("alternateplot"), width = NULL)
          ),
          column(width = 4,
              box("Top 10 enriched genes in cluster:", uiOutput("altclust"), width = NULL),
              box(textInput("querygene",label = "Query a gene"), width = NULL),
              box(textOutput("geneoutput"), width = NULL)
          )
        )
)
)
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
    testing %>%
      highlight_key(~mclust) %>%
      plot_ly(x = ~tsne1, y = ~tsne2, 
              color = ~Cell_type,
              colors = "Set3",
              type = "scatter",
              mode = "markers",
              source = "testing",
              customdata = ~mclust,
              text = ~paste("Cluster ID:", mclust, 
                            "<br>Cell Type:", Cell_type),
              hoverinfo = "text")%>%
      highlight(on = "plotly_hover", off = "plotly_doubleclick") 
  })
  
  observe({
    click_data <- event_data(event = "plotly_click", source = "altdata")
    clickval <- as.integer(click_data$key) 
    val(clickval)

  })
  
  output$altclust <- renderUI({
    altclust[val(),]
  })
  
  output$alternateplot <- renderPlotly({
    altdata %>%
      highlight_key(~clusters) %>%
      plot_ly(x = ~UMAP_1, y = ~UMAP_2, 
              type = 'scatter', 
              color = ~Cell_type, 
              colors = "Set3",
              source = "altdata",
              mode = 'markers',  
              text = ~paste("Cluster ID:", clusters,"<br>Cell Type:", Cell_type), 
              hoverinfo = "text") %>%
      highlight(on = "plotly_hover", off = "plotly_doubleclick")
    
  })
  output$PLS <- renderUI({
    testdata[val(),]
  })
  
  querytest <- reactive({
    req(input$querygene)
    grepstuff = grep(input$querygene, altclust$genes)
    if (identical(grepstuff, integer(0))){
      print("No genes found")
      
    }else{
      results = paste(grep(input$querygene, altclust$genes), collapse = " ")
      paste("Genes found in clustert(s):",results, collapse = " ")
    }
  })
  
  output$geneoutput <- renderText({
    querytest()
  })
} 

shinyApp(ui, server)