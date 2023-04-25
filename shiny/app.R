library(shiny)
library(shinydashboard)
library(plotly)
library(ggplot2)
library(dplyr)
library(mclust)
library(Rtsne)


origclust <- read.csv("data/origtop.csv")
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
  dashboardHeader(title = "Group A submission website"),
  sidebar <- dashboardSidebar(
    sidebarMenu(
      menuItem("Welcome!", tabName = "welcome", icon = icon("dashboard")),
      menuItem("Original pipeline", icon = icon("th"), tabName = "originalpipeline"),
      menuItem("Alternate pipeline", icon = icon("th"), tabName = "alternatepipeline")
    )
  ),
  
  dashboardBody(
    tabItems(
      tabItem(tabName = "welcome",
              fluidRow(
                column(width = 12,align = "center",
                       titlePanel("Group A submission")
                ),
              box("In 2015 a study was carried out by Damian in which an attempt was carried out to build a brain cell atlas. Cell samples were obtained from the brain and submitted to a scRNAseq annalysis to study their transcriptomes. 

The following is an interactive site that allows for the visualization of the analysis of the data gathered in PAPER. Two analysis were carried out. The first was the recreation of the author's original pipeline (EXPLAIn), followed by the outcome of implementing an alternative pipeline (EXPLAIN). Through the annalysis, the cells were clustered according to their molecular signatures. Both versions display the found clusters in the data as well as the top ten genes in those clusters and their tissue of origin.", width = 12)
      )
              ),
      tabItem(tabName = "originalpipeline",
    fluidRow(
      column(h3("Results for the original pipeline"), width = 8,
             box(h3("Clusters of genes for original pipelines"), h5("Hover over the clusters in the graph to see their cluster ID and cell type. On the pannel on the right, there is a list of the different cell types which can be clicked on to hide or show them on the graph. If a dot on the graph is selected, their cell types will appear at the end of the list of the pannel on the right. "), plotlyOutput("mclustplot"), width = NULL)
      ),
      column(width = 4,
             box("Select a cluster in the graph to see the top 10 enriched genes.",uiOutput("origclust"), width = NULL),
             box(textInput("origquerygene",label = "Type the name of a gene to see in which clusters it appears."), width = NULL),
             box(textOutput("origgeneoutput"), width = NULL)
      )
    )
  ),
      tabItem(tabName = "alternatepipeline",
        fluidRow(
          column(h3("Results for the alternative pipeline"), width = 8,
              box(h3("Clusters of genes for alternative pipeline"), h5("Hover over the clusters in the graph to see their cluster ID and cell type. On the pannel on the right, there is a list of the different cell types which can be clicked on to hide or show them on the graph. If a dot on the graph is selected, their cell types will appear at the end of the list of the pannel on the right. "), plotlyOutput("alternateplot"), width = NULL)
          ),
          column(width = 4,
              box("Select a cluster in the graph to see the top 10 enriched genes.", uiOutput("altclust"), width = NULL),
              box(textInput("querygene",label = "Type the name of a gene to see in which clusters it appears"), width = NULL),
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
  output$origclust <- renderUI({
    origclust[val(),]
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
  
  origquerytest <- reactive({
    req(input$origquerygene)
    grepstuff = grep(input$origquerygene, origclust$Genes)
    print(grepstuff)
    if (identical(grepstuff, integer(0))){
      print("No genes found")
      
    }else{
      results = paste(grep(input$origquerygene, origclust$Genes), collapse = " ")
      paste("Genes found in clustert(s):",results, collapse = " ")
    }
  })
  output$geneoutput <- renderText({
    querytest()
  })
  output$origgeneoutput <- renderText({
    origquerytest()
  })
} 

shinyApp(ui, server)