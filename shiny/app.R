library(shiny)
library(shinydashboard)
library(plotly)
library(ggplot2)
library(dplyr)
library(mclust)
library(Rtsne)

#To whoever reads this
#I am aware that this code could be a lot better but
#1) I'm not super comfortable in R 
#2) I essentially developed the framework for this by myself alongside
##  the original pipeline
#I am sorry that it is a mess, but it works 

#Reading in necessary files
origclust <- read.csv("data/origtop.csv")
altclust <- read.csv("data/cluster_genes.csv")
geneinfo <- read.csv("data/geneinfo.csv")
metadata <- read.csv("data/SraRunTableMod.csv", row.names = 1)
#Subsetting the dataframe to extract cell types 
features <-subset(metadata, select = c(Cell_type))
distance_matrix <- read.csv("data/distance_matrix.csv", row.names = "X")
altdata <- read.csv("data/alternate.csv")


#replace all nans with 0
#Performing TSNE, this code is near identical to that in RNAseq.r
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

#UI consists of an shinydashboard page
#Up to line 44 sets up the tabs 
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
                       #Simple title pages
                       titlePanel("Group A submission"),
                       verbatimTextOutput("welcometext"),
                       titlePanel("Instructions for site use:")
                ),
                verbatimTextOutput("instructiontext")
      )
              ),
      tabItem(tabName = "originalpipeline",
    fluidRow(
      column(h3("Results for the original pipeline"), align = "center", width = 8,
             #Plot output, alongside checkbox input 
             box(h3("Clusters of genes for original pipelines"), plotlyOutput("mclustplot"), checkboxGroupInput("origcheckcluster", label = "Select clusters:", choices = c(1,2,3,4,5,6,7,8,9), inline = TRUE,
                                                                                                               selected = c(1,2,3,4,5,6,7,8,9)), width = NULL)
             
      ),
      column(width = 4,
             box("Select a cluster in the graph to see the top 10 enriched genes.",uiOutput("origclust"), width = NULL),
             #Where the cluster numbers are printed, along with the description
             box(textInput("origquerygene",label = "Type the name of a gene to see in which clusters it appears."), textOutput("origgeneoutput"), textOutput("origgeneinfo"), width = NULL, height = "100vh")
      )
    )
  ),
      tabItem(tabName = "alternatepipeline",
        fluidRow(
          column(h3("Results for the alternative pipeline"), align = "center", width = 8,
                 #Alternate plot output
              box(h3("Clusters of genes for alternative pipeline"), plotlyOutput("alternateplot"), checkboxGroupInput("altcheckcluster", label = "Select clusters:", choices = c(1,2,3,4,5,6,7,8,9), inline = TRUE,
                                                                                                                     selected = c(1,2,3,4,5,6,7,8,9)), width = NULL)
          ),
          column(width = 4,
              box("Select a cluster in the graph to see the top 10 enriched genes.", uiOutput("altclust"), width = NULL),
              #Alternate plot cluster and cell info 
              box(textInput("querygene",label = "Type the name of a gene to see in which clusters it appears."), textOutput("geneoutput"), textOutput("alternategeneinfo"), width = NULL, height = "100vh")
          )
          )
        )
)
)
)


server <- function(input, output) {
  #Writing this at midnight before the deadline
  #I am too scared to change this away from a testing variable in case everything breaks
  #This code block detects clicks from the alternate plot and stores what was clicked
  #I know this is a weird order but again, im too scared to move it
  val = reactiveVal(1)
  observe({
    testingclick <- event_data(event = "plotly_click", source = "testing")
    clickval <- testingclick$customdata
    val(clickval)
    
  })
  #THis renders the original pipeline's plot
  output$mclustplot <- renderPlotly({
    testing %>%
      #this checks the values of the checkboxes to know which clusters to include
      subset(mclust %in% input$origcheckcluster) %>%
      #I have no idea how this works but this somehow allows for cluster highlighting
      highlight_key(~mclust) %>%
      plot_ly(x = ~tsne1, y = ~tsne2, 
              color = ~Cell_type,
              colors = "Set3",
              type = "scatter",
              mode = "markers",
              source = "testing",
              #Following two lines are necessary for the customised hover text
              customdata = ~mclust,
              text = ~paste("Cluster ID:", mclust, 
                            "<br>Cell Type:", Cell_type),
              hoverinfo = "text")%>%
      #This also somehow makes the hovering work 
      highlight(on = "plotly_hover", off = "plotly_doubleclick") 
  })
  #This block is essentially identical to the code block above for the
  #Clicking of a plot, except for altdata this time 
  observe({
    click_data <- event_data(event = "plotly_click", source = "altdata")
    clickval <- as.integer(click_data$key) 
    val(clickval)

  })
  #This outputs the top 10 enriched genes for the alternate plot
  output$altclust <- renderUI({
    altclust[val(),]
  })
  #Same code block as above, but for the alternate plot
  output$alternateplot <- renderPlotly({
    altdata %>%
      subset(clusters %in% input$altcheckcluster) %>%
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
  #Outputs top 10 genes for original plot
  output$origclust <- renderUI({
    origclust[val(),]
  })
  #This block allows for the searching of genes for the cluster output
  querytest <- reactive({
    req(input$querygene)
    uppquery <- toupper(input$querygene)
    #simple grepping to search for the genes
    grepstuff = grep(uppquery, altclust$genes)
    if (identical(grepstuff, integer(0))){
      print("No genes found")
      
    }else{
      results = paste(grep(uppquery, altclust$genes), collapse = " ")
      paste("Genes found in cluster(s):",results, collapse = " ")
    }
  })
  #Does the same as the code block above
  origquerytest <- reactive({
    req(input$origquerygene)
    uppquery <- toupper(input$origquerygene)
    grepstuff = grep(uppquery, origclust$Genes)
    print(grepstuff)
    if (identical(grepstuff, integer(0))){
      print("No genes found")
      
    }else{
      results = paste(grep(uppquery, origclust$Genes), collapse = " ")
      paste("Genes found in cluster(s):",results, collapse = " ")
    }
  })
  #generates the output from above on the website
  output$geneoutput <- renderText({
    querytest()
  })
  output$origgeneoutput <- renderText({
    origquerytest()
  })
  #Simple instructions for the titlepage
  output$instructiontext <- renderText({
    paste("All plots on the site are fully interactive:",
          "- Clicking cell types within the plot legend removes them from the plot",
          "- Double clicking on a cell type within the plot legend isolates the cell type",
          "- Number of displayed clusters can be altered through the use of the checkboxes underneath the graph.",
          "- Hovering over a cluster highlights all points within the cluster - double cick to clear",
          "- Clicking on a cluster will show the top10 enriched genes for that cluster",
          "- Querying a gene is possible by entering it into the search box, this will result in the clusters it ",
          "  is present in as well as additional information", sep = "\n")
  })
  #Simple output for the welcome page
  output$welcometext <- renderText({
    paste("In 2015 a study was carried out by Darmanis et al. to attempt to build a brain cell atlas.
Cell samples were obtained from the brain, from which single cell transcriptomic analysis was performed.
The following is an interactive site that allows for the visualization of the reanalysis. 
Two analyses were carried out. The first was the recreation of the author's original pipeline, using SCDE,
Rtsne, Mclust and FactoMineR, followed by the outcome of implementing an alternative pipeline using Seurat.
Through the analysis, the cells were clustered according to their molecular signatures. Both versions display 
the found clusters in the data as well as the top ten genes in those clusters and their predicted cell type")
  })
  #This is similar to the code blocks above, except it searches for the gene and gets 
  #its description
  origreactivegeneinfo <- reactive({
    req(input$origquerygene)
    uppquery <- toupper(input$origquerygene)
    grepinfo <- grep(uppquery, geneinfo$Gene)[1]
    if (is.na(grepinfo)){
      print("No genes found")
      
    }else{
      genename <- geneinfo[grepinfo,2]
      genedesc <- geneinfo[grepinfo,3]
      paste(genename,":", genedesc, sep ="\n")
    }
  })
  #outputs gene info
  output$origgeneinfo <- renderText({
    origreactivegeneinfo()
  })
  #Gets gene description again
  alternatereactivegeneinfo <- reactive({
    req(input$querygene)
    uppquery <- toupper(input$querygene)
    grepinfo <- grep(uppquery, geneinfo$Gene)[1]
    
    if (is.na(grepinfo)){
      print("No genes found")
      
    }else{
      genename <- geneinfo[grepinfo,2]
      genedesc <- geneinfo[grepinfo,3]
      paste(genename,":", genedesc, sep ="\n")
    }
  })
  #outputs gene info 
  output$alternategeneinfo <- renderText({
    alternatereactivegeneinfo()
  })
} 

shinyApp(ui, server)