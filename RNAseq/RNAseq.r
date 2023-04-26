#PUT THIS INTO R TERMINAL TO INSTALL SCDE

#require(devtools)
#devtools::install_version('flexmix', '2.3-13')
#devtools::install_github('hms-dbmi/scde', build_vignettes = FALSE)

#Install packages if you don't have them, might automate this t some point
#install.packages("mclust")
#install.packages("parallel")
#install.packages("Rtsne")

library(scde)
library(parallel)
library(Rtsne)
library(mclust)
library(FactoMineR)
library(factoextra)
library(plotly)
#Following the tutorial for SCDE

#Reads in data
metadata <- read.csv("data/SraRunTableMod.csv", row.names = 1)
features_cell <-subset(metadata, select = c(Cell_type))
dir <- "data/transformed_data.csv"

seq_data <- read.csv(dir, row.names = "Gene_Id")
test_data <- read.csv(dir)

#Convert all non gene columns to numeric
seq_data[,1:466] <- sapply(seq_data[,1:466], as.integer)

#Clean the data, at least 1000 reads per cell
cleaned <- clean.counts(seq_data, min.lib.size=1000, min.reads = 1, min.detected = 1)


#do NOT run this with 8 cores - crashes laptop
#Will check to see if the model.csv exists to skip the running in this block
#Much of this code is repurposed from the SCDE tutorial
if (file.exists("data/model.csv"))
  {
  orig_model <- read.csv("data/model.csv", row.names = "X")
} else{
  orig_model <- scde.error.models(counts = cleaned, n.cores = 2, 
                                  threshold.segmentation = TRUE, 
                                  save.crossfit.plots = FALSE, 
                                  save.model.plots = FALSE, 
                                  verbose = 1)
  
  valid.cells <- orig_model$corr.a > 0
  orig_model <- orig_model[valid.cells,]
  write.csv(as.data.frame(orig_model),"data/model.csv")
}

#Do NOT run this past one core ,my laptop nearly exploded 
#Avoid calculation if the distance matrix already exists
#If you want to run this block, delete the distance_matrix.csv in the data folder
if (file.exists("data/distance_matrix.csv"))
{
  distance_matrix <- read.csv("data/distance_matrix.csv", row.names = "X")
  #replace all nans with 0
  distance_matrix[is.na(distance_matrix)] <- 0
  distance_matrix <- as.matrix(distance_matrix)
} else 
{
  p.self.fail <- scde.failure.probability(models = orig_model, counts = cleaned)
  n.simulations <- 500; k <- 0.9;
  cell.names <- colnames(cleaned); names(cell.names) <- cell.names;
  dl <- mclapply(1:n.simulations,function(i) {
    scd1 <- do.call(cbind,lapply(cell.names,function(nam) {
      x <- cleaned[,nam]
      # replace predicted drop outs with NAs
      #Not entirely sure why, but I had to replace the value that was in self.fail
      #In the tutorial with the number of columns - 1
      #Fixes it, not sure why
      #No-one is helping me with this so I am limited with what I can look into 
      x[!as.logical(rbinom(length(x),1,1-p.self.fail[,464]*k))] <- NA;
      x;
    }))
    rownames(scd1) <- rownames(cleaned); 
    # calculate correlation on the complete observation pairs
    cor(log10(scd1+1),use="pairwise.complete.obs");
  }, mc.cores = 1)
  #Get distance matrix
  direct.dist <- as.dist(1-Reduce("+",dl)/length(dl))
  distance_matrix<-as.data.frame(as.matrix(direct.dist))
  #replace all nans with 0
  distance_matrix[is.na(distance_matrix)] <- 0
  write.csv(distance_matrix,file = "data/distance_matrix.csv")
  distance_matrix <- as.matrix(direct.dist)
}


#If you don't have this as a matrix it gets EXTREMELY upset and gives you
#a nonsensical error that isnt actually the problem
#Normalise the matrix and put it into TSNE
normalised_distance = normalize_input(distance_matrix)
output_tsne <- Rtsne(normalised_distance, theta = 0.0)
#Because the LAMP server doesn't take plotly, have to extract the data manually
needed_output <- as.data.frame(output_tsne$Y)
tsne1 <- output_tsne$Y[,1]
tsne2 <- output_tsne$Y[,2]
#Binds the extracted data together, second line binds it alongside the cell names
clust <- data.frame(tsne1, tsne2)
testing <- cbind(clust, features_cell)


#Running clustering using mclust
clustering <- Mclust(needed_output)
testing$mclust <- clustering$classification
#Have to use plotly because of rendering issues on LAMP
plot_ly(testing,
        x = ~tsne1, y = ~tsne2, 
        color = ~Cell_type,
        colors = "Set3",
        type = "scatter",
        mode = "markers",
        customdata = ~mclust,
        text = ~paste("Cluster ID:", mclust, 
                      "<br>Cell Type:", Cell_type),
        hoverinfo = "text")

#Mclust clusters
plot(clustering, what= "classification")
#Optimal number of components
plot(clustering, what ="BIC", xlab = "Number of Components")

#FactoMineR
#Not sure if we actually need this but its here
res.pca <- PCA(needed_output, ncp = 3, graph = TRUE)
fviz_pca_ind(res.pca, axes = c(1,2), geom.ind = "point", palette = c("green", "purple", "black"))