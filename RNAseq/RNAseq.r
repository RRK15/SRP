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
#Following the tutorial for SCDE


dir <- "data/transformed_data.csv"
 #Needs to have it read the row names in or it doesnt fcking work
seq_data <- read.csv(dir, row.names = "Gene_Id")
#Convert all non gene columns to numeric
seq_data[,1:466] <- sapply(seq_data[,1:466], as.integer)

#Cleaning counts - no idea what this does tbh 
#get someone to read into these parameters
cleaned <- clean.counts(seq_data, min.lib.size=1000, min.reads = 1, min.detected = 1)

#This takes a fucking long time to run 
#do NOT run this with 8 cores - crashes laptop
#May be worth looking into saving the model when ran
#Also I have no idea what this does
#SOMETIMES THE MODEL JUST RANDOMLY FUCKS UP WHY
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

#Again, no idea what the fuck this is doing
#Also do not run this past 1 core

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
      x[!as.logical(rbinom(length(x),1,1-p.self.fail[,464]*k))] <- NA;
      x;
    }))
    rownames(scd1) <- rownames(cleaned); 
    # calculate correlation on the complete observation pairs
    cor(log10(scd1+1),use="pairwise.complete.obs");
  }, mc.cores = 1)
  #I think these are the distances we need for the next step
  direct.dist <- as.dist(1-Reduce("+",dl)/length(dl))
  distance_matrix<-as.data.frame(as.matrix(direct.dist))
  #replace all nans with 0
  distance_matrix[is.na(distance_matrix)] <- 0
  write.csv(distance_matrix,file = "data/distance_matrix.csv")
  distance_matrix <- as.matrix(direct.dist)
}


#WHO KNOWS WHAT THE FUCK THIS IS DOING
#If you don't have this as a matrix it gets EXTREMELy upset and gives you
#a nonsensical error that isnt actually the problem
normalised_distance = normalize_input(distance_matrix)
output_tsne <- Rtsne(normalised_distance, theta = 0.0)
needed_output <- as.data.frame(output_tsne$Y)

#Holy FUCK why did they have some options that were capitalised and some not
#It took me like 10 minutes to figur eout where i was going wrong
clustering <- Mclust(needed_output)

#Neat clusters thank you mclust
plot(clustering, what= "classification")
#Something about optimal clusters idk
plot(clustering, what ="BIC", xlab = "Number of Components")