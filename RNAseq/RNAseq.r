library(scde)
library(parallel)
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
orig_model <- scde.error.models(counts = cleaned, n.cores = 4, 
                                threshold.segmentation = TRUE, 
                                save.crossfit.plots = FALSE, 
                                save.model.plots = FALSE, 
                                verbose = 1)

valid.cells <- orig_model$corr.a > 0
orig_model <- orig_model[valid.cells,]
#Again, no idea what the fuck this is doing
#Also do not run this past 1 core
p.self.fail <- scde.failure.probability(models = orig_model, counts = cleaned)
n.simulations <- 10; k <- 0.9;
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
mydf<-as.data.frame(as.matrix(direct.dist))
write.csv(mydf,file = "filename.csv")