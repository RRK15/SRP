# Load required libraries
library(Seurat)
library(dplyr)
library(patchwork)
library(Matrix)

# Read in the raw count matrix from a CSV file, with the first column as row names
raw_counts<-read.csv("data/orign_matrix_clean.csv", row.names = "Gene_Id")

# Convert the count data from character to integer
raw_counts[,1:466] <- sapply(raw_counts[,1:466], as.integer)

# Create a Seurat object with the count data
mydata <- CreateSeuratObject(raw_counts)

# Plot the relationship between the number of RNA molecules detected and the number of genes detected in each cell
plot1 <- FeatureScatter(mydata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1

# Normalize the count data to log base 10 and counts per million (CPM)
mydata <- NormalizeData(mydata, normalization.method = "LogNormalize", scale.factor = 1000000)

# Identify highly variable genes using the variance stabilizing transformation (VST) and select the top 2000 variable genes
mydata <- FindVariableFeatures(mydata, selection.method = "vst", nfeatures = 2000)

# Plot the top 20 variable genes and label them on the plot
top20 <- head(VariableFeatures(mydata), 20)
plot1 <- VariableFeaturePlot(mydata)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
plot1 + plot2

# Get the names of all genes
all.genes <- rownames(mydata)

# Scale the data so that each gene has mean 0 and standard deviation 1
mydata <- ScaleData(mydata, features = all.genes)

# Run principal component analysis (PCA) on the scaled data, using the top 2000 variable genes
mydata <- RunPCA(mydata, features = VariableFeatures(object = mydata))

# Visualize the loadings of the top two principal components
VizDimLoadings(mydata, dims = 1:2, reduction = "pca")

# Visualize the cells in two dimensions, based on the first two principal components
DimPlot(mydata, reduction = "pca")

# Visualize a heatmap of the first 15 principal components, for the first 500 cells
DimHeatmap(mydata, dims = 1:15, cells = 500, balanced = TRUE)
