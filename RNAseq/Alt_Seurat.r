# Load required packages
library(Seurat)
library(dplyr)
library(patchwork)
library(Matrix)

# Read in gene expression data
raw_counts <- read.csv("data/transformed_data.csv", row.names = "Gene_Id")

# Convert expression data to integers
raw_counts[, 1:466] <- sapply(raw_counts[, 1:466], as.integer)

# Create Seurat object
mydata <- CreateSeuratObject(raw_counts)

# Visualize the distribution of number of RNA molecules and number of detected features
plot1 <- FeatureScatter(mydata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1

# Normalize the data using logarithm transformation (don't run the below line if using transformed_data.csv)
# mydata <- NormalizeData(mydata, normalization.method = "LogNormalize", scale.factor = 1000000)

# Find variable features using the variance-stabilizing transformation method
mydata <- FindVariableFeatures(mydata, selection.method = "vst", nfeatures = 2000)

# Plot the top 20 variable features
top20 <- head(VariableFeatures(mydata), 20)
plot1 <- VariableFeaturePlot(mydata)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
plot1 + plot2

# Scale the data and perform PCA
all.genes <- rownames(mydata)
mydata <- ScaleData(mydata, features = all.genes)
mydata <- RunPCA(mydata, features = VariableFeatures(object = mydata))

# Visualize the PCA results
VizDimLoadings(mydata, dims = 1:2, reduction = "pca")
DimPlot(mydata, reduction = "pca")
DimHeatmap(mydata, dims = 1:15, cells = 500, balanced = TRUE)

# Identify significant principal components using JackStraw
mydata <- JackStraw(mydata, num.replicate = 100)
mydata <- ScoreJackStraw(mydata, dims = 1:15)
JackStrawPlot(mydata, dims = 1:15)
ElbowPlot(mydata)

# Perform clustering and visualize the results
mydata <- FindNeighbors(mydata, dims = 1:10)
mydata <- FindClusters(mydata, resolution = 0.5)
head(Idents(mydata), 5)

# You may get error like no umap pakage is installed for the below line if it shows error the run next line
# reticulate::py_install("umap-learn")

# Perform UMAP dimensionality reduction and visualize the results
mydata <- RunUMAP(mydata, dims = 1:10)
DimPlot(mydata, reduction = "umap")

# Find differentially expressed genes between clusters and generate heatmap
mydata.markers <- FindAllMarkers(mydata, only.pos = TRUE, min.pct = 0.25)

mydata.markers %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC)

mydata.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

DoHeatmap(mydata, features = top10$gene)

