# Load required packages
library(Seurat)
library(dplyr)
library(plotly)
library(ggplot2)

# Read in gene expression data
raw_counts <- read.csv("data/transformed_data.csv", row.names = "Gene_Id")

# Convert expression data to integers
raw_counts[, 1:466] <- sapply(raw_counts[, 1:466], as.integer)

# Create Seurat object
mydata <- CreateSeuratObject(raw_counts)

# adding metadata to seurat object
metadata <- read.csv("data/SraRunTable.csv", row.names = 1)
mydata <- AddMetaData(mydata, metadata)

# Please ignore this
# common_cells <- intersect(colnames(mydata), rownames(metadata))
# metadata <- metadata[common_cells,]
# mydata <- mydata[,common_cells]
# mydata <- AddMetaData(mydata, metadata, subset = colnames(mydata))

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
plot1
plot2

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

# You may get error like no umap package is installed for the below line if it shows error then run next line
# reticulate::py_install("umap-learn")

# Perform UMAP dimensionality reduction and visualize the results
mydata <- RunUMAP(mydata, dims = 1:10, n.components = 3L)
umap_df <- as.data.frame(mydata@reductions$umap@cell.embeddings)
umap_df$clusters <- Idents(mydata)

set.seed(0)

########### unbiased approach coloured based on umap clustering ############
# normal umap #
fig1 <- plot_ly(data = umap_df, x = ~UMAP_1, y = ~UMAP_2, type = 'scatter', color = ~clusters, mode = 'markers')
###3d plot#####
fig2 <- plot_ly(data = umap_df, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, color = ~clusters, mode = 'markers') %>% add_markers(size = 8)

fig1
fig2
############################################

############ biased approach coloured based on given metadata celltype ###########
features <-subset(metadata, select = c(Cell_type))
pbd <- cbind(umap_df, features)

fig3 <- plot_ly(data = pbd, x = ~UMAP_1, y = ~UMAP_2, type = 'scatter', color = ~Cell_type, mode = 'markers')
fig4 <- plot_ly(data = pbd, x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, color = ~Cell_type, mode = 'markers') %>% add_markers(size = 8)

fig3
fig4

############################################

# Seurat Dimention plot
DimPlot(mydata, reduction = "umap", label = TRUE)

##### tried using plotly ######
plot1 <- ggplotly(
  DimPlot(
    mydata, 
    reduction = "umap", 
    label = TRUE
  ) + 
    ggtitle("UMAP Plot") +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10)
    )
)
plot1


# Same thing if want to try tsne
# mydata <- RunTSNE(mydata)
# DimPlot(mydata, reduction = "tsne", label = TRUE)

# Find differentially expressed genes between clusters and generate heatmap
mydata.markers <- FindAllMarkers(mydata, only.pos = TRUE, min.pct = 0.25)

# creating a marker list of top 20 genes with data
marker_list <- mydata.markers %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC)

# samething but for top10
mydata.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

# just markergene names of each cluster
cluster_genes <- marker_list %>% 
  group_by(cluster) %>% 
  summarize(genes = paste(gene, collapse = ", "))

# View the list of genes for each cluster
cluster_genes

# heatmap for top20 and top10 genes in all clusters(top20 looks a bit confusing)
DoHeatmap(mydata, features = marker_list$gene)
DoHeatmap(mydata, features = top10$gene)

# Plot the single marker in UMAP
FeaturePlot(object = mydata, features = "GAD1", pt.size = 0.5, combine = TRUE)

###### plotly version ######
ggplotly(
  FeaturePlot(mydata, features = "GAD1", pt.size = 0.5) +
  ggtitle("GAD1 expression in UMAP") +
  theme(plot.title = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10))
)

# creating a violin plot for single marker
VlnPlot(mydata, features = "GAD1")

####### plotly version #########
ggplotly(
  VlnPlot(mydata, features = "GAD1") + 
  ggtitle("GAD1 expression")
)

# same thing but a bit different(it's a feature plot in seurat can convert into ggplot but not into plotly)
RidgePlot(mydata, features = "GAD1")


###### just ignore this #######

# Define a list of known marker genes for different cell types
marker_genes <- list(
  Astrocytes = c("ABCA1", "ACSL6", "ADCYAP1R1", "AGT", "ALDOC", "ANLN", "APOE", "ARHGEF26", "ATF3", "ATP13A4", "ATP1B2", "BBOX1", "BTG2", "C1orf61", "CA12", "CASQ1", "CBS", "CDC42EP4", "CHST9", "CLU", "CPE", "CPNE5", "CPVL", "CRYAB", "CST3", "DAND5", "DGKG", "DKK3", "DNAJB1", "DNASE2", "DOK5", "DTNA", "EDNRB", "EEPD1", "EFEMP1", "EGLN3", "EPAS1", "EZR", "F3", "FOS", "FOSB", "FTH1", "FXYD7", "GABRB1", "GADD45B", "GATSL3", "GLUL", "GPR75", "GRAMD3", "GRIA1", "HEPACAM", "HEPN1", "HEY1", "HIF1A", "HIF3A", "HLA-E", "HNMT", "HSPB8", "ID1", "ID4", "IER2", "IGFBP7", "IL11RA", "IL33", "JUN", "JUNB", "JUND", "KCNH7", "KCNIP2", "KIF21A", "L1CAM", "LCAT", "LHFP", "LIX1", "LPL", "LRIG1", "LRRC8A", "LYPD1", "MAFB", "METTL7A", "MLC1", "MMD2", "MT1X", "MT2A", "MTHFD2", "NDRG2", "NFKBIA", "NHSL1", "NOG", "NRP1", "NTRK2", "P2RY1", "PAPLN", "PEA15", "PER1", "PFKFB3", "PLTP", "PON2", "PRLHR", "PRRT2", "PSAP", "RASL10A", "RASSF4", "RFX4", "RGMA", "RHOB", "RND3", "SCG2", "SCG3", "SEMA6A", "SERPINA3", "SF3A1", "SLC1A2", "SLC1A3", "SLC39A11", "SORL1", "SOX9", "SPARCL1", "SPOCK1", "SPON1", "SRPX", "SSTR2", "ST6GAL2", "TACR1", "TBC1D10A", "TIMP3", "TMEM151B", "TNS1", "TOB2", "TPCN1", "TRIL", "TSC22D4", "TSPAN12", "ZFAND5", "ZFP36", "ZFP36L1", "ZFP36L2"),
  Endothelial = c("A2M", "APOLD1", "FLT1", "TM4SF1"),
  Fetal_quiescent = c("STMN2", "NEFM", "SOX11"),
  Fetal_replicating = c("TOP2A", "MKI67", "RFC2"),
  Hybrid = c("DCX", "PAX6", "NEUROD1"),
  Microglia = c("CD68", "CD163", "P2RY12"),
  Neurons = c("NEFL", "SYP", "MAP2"),
  Oligodendrocytes = c("MOBP", "CNP", "QKI"),
  OPC = c("PDGFRA", "SOX10", "OLIG2")
)
