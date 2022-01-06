library(dplyr)
library(Seurat)
library(patchwork)

SAN.data <- Read10X(data.dir = "./SAN_new/")

SAN <- CreateSeuratObject(counts = SAN.data, project = "SAN", min.cells = 3, min.features = 200)

SAN

# Lets examine a few genes in the first thirty cells
SAN.data[c("SHOX2", "TBX3", "TBX5"), 1:30]

dense.size <- object.size(as.matrix(SAN.data))

dense.size

sparse.size <- object.size(SAN.data)

sparse.size

dense.size/sparse.size



# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
SAN[["percent.mt"]] <- PercentageFeatureSet(SAN, pattern = "^MT-")


# Show QC metrics for the first 5 cells
head(SAN@meta.data, 5)



# Visualize QC metrics as a violin plot
VlnPlot(SAN, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(SAN, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot1
plot2 <- FeatureScatter(SAN, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2


SAN <- subset(SAN, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 5)



SAN <- NormalizeData(SAN, normalization.method = "LogNormalize", scale.factor = 10000)

SAN <- FindVariableFeatures(SAN, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(SAN), 10)


plot3 <- VariableFeaturePlot(SAN)
plot3
plot4 <- LabelPoints(plot = plot3, points = top10, repel = TRUE)
plot4

#Scaling the data
all.genes <- rownames(SAN)
SAN <- ScaleData(SAN, features = all.genes)

#'regress out' heterogeneity associated with mitochondrial contamination
SAN <- ScaleData(SAN, vars.to.regress = "percent.mt")


#Perform linear dimensional reduction
SAN <- RunPCA(SAN, features = VariableFeatures(object = SAN))


# Examine and visualize PCA results a few different ways
print(SAN[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(SAN, dims = 1:2, reduction = "pca")

DimPlot(SAN, reduction = "pca")

DimHeatmap(SAN, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(SAN, dims = 1:15, cells = 500, balanced = TRUE)


# Determine the 'dimensionality' of the dataset
SAN <- JackStraw(SAN, num.replicate = 100)
SAN <- ScoreJackStraw(SAN, dims = 1:20)

JackStrawPlot(SAN, dims = 1:20,ymax = 0.8)

ElbowPlot(SAN)


#Cluster the cells
SAN <- FindNeighbors(SAN, dims = 1:15)
SAN <- FindClusters(SAN, resolution = 0.1)
# resolution can adjust from 0.2 to 1.5

head(Idents(SAN), 5)



#Run non-linear dimensional reduction (UMAP/tSNE)
SAN <- RunUMAP(SAN, dims = 1:15)

DimPlot(SAN, reduction = "umap")


SAN <- RunTSNE(SAN, dims = 1:15)

DimPlot(SAN, reduction = "tsne")


UMAPPlot(SAN,label=TRUE)

TSNEPlot(SAN,label=TRUE)

# find all markers of cluster 0
cluster0.markers <- FindMarkers(SAN, ident.1 = 0, min.pct = 0.25)
head(cluster0.markers, n = 20)

# find markers distinguishing cluster 1 from cluster 0
cluster1_vs_cluster0<-FindMarkers(SAN,ident.1 = 1,ident.2 = 0,min.pct = 0.25)
head(cluster1_vs_cluster0,n=10)
write.csv(cluster1_vs_cluster0,file = "cluster1_versus_cluster2_markers.csv")

# find markers for every cluster compared to all remaining cells, report only the positive ones
SAN.markers <- FindAllMarkers(SAN, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
SAN.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)


top10 <- SAN.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top15 <- SAN.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_logFC)
DoHeatmap(SAN, features = top10$gene) + NoLegend()

#find each cluster characters
VlnPlot(SAN, features = c("SHOX2"))
FeaturePlot(SAN,features = c("SHOX2"))

VlnPlot(SAN, features = c("EGFP"))
FeaturePlot(SAN,features = c("EGFP"))

VlnPlot(SAN, features = c("MYH6"))
FeaturePlot(SAN,features = c("MYH6"))

VlnPlot(SAN, features = c("mcherry"))
FeaturePlot(SAN,features = c("mcherry"))


saveRDS(SAN, file = "./SAN_new.rds")



