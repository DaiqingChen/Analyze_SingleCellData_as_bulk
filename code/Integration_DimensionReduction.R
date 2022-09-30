################################################################################
# the steps we are going to go over in this pipeline:
#(1) quality control
#(2) data normalization
#(3) integration
#(4) dimensional reduction
#(5) clustering 
#(6) information extraction

################## (3)&(4) Integration and Dimension reduction #################
# Let's take a look of how dimension reduction would look like without integration
# combine the 6 samples into one seurat object
atvb <- merge(week4_rep1,week4_rep2,week26_rep1,week26_rep2,week86_rep1,week86_rep2, add.cell.ids = c("S1", "S2", "S3", "S4", "S5", "S6"), 
                       project = "atvb")
# scale the data 
atvb <- ScaleData(atvb, features = all.genes)
# find variable features
atvb <- FindVariableFeatures(atvb, selection.method = "vst", nfeatures = 2000)
features <- SelectIntegrationFeatures(atvb)
# dimension reduction by PCA and UMAP
atvb <- FindVariableFeatures(object = atvb)
atvb <- RunPCA(atvb, features = VariableFeatures(object = atvb))
# take a look at
VizDimLoadings(atvb, dims = 1:2, reduction = "pca")
DimPlot(atvb, reduction = "pca")
atvb <- RunUMAP(atvb, dims = 1:10)
DimPlot(atvb, reduction = "umap")

# we can see that the cluster is seperated by sample batches.


############################ (3) Perform integration ###########################
# find features among the 6 samples before merging them 
raw.list <- lapply(X = c(week4_rep1,week4_rep2,week26_rep1,week26_rep2,week86_rep1,week86_rep2), FUN = function(x) {
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = raw.list)
# find anchors
atvb.anchors <- FindIntegrationAnchors(object.list = raw.list, anchor.features = features) #took several hours to run
# this command creates an 'integrated' data assay
atvb.combined <- IntegrateData(anchorset = atvb.anchors) 

# save the object as RDS file
saveRDS(raw.list, file = "raw.list.rds") 
saveRDS(features, file = "features.rds")
saveRDS(atvb.anchors, file = "atvb.anchors.rds")

# specify that we will perform downstream analysis on the corrected data 
# note that the original unmodified data still resides in the 'RNA' assay
DefaultAssay(atvb.combined) <- "integrated"


######################### (4) dimensional reduction ############################
# Run the standard workflow again
# to compare the result of dimension reduction
atvb.combined <- ScaleData(atvb.combined, verbose = FALSE) 
atvb.combined <- RunPCA(atvb.combined, npcs = 30, verbose = FALSE)
atvb.combined <- RunUMAP(atvb.combined, reduction = "pca", dims = 1:30)
atvb.combined <- RunTSNE(atvb.combined)

# Need to decide how many features to use for dimension reduction before plotting it
# Doing this by feature selection
ElbowPlot(atvb.combined)
# 10-15 PCs is fine
#The signal is captured in the first 10 PCs.
DimPlot(atvb.combined, reduction = "pca")
DimPlot(atvb.combined, reduction = "umap")

# save the object as RDS file
saveRDS(atvb.combined, file = "atvb.combined.rds")

##### end of this session 
##### use the next R code file to run the following steps


