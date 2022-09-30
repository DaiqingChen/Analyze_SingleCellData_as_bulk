################################################################################
# the steps we are going to go over in this pipeline:
#(1) quality control
#(2) data normalization
#(3) integration
#(4) dimensional reduction
#(5) clustering 
#(6) information extraction

################ (5)&(6) Clustering and Info extraction ########################

########################## (5) Cluster the cell ##############################
atvb.combined <- FindNeighbors(atvb.combined, reduction = "pca", compute.SNN = TRUE)

# There should be cell cluster number and biomarkers of each cluster of cells in the publication.
# for this paper, there are 16 cell clusters. 
# We can test out the parameter "resolution" here to find the best value that give same number of clusters after dimension reduction.
# trying to find a number that gives 16 clusters
# usually setting this parameter between 0.4-1.2 returns good results for single-cell datasets of around 3K cells. 
atvb.combined <- FindClusters(atvb.combined, resolution = 0.2)
DimPlot(atvb.combined, label = TRUE,reduction="tsne")

# then find the markers of each cluster
atvb.markers <- FindAllMarkers(atvb.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use="wilcox" )
# save markers as RDS file
saveRDS(atvb.markers, file = "atvb.markers.rds")

# find the corresponding cluster from the publication.
# assign the name according to their cell number
new.cluster.ids <- c("B cells","fibroblast","monocytes","EC1","CD4+ T cell","CD8+ Tcell","VSMC","EC2","neutrophils","plasmocytes","lymphoEC","mesothelial cell","pericytes","neuron","dendritic cell","mastocytes")
names(new.cluster.ids) <- levels(atvb.combined)
atvb.combined.cell <- RenameIdents(atvb.combined, new.cluster.ids)
DimPlot(atvb.combined.cell, reduction = "tsne", label = TRUE, pt.size = 0.5) 

# save combined object with cell name in each cluster as RDS file
saveRDS(atvb.combined.cell, file = "atvb.combined.cell.rds")

############################# (6) extract the information ######################
# in my project, we are interested in the Rock2 expression level in endothelial cells
# take a look at the expression level across all cell types
VlnPlot(atvb.combined.cell,features="Rock2") 

# extracting endothelial cells
EC1<-subset(atvb.combined.cell, idents = "EC1")
EC2<-subset(atvb.combined.cell, idents = "EC2")

# save as RDS file
saveRDS(EC1,"EC1.rds")
saveRDS(EC2,"EC2.rds")

#### Then we are going to extract the gene expression level from the object
# We can take a look in the meta-data
# nFeature_RNA is the number of genes detected in each cell. 
# nCount_RNA is the total number of molecules detected within a cell.
# We use nFeature_RNA, as our gene expression value.

# I want to compare the Rock2 expression level in EC1 at different age group: week4,26 and 86
# but there is no column having this information
# we need to set the column by ourself

# the cell type information is saved as "ident" in our object
# we need to write that into the metadata slot to connect the condition information and cell type information

# First, create a column in the meta.data slot to hold both the cell type and stimulation information 
age<-as.array(atvb.combined.cell@meta.data$orig.ident)
# clean up the name of each condition 
age<-strsplit(age, "_", fixed = TRUE)
for (row in (1:length(age))){
  age_new[row]<-age[[row]][[1]]
}
# Then, save the current ident to that column. 
atvb.combined.cell$age <- age_new
atvb.combined.cell$celltype.age <- paste(Idents(atvb.combined.cell), age_new, sep = "_")
atvb.combined.cell$celltype <- Idents(atvb.combined.cell)

# rename our Ident
Idents(atvb.combined.cell) <- "celltype.age"

# Finally, we can use FindMarkers to find the genes that are different between week4 and week26 in EC1 cells
# kind of like the differential expression gene list
EC1.4_26 <- FindMarkers(atvb.combined.cell, ident.1 = "EC1_week4", ident.2 = "EC1_week26", verbose = TRUE)
head(EC1.4_26, n = 15)

# write to file
# you can start analyze it just as bulk RNA-seq now!
EC1.4_26$gene <- rownames(EC1.4_26)
write.csv(EC1.4_26,"result/EC1.4_26.csv")

##### end of this session 
##### Coungrats !!


