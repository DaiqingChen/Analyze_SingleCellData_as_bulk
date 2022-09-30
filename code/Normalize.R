################################################################################
# the steps we are going to go over in this pipeline:
#(1) quality control
#(2) data normalization
#(3) integration
#(4) dimensional reduction
#(5) clustering 
#(6) information extraction

################ (1)&(2) Quality control and data normalization ################
### import package
library(dplyr)
library(Seurat)

### read in data
# There are 6 samples need to be downloaded from GSE164585, then up-zip and put into your prefered folder.
# I have already downloaded and upziped the folder
# I am going to read in by the sample name and create Seurat object for each sample
# S1-6 means the 6 samples.
#S1: week4  rep1
#S2: week4  rep2
#S3: week26 rep1
#S4: week26 rep2
#S5: week86 rep1
#S6: week86 rep2

## S1
week4_rep1.data <- Read10X(data.dir = "/Users/lasia/Desktop/LiaoLab/HawChih/Aging aorta paper/ATVB/GSM5014406")
# Initialize the Seurat object with the raw (non-normalized data).
week4_rep1 <- CreateSeuratObject(counts = week4_rep1.data, project = "week4_rep1", min.cells = 5, min.features = 200)
week4_rep1 # take a look at the object
# repeat for the rest
## S2
week4_rep2.data <- Read10X(data.dir = "/Users/lasia/Desktop/LiaoLab/HawChih/Aging aorta paper/ATVB/GSM5014407")
week4_rep2 <- CreateSeuratObject(counts = week4_rep2.data, project = "week4_rep2", min.cells = 5, min.features = 200)
## S3
week26_rep1.data <- Read10X(data.dir = "/Users/lasia/Desktop/LiaoLab/HawChih/Aging aorta paper/ATVB/GSM5014408")
week26_rep1 <- CreateSeuratObject(counts = week26_rep1.data, project = "week26_rep1", min.cells = 5, min.features = 200)
## S4
week26_rep2.data <- Read10X(data.dir = "/Users/lasia/Desktop/LiaoLab/HawChih/Aging aorta paper/ATVB/GSM5014409")
week26_rep2 <- CreateSeuratObject(counts = week26_rep2.data, project = "week26_rep2", min.cells = 5, min.features = 200)
## S5
week86_rep1.data <- Read10X(data.dir = "/Users/lasia/Desktop/LiaoLab/HawChih/Aging aorta paper/ATVB/GSM5014410")
week86_rep1 <- CreateSeuratObject(counts = week86_rep1.data, project = "week86_rep1", min.cells = 5, min.features = 200)
## S6
week86_rep2.data <- Read10X(data.dir = "/Users/lasia/Desktop/LiaoLab/HawChih/Aging aorta paper/ATVB/GSM5014411")
week86_rep2 <- CreateSeuratObject(counts = week86_rep2.data, project = "week86_rep2", min.cells = 5, min.features = 200)


########################### (1) quality control ################################
# check the mitochondria rate
week4_rep1[["percent.mt"]] <- PercentageFeatureSet(week4_rep1, pattern = "^MT-") 
# Visualize QC metrics as a violin plot
VlnPlot(week4_rep1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# we can see that there's no mt count
# the authors have already cleaned up this data-set
# no need to repeat for the rest of samples

########################### (2) Normalization ##################################
# testing the fuction by the first sample
week4_rep1 <- NormalizeData(week4_rep1, normalization.method = "LogNormalize", scale.factor = 10000)
week4_rep1 <- NormalizeData(week4_rep1) 
# take a look at the result
week4_rep1[["RNA"]]@data

# Put all 6 samples into a list, then run the normalization together
raw.list <- lapply(X = c(week4_rep1,week4_rep2,week26_rep1,week26_rep2,week86_rep1,week86_rep2), FUN = function(x) {
  x <- NormalizeData(x)
})

# save the raw.list as a RDS file 
saveRDS(raw.list, file = "raw.list.rds") 

##### end of this session 
##### use the next R code file to run the following steps



