# Analyze_SingleCellData_as_bulk
When we are trying to get our interested gene's expression level from online dataset, my favorite case is when the authors attach a normalized gene count table, so we can run a simple RNA-seq pipeline to get the foldchange of our interested gene (or differential expressed gene list) in each condition. But there's a high chance we running into single cell RNA-seq, while most of the tutorials of single cell RNA-seq analysis are talking about how to get marker genes from each cell cluster. It will take some time to figure out how to get the similar information as bulk RNA-seq: a normalized gene count table in different conditions of one cell type. 
So I'm going to share my pipeline of how to extract the information of gene expression level in one cell type.

## Data set used
The dataset I'm using is [GSE164585](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164585). There are 6 set of samples in their experiment, 2 replicates in 4-week, 26-week and  86-week mice. The publication discribing their experiment is [Xie et al.](https://pubmed.ncbi.nlm.nih.gov/34879708/).

## Outline
When I first trying to analyze this dataset following the standard workflow, after the step of dimension reduction, I got 6 clusters of cell, which represents the 6 set of samples we downloaded instead of cell type. To overcome this problem, we need to integrate the 6 samples together, then conduct another dimension reduction. Another problem I encountered was it's quite confusing of how to extract iformation from a Seurat object, I also included this part in my outline.

The steps we are going to go over in this pipeline are:
1. quality control
2. data normalization
3. integration
4. dimensional reduction
5. clustering 
6. information extraction


## Code and intermediate files
I put my code in 3 different R files under the 'code' folder, and saved those important intermediate objects as RDS files.
Feel free to download and run the code!
