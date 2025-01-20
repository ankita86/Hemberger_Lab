## R script to perform deconvolution analysis -

library(Biobase)
library(BisqueRNA)
 
library(Seurat)
library(dplyr)
library(DeconvoBuddies)
library(SummarizedExperiment)
library(ggplot2)
library(RColorBrewer)

source("Deconvo/functions.R")

# read single cell data created using scRNA.R 
x <-readRDS("UT_Winkler_aging.rds")

#extract relevant cell types based on the expression profile of classical markers and convert seurat object to SingleCellExperiment object -
y <-subset(x, idents=c(6,8,0,1,2,5,7,9,3,4,10,11,12,13,14,19,21,24))

metadata <-y@meta.data

metadata$cluster_name <- NA

## Annotate epithelial, endothelial, immune and stromal clusters -

metadata$cluster_name[which(y@meta.data$seurat_clusters =="6" | y@meta.data$seurat_clusters =="8")] <- "Epithelial"

#metadata$cluster_name[which(y@meta.data$seurat_clusters =="6")] <- "LE"
#metadata$cluster_name[which(y@meta.data$seurat_clusters =="8")] <- "GE"

metadata$cluster_name[which(y@meta.data$seurat_clusters =="0" |y@meta.data$seurat_clusters =="1" |  y@meta.data$seurat_clusters =="2" | y@meta.data$seurat_clusters =="5" |y@meta.data$seurat_clusters =="7"| y@meta.data$seurat_clusters =="9" )] <- "Stromal"

### immune cells further sub-classified
metadata$cluster_name[which(y@meta.data$seurat_clusters =="12" |y@meta.data$seurat_clusters =="13")] <- "T cell"
metadata$cluster_name[which(y@meta.data$seurat_clusters =="14" |y@meta.data$seurat_clusters =="11")] <- "DC"
metadata$cluster_name[which(y@meta.data$seurat_clusters =="3")] <- "NKC"
metadata$cluster_name[which(y@meta.data$seurat_clusters =="4")] <- "Macrophage"
metadata$cluster_name[which(y@meta.data$seurat_clusters =="19")] <- "B cell"
metadata$cluster_name[which(y@meta.data$seurat_clusters =="10")] <- "Neutrophil"


metadata$cluster_name[which(y@meta.data$seurat_clusters =="21" |y@meta.data$seurat_clusters =="24" )] <- "Endothelial"
y@meta.data<-metadata
y.sce <- as.SingleCellExperiment(y)



## calling function - full_dataset - This function is written based on the analysis script of a tool - DeconvoBuddies. This tool allows implementation of Mean Ratio method to select cell-type specific marker genes which is used as an input for further downstream analysis by bisque.

#Link : https://research.libd.org/DeconvoBuddies/articles/Deconvolution_Benchmark_DLPFC.html

full_sketched<-full_dataset(sce=y.sce)
#subset_8_sketched<-subset_dataset(sce=y.sce)


## plotting
plot_data<-full_sketched$plot_data

#V1
plot_data$cluster_name <- factor(plot_data$cluster_name , levels = c("LE", "GE", "Endothelial", "Stromal", "B cell", "T cell", "DC", "NKC", "Neutrophil", "Macrophage"))

nb.cols <- 10
mycolors <- colorRampPalette(brewer.pal(8, "RdBu"))(nb.cols)

#pdf("Full_dataset.pdf")
 
plot_composition_bar(
     plot_data,
     sample_col = "Sample_ID",
     x_col = "Sample_ID",
     ct_col ="cluster_name",
     add_text = FALSE
 ) + scale_fill_manual(values = mycolors)
 

 
