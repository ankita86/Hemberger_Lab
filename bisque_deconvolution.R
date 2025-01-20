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

#Annotate and extract relevant cell types based on the expression profile of classical markers 
#Convert seurat object to SingleCellExperiment object -

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



## calling function - deconvo_data - This function is written based on the analysis script of a tool - DeconvoBuddies. This tool allows implementation of Mean Ratio method to select cell-type specific marker genes which is used as an input for further downstream analysis by bisque.

#Link : https://research.libd.org/DeconvoBuddies/articles/Deconvolution_Benchmark_DLPFC.html

results<-deconvo_data(sce=y.sce)



## plotting
plot_data<-results$plot_data

#V1
plot_data$cluster_name <- factor(plot_data$cluster_name , levels = c("LE", "GE", "Endothelial", "Stromal", "B cell", "T cell", "DC", "NKC", "Neutrophil", "Macrophage"))

nb.cols <- 10
mycolors <- colorRampPalette(brewer.pal(8, "RdBu"))(nb.cols)


 
plot_composition_bar(
     plot_data,
     sample_col = "Sample_ID",
     x_col = "Sample_ID",
     ct_col ="cluster_name",
     add_text = FALSE
 ) + scale_fill_manual(values = mycolors)
 


########## function - deconvo_data -##########################################################

## this function can be customised for the number of markers used 
deconvo_data <- function (sce) {
 
 marker_stats <-get_mean_ratio(sce,
    cellType_col = "cluster_name")
   
  
bulk <- read.table("Atp11a_raw_counts.txt", sep="\t", header=T)
table(duplicated(bulk$Gene))
rownames(bulk) <- make.names(bulk$Gene, unique = TRUE)
bulk <-bulk[,-1]


marker_genes <- marker_stats |>
    filter(MeanRatio.rank <= 35 & gene %in% rownames(bulk))
	
top_35_markers <- marker_stats |> filter(MeanRatio.rank <= 35)
	
#marker_genes |> count(cellType.target)
	
marker_genes <- marker_genes |> pull(gene)


 samples <-colnames(bulk)
 genes <-rownames(bulk)
 
 
se <- SummarizedExperiment(assays = list(counts = as.matrix(bulk)),
                           colData = samples,
                           rowData = genes)
 
 colnames(colData(se))[colnames(colData(se)) == "X"] <- "Sample_ID"
 
 

 
 ### make objects for bisque
 ## for bulk
 
 exp_set_bulk <- Biobase::ExpressionSet(
    assayData = assays(se[marker_genes, ])$counts,
    phenoData = AnnotatedDataFrame(
        as.data.frame(colData(se))[c("Sample_ID")]
    )
)
 

 ## for sce
 exp_set_sce <- Biobase::ExpressionSet(
     assayData = as.matrix(assays(sce[marker_genes, ])$counts),
     phenoData = AnnotatedDataFrame(
         as.data.frame(colData(sce))[, c("cluster_name", "orig.ident")]
     ) )
	 
zero_cell_filter <- colSums(exprs(exp_set_sce)) != 0
message("Exclude ", sum(!zero_cell_filter), " cells")
#> Exclude 0 cells

exp_set_sce <- exp_set_sce[, zero_cell_filter] 
	 
	 
est_prop <- ReferenceBasedDecomposition(
    bulk.eset = exp_set_bulk,
    sc.eset = exp_set_sce,
    cell.types = "cluster_name",
    subject.names = "orig.ident",
    use.overlap = FALSE)


est_prop$bulk.props <- t(est_prop$bulk.props)


write.table(est_prop$bulk.props, file="Estimated_proportions_top35markers.txt", sep="\t", quote=FALSE)

### plotting
pd <- colData(se) |>
    as.data.frame() |>
    select(Sample_ID)
	
	prop_long <- est_prop$bulk.props |>
    as.data.frame() |>
    tibble::rownames_to_column("Sample_ID") |>
    tidyr::pivot_longer(!Sample_ID, names_to = "cluster_name", values_to = "prop") |>
    left_join(pd)
	
prop_long$Sample_ID <- factor(prop_long$Sample_ID , levels = c("WT5", "WT7", "WT8", "WT9", "KO5", "KO6", "KO8", "KO9"))
 return(list(all=est_prop,prop=est_prop$bulk.props, markers=marker_stats, plot_data=prop_long, top_35_Markers_deconvo=top_35_markers))

}
 

