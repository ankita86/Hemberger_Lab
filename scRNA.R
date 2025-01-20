## Script to create aging-based scRNA reference dataset of mouse uterus

## packages to install
library(Seurat)
library(BPCells)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
# set this option when analyzing large datasets
options(future.globals.maxSize = 3e+09)


##Study accession numbers: Arrayexpression accession numbers: E-MTAB-11491, E-MTAB-12889. Processing and QC of expression matrices for 27 samples was performed using code adapted from Winkler et al. For details, #please visit - https://zenodo.org/records/10259662
## Adding age as a metadata "Condition" to scRNA expression profiles using Table S1 from Winkler et al. 2024 and post-qced expression matrices were merged.
## Link: #https://www.sciencedirect.com/science/article/pii/S0092867424000588?via%3Dihub 


## read merged dataset and split it into layers based on the metadata column - Condition (which is indicative of age) to facilitate integration
object <- readRDS("Merged.rds") 
DefaultAssay(object=object) <- "RNA"
object$Condition <- factor(object$Condition , levels = c("Young", "Age_9M", "Age_12M", "Age_15M", "Age_18M"))
object <- NormalizeData(object)

#Split the seurat object 
object[["RNA"]] <- split(object[["RNA"]], f = object$Condition)
#object <-SCTransform(object, vst.flavor = "v2")
object <- FindVariableFeatures(object, verbose = T)


## sample cells from each dataset. 
object <- SketchData(object = object, ncells = 6500, method = "LeverageScore", sketched.assay = "sketch", assay = "RNA")
DefaultAssay(object) <- "sketch"
object <- FindVariableFeatures(object, verbose = F)
object <- ScaleData(object, verbose = F)
object <- RunPCA(object, verbose = F)

# Perform integration on the sketched cells across samples splitted based on Condition
object <- IntegrateLayers(object, method = HarmonyIntegration, orig = "pca", new.reduction = "harmony",
    dims = 1:30, k.anchor = 20, reference = which(Layers(object, search = "data") %in% c("data.Young")),
    verbose = F)
# cluster the integrated data
object <- FindNeighbors(object, reduction = "harmony", dims = 1:30)
object <- FindClusters(object, resolution = 0.7)

object <- RunUMAP(object, reduction = "harmony", dims = 1:30, return.model = T, verbose = F)

print (object)

## save sketched object

sc = JoinLayers(object)
sc_2 = DietSeurat(sc, assays="sketch", dimreducs = c("umap", "harmony", "pca"))

saveRDS(sc_2, file="UT_Winkler_aging.rds")
