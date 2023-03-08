library(ggplot2)
library(Seurat)
library(Matrix)
library(stringr)
library(readr)
library(here)
library(fitdistrplus)
library(dplyr)
library(data.table)
options(future.globals.maxSize = 10000 * 1024^2)

##### restart R to run this command first
library(reticulate)
path_to_python <- "/Users/naveenjain/opt/anaconda3/"
use_python(path_to_python)
reticulate::use_condaenv(condaenv = "py36", required = TRUE)
scanorama = reticulate::import(module = "scanorama")
##### restart R to run this command first

homeDirectory <-("/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/extractedData/fateMap10X/FM3/")

# Load the booster reprogramming dataset
s1_data <- Read10X(data.dir = "/Users/naveenjain/Downloads/FM3_boosterReprogramming/10XCellRangerOuts/1_DMSO_A/filtered_feature_bc_matrix/")
s2_data <- Read10X(data.dir = "/Users/naveenjain/Downloads/FM3_boosterReprogramming/10XCellRangerOuts/2_DMSO_B/filtered_feature_bc_matrix/")
s3_data <- Read10X(data.dir = "/Users/naveenjain/Downloads/FM3_boosterReprogramming/10XCellRangerOuts/3_LSD1i_A/filtered_feature_bc_matrix/")
s4_data <- Read10X(data.dir = "/Users/naveenjain/Downloads/FM3_boosterReprogramming/10XCellRangerOuts/4_LSD1i_B/filtered_feature_bc_matrix/")
s5_data <- Read10X(data.dir = "/Users/naveenjain/Downloads/FM3_boosterReprogramming/10XCellRangerOuts/5_DOT1Li_A/filtered_feature_bc_matrix/")
s6_data <- Read10X(data.dir = "/Users/naveenjain/Downloads/FM3_boosterReprogramming/10XCellRangerOuts/6_DOT1Li_B/filtered_feature_bc_matrix/")

s1 <- CreateSeuratObject(counts = s1_data, project = "DMSO_A", min.cells = 3, min.features = 200)
s2 <- CreateSeuratObject(counts = s2_data, project = "DMSO_B", min.cells = 3, min.features = 200)
s3 <- CreateSeuratObject(counts = s3_data, project = "LSD1i_A", min.cells = 3, min.features = 200)
s4 <- CreateSeuratObject(counts = s4_data, project = "LSD1i_B", min.cells = 3, min.features = 200)
s5 <- CreateSeuratObject(counts = s5_data, project = "DOT1Li_A", min.cells = 3, min.features = 200)
s6 <- CreateSeuratObject(counts = s6_data, project = "DOT1Li_B", min.cells = 3, min.features = 200)

s1[["percent.mt"]] <- PercentageFeatureSet(object = s1, pattern = "^MT-")
s2[["percent.mt"]] <- PercentageFeatureSet(object = s2, pattern = "^MT-")
s3[["percent.mt"]] <- PercentageFeatureSet(object = s3, pattern = "^MT-")
s4[["percent.mt"]] <- PercentageFeatureSet(object = s4, pattern = "^MT-")
s5[["percent.mt"]] <- PercentageFeatureSet(object = s5, pattern = "^MT-")
s6[["percent.mt"]] <- PercentageFeatureSet(object = s6, pattern = "^MT-")

VlnPlot(s1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(s2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(s3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(s4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(s5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(s6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

s1 <- subset(x = s1, subset = nFeature_RNA > 200 & percent.mt < 30)
s2 <- subset(x = s2, subset = nFeature_RNA > 200 & percent.mt < 30)
s3 <- subset(x = s3, subset = nFeature_RNA > 200 & percent.mt < 30)
s4 <- subset(x = s4, subset = nFeature_RNA > 200 & percent.mt < 30)
s5 <- subset(x = s5, subset = nFeature_RNA > 200 & percent.mt < 30)
s6 <- subset(x = s6, subset = nFeature_RNA > 200 & percent.mt < 30)

scTransform <- merge(s1, y = c(s2, s3, s4, s5, s6), add.cell.ids = c("S1", "S2", "S3", "S4", "S5", "S6"))

extractDataScanorama <- function(seurat.object, assay = "RNA", slot = "counts", groupingVar = "orig.ident", group_name){
  return(t(as.matrix(GetAssayData(seurat.object, assay = assay, slot = slot)))[colnames(seurat.object)[seurat.object@meta.data[,groupingVar] == group_name],])
}

samples = c("DMSO_A","DMSO_B", "LSD1i_A", "LSD1i_B", "DOT1Li_A", "DOT1Li_B")

datasets = list()
gene_list = list()

for (i in 1:length(samples)){
  datasets[[i]] <- extractDataScanorama(scTransform, group_name = samples[[i]])
  gene_list[[i]] <- rownames(scTransform)
}
rm(s1)
rm(s1_data)
rm(s2)
rm(s2_data)
rm(s3)
rm(s3_data)
rm(s4)
rm(s4_data)
rm(s5)
rm(s5_data)
rm(s6)
rm(s6_data)

integrated_corrected_data = scanorama$correct(datasets, gene_list, return_dimred = TRUE, return_dense = TRUE, ds_names = samples, verbose = TRUE)

corrected_scanorama <- t(do.call(rbind, integrated_corrected_data[[2]]))
colnames(corrected_scanorama) <- colnames(scTransform)
rownames(corrected_scanorama) <- integrated_corrected_data[[3]]
corrected_scanorama_pca <- t(do.call(rbind, integrated_corrected_data[[1]]))
colnames(corrected_scanorama_pca) <- colnames(scTransform)

# add in assay and format as  seurat object
scanorama_assay <- CreateAssayObject(data = corrected_scanorama)
scTransform[["scanorama"]] <- scanorama_assay
DefaultAssay(scTransform) <- "scanorama"

# preprocess scanorama values and perform PCA
scTransform <- FindVariableFeatures(scTransform, assay = "scanorama", selection.method = "vst", nfeatures = 7000)
all.genes <- rownames(scTransform)
scTransform <- ScaleData(scTransform, features = all.genes)
scTransform <- RunPCA(object = scTransform, assay = "scanorama", reduction.name = "pca_scanorama")
scTransform <- FindNeighbors(object=scTransform, dims = 1:50, reduction = "pca_scanorama", force.recalc = TRUE, graph.name = "scanorama_snn")
scTransform <- FindClusters(object=scTransform,graph.name = "scanorama_snn", resolution = 0.8)
scTransform <- RunUMAP(object = scTransform, reduction = "pca_scanorama", dims = 1:50, reduction.name = "umap_scanorama")
clusterUMAP = DimPlot(scTransform, reduction = "umap_scanorama", label = TRUE, group.by = "scanorama_snn_res.0.8")
sampleUMAP = DimPlot(scTransform, reduction = "umap_scanorama", group.by = "orig.ident")
sampleSplitUMAP = DimPlot(scTransform, reduction = "umap_scanorama", split.by = "orig.ident")

scTransform@assays$scanorama@scale.data

saveRDS(scTransform, file = paste0(homeDirectory, "scanorma_filter.rds"))
rm(scTransform)

#############################################################################
#DON'T NEED TO RUN BLOCK THIS AGAIN
#############################################################################
scanorama_filter <- readRDS(paste0(homeDirectory, "scanorama_filter.rds"))
logNormalizedCounts = scanorama_filter@assays$scanorama@scale.data
logNormalizedCountsRound = round(logNormalizedCounts, 4)

cells_count = sub("-1", "", colnames(logNormalizedCountsRound))
cells_count_cellID = sub("S\\d_", "", cells_count)
cells_count_sampleNum = gsub("[^S123456]", "", cells_count)
logNormalizedCountsRound = as_tibble(as.data.frame((t(as.matrix(logNormalizedCountsRound)))))
logNormalizedCountsRound = logNormalizedCountsRound %>% mutate(cellID = cells_count_cellID,
                                                               sampleNum = cells_count_sampleNum)

umapCoordinates = (scanorama_filter[['umap_scanorama']])@cell.embeddings
cells_UMAP = sub("-1", "", rownames(umapCoordinates))
cells_UMAP_cellID = sub("S\\d_", "", cells_UMAP)
cells_UMAP_sampleNum = gsub("[^S123456]", "", cells_UMAP)
umapCoordinates = as_tibble(umapCoordinates)
umapCoordinates = umapCoordinates %>% mutate(cellID = cells_UMAP_cellID,
                                             sampleNum = cells_UMAP_sampleNum)

umapClusters = (scanorama_filter[['scanorama_snn_res.0.3']])
cells_Clusters = sub("-1", "", rownames(umapClusters))
cells_Clusters_cellID = sub("S\\d_", "", cells_Clusters)
cells_Clusters_sampleNum = gsub("[^S123456]", "", cells_Clusters)
umapClusters = as_tibble(umapClusters)
umapClusters = umapClusters %>% mutate(cellID = cells_Clusters_cellID,
                                       sampleNum = cells_Clusters_sampleNum)

pcaCoordinates = (scanorama_filter[['pca_scanorama']])@cell.embeddings
cells_PCA = sub("-1", "", rownames(pcaCoordinates))
cells_PCA_cellID = sub("S\\d_", "", cells_PCA)
cells_PCA_sampleNum = gsub("[^S123456]", "", cells_PCA)
pcaCoordinates = as_tibble(pcaCoordinates)
pcaCoordinates = pcaCoordinates %>% mutate(cellID = cells_PCA_cellID,
                                             sampleNum = cells_PCA_sampleNum)

write.table(pcaCoordinates, file=paste0(homeDirectory, 'pcaCoordinates_Scanorama_50pcs_filter.tsv'), col.names = TRUE, sep='\t')
write.table(umapCoordinates, file=paste0(homeDirectory, 'umapCoordinates_Scanorama_50pcs_filter.tsv'), col.names = TRUE, sep='\t')
write.table(umapClusters, file=paste0(homeDirectory, 'umapClusters_Scanorama_50pcs_filter.tsv'), col.names = TRUE, sep='\t')
write.table(logNormalizedCountsRound, file=paste0(homeDirectory, 'logNormalizedCounts_Scanorama_50pcs_filterRound.tsv'), col.names = TRUE, sep='\t')

rm(scanorama_filter)
rm(umapCoordinates)
rm(umapClusters)
rm(logNormalizedCounts)
rm(logNormalizedCountsRound)

#############################################################################
#DON'T NEED TO RUN BLOCK THIS AGAIN
#############################################################################

library(cowplot)
plots <- FeaturePlot(scanorama_filter, features = c("NANOG", "PODXL", "ALPL", "TDGF1",
                                                    "SNAI2", "LUM", "COL1A1", "FOSL1",
                                                    "POU5F1", "SOX2", "KLF4", "MYC",
                                                    "TOP2A", "MKI67", "PCNA", "TP53"), reduction = "umap_scanorama", min.cutoff = 'q5', slot = 'scale.data', combine = FALSE)

for(i in 1:length(plots)) {
  plots[[i]] <- plots[[i]] + NoLegend() + NoAxes()
}
plot_grid(plotlist = plots)

scanorama_filter <- FindClusters(object = scanorama_filter, resolution = 0)
scanorama_filter <- FindClusters(object = scanorama_filter, resolution = 0.1)
scanorama_filter <- FindClusters(object = scanorama_filter, resolution = 0.2)
scanorama_filter <- FindClusters(object = scanorama_filter, resolution = 0.3)
scanorama_filter <- FindClusters(object = scanorama_filter, resolution = 0.4)
scanorama_filter <- FindClusters(object = scanorama_filter, resolution = 0.5)
scanorama_filter <- FindClusters(object = scanorama_filter, resolution = 0.6)
scanorama_filter <- FindClusters(object = scanorama_filter, resolution = 0.7)
scanorama_filter <- FindClusters(object = scanorama_filter, resolution = 0.8)
scanorama_filter <- FindClusters(object = scanorama_filter, resolution = 0.9)
scanorama_filter <- FindClusters(object = scanorama_filter, resolution = 1)

library(clustree)
clustree(scanorama_filter, prefix = "scanorama_snn_res.")

DimPlot(scanorama_filter, group.by = "scanorama_snn_res.0.3", reduction = "umap_scanorama")
Idents(scanorama_filter) <- scanorama_filter$scanorama_snn_res.0.3
scanorama_filter[["percent.ribo"]] <- PercentageFeatureSet(object = scanorama_filter, pattern = "^RP[SL]", assay = "RNA")
FeaturePlot(scanorama_filter, label = TRUE, features = c("percent.ribo"), min.cutoff = "q5", max.cutoff = "q95")
VlnPlot(object = scanorama_filter, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 3)

### for making subsets used in different analyses ###
scanorama_filter <- readRDS(paste0(homeDirectory, "scanorama_filter.rds"))
DimPlot(scanorama_filter)

scanorama_filter <- FindClusters(object = scanorama_filter, resolution = 0)
scanorama_filter <- FindClusters(object = scanorama_filter, resolution = 0.1)
scanorama_filter <- FindClusters(object = scanorama_filter, resolution = 0.2)
scanorama_filter <- FindClusters(object = scanorama_filter, resolution = 0.3)
scanorama_filter <- FindClusters(object = scanorama_filter, resolution = 0.4)
scanorama_filter <- FindClusters(object = scanorama_filter, resolution = 0.5)
scanorama_filter <- FindClusters(object = scanorama_filter, resolution = 0.6)
scanorama_filter <- FindClusters(object = scanorama_filter, resolution = 0.7)
scanorama_filter <- FindClusters(object = scanorama_filter, resolution = 0.8)
scanorama_filter <- FindClusters(object = scanorama_filter, resolution = 0.9)
scanorama_filter <- FindClusters(object = scanorama_filter, resolution = 1)

library(clustree)
clustree(scanorama_filter, prefix = "scanorama_snn_res.")
ggsave(filename = paste0(plotDirectory, "FM3_clustreePlot.pdf"), unit = "in", height = 5, width = 5)

scanorama_filter_subset <- subset(scanorama_filter, subset = seurat_clusters %in% c(0, 1, 2, 3, 4, 5))
scanorama_filter_subset <- AddMetaData(object = scanorama_filter_subset, metadata = scanorama_filter_subset@reductions$umap_scanorama@cell.embeddings[,2], col.name = "UMAP_2")
scanorama_filter_subset <- subset(scanorama_filter_subset, subset = UMAP_2 < 10)
DimPlot(scanorama_filter_subset)

saveRDS(object = scanorama_filter_subset, file = paste0(homeDirectory, "scanorama_filter_subset.rds"))

scanorama_filter_iPSC <- subset(scanorama_filter_subset, subset = seurat_clusters %in% c(0, 1, 2, 3, 5))
saveRDS(object = scanorama_filter_iPSC, file = paste0(homeDirectory, "scanorama_filter_iPSC.rds"))