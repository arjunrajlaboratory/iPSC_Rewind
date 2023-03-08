rm(list=ls())
gc()

library(tidyverse)
library(reshape2)
library(ggrepel)
library(ggforce)
library(ggsignif)
library(ggridges)
library(ggrastr)
library(RColorBrewer)
library(viridis)
library(egg)
library(Seurat)
library(biomaRt)
library(spgs)
library(matrixStats)

theme_set(theme_classic())

homeDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/extractedData/rewind10X/R3/"
plotDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/plots/rewind10X/R3/"

scanorama_filter <- readRDS(paste0(homeDirectory, "scanorma_filter.rds"))
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
scanorama_filter <- CellCycleScoring(scanorama_filter, s.features = s.genes, g2m.features = g2m.genes)

umapCoordinates = as_tibble(read.table(file = paste0(homeDirectory, "umapCoordinates_Scanorama_50pcs_filter.tsv"), header = TRUE, stringsAsFactors = F, sep = "\t"))
umapClusters = as_tibble(read.table(file = paste0(homeDirectory, "umapClusters_Scanorama_50pcs_filter.tsv"), header = TRUE, stringsAsFactors = F, sep = "\t"))
umapClusters <- umapClusters %>% dplyr::rename(cluster = scanorama_snn_res.0.3)

primedCellIDList <- readRDS(file = paste0(homeDirectory, "primedCellIDList.rds"))
cutoffList <- c(10, 25, 50, 100, 150, 200, 250, 500, 1000)
i = 6
primedAllUMAP <- filter(umapCoordinates, cellID %in% c(unlist(primedCellIDList[[i]]), unlist(primedCellIDList[[i+length(cutoffList)]]), unlist(primedCellIDList[[i+2*length(cutoffList)]]))) %>%
  dplyr::filter(sampleNum %in% c("S1", "S2", "S3"))

overlapTableList <- readRDS(paste0(homeDirectory, "overalapTableList.rds"))
linCountToOverlaps <- read.table(file = paste0(homeDirectory, "filtered10XCells.txt"), sep = "\t")

#### check if fast versus slow primed cells have different fates ####
#######################################################################################################################################################
overlapTableAll <- bind_rows(overlapTableList[[1]], overlapTableList[[2]], overlapTableList[[3]])
overlapTableAll <- overlapTableAll %>% rowwise() %>% mutate(nUMIMax = max(nUMINorm.x, nUMINorm.y)) %>% dplyr::select(BC50StarcodeD8, nUMIMax) %>% ungroup()
overlapTableAll <- overlapTableAll %>% group_by(BC50StarcodeD8) %>% filter(nUMIMax == max(nUMIMax))

primedFateTable <- inner_join(linCountToOverlaps %>% dplyr::select(1:2), overlapTableAll, by = c("barcode" = "BC50StarcodeD8"))
primedFateTable <- inner_join(primedFateTable, umapCoordinates, by = c("cellID"))
primedFateTable$sampleNum <- factor(primedFateTable$sampleNum, levels = c("S3", "S1", "S2"), labels = c("slow", "control", "fast"))

ggplot(primedFateTable %>% filter(nUMIMax > 15), aes(x = sampleNum, y = log(nUMIMax))) +
  geom_boxplot() +
  geom_jitter(seed = 1234) +
  geom_signif(comparisons = list(c("control", "fast"), c("control", "slow"))) +
  geom_signif(comparisons = list(c("fast", "slow")), position = position_nudge(x = 0, y = 0.5)) +
  ylab("log(number of cells in iPSC colonies)") +
  theme(axis.title.x = element_blank())
ggsave(filename = paste0(plotDirectory, "R3_fateDifference_largeColonies.pdf"), unit = "in", height = 5, width = 5)

ggplot(primedFateTable, aes(x = log(nUMIMax), y = sampleNum, height = stat(density))) +
  geom_density_ridges2(stat = "binline", bins = 40, scale = 0.95, draw_baseline = TRUE) +
  xlab("log(number of cells in iPSC colonies)") +
  geom_vline(xintercept = log(5), linetype = "dashed") +
  theme(axis.title.y = element_blank()) + NoLegend() + xlim(0, 11)
ggsave(filename = paste0(plotDirectory, "R3_fateDifference_allColonies.pdf"), unit = "in", height = 5, width = 3)

#### check PCA of primed cells ####
#######################################################################################################################################################
labelsToAdd <- umapCoordinates %>% mutate(label = "none")
labelsToAdd$label <- ifelse(labelsToAdd$cellID %in% linCountToOverlaps$cellID, "barcoded", labelsToAdd$label)
labelsToAdd$label <- ifelse(labelsToAdd$cellID %in% primedAllUMAP$cellID, "primed", labelsToAdd$label)
labelsToAdd <- labelsToAdd %>% dplyr::select(label)

scanorama_filter_labels <- AddMetaData(scanorama_filter, metadata = labelsToAdd)
scanorama_filter_labels$label <- labelsToAdd$label
Idents(scanorama_filter_labels) <- scanorama_filter_labels$label
DimPlot(scanorama_filter_labels)

scanorama_subset <- subset(scanorama_filter_labels, idents = c("primed"))
scTransform <- FindVariableFeatures(scanorama_subset, assay = "scanorama", selection.method = "vst", nfeatures = 7000)
all.genes <- rownames(scTransform)

scTransform <- ScaleData(scTransform, features = all.genes)
scTransform <- RunPCA(object = scTransform, assay = "scanorama", reduction.name = "pca_scanorama")
scTransform <- FindNeighbors(object=scTransform, dims = 1:50, reduction = "pca_scanorama", force.recalc = TRUE, graph.name = "scanorama_snn")
DimPlot(scTransform, reduction = "pca_scanorama", label = TRUE)

ElbowPlot(object = scTransform, reduction = "pca_scanorama", ndims = 20)

# pcLoadings <- as_tibble(scTransform[["pca_scanorama"]][, 1:2])
# pcLoadings$gene <- rownames(scTransform[["pca_scanorama"]][, 1:2])
# pcLoadingsMelt <- melt(pcLoadings, id.vars = "gene")
# mat <- Seurat::GetAssayData(scTransform, assay = "scanorama", slot = "scale.data")
# pca <- scTransform[["pca_scanorama"]]
# total_variance <- sum(matrixStats::rowVars(mat))
# eigValues = (pca@stdev)^2 
# varExplained = eigValues / total_variance

scTransform <- FindClusters(object = scTransform, resolution = 0.6)
scTransform <- RunUMAP(object = scTransform, reduction = "pca_scanorama", dims = 1:50, reduction.name = "umap_scanorama")
DimPlot(scTransform, reduction = "umap_scanorama", label = TRUE, group.by = "scanorama_snn_res.0.6")

umapCoordinates = (scTransform[['umap_scanorama']])@cell.embeddings
cells_UMAP = sub("-1", "", rownames(umapCoordinates))
cells_UMAP_cellID = sub("S\\d_", "", cells_UMAP)
cells_UMAP_sampleNum = gsub("[^S123456]", "", cells_UMAP)
umapCoordinates = as_tibble(umapCoordinates)
umapCoordinates = umapCoordinates %>% mutate(cellID = cells_UMAP_cellID,
                                             sampleNum = cells_UMAP_sampleNum)

umapClusters = (scTransform[['scanorama_snn_res.0.6']])
cells_Clusters = sub("-1", "", rownames(umapClusters))
cells_Clusters_cellID = sub("S\\d_", "", cells_Clusters)
cells_Clusters_sampleNum = gsub("[^S123456]", "", cells_Clusters)
umapClusters = as_tibble(umapClusters)
umapClusters = umapClusters %>% mutate(cellID = cells_Clusters_cellID,
                                       sampleNum = cells_Clusters_sampleNum)

pcaCoordinates = (scTransform[['pca_scanorama']])@cell.embeddings
cells_PCA = sub("-1", "", rownames(pcaCoordinates))
cells_PCA_cellID = sub("S\\d_", "", cells_PCA)
cells_PCA_sampleNum = gsub("[^S123456]", "", cells_PCA)
pcaCoordinates = as_tibble(pcaCoordinates)
pcaCoordinates = pcaCoordinates %>% mutate(cellID = cells_PCA_cellID, sampleNum = cells_PCA_sampleNum)
pcaCoordinates <- pcaCoordinates %>% dplyr::filter(PC_2 > -200)

umapCoordClst <- inner_join(umapCoordinates, umapClusters)
ggplot(umapCoordClst, aes(x = UMAP_1, y = UMAP_2, color = sampleNum)) +
  geom_point() + NoLegend()
ggsave(filename = paste0(plotDirectory, "R3_primedUMAP_samples.pdf"), unit = "in", height = 2, width = 2, useDingbats = FALSE)
ggplot(umapCoordClst, aes(x = UMAP_1, y = UMAP_2, color = scanorama_snn_res.0.6)) +
  geom_point() + NoLegend()
ggsave(filename = paste0(plotDirectory, "R3_primedUMAP_clusters.pdf"), unit = "in", height = 2, width = 2, useDingbats = FALSE)

Idents(scTransform) <- scTransform$orig.ident
DefaultAssay(scTransform) <- "scanorama"
FeaturePlot(scTransform, features = c("SPP1", "GDF15", "CDKN1A", "FTH1", "CENPF", "MKI67", "TOP2A", "SOX21"), reduction = "umap_scanorama", slot = 'scale.data', ncol = 4) &
  NoAxes() & scale_color_gradient2(low = "#0039A6", mid = "lightgray", high = "#EE352E", midpoint = 0.5)

DefaultAssay(scTransform) <- "RNA"
scTransform <- NormalizeData(scTransform)
markers <- FindAllMarkers(scTransform, logfc.threshold = 0)
markersFSpecific <- markers %>% filter(cluster == "2_fast")
ggplot(markersFSpecific, aes(x = avg_log2FC, y = -log(p_val_adj))) +
  rasterise(geom_point(), dpi = 300) +
  geom_text_repel(aes(label = gene))
ggsave(filename = paste0(plotDirectory, "R3_volcanoPlot_fastPrimedVsSlowPrimed.pdf"), unit = "in", height = 4, width = 4, useDingbats = FALSE)

pcaLoadings <- Loadings(scTransform, reduction = "pca_scanorama")[, 1:5]
genes <- row.names(pcaLoadings)
pcaLoadings <- as_tibble(pcaLoadings) %>% mutate(genes = genes)

# pcaLoadings <- pcaLoadings %>% rowwise() %>% mutate(max = max(abs(PC_1), abs(PC_2), abs(PC_3), abs(PC_4), abs(PC_5))) %>%
#   ungroup() %>% mutate(mean = rowMeans2(abs(as.matrix(pcaLoadings[,c(1:5)]))), sd = rowSds(abs(as.matrix(pcaLoadings[,c(1:5)]))))

pcaLoadings <- pcaLoadings %>% rowwise() %>% mutate(max = max(abs(PC_4), abs(PC_5))) %>%
  ungroup() %>% mutate(mean = rowMeans2(abs(as.matrix(pcaLoadings[,c(1:5)]))), sd = rowSds(abs(as.matrix(pcaLoadings[,c(1:5)]))))
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = 'https://www.ensembl.org')
pcGenes <- biomaRt::getBM(attributes = c("external_gene_name", "transcript_biotype"), filters = c("transcript_biotype"), values = list("protein_coding"), mart = mart)
pcaLoadingsFilter <- pcaLoadings %>% dplyr::filter(genes %in% pcGenes$external_gene_name)

ggplot(pcaLoadingsFilter, aes(x = mean, y = max)) +
  geom_point() +
  geom_text_repel(data = pcaLoadingsFilter %>% dplyr::slice_max(., n = 50, order_by = max), aes(label = genes), max.overlaps = 20, color = "blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")

# diffGeneList <- pcaLoadingsFilter %>% dplyr::slice_max(., n = 50, order_by = max) %>% .$genes
# diffGeneList <- c("CEBPA", "RGCC", "TSNARE1", "PLEKHD1", "DGCR6", "TSPAN19", "LVRN", "VASH2", "MN1", "TACC3", "IGFBP2", "MT2A", "COX7C", "FTH1", "POLR2L", "TMA7", "RPS27", "SPP1", "TGFB1", "ASH2L", "TGFBR1", "SMAD2", "SMAD3")
# diffGeneList <- c("RGCC", "MT2A", "MN1", "FTH1", "NQO1", "IGFBP4")
# FeaturePlot(scTransform, features = diffGeneList, reduction = "umap_scanorama", slot = 'scale.data', ncol = 3, max.cutoff = "q90") &
#   NoAxes() & NoLegend() & scale_color_gradient2(low = "red", mid = "lightgray", high = "blue", midpoint = 0.5)
# Idents(scTransform) <- scTransform$orig.ident
# VlnPlot(scTransform, features = diffGeneList, slot = "scale.data")

logNormalizedCounts = scTransform@assays$RNA@data
logNormalizedCountsRound = round(logNormalizedCounts, 4)
cells_count = sub("-1", "", colnames(logNormalizedCountsRound))
cells_count_cellID = sub("S\\d_", "", cells_count)
cells_count_sampleNum = gsub("[^S123456]", "", cells_count)
logNormalizedCountsRound = as_tibble(as.data.frame((t(as.matrix(logNormalizedCountsRound)))))
logNormalizedCountsRound = logNormalizedCountsRound %>% mutate(cellID = cells_count_cellID, sampleNum = cells_count_sampleNum)
logNormalizedCountsFilter <- logNormalizedCountsRound %>% dplyr::select(., POLR2L, MT2A, S100A6, FTH1, NQO1, IGFBP4, cellID, sampleNum)
logNormalizedCountsFilterMelt <- melt(logNormalizedCountsFilter, id.vars = c("cellID", "sampleNum"))
logNormalizedCountsFilterMelt$sampleNum <- factor(logNormalizedCountsFilterMelt$sampleNum, levels = c("S3", "S1", "S2"), labels = c("slow", "control", "fast"))

ggplot(logNormalizedCountsFilterMelt %>% dplyr::filter(sampleNum != "control"), aes(x = sampleNum, y = value)) +
  geom_jitter(size = 0.25) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  facet_wrap(~variable, ncol = 3) +
  geom_signif(comparisons = list(c("fast", "slow")))

diffGeneList <- c("POLR2L", "MT2A", "S100A6", "FTH1", "NQO1", "IGFBP7")
DefaultAssay(scTransform) <- "RNA"
scTransform <- ScaleData(scTransform)
FeaturePlot(scTransform, features = diffGeneList, reduction = "umap_scanorama", slot = "scale.data", ncol = 3, , min.cutoff = "q10", max.cutoff = "q90", raster = TRUE, pt.size = 5, raster.dpi = c(300, 300)) &
  NoAxes() & scale_color_gradient2(low = "#0039A6", mid = "lightgray", high = "#EE352E", midpoint = 0.6) & NoLegend() 
ggsave(filename = paste0(plotDirectory, "markerUMAPSubsetPCA.pdf"), units = "in", height = 4, width = 5)

markersProlif <- readRDS(file = "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/extractedData/rewind10X/R3/DEProlifSpeed.rds")
markersProlifComp <- inner_join(markersFSpecific, markersProlif, by = c("gene" = "gene_name"))
ggplot(markersProlifComp, aes(x = avg_log2FC, y = log2FoldChange)) +
  geom_point() +
  geom_text_repel(aes(label = gene)) +
  geom_text_repel(data = markersProlifComp %>% dplyr::filter(gene %in% g2m.genes),
                  aes(label = gene, x = avg_log2FC, y = log2FoldChange), color = "red") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")

#### compare distribution of cells in each cell cycle phase across prolif speed and across primed cells ####
#######################################################################################################################################################
Idents(scTransform) <- scTransform$Phase
DimPlot(scTransform, reduction = "umap_scanorama")
speed = (scTransform[['orig.ident']]) %>% as_tibble()
phase = (scTransform[['Phase']]) %>% as_tibble()
speedPhasePlot <- bind_cols(speed, phase)

speedPhasePlot$Phase <- as.factor(speedPhasePlot$Phase)
speedPhasePlot$orig.ident <- factor(speedPhasePlot$orig.ident, levels = c("3_slow", "1_control", "2_fast"), labels = c("slow", "control", "fast"))
ggplot(speedPhasePlot, aes(x = orig.ident, fill = Phase)) +
  geom_bar(aes(fill = as.factor(Phase)), position = "fill")

Idents(scanorama_filter) <- scanorama_filter$Phase
DimPlot(scanorama_filter, reduction = "umap_scanorama")
speed = (scanorama_filter[['orig.ident']]) %>% as_tibble()
phase = (scanorama_filter[['Phase']]) %>% as_tibble()
speedPhasePlotAll <- bind_cols(speed, phase)

speedPhasePlotAll$Phase <- as.factor(speedPhasePlotAll$Phase)
speedPhasePlotAll$orig.ident <- factor(speedPhasePlotAll$orig.ident, levels = c("3_slow", "1_control", "2_fast"), labels = c("slow_all", "control_all", "fast_all"))
ggplot(speedPhasePlotAll, aes(x = orig.ident, fill = Phase)) +
  geom_bar(aes(fill = as.factor(Phase)), position = "fill")

speedPhasePlotComb <- bind_rows(speedPhasePlot, speedPhasePlotAll)
speedPhasePlotCombSum <-  speedPhasePlotComb %>% dplyr::group_by(orig.ident, Phase) %>% dplyr::summarize(nCells = n()) %>% ungroup() %>%
  group_by(orig.ident) %>% mutate(prop = nCells/sum(nCells))
speedPhasePlotCombSum$orig.ident <-factor(speedPhasePlotCombSum$orig.ident, levels = c("slow_all", "slow", "control_all", "control", "fast_all", "fast"))
ggplot(speedPhasePlotCombSum, aes(x = orig.ident, y = prop, fill = factor(Phase))) +
  geom_col(position = "fill") +
  geom_text(aes(label = paste0(round(prop, 3)*100, "%")), position = position_fill(vjust = 0.5)) + NoLegend()
ggsave(filename = paste0(plotDirectory, "R3_cellPhaseDistPerCond.pdf"), height = 4, width = 4)

speedPhasePlotCombSum$label <- "G1"
speedPhasePlotCombSum$label <- ifelse(speedPhasePlotCombSum$Phase == "G1", speedPhasePlotCombSum$label, "non-G1")
speedPhasePlotCombSumFact <- speedPhasePlotCombSum %>% dplyr::filter(label == "G1")

ggplot(speedPhasePlotCombSumFact, aes(x = orig.ident, y = 1-prop)) +
  geom_col() + NoLegend() + ylim(0, 1)
ggsave(filename = paste0(plotDirectory, "R3_cellPhaseDistPerCond_fracNonG1.pdf"), height = 4, width = 4)

DefaultAssay(scanorama_filter) <- "RNA"
scanorama_filter <- NormalizeData(scanorama_filter)

logNormalizedCountsRound = as_tibble(read.table(file = paste0(homeDirectory, "logNormalizedCounts_Scanorama_50pcs_filterRound.tsv"), header = TRUE, stringsAsFactors = F, sep = "\t"))

cellPhaseVsSPP1 <- logNormalizedCountsRound %>% dplyr::select(cellID, sampleNum, SPP1) %>% bind_cols(phase)
hist(cellPhaseVsSPP1$SPP1, breaks = 100)
sum(cellPhaseVsSPP1$SPP1 < -0.375)/length(cellPhaseVsSPP1$SPP1)
cellPhaseVsSPP1$label <- "all"
SPP12_5PerCells <- cellPhaseVsSPP1 %>% slice_min(., order_by = SPP1, prop = 0.025) %>% mutate(label = "SPP12.5Per")
SPP15PerCells <- cellPhaseVsSPP1 %>% slice_min(., order_by = SPP1, prop = 0.05) %>% mutate(label = "SPP15Per")
SPP110PerCells <- cellPhaseVsSPP1 %>% slice_min(., order_by = SPP1, prop = 0.10) %>% mutate(label = "SPP110Per")
cellPhaseVsSPP1 <- bind_rows(cellPhaseVsSPP1, SPP12_5PerCells, SPP15PerCells, SPP110PerCells)

cellPhaseVsSPP1Plot <- cellPhaseVsSPP1 %>% group_by(label, Phase) %>% dplyr::summarize(nCells = n()) %>% ungroup() %>%
  group_by(label) %>% mutate(prop = nCells/sum(nCells))

cellPhaseVsSPP1Plot$label <- factor(cellPhaseVsSPP1Plot$label, levels = c("all", "SPP110Per", "SPP15Per", "SPP12.5Per"))

ggplot(cellPhaseVsSPP1Plot, aes(x = label, y = prop, fill = factor(Phase))) +
  geom_col(position = "fill") +
  geom_text(aes(label = paste0(round(prop, 3)*100, "%")), position = position_fill(vjust = 0.5)) + NoLegend()
ggsave(filename = paste0(plotDirectory, "R3_cellPhaseDistPerSPP1Level.pdf"), height = 4, width = 4)

ggplot(cellPhaseVsSPP1Plot %>% dplyr::filter(Phase == "G1"), aes(x = label, y = 1-prop)) +
  geom_col() + NoLegend() + ylim(0, 1)
ggsave(filename = paste0(plotDirectory, "R3_cellPhaseDistPerSPP1Level_fracNonG1.pdf"), height = 4, width = 4)

speedPhasePlotAllPlot <- speedPhasePlotAll %>% group_by(Phase) %>% summarise(nCells = n()) %>% ungroup() %>% mutate(prop = nCells/sum(nCells))

ggplot(speedPhasePlotAllPlot, aes(x = "", y = prop, fill = factor(Phase))) +
  geom_col(position = "fill") +
  geom_text(aes(label = paste0(round(prop, 3)*100, "%")), position = position_fill(vjust = 0.5)) + NoLegend()
ggsave(filename = paste0(plotDirectory, "R3_cellPhaseDistPerSPP1Level.pdf"), height = 4, width = 4)

#### compare overlap across samples as well as fraction in G1 per barcode ####
#######################################################################################################################################################
linCountToOverlaps <- read.table(file = paste0(homeDirectory, "filtered10XCells.txt"), sep = "\t")
linCountToOverlapsCycling <- linCountToOverlaps %>% dplyr::filter(SampleNum %in% c("S1", "S2", "S3"))

umapCoordinates = as_tibble(read.table(file = paste0(homeDirectory, "umapCoordinates_Scanorama_50pcs_filter.tsv"), header = TRUE, stringsAsFactors = F, sep = "\t"))

jointUMAP <- inner_join(umapCoordinates, linCountToOverlaps, by = c("cellID", "sampleNum" = "SampleNum"))

phase = as_tibble(data.frame(cellID = rownames(scanorama_filter[['Phase']]) %>% substr(., 4, 19),
                  phase = scanorama_filter[['Phase']]))
speedPhasePlotAll <- inner_join(jointUMAP, phase, by = "cellID") %>% dplyr::select(-nUMI, -fracUMI, -nLineages)
speedPhasePlotAll$cyclingPhase <- "G1"
speedPhasePlotAll$cyclingPhase <- ifelse(speedPhasePlotAll$Phase %in% c("G2M", "S"), "non-G1", speedPhasePlotAll$cyclingPhase)

nCellsPerSpeedPerBC <- speedPhasePlotAll %>% group_by(barcode, sampleNum) %>% summarise(nCells = n()) %>% ungroup() %>%
  complete(barcode, sampleNum) %>% replace(., is.na(.), 0) %>%
  group_by(barcode) %>% mutate(nCellsAll = sum(nCells), prop = nCells / sum(nCells))
nCellsPerSpeedPerBCCast <- dcast(nCellsPerSpeedPerBC, formula = barcode + nCellsAll ~ sampleNum)

nBCAll <- nCellsPerSpeedPerBCCast %>% .$barcode %>% unique() %>% length()
nBCControlAll <- nCellsPerSpeedPerBCCast %>% dplyr::filter(S1 > 0) %>% .$barcode %>% unique %>% length()
nBCFastAll <- nCellsPerSpeedPerBCCast %>% dplyr::filter(S2 > 0) %>% .$barcode %>% unique %>% length()
nBCSlowAll <- nCellsPerSpeedPerBCCast %>% dplyr::filter(S3 > 0) %>% .$barcode %>% unique %>% length()
nBCControlOnly <- nCellsPerSpeedPerBCCast %>% dplyr::filter(S1 == 1) %>% .$barcode %>% unique %>% length()
nBCFastOnly <- nCellsPerSpeedPerBCCast %>% dplyr::filter(S2 == 1) %>% .$barcode %>% unique %>% length()
nBCSlowOnly <- nCellsPerSpeedPerBCCast %>% dplyr::filter(S3 == 1) %>% .$barcode %>% unique %>% length()
nBCControlFast <- nCellsPerSpeedPerBCCast %>% dplyr::filter(S1 > 0 & S2 > 0) %>% .$barcode %>% unique %>% length()
nBCControlSlow <- nCellsPerSpeedPerBCCast %>% dplyr::filter(S1 > 0 & S3 > 0) %>% .$barcode %>% unique %>% length()
nBCFastSlow <- nCellsPerSpeedPerBCCast %>% dplyr::filter(S2 > 0 & S3 > 0) %>% .$barcode %>% unique %>% length()
nBCAllCond <- nCellsPerSpeedPerBCCast %>% dplyr::filter(S1 > 0 & S2 > 0 & S3 > 0) %>% .$barcode %>% unique %>% length()

library(VennDiagram)
grid.newpage()
venndiagram <- draw.triple.venn(area1 = nBCControlAll, area2 = nBCFastAll, area3 = nBCSlowAll,
                                n12 = nBCControlFast, n23 = nBCFastSlow, n13 = nBCControlSlow, n123 = nBCAllCond,
                                euler.d = TRUE, fill = c(alpha("gray", 0.5), alpha("gray", 0.5), alpha("gray", 0.5)),
                                fontfamily = "sans", cex = 2.5, print.mode = "percent", sigdigs = 3)
ggsave(plot = venndiagram, filename = paste0(plotDirectory, 'overlapVennDiagram_allSlowControlFast.pdf'), width = 5, height = 5, useDingbats = F)

nCellsPerSpeedPerBCCast <- dcast(nCellsPerSpeedPerBC, formula = barcode + nCellsAll ~ sampleNum, value.var = "nCells")
overlapCond <- nCellsPerSpeedPerBCCast %>% dplyr::filter(S2 > 2 & S3 > 2)

fracG1PerSpeedPerBarcode <- speedPhasePlotAll %>% dplyr::filter(barcode %in% overlapCond$barcode)
fracG1PerSpeedPerBarcode <- fracG1PerSpeedPerBarcode %>% group_by(barcode, sampleNum, cyclingPhase) %>% summarise(nCells = n()) %>% ungroup() %>%
  complete(barcode, sampleNum, cyclingPhase) %>% replace(., is.na(.), 0) %>%
  group_by(barcode) %>% mutate(nCellsPerBC = sum(nCells)) %>% ungroup() %>% 
  group_by(barcode, nCellsPerBC, sampleNum) %>% mutate(nCellsPerSpeed = n(), prop = nCells/sum(nCells))
fracG1PerSpeedPerBarcode$sampleNum <- factor(fracG1PerSpeedPerBarcode$sampleNum, levels = c("S3", "S1", "S2"), labels = c("slow", "ungated", "fast"))

ggplot(fracG1PerSpeedPerBarcode %>% dplyr::filter(cyclingPhase == "G1"), aes(x = sampleNum, y = prop)) +
  geom_point() + geom_line(aes(group = barcode)) + facet_wrap(~barcode)

ggplot(fracG1PerSpeedPerBarcode %>% dplyr::filter(cyclingPhase == "G1", sampleNum %in% c("slow", "fast")), aes(x = sampleNum, y = prop)) +
  geom_point() + geom_line(aes(group = barcode)) + facet_wrap(~barcode)

popAverage <- 1 - speedPhasePlotAllPlot %>% dplyr::filter(Phase == "G1") %>% .$prop

fracG1PerSpeedPerBarcodeCast <- fracG1PerSpeedPerBarcode %>% dplyr::filter(cyclingPhase == "non-G1", sampleNum %in% c("slow", "fast")) %>%
  dcast(., formula = barcode ~ sampleNum)
ggplot(fracG1PerSpeedPerBarcodeCast, aes(x = fast, y = slow)) +
  geom_point() + geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(aes(y = (26.7 + 12.8)/100, x = (18.2+43.9)/100)) +
  xlim(0, 1) + ylim(0, 1) +
  geom_hline(yintercept = popAverage, linetype = "dashed") + geom_vline(xintercept = popAverage, linetype = "dashed")
ggsave(filename = paste0(plotDirectory, "R3_fracNonCycAcrossTwins.pdf"), height = 4, width = 4)

#### what fraction of primed barcodes in each condition are found in other conditions ####
#######################################################################################################################################################
primedCellsTable <- linCountToOverlaps %>% dplyr::filter(cellID %in% primedAllUMAP$cellID, SampleNum == "S1")
primedCellsOverlap <- nCellsPerSpeedPerBCCast %>% dplyr::filter(barcode %in% primedCellsTable$barcode)

nBCAll <- primedCellsOverlap %>% .$barcode %>% unique() %>% length()
nBCControlAll <- primedCellsOverlap %>% dplyr::filter(S1 > 0) %>% .$barcode %>% unique %>% length()
nBCFastAll <- primedCellsOverlap %>% dplyr::filter(S2 > 0) %>% .$barcode %>% unique %>% length()
nBCSlowAll <- primedCellsOverlap %>% dplyr::filter(S3 > 0) %>% .$barcode %>% unique %>% length()
nBCControlFast <- primedCellsOverlap %>% dplyr::filter(S1 > 0 & S2 > 0) %>% .$barcode %>% unique %>% length()
nBCControlSlow <- primedCellsOverlap %>% dplyr::filter(S1 > 0 & S3 > 0) %>% .$barcode %>% unique %>% length()
nBCFastSlow <- primedCellsOverlap %>% dplyr::filter(S2 > 0 & S3 > 0) %>% .$barcode %>% unique %>% length()
nBCAllCond <- primedCellsOverlap %>% dplyr::filter(S1 > 0 & S2 > 0 & S3 > 0) %>% .$barcode %>% unique %>% length()

library(VennDiagram)
grid.newpage()
venndiagram <- draw.triple.venn(area1 = nBCControlAll, area2 = nBCFastAll, area3 = nBCSlowAll,
                                n12 = nBCControlFast, n23 = nBCFastSlow, n13 = nBCControlSlow, n123 = nBCAllCond,
                                euler.d = TRUE, fill = c(alpha("gray", 0.5), alpha("gray", 0.5), alpha("gray", 0.5)),
                                fontfamily = "sans", cex = 2.5, print.mode = "percent", sigdigs = 3)
ggsave(plot = venndiagram, filename = paste0(plotDirectory, 'overlapVennDiagramPrimed_ungated.pdf'), width = 5, height = 5, useDingbats = F)

primedCellsTable <- linCountToOverlaps %>% dplyr::filter(cellID %in% primedAllUMAP$cellID, SampleNum == "S2")
primedCellsOverlap <- nCellsPerSpeedPerBCCast %>% dplyr::filter(barcode %in% primedCellsTable$barcode)

nBCAll <- primedCellsOverlap %>% .$barcode %>% unique() %>% length()
nBCControlAll <- primedCellsOverlap %>% dplyr::filter(S1 > 0) %>% .$barcode %>% unique %>% length()
nBCFastAll <- primedCellsOverlap %>% dplyr::filter(S2 > 0) %>% .$barcode %>% unique %>% length()
nBCSlowAll <- primedCellsOverlap %>% dplyr::filter(S3 > 0) %>% .$barcode %>% unique %>% length()
nBCControlFast <- primedCellsOverlap %>% dplyr::filter(S1 > 0 & S2 > 0) %>% .$barcode %>% unique %>% length()
nBCControlSlow <- primedCellsOverlap %>% dplyr::filter(S1 > 0 & S3 > 0) %>% .$barcode %>% unique %>% length()
nBCFastSlow <- primedCellsOverlap %>% dplyr::filter(S2 > 0 & S3 > 0) %>% .$barcode %>% unique %>% length()
nBCAllCond <- primedCellsOverlap %>% dplyr::filter(S1 > 0 & S2 > 0 & S3 > 0) %>% .$barcode %>% unique %>% length()

library(VennDiagram)
grid.newpage()
venndiagram <- draw.triple.venn(area1 = nBCControlAll, area2 = nBCFastAll, area3 = nBCSlowAll,
                                n12 = nBCControlFast, n23 = nBCFastSlow, n13 = nBCControlSlow, n123 = nBCAllCond,
                                euler.d = TRUE, fill = c(alpha("gray", 0.5), alpha("gray", 0.5), alpha("gray", 0.5)),
                                fontfamily = "sans", cex = 2.5, print.mode = "percent", sigdigs = 3)
ggsave(plot = venndiagram, filename = paste0(plotDirectory, 'overlapVennDiagramPrimed_fast.pdf'), width = 5, height = 5, useDingbats = F)

primedCellsTable <- linCountToOverlaps %>% dplyr::filter(cellID %in% primedAllUMAP$cellID, SampleNum == "S3")
primedCellsOverlap <- nCellsPerSpeedPerBCCast %>% dplyr::filter(barcode %in% primedCellsTable$barcode)

nBCAll <- primedCellsOverlap %>% .$barcode %>% unique() %>% length()
nBCControlAll <- primedCellsOverlap %>% dplyr::filter(S1 > 0) %>% .$barcode %>% unique %>% length()
nBCFastAll <- primedCellsOverlap %>% dplyr::filter(S2 > 0) %>% .$barcode %>% unique %>% length()
nBCSlowAll <- primedCellsOverlap %>% dplyr::filter(S3 > 0) %>% .$barcode %>% unique %>% length()
nBCControlFast <- primedCellsOverlap %>% dplyr::filter(S1 > 0 & S2 > 0) %>% .$barcode %>% unique %>% length()
nBCControlSlow <- primedCellsOverlap %>% dplyr::filter(S1 > 0 & S3 > 0) %>% .$barcode %>% unique %>% length()
nBCFastSlow <- primedCellsOverlap %>% dplyr::filter(S2 > 0 & S3 > 0) %>% .$barcode %>% unique %>% length()
nBCAllCond <- primedCellsOverlap %>% dplyr::filter(S1 > 0 & S2 > 0 & S3 > 0) %>% .$barcode %>% unique %>% length()

# col = c("red", "green", "blue") #insert for determining which circle is which

library(VennDiagram)
grid.newpage()
venndiagram <- draw.triple.venn(area1 = nBCControlAll, area2 = nBCFastAll, area3 = nBCSlowAll,
                                n12 = nBCControlFast, n23 = nBCFastSlow, n13 = nBCControlSlow, n123 = nBCAllCond,
                                euler.d = TRUE, fill = c(alpha("gray", 0.5), alpha("gray", 0.5), alpha("gray", 0.5)),
                                fontfamily = "sans", cex = 2.5, print.mode = "percent", sigdigs = 3)
ggsave(plot = venndiagram, filename = paste0(plotDirectory, 'overlapVennDiagramPrimed_slow.pdf'), width = 5, height = 5, useDingbats = F)

#### analysis of what markers are important for fast primed cells versus fast nonprimed cells ####
#######################################################################################################################################################
labelsToAdd <- umapCoordinates %>% mutate(label = "none")
labelsToAdd$label <- ifelse(labelsToAdd$cellID %in% linCountToOverlaps$cellID, "barcoded", labelsToAdd$label)
labelsToAdd$label <- ifelse(labelsToAdd$cellID %in% primedAllUMAP$cellID, "primed", labelsToAdd$label)
labelsToAdd <- labelsToAdd %>% dplyr::select(label)

scanorama_filter_labels <- AddMetaData(scanorama_filter, metadata = labelsToAdd)
scanorama_filter_labels$label <- labelsToAdd$label
Idents(scanorama_filter_labels) <- scanorama_filter_labels$label
DimPlot(scanorama_filter_labels)

DefaultAssay(scanorama_filter_labels) <- "RNA"
scanorama_filter_labels <- NormalizeData(scanorama_filter_labels)

Idents(scanorama_filter_labels) <- scanorama_filter_labels$orig.ident
scanorama_filter_labels$label.combined <- paste(Idents(scanorama_filter_labels), scanorama_filter_labels$label, sep = "_")
Idents(scanorama_filter_labels) <- scanorama_filter_labels$label.combined
DimPlot(scanorama_filter_labels)

#### detour for aggregating into pseudobulk for generating PCA
scanorama_aggregate <- scanorama_filter_labels %>% AggregateExpression(return.seurat = TRUE)
scanorama_aggregate <- NormalizeData(scanorama_aggregate)
scanorama_aggregate <- FindVariableFeatures(scanorama_aggregate, nfeatures = 8000)
all.genes <- rownames(scanorama_aggregate)
scanorama_aggregate <- ScaleData(scanorama_aggregate, features = all.genes)
scanorama_aggregate <- RunPCA(scanorama_aggregate, npcs = 8, approx = FALSE)
DimPlot(scanorama_aggregate, reduction = "pca", dims = c(1, 2), pt.size = 5)

scanorama_aggregate_subset <- subset(scanorama_aggregate, idents = c("1_control_barcoded", "2_fast_barcoded", "3_slow_barcoded",
                                                                     "1_control_primed", "2_fast_primed", "3_slow_primed"))
DimPlot(scanorama_aggregate_subset, reduction = "pca", pt.size = 5)
scanorama_aggregate_subset <- RunPCA(scanorama_aggregate_subset, npcs = 5, approx = FALSE)
DimPlot(scanorama_aggregate_subset, reduction = "pca", pt.size = 5)
DimPlot(scanorama_aggregate_subset, reduction = "pca", pt.size = 5) & NoLegend()
ggsave(filename = paste0(plotDirectory, "R3_PCAforCyclingSpeedAndPriming.pdf"), height = 4, width = 4)
pcLoadings <- as_tibble(scanorama_aggregate_subset[["pca"]][, 1:2])
pcLoadings$gene <- rownames(scanorama_aggregate_subset[["pca"]][, 1:2])
pcLoadingsMelt <- melt(pcLoadings, id.vars = "gene")

ggplot(pcLoadingsMelt, aes(x = variable, y = value)) +
  geom_jitter()

markersPrevious <- readRDS("/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/extractedData/rewind10X/markersComb.rds")
negMarkers <- slice_min(markersPrevious, order_by = mean, n = 25)
posMarkers <- slice_max(markersPrevious, order_by = mean, n = 25)

pcLoadingsMeltPlot <- pcLoadingsMelt %>% dplyr::filter(variable == "PC_1")
ggplot() +
  rasterise(geom_jitter(data = pcLoadingsMeltPlot, aes(x = variable, y = value), color = "lightgray"), dpi = 300) +
  geom_point(data = pcLoadingsMeltPlot %>% dplyr::filter(gene %in% negMarkers$gene),
             aes(x = variable, y = value), color = "blue", position = position_jitter(seed = 1234)) +
  geom_text_repel(data = pcLoadingsMeltPlot %>% dplyr::filter(gene %in% negMarkers$gene),
                  aes(x = variable, y = value, label = gene), color = "blue", position = position_jitter(seed = 1234)) +
  geom_point(data = pcLoadingsMeltPlot %>% dplyr::filter(gene %in% posMarkers$gene),
             aes(x = variable, y = value), color = "red", position = position_jitter(seed = 1234)) +
  geom_text_repel(data = pcLoadingsMeltPlot %>% dplyr::filter(gene %in% posMarkers$gene),
                  aes(x = variable, y = value, label = gene), color = "red", position = position_jitter(seed = 1234)) +
  geom_hline(yintercept = 0, linetype = "dashed") + coord_flip()
ggsave(filename = paste0(plotDirectory, "R3_PCAforCyclingSpeedAndPriming_Loadings.pdf"), height = 4, width = 6)

pcLoadingsMeltPlot <- pcLoadingsMelt %>% dplyr::filter(variable == "PC_2")
ggplot() +
  rasterise(geom_jitter(data = pcLoadingsMeltPlot, aes(x = variable, y = value), color = "lightgray"), dpi = 300) +
  geom_point(data = pcLoadingsMeltPlot %>% dplyr::filter(gene %in% negMarkers$gene),
             aes(x = variable, y = value), color = "blue", position = position_jitter(seed = 1234)) +
  geom_text_repel(data = pcLoadingsMeltPlot %>% dplyr::filter(gene %in% negMarkers$gene),
                  aes(x = variable, y = value, label = gene), color = "blue", position = position_jitter(seed = 1234)) +
  geom_point(data = pcLoadingsMeltPlot %>% dplyr::filter(gene %in% posMarkers$gene),
             aes(x = variable, y = value), color = "red", position = position_jitter(seed = 1234)) +
  geom_text_repel(data = pcLoadingsMeltPlot %>% dplyr::filter(gene %in% posMarkers$gene),
                  aes(x = variable, y = value, label = gene), color = "red", position = position_jitter(seed = 1234)) +
  geom_hline(yintercept = 0, linetype = "dashed") + coord_flip()
ggsave(filename = paste0(plotDirectory, "R3_PCAforCyclingSpeedAndPriming_Loadings_PC2.pdf"), height = 4, width = 6)

mat <- Seurat::GetAssayData(scanorama_aggregate_subset, assay = "RNA", slot = "scale.data")
pca <- scanorama_aggregate_subset[["pca"]]
total_variance <- sum(matrixStats::rowVars(mat))
eigValues = (pca@stdev)^2 
varExplained = eigValues / total_variance

# DimPlot(scanorama_aggregate_subset, reduction = "pca", pt.size = 5, dims = c(1, 5))
# pcaPlotAggregate <- data.frame((scanorama_aggregate_subset[["pca"]])@cell.embeddings)
# colnames(pcaPlotAggregate) <- colnames((scanorama_aggregate_subset[["pca"]])@cell.embeddings)
# pcaPlotAggregate$sample <- rownames((scanorama_aggregate_subset[["pca"]])@cell.embeddings)
# 
# ggplot(pcaPlotAggregate, aes(x = PC_1, y = "", color = sample)) +
#   geom_point(size = 5) +
#   geom_text(aes(label = sample), angle = 90, vjust = 0, hjust = -0.25) + NoLegend()
# 
# scanorama_subset <- scanorama_filter_labels %>% subset(idents = c("1_control_barcoded", "2_fast_barcoded", "3_slow_barcoded",
#                                                                   "1_control_primed", "2_fast_primed", "3_slow_primed"))
# scanorama_subset <- NormalizeData(scanorama_subset)
# scanorama_subset <- FindVariableFeatures(scanorama_subset, nfeatures = 8000)
# all.genes <- rownames(scanorama_subset)
# scanorama_subset <- ScaleData(scanorama_subset, features = all.genes)
# scanorama_subset <- RunPCA(scanorama_subset)
# DimPlot(scanorama_subset, reduction = "pca", pt.size = 2, split.by = "label.combined", group.by = "Phase")
# 
# ElbowPlot(scanorama_subset)
# scanorama_subset <- RunUMAP(scanorama_subset, dims = 1:50)
# DimPlot(scanorama_subset, reduction = "umap", pt.size = 2, split.by = "label.combined", group.by = "Phase")
# 
# scanorama_aggregate_subset <- subset(scanorama_aggregate, idents = c("1_control_barcoded", "2_fast_barcoded", "3_slow_barcoded",
#                                                                      "1_control_primed", "2_fast_primed", "3_slow_primed"))
# DimPlot(scanorama_aggregate_subset, reduction = "pca", pt.size = 5)
# scanorama_aggregate_subset <- RunPCA(scanorama_aggregate_subset, npcs = 5, approx = FALSE)
# DimPlot(scanorama_aggregate_subset, reduction = "pca", pt.size = 5)
# pcLoadings <- as_tibble(scanorama_aggregate_subset[["pca"]][, 1:2])
# pcLoadings$gene <- rownames(scanorama_aggregate_subset[["pca"]][, 1:2])
# pcLoadingsMelt <- melt(pcLoadings, id.vars = "gene")
# 
# ggplot(pcLoadingsMelt, aes(x = variable, y = value)) +
#   geom_jitter()

################################################################

# s1s2_scTransform.features <- SelectIntegrationFeatures(object.list = s1s2_scTransform.list, nfeatures = 8000)
# s1s2_scTransform.list <- PrepSCTIntegration(object.list = s1s2_scTransform.list, anchor.features = s1s2_scTransform.features, verbose = FALSE)
# 
# s1s2_scTransform.anchors <- FindIntegrationAnchors(object.list = s1s2_scTransform.list, normalization.method = "SCT", anchor.features = s1s2_scTransform.features, verbose = FALSE)
# s1s2_scTransform.integrated <- IntegrateData(anchorset = s1s2_scTransform.anchors, normalization.method = "SCT", verbose = FALSE)
# 
# rm(sample1)
# rm(sample2)
# rm(sample1.data)
# rm(sample2.data)
# rm(s1s2_scTransform)
# rm(s1s2_scTransform.features)
# rm(s1s2_scTransform.list)
# rm(s1s2_scTransform.anchors)
# 
# s1s2_scTransform.integrated <- ScaleData(s1s2_scTransform.integrated, verbose = FALSE)
# s1s2_scTransform.integrated <- RunPCA(s1s2_scTransform.integrated, verbose = FALSE, nfeatures = 50)

markersPrevious <- readRDS("/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/extractedData/rewind10X/R1/primedMarkersAll.rds")

markersFast <- FindMarkers(object = scanorama_filter_labels, ident.1 = "2_fast_primed", ident.2 = "2_fast_barcoded", logfc.threshold = 0)
markersFast <- markersFast %>% mutate(gene = rownames(.))
markersFast <- markersFast %>% arrange(., avg_log2FC) %>% mutate(number = seq.int(nrow(markersFast)))
ggplot(markersFast, aes(x = number, y = avg_log2FC, label = gene)) +
  geom_point(alpha = 0.1) +
  geom_point(data = filter(markersFast, gene %in% slice_max(markersPrevious, avg_log2FC, n = 25)$gene), aes(x = number, y = avg_log2FC), size = 2.5, color = "red") +
  geom_point(data = filter(markersFast, gene %in% slice_min(markersPrevious, avg_log2FC, n = 25)$gene), aes(x = number, y = avg_log2FC), size = 2.5, color = "blue") +
  geom_text_repel(data = filter(markersFast, gene %in% slice_max(markersPrevious, avg_log2FC, n = 25)$gene), aes(x = number, y = avg_log2FC, label = gene), size = 2.5, color = "red") +
  geom_text_repel(data = filter(markersFast, gene %in% slice_min(markersPrevious, avg_log2FC, n = 25)$gene), aes(x = number, y = avg_log2FC, label = gene), size = 2.5, color = "blue")

scanorama_filter_labels_subset <- subset(scanorama_filter_labels, subset = label.combined %in% c("2_fast_primed", "2_fast_barcoded"))
VlnPlot(scanorama_filter_labels_subset, features = c("CCN2", "ACTA2", "SERPINE2", "PTX3", "CDC20", "CCNB1"), same.y.lims = TRUE, pt.size = 0) &
  stat_summary(fun = mean, geom = "point", size = 1) &
  stat_summary(fun = mean, geom = "line", aes(group = 1)) &
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.ticks.x = element_blank())
ggsave(filename = paste0(plotDirectory, 'nonprimedVsPrimedFastViolinPlots.pdf'), width = 4, height = 3, useDingbats = F)

markersSlow <- FindMarkers(object = scanorama_filter_labels, ident.1 = "3_slow_primed", ident.2 = "3_slow_barcoded", logfc.threshold = 0)
markersSlow <- markersSlow %>% mutate(gene = rownames(.))
markersSlow <- markersSlow %>% arrange(., avg_log2FC) %>% mutate(number = seq.int(nrow(markersSlow)))
ggplot(markersSlow, aes(x = number, y = avg_log2FC, label = gene)) +
  geom_point() +
  geom_point(data = filter(markersSlow, gene %in% slice_max(markersPrevious, avg_log2FC, n = 25)$gene), aes(x = number, y = avg_log2FC), size = 2.5, color = "red") +
  geom_point(data = filter(markersSlow, gene %in% slice_min(markersPrevious, avg_log2FC, n = 25)$gene), aes(x = number, y = avg_log2FC), size = 2.5, color = "blue") +
  geom_text_repel(data = filter(markersSlow, gene %in% slice_max(markersPrevious, avg_log2FC, n = 25)$gene), aes(x = number, y = avg_log2FC, label = gene), size = 2.5, color = "red") +
  geom_text_repel(data = filter(markersSlow, gene %in% slice_min(markersPrevious, avg_log2FC, n = 25)$gene), aes(x = number, y = avg_log2FC, label = gene), size = 2.5, color = "blue")

markersControl <- FindMarkers(object = scanorama_filter_labels, ident.1 = "1_control_primed", ident.2 = "1_control_barcoded", logfc.threshold = 0)
markersControl <- markersControl %>% mutate(gene = rownames(.))
markersControl <- markersControl %>% arrange(., avg_log2FC) %>% mutate(number = seq.int(nrow(markersControl)))
ggplot(markersControl, aes(x = number, y = avg_log2FC, label = gene)) +
  geom_point() +
  geom_point(data = filter(markersControl, gene %in% slice_max(markersPrevious, avg_log2FC, n = 25)$gene), aes(x = number, y = avg_log2FC), size = 2.5, color = "red") +
  geom_point(data = filter(markersControl, gene %in% slice_min(markersPrevious, avg_log2FC, n = 25)$gene), aes(x = number, y = avg_log2FC), size = 2.5, color = "blue") +
  geom_text_repel(data = filter(markersControl, gene %in% slice_max(markersPrevious, avg_log2FC, n = 25)$gene), aes(x = number, y = avg_log2FC, label = gene), size = 2.5, color = "red") +
  geom_text_repel(data = filter(markersControl, gene %in% slice_min(markersPrevious, avg_log2FC, n = 25)$gene), aes(x = number, y = avg_log2FC, label = gene), size = 2.5, color = "blue")

# markersOverlap <- c(intersect(markersFast$gene, markersControl$gene),
#                     intersect(markersSlow$gene, markersControl$gene),
#                     intersect(markersFast$gene, markersSlow$gene)) %>% unique()

markersOverlapInt <- c(intersect(intersect(markersFast$gene, markersControl$gene), markersSlow$gene))

markersAll <- bind_rows(markersControl %>% mutate(label = "control"),
                        markersFast %>% mutate(label = "fast"),
                        markersSlow %>% mutate(label = "slow")) %>% dplyr::filter(gene %in% markersOverlapInt)
markersAllSelect <- markersAll %>% group_by(gene) %>% summarise(mean = mean(avg_log2FC)) %>% dplyr::filter(abs(mean) > 0.2) %>% .$gene
markersAllPlot <- markersAll %>% dplyr::filter(gene %in% markersAllSelect)
markersAllPlot$gene = with(markersAllPlot, reorder(gene, avg_log2FC, mean))
ggplot(markersAllPlot, aes(x = gene, y = avg_log2FC)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.25) +
  stat_summary(fun = mean, geom = "text", aes(label = gene), angle = -45, hjust = 0, vjust = 0, size = 3, position = position_nudge(x = 0.1, y = -0.05)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

numberCut <- 30
plot1 <- ggplot(markersFast, aes(x = 1, y = avg_log2FC, label = gene)) +
  rasterize(geom_jitter(position = position_jitter(seed = 1234), alpha = 0.1), dpi = 100) +
  geom_jitter(data = filter(markersFast, gene %in% slice_max(markersPrevious, avg_log2FC, n = numberCut)$gene), aes(x = 1, y = avg_log2FC), size = 2.5, color = "#ee352e", position = position_jitter(seed = 1234)) +
  geom_jitter(data = filter(markersFast, gene %in% slice_min(markersPrevious, avg_log2FC, n = numberCut)$gene), aes(x = 1, y = avg_log2FC), size = 2.5, color = "#a7a9ac", position = position_jitter(seed = 1234)) +
  geom_text_repel(data = filter(markersFast, gene %in% slice_max(markersPrevious, avg_log2FC, n = numberCut)$gene), aes(x = 1, y = avg_log2FC, label = gene), size = 2.5, color = "#ee352e", position = position_jitter(seed = 1234)) +
  geom_text_repel(data = filter(markersFast, gene %in% slice_min(markersPrevious, avg_log2FC, n = numberCut)$gene), aes(x = 1, y = avg_log2FC, label = gene), size = 2.5, color = "#a7a9ac", position = position_jitter(seed = 1234)) +
  ylim(-2, 0.5) + ggtitle("fast primed / fast nonprimed") +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
plot2 <- ggplot(markersSlow, aes(x = 1, y = avg_log2FC, label = gene)) +
  rasterize(geom_jitter(position = position_jitter(seed = 1234), alpha = 0.1), dpi = 100) +
  geom_jitter(data = filter(markersSlow, gene %in% slice_max(markersPrevious, avg_log2FC, n = numberCut)$gene), aes(x = 1, y = avg_log2FC), size = 2.5, color = "#ee352e", position = position_jitter(seed = 1234)) +
  geom_jitter(data = filter(markersSlow, gene %in% slice_min(markersPrevious, avg_log2FC, n = numberCut)$gene), aes(x = 1, y = avg_log2FC), size = 2.5, color = "#a7a9ac", position = position_jitter(seed = 1234)) +
  geom_text_repel(data = filter(markersSlow, gene %in% slice_max(markersPrevious, avg_log2FC, n = numberCut)$gene), aes(x = 1, y = avg_log2FC, label = gene), size = 2.5, color = "#ee352e", position = position_jitter(seed = 1234)) +
  geom_text_repel(data = filter(markersSlow, gene %in% slice_min(markersPrevious, avg_log2FC, n = numberCut)$gene), aes(x = 1, y = avg_log2FC, label = gene), size = 2.5, color = "#a7a9ac", position = position_jitter(seed = 1234)) +
  ylim(-2, 0.5) + ggtitle("slow primed / slow nonprimed") +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
plot3 <- ggplot(markersControl, aes(x = 1, y = avg_log2FC, label = gene)) +
  rasterize(geom_jitter(position = position_jitter(seed = 1234), alpha = 0.1), dpi = 100) +
  geom_jitter(data = filter(markersControl, gene %in% slice_max(markersPrevious, avg_log2FC, n = numberCut)$gene), aes(x = 1, y = avg_log2FC), size = 2.5, color = "#ee352e", position = position_jitter(seed = 1234)) +
  geom_jitter(data = filter(markersControl, gene %in% slice_min(markersPrevious, avg_log2FC, n = numberCut)$gene), aes(x = 1, y = avg_log2FC), size = 2.5, color = "#a7a9ac", position = position_jitter(seed = 1234)) +
  geom_text_repel(data = filter(markersControl, gene %in% slice_max(markersPrevious, avg_log2FC, n = numberCut)$gene), aes(x = 1, y = avg_log2FC, label = gene), size = 2.5, color = "#ee352e", position = position_jitter(seed = 1234)) +
  geom_text_repel(data = filter(markersControl, gene %in% slice_min(markersPrevious, avg_log2FC, n = numberCut)$gene), aes(x = 1, y = avg_log2FC, label = gene), size = 2.5, color = "#a7a9ac", position = position_jitter(seed = 1234)) +
  ylim(-2, 0.5) + ggtitle("control primed / control nonprimed") +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
plot <- egg::ggarrange(plots = list(plot2, plot3, plot1), ncol = 3)
ggsave(plot, filename = paste0(plotDirectory, "R3_primedMarkersAcrossProlifCondition.pdf"), units = "in", height = 3, width = 5.5, useDingbats = FALSE)

numberCut <- 50

markersPrevHigh <- slice_max(markersPrevious, avg_log2FC, n = numberCut)$gene
markersPrevLow <- slice_min(markersPrevious, avg_log2FC, n = numberCut)$gene

markersFastTop <- markersFast %>% mutate(marker = "none")
markersFastTop$marker <- ifelse(markersFastTop$gene %in% markersPrevHigh, "high", markersFastTop$marker)
markersFastTop$marker <- ifelse(markersFastTop$gene %in% markersPrevLow, "low", markersFastTop$marker)
markersFastTop <- markersFastTop %>% mutate(cond = "fast")

markersSlowTop <- markersSlow %>% mutate(marker = "none")
markersSlowTop$marker <- ifelse(markersSlowTop$gene %in% markersPrevHigh, "high", markersSlowTop$marker)
markersSlowTop$marker <- ifelse(markersSlowTop$gene %in% markersPrevLow, "low", markersSlowTop$marker)
markersSlowTop <- markersSlowTop %>% mutate(cond = "slow")

markersControlTop <- markersControl %>% mutate(marker = "none")
markersControlTop$marker <- ifelse(markersControlTop$gene %in% markersPrevHigh, "high", markersControlTop$marker)
markersControlTop$marker <- ifelse(markersControlTop$gene %in% markersPrevLow, "low", markersControlTop$marker)
markersControlTop <- markersControlTop %>% mutate(cond = "control")

markersTopAll <- bind_rows(markersFastTop, markersSlowTop, markersControlTop)
markersTopAll$marker <- factor(markersTopAll$marker, levels = c("low", "none", "high"), labels = c("negative", "none", "positive"))
markersTopAll$cond <- factor(markersTopAll$cond, levels = c("slow", "control", "fast"))
ggplot(markersTopAll %>% dplyr::filter(marker != "none"), aes(x = cond, y = avg_log2FC)) +
  geom_jitter() +
  geom_boxplot(size = 0.5) +
  geom_signif(comparisons = list(c("slow", "control"), c("control", "fast"))) +
  geom_signif(comparisons = list(c("slow", "fast")), position = position_nudge(x = 0, y = 0.1)) +
  facet_wrap(~marker)
ggsave(filename = paste0(plotDirectory, "R3_primedMarkersAcrossProlifConditionsBoxplot.pdf"), units = "in", height = 3, width = 3, useDingbats = FALSE)

#### analysis of differences for score modules across proliferation and priming ####
#######################################################################################################################################################
genesPluri <- list(c("NANOG", "ALPL", "POU5F1", "SOX2", "PODXL", "LIN28A", "UTF1", "TERT", "ZFP42", "DNMT3B"))
genesFibro <- list(c("LUM", "S100A4", "THY1", "PDGFRA", "COL1A1", "COL5A1", "LOXL1", "FBLN1", "FBLN2", "VTN"))
genesEpi <- list(c("CDH1", "CLDN3", "KRT3", "OCLN", "EPCAM", "ANPEP", "MUC1", "CD24"))
genesMes <- list(c("CDH2", "VIM", "FN1", "ZEB1", "SNAI2", "TWIST1", "TWIST2", "TGFB1"))
genesMyo <- list(c("ACTA2", "TPX2", "TAGLN", "MYL9", "CDKN1A", "CDKN2A", "CCN2", "GLI2"))
scanorama_filter_labels <- AddModuleScore(scanorama_filter_labels, features = genesFibro, name = 'fibroMarkers')
scanorama_filter_labels <- AddModuleScore(scanorama_filter_labels, features = genesPluri, name = 'pluriMarkers')
scanorama_filter_labels <- AddModuleScore(scanorama_filter_labels, features = genesEpi, name = 'epiMarkers')
scanorama_filter_labels <- AddModuleScore(scanorama_filter_labels, features = genesMes, name = 'mesMarkers')
scanorama_filter_labels <- AddModuleScore(scanorama_filter_labels, features = genesMyo, name = 'myoMarkers')

FeaturePlot(scanorama_filter_labels, features = c('fibroMarkers1', 'pluriMarkers1', 'epiMarkers1', 'mesMarkers1', 'myoMarkers1'), min.cutoff = 'q10', max.cutoff = 'q90') & NoAxes()

fibroScore = (scanorama_filter_labels[['fibroMarkers1']]) %>% as_tibble()
pluriScore = (scanorama_filter_labels[['pluriMarkers1']]) %>% as_tibble()
mesScore = (scanorama_filter_labels[['mesMarkers1']]) %>% as_tibble()
epiScore = (scanorama_filter_labels[['epiMarkers1']]) %>% as_tibble()
myoScore = (scanorama_filter_labels[['myoMarkers1']]) %>% as_tibble()
cells_Clusters = sub("-1", "", colnames(scanorama_filter_labels))
cells_Clusters_cellID = sub("S\\d_", "", cells_Clusters)
cells_Clusters_sampleNum = gsub("[^S123456]", "", cells_Clusters)
scoreModules = fibroScore %>% bind_cols(., pluriScore, mesScore, epiScore, myoScore) %>% mutate(cellID = cells_Clusters_cellID, sampleNum = cells_Clusters_sampleNum)

umapCoordinates = as_tibble(read.table(file = paste0(homeDirectory, "umapCoordinates_Scanorama_50pcs_filter.tsv"), header = TRUE, stringsAsFactors = F, sep = "\t"))
labelsToAdd <- umapCoordinates %>% mutate(label = "none")
labelsToAdd$label <- ifelse(labelsToAdd$cellID %in% linCountToOverlaps$cellID, "barcoded", labelsToAdd$label)
labelsToAdd$label <- ifelse(labelsToAdd$cellID %in% primedAllUMAP$cellID, "primed", labelsToAdd$label)
scoreModulesLabel <- inner_join(scoreModules, labelsToAdd, by = c("cellID", "sampleNum"))

scoreModulesLabelMelt <- scoreModulesLabel %>% melt(., id.vars = c("cellID", "sampleNum", "UMAP_1", "UMAP_2", "label"))
scoreModulesLabelMelt$sampleNum <- factor(scoreModulesLabelMelt$sampleNum, levels = c("S3", "S1", "S2"), labels = c("slow", "control", "fast"))

ggplot(scoreModulesLabelMelt %>% filter(label != "none"), aes(x = label, y = value)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + facet_wrap(sampleNum~variable, ncol = 5) +
  stat_summary(fun = mean, geom = "point") +
  geom_signif(comparisons = list(c("barcoded", "primed")), test = "t.test")

ggplot(scoreModulesLabelMelt %>% filter(label != "none"), aes(x = sampleNum, y = value)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + facet_wrap(~variable, ncol = 5) +
  stat_summary(fun = mean, geom = "point") +
  geom_signif(comparisons = list(c("slow", "control"), c("control", "fast")), test = "t.test") +
  geom_signif(comparisons = list(c("slow", "fast")), position = position_nudge(x = 0, y = 0.2), test = "t.test")

ggplot(scoreModulesLabelMelt %>% filter(label == "primed") %>% filter(sampleNum != "control"), aes(x = sampleNum, y = value)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + facet_wrap(~variable, ncol = 5) +
  stat_summary(fun = mean, geom = "point") +
  geom_signif(comparisons = list(c("slow", "fast")), test = "t.test") 

ggplot(scoreModulesLabelMelt %>% filter(label != "none"), aes(x = label, y = value)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + facet_wrap(~variable, ncol = 5) +
  stat_summary(fun = mean, geom = "point") +
  geom_signif(comparisons = list(c("barcoded", "primed")), test = "t.test")

#### identification of and measurement of overlap in primed state markers across replicates ####
#######################################################################################################################################################
linCountToOverlaps <- read.table(file = paste0(homeDirectory, "filtered10XCells.txt"), sep = "\t")

labelsToAdd <- umapCoordinates %>% mutate(label = "none")
labelsToAdd$label <- ifelse(labelsToAdd$cellID %in% linCountToOverlaps$cellID, "barcoded", labelsToAdd$label)
labelsToAdd$label <- ifelse(labelsToAdd$cellID %in% primedAllUMAP$cellID, "primed", labelsToAdd$label)
labelsToAdd <- labelsToAdd %>% dplyr::select(label)

scanorama_filter_labels <- AddMetaData(scanorama_filter, metadata = labelsToAdd)
scanorama_filter_labels$label <- labelsToAdd$label
Idents(scanorama_filter_labels) <- scanorama_filter_labels$label
DimPlot(scanorama_filter_labels)

DefaultAssay(scanorama_filter_labels) <- "RNA"
scanorama_filter_labels <- NormalizeData(scanorama_filter_labels)
markers <- FindMarkers(object = scanorama_filter_labels, ident.1 = "primed", ident.2 = "barcoded", logfc.threshold = 0)
markers <- markers %>% mutate(gene = rownames(.))
saveRDS(object = markers, file = paste0(homeDirectory, "primedMarkersAll.rds"))