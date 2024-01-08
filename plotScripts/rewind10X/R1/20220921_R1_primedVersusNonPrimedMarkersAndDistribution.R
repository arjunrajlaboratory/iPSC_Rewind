library(tidyverse)
library(reshape2)
library(ggsignif)
library(ggrepel)
library(ggrastr)
library(egg)
library(spgs)
library(Seurat)
library(viridis)
library(matrixStats)
options(future.globals.maxSize = 8000 * 1024 ^ 2)
theme_set(theme_classic())

homeDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/extractedData/rewind10X/R1/"
plotDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/plots/rewind10X/R1/"

scanorama <- readRDS(paste0(homeDirectory, "scTransform.integrated"))
scanorama <- FindClusters(object = scanorama, resolution = 0.45)
Idents(scanorama) <- scanorama$seurat_clusters
DimPlot(scanorama)

#### look for cellIDs corresponding to primed cells based on iPSC reads ####
barcodesR1 <- read.table("/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/extractedData/rewind10X/R1/stepThreeStarcodeShavedReads_BC_10XAndGDNA.txt", header = TRUE, sep = '\t')
linCountToOverlaps <- read.table(file = paste0(homeDirectory, "filtered10XCells.txt"), sep = "\t")
probedBarcodesR1 <- barcodesR1 %>% filter(cellID == "dummy") %>% dplyr::select(UMI, BC50StarcodeD8, SampleNum)
probedBarcodesR1$UMI <- as.numeric(probedBarcodesR1$UMI)
probedBarcodesR1 <- probedBarcodesR1 %>% group_by(BC50StarcodeD8, SampleNum) %>% summarise(nUMI = sum(UMI))

umapCoordinates <- read.table(paste0(homeDirectory, "umapCoordinates_Scanorama_50pcs_filter.tsv"), header = T, sep = "\t")
umapClusters <- read.table(paste0(homeDirectory, "umapClusters_Scanorama_50pcs_filter.tsv"), header = T, sep = "\t")
umapCoordClust <- inner_join(umapCoordinates, umapClusters, by = c("cellID", "sampleNum")) %>% rename(cluster = 5)

probedBarcodesTopR1 <- probedBarcodesR1 %>% ungroup() %>% slice_max(., order_by = nUMI, n = 100)
primedCells <- inner_join(probedBarcodesTopR1, linCountToOverlaps, by = "BC50StarcodeD8")
seuratLabelPriming <- left_join(umapCoordClust, primedCells, by = "cellID") %>% mutate(label = "none")
seuratLabelPriming$label <- ifelse(seuratLabelPriming$cellID %in% linCountToOverlaps$cellID, "nonprimed", seuratLabelPriming$label)
seuratLabelPriming$label <- ifelse(seuratLabelPriming$cellID %in% primedCells$cellID, "primed", seuratLabelPriming$label)

#### plot UMAPs with different clusters and primed cells indicated ####
seuratLabelPrimingFilter <- seuratLabelPriming %>% filter(!(cluster %in% c(4, 9, 10))) %>% filter(UMAP_1 > -5)
seuratLabelPrimingFilter$cluster <- factor(seuratLabelPrimingFilter$cluster, levels = seuratLabelPrimingFilter$cluster %>% unique() %>% sort())
seuratLabelPrimingFilterPlot <- seuratLabelPrimingFilter %>% dplyr::select(cluster, UMAP_1, UMAP_2, label)
seuratLabelPrimingFilterPlot$text <- ""
centroids <- aggregate(cbind(UMAP_1, UMAP_2) ~ cluster, seuratLabelPrimingFilterPlot, mean) %>%
  mutate(label = "centroid") %>% mutate(text = cluster)
seuratLabelPrimingFilterPlot <- bind_rows(seuratLabelPrimingFilterPlot, centroids)

ggplot(data = seuratLabelPrimingFilterPlot %>% filter(label != "none"), aes(x = UMAP_1, y = UMAP_2, color = cluster, label = text)) +
  rasterise(geom_point(), dpi = 100) +
  # scale_color_viridis_d(option = "viridis") +
  geom_text_repel(color = "black", max.overlaps = Inf, size = 7.5, seed = 2468, point.size = 12.5, box.padding = 0.75) +
  theme(legend.position = "none")
ggsave(filename = paste0(plotDirectory, "R1_UMAP_clusters.pdf"), unit = "in", height = 5, width = 5)

ggplot() +
  rasterise(geom_point(data = seuratLabelPrimingFilter %>% filter(label == "nonprimed"), aes(x = UMAP_1, y = UMAP_2), color = "#D8D8D8"), dpi = 100) +
  geom_point(data = seuratLabelPrimingFilter %>% filter(label == "primed"), aes(x = UMAP_1, y = UMAP_2), color = "#ED1C24", size = 2.5) +
  geom_text_repel(data = seuratLabelPrimingFilterPlot %>% filter(label != "none"), aes(x = UMAP_1, y = UMAP_2, label = text),
                  color = "black", max.overlaps = Inf, size = 7.5, seed = 2468, point.size = 12.5, box.padding = 0.75)
ggsave(filename = paste0(plotDirectory, "R1_UMAP_primedCells.pdf"), unit = "in", height = 5, width = 5)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
scanorama <- CellCycleScoring(scanorama, s.features = s.genes, g2m.features = g2m.genes, search = TRUE)

ccScoring = (scanorama[['Phase']])
cells_ccScore = sub("-1", "", rownames(ccScoring))
cells_ccScore_cellID = sub("S\\d_", "", cells_ccScore)
cells_ccScore_sampleNum = gsub("[^S123456]", "", cells_ccScore)
ccScoring = as_tibble(ccScoring)
ccScoring = ccScoring %>% mutate(cellID = cells_ccScore_cellID, sampleNum = cells_ccScore_sampleNum)

ccScoringPlot <- left_join(seuratLabelPrimingFilterPlot, umapCoordinates %>% dplyr::select(1:3), by = c("UMAP_1", "UMAP_2")) %>% left_join(., ccScoring, by = "cellID")
ggplot(data = ccScoringPlot %>% filter(label != "none"), aes(x = UMAP_1, y = UMAP_2, color = Phase, label = text)) +
  rasterise(geom_point(data = ccScoringPlot %>% filter(label != "none") %>% filter(label != "centroid")), dpi = 100) +
  # scale_color_viridis_d(option = "viridis") +
  geom_text_repel(color = "black", max.overlaps = Inf, size = 7.5, seed = 2468, point.size = 12.5, box.padding = 0.75) +
  theme(legend.position = "none")
ggsave(filename = paste0(plotDirectory, "R1_UMAP_Phase.pdf"), unit = "in", height = 5, width = 5)

ccScoringPieNonprimed <- ccScoringPlot %>% dplyr::filter(label == "nonprimed") %>% group_by(Phase) %>% summarise(count = n()) %>%
  arrange(desc(Phase)) %>%
  mutate(prop = count / sum(ccScoringPieNonprimed$count) * 100) %>%
  mutate(ypos = cumsum(prop) - 0.5*prop)
ggplot(ccScoringPieNonprimed, aes(x = "", y = prop, fill = Phase)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  coord_polar("y", start = 0) +
  theme_void() + theme(legend.position = "none") +
  geom_text(aes(y = ypos, label = paste0(Phase, "\n", round(prop, 1), "%")), color = "black", size = 6)
ggsave(filename = paste0(plotDirectory, "R1_nonprimedPhasePieChart.pdf"), height = 3, width = 3)

ccScoringPiePrimed <- ccScoringPlot %>% dplyr::filter(label == "primed") %>% group_by(Phase) %>% summarise(count = n()) %>%
  arrange(desc(Phase)) %>%
  mutate(prop = count / sum(ccScoringPiePrimed$count) * 100) %>%
  mutate(ypos = cumsum(prop) - 0.5*prop)
ggplot(ccScoringPiePrimed, aes(x = "", y = prop, fill = Phase)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  coord_polar("y", start = 0) +
  theme_void() + theme(legend.position = "none") +
  geom_text(aes(y = ypos, label = paste0(Phase, "\n", round(prop, 1), "%")), color = "black", size = 6)
ggsave(filename = paste0(plotDirectory, "R1_primedPhasePieChart.pdf"), height = 3, width = 3)

#### plot normalized distribution of primed cells in different clusters ####
colonyCutoff <- c(25, 50, 100, 250, 500, 1000)
plotList <- list()
for(j in 1:length(colonyCutoff)) {
probedBarcodesTopR1 <- probedBarcodesR1 %>% ungroup() %>% slice_max(., order_by = nUMI, n = colonyCutoff[j])
primedCells <- inner_join(probedBarcodesTopR1, linCountToOverlaps, by = "BC50StarcodeD8")
seuratLabelPriming <- left_join(umapCoordClust, primedCells, by = "cellID") %>% mutate(label = "none")
seuratLabelPriming$label <- ifelse(seuratLabelPriming$cellID %in% linCountToOverlaps$cellID, "nonprimed", seuratLabelPriming$label)
seuratLabelPriming$label <- ifelse(seuratLabelPriming$cellID %in% primedCells$cellID, "primed", seuratLabelPriming$label)
seuratLabelPrimingFilter <- seuratLabelPriming %>% filter(!(cluster %in% c(4, 9, 10)))
seuratLabelPrimingFilter$cluster <- factor(seuratLabelPrimingFilter$cluster, levels = seuratLabelPrimingFilter$cluster %>% unique() %>% sort())

clusterTable <- seuratLabelPrimingFilter %>% dplyr::select(UMAP_1, UMAP_2, cellID, cluster, label)
clusterTableAll <- clusterTable %>% group_by(cluster) %>% summarise(cellsperclusterall = n()) %>% ungroup()
clusterTableBarcoded <- clusterTable %>% filter(label != "none") %>% group_by(cluster) %>% summarise(cellsperclusterbc = n()) %>% ungroup()
clusterTablePrimed <- clusterTable %>% filter(label == "primed") %>% group_by(cluster) %>% summarise(cellsperclusterprimed = n()) %>% ungroup()

clusterTableRandomValueList <- list()
for(i in 1:1000) {
  clusterTableRandomTemp <- clusterTable %>% filter(label != "none") %>% sample_n(size = nrow(clusterTable %>% filter(label == "primed"))) %>% group_by(cluster) %>% summarise(cellsperclusterrandom = n()) %>% ungroup()
  clusterTableRandomTemp <- left_join(clusterTableBarcoded, clusterTableRandomTemp, by = "cluster")
  clusterTableRandomTemp[is.na(clusterTableRandomTemp)] <- 0
  clusterTableRandomTemp <-  clusterTableRandomTemp %>% rowwise() %>%
    mutate(propPerClust = cellsperclusterrandom/cellsperclusterbc) %>% ungroup() %>%
    mutate(propPerClustNorm = propPerClust/sum(propPerClust))
  clusterTableRandomValueList[[i]] <- clusterTableRandomTemp$propPerClustNorm
}
rowMeans(sapply(clusterTableRandomValueList, unlist))
rowSds(sapply(clusterTableRandomValueList, unlist))
clusterMatrixRandom <- sapply(clusterTableRandomValueList, unlist)

clusterTableMerge <- left_join(clusterTableAll, clusterTableBarcoded, by = "cluster") %>% left_join(., clusterTablePrimed, by = "cluster") %>% mutate(cellsperclusterrandom = rowMeans(sapply(clusterTableRandomValueList, unlist)))
clusterTableMerge[is.na(clusterTableMerge)] <- 0
clusterTableMerge <- clusterTableMerge %>% rowwise() %>%
  mutate(propPerClust = cellsperclusterprimed/cellsperclusterbc) %>% ungroup() %>%
  mutate(propPerBarClustNorm = propPerClust/sum(propPerClust)) %>%
  mutate(propPerClustAll = cellsperclusterprimed/cellsperclusterall) %>% ungroup() %>%
  mutate(propPerBarClustNormAll = propPerClustAll/sum(propPerClustAll)) %>%
  mutate(propPerBarClustProp = (propPerBarClustNormAll + 0.001)/(propPerBarClustNorm + 0.001)) %>% ungroup()

clusterPValues <- data.frame(cluster = c(0:(nrow(clusterMatrixRandom)-1)), pvalue = 0)
for(i in 1:nrow(clusterMatrixRandom)) {
  clusterPValues$pvalue[i] <- sum(abs(clusterMatrixRandom[i, ] - mean(clusterMatrixRandom[i, ])) > abs(clusterTableMerge$propPerBarClustNorm[i] - mean(clusterMatrixRandom[i, ])))/ncol(clusterMatrixRandom)
}

plotList[[j]] <- ggplot(clusterTableMerge, aes(x = cluster, y = propPerBarClustNorm)) +
  geom_col(fill = "red") + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  geom_hline(yintercept = mean(clusterTableMerge$cellsperclusterrandom), linetype = "dashed")
}
egg::ggarrange(plots = plotList, ncol = 3)

probedBarcodesTopR1 <- probedBarcodesR1 %>% ungroup() %>% slice_max(., order_by = nUMI, n = 100)
primedCells <- inner_join(probedBarcodesTopR1, linCountToOverlaps, by = "BC50StarcodeD8")
seuratLabelPriming <- left_join(umapCoordClust, primedCells, by = "cellID") %>% mutate(label = "none")
seuratLabelPriming$label <- ifelse(seuratLabelPriming$cellID %in% linCountToOverlaps$cellID, "nonprimed", seuratLabelPriming$label)
seuratLabelPriming$label <- ifelse(seuratLabelPriming$cellID %in% primedCells$cellID, "primed", seuratLabelPriming$label)
seuratLabelPrimingFilter <- seuratLabelPriming %>% filter(!(cluster %in% c(4, 9, 10)))
seuratLabelPrimingFilter$cluster <- factor(seuratLabelPrimingFilter$cluster, levels = seuratLabelPrimingFilter$cluster %>% unique() %>% sort())

clusterTable <- seuratLabelPrimingFilter %>% dplyr::select(UMAP_1, UMAP_2, cellID, cluster, label)
clusterTableAll <- clusterTable %>% group_by(cluster) %>% summarise(cellsperclusterall = n()) %>% ungroup()
clusterTableBarcoded <- clusterTable %>% filter(label != "none") %>% group_by(cluster) %>% summarise(cellsperclusterbc = n()) %>% ungroup()
clusterTablePrimed <- clusterTable %>% filter(label == "primed") %>% group_by(cluster) %>% summarise(cellsperclusterprimed = n()) %>% ungroup()

set.seed(1234)
clusterTableRandomValueList <- list()
for(i in 1:10000) {
  clusterTableRandomTemp <- clusterTable %>% filter(label != "none") %>% sample_n(size = nrow(clusterTable %>% filter(label == "primed"))) %>% group_by(cluster) %>% summarise(cellsperclusterrandom = n()) %>% ungroup()
  clusterTableRandomTemp <- left_join(clusterTableBarcoded, clusterTableRandomTemp, by = "cluster")
  clusterTableRandomTemp[is.na(clusterTableRandomTemp)] <- 0
  clusterTableRandomTemp <-  clusterTableRandomTemp %>% rowwise() %>%
    mutate(propPerClust = cellsperclusterrandom/cellsperclusterbc) %>% ungroup() %>%
    mutate(propPerClustNorm = propPerClust/sum(propPerClust))
  clusterTableRandomValueList[[i]] <- clusterTableRandomTemp$propPerClustNorm
}
rowMeans(sapply(clusterTableRandomValueList, unlist))
rowSds(sapply(clusterTableRandomValueList, unlist))
clusterMatrixRandom <- sapply(clusterTableRandomValueList, unlist)

clusterTableMerge <- left_join(clusterTableAll, clusterTableBarcoded, by = "cluster") %>% left_join(., clusterTablePrimed, by = "cluster") %>% mutate(cellsperclusterrandom = rowMeans(sapply(clusterTableRandomValueList, unlist)))
clusterTableMerge[is.na(clusterTableMerge)] <- 0
clusterTableMerge <- clusterTableMerge %>% rowwise() %>%
  mutate(propPerClust = cellsperclusterprimed/cellsperclusterbc) %>% ungroup() %>%
  mutate(propPerBarClustNorm = propPerClust/sum(propPerClust)) %>%
  mutate(propPerClustAll = cellsperclusterprimed/cellsperclusterall) %>% ungroup() %>%
  mutate(propPerBarClustNormAll = propPerClustAll/sum(propPerClustAll)) %>%
  mutate(propPerBarClustProp = (propPerBarClustNormAll + 0.001)/(propPerBarClustNorm + 0.001)) %>% ungroup()

clusterPValues <- data.frame(cluster = c(0:(nrow(clusterMatrixRandom)-1)), pvalue = 0)
for(i in 1:nrow(clusterMatrixRandom)) {
  clusterPValues$pvalue[i] <- sum(abs(clusterMatrixRandom[i, ] - mean(clusterMatrixRandom[i, ])) > abs(clusterTableMerge$propPerBarClustNorm[i] - mean(clusterMatrixRandom[i, ])))/ncol(clusterMatrixRandom)
}

sigLabels <- c("", "", "*", "", "**", "***", "", "")
plot1 <- ggplot(clusterTableMerge, aes(x = cluster, y = propPerBarClustNorm)) +
  geom_col(fill = "#ED1C24") + ylim(0, 0.4) + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  geom_hline(yintercept = mean(clusterTableMerge$cellsperclusterrandom), linetype = "dashed") +
  geom_text(aes(label = sigLabels), vjust = -0.5, size = 5, y = 0.3)
plot2 <- ggplot(clusterTableMerge, aes(x = cluster, y = cellsperclusterrandom)) +
  geom_col(fill = "#D8D8D8") + ylim(0, 0.4) + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  geom_hline(yintercept = mean(clusterTableMerge$cellsperclusterrandom), linetype = "dashed")
plot3 <- ggarrange(plot1, plot2, ncol = 1, nrow = 2)
ggsave(plot3, filename = paste0(plotDirectory, "clusterDistributionTables.pdf"), unit = "in", width = 4, height = 2.5)

#### check for bias in distribution of barcoded cells in UMAP space ####
ggplot() +
  rasterise(geom_point(data = seuratLabelPrimingFilter %>% filter(label == "none"), aes(x = UMAP_1, y = UMAP_2), color = "#D8D8D8"), dpi = 100) +
  geom_point(data = seuratLabelPrimingFilter %>% filter(label != "none"), aes(x = UMAP_1, y = UMAP_2), color = "#06BF06") +
  geom_text_repel(data = seuratLabelPrimingFilterPlot %>% filter(label != "none"), aes(x = UMAP_1, y = UMAP_2, label = text),
                  color = "black", max.overlaps = Inf, size = 7.5, seed = 2468, point.size = 12.5, box.padding = 0.75)
ggsave(filename = paste0(plotDirectory, "R1_UMAP_barcodeddCells.pdf"), unit = "in", height = 5, width = 5)

clusterTable <- seuratLabelPrimingFilter %>% dplyr::select(UMAP_1, UMAP_2, cellID, cluster, label)
clusterTableAll <- clusterTable %>% group_by(cluster) %>% summarise(cellsperclusterall = n()) %>% ungroup()
clusterTableBarcoded <- clusterTable %>% filter(label != "none") %>% group_by(cluster) %>% summarise(cellsperclusterbc = n()) %>% ungroup()

set.seed(1234)
clusterTableRandomValueList <- list()
for(i in 1:10000) {
  clusterTableRandomTemp <- clusterTable %>% sample_n(size = nrow(clusterTable %>% filter(label != "none"))) %>% group_by(cluster) %>% summarise(cellsperclusterrandom = n()) %>% ungroup()
  clusterTableRandomTemp <- left_join(clusterTableAll, clusterTableRandomTemp, by = "cluster")
  clusterTableRandomTemp[is.na(clusterTableRandomTemp)] <- 0
  clusterTableRandomTemp <-  clusterTableRandomTemp %>% rowwise() %>%
    mutate(propPerClust = cellsperclusterrandom/cellsperclusterall) %>% ungroup() %>%
    mutate(propPerClustNorm = propPerClust/sum(propPerClust))
  clusterTableRandomValueList[[i]] <- clusterTableRandomTemp$propPerClustNorm
}
rowMeans(sapply(clusterTableRandomValueList, unlist))
rowSds(sapply(clusterTableRandomValueList, unlist))
clusterMatrixRandom <- sapply(clusterTableRandomValueList, unlist)

clusterTableMerge <- left_join(clusterTableAll, clusterTableBarcoded, by = "cluster") %>% mutate(cellsperclusterrandom = rowMeans(sapply(clusterTableRandomValueList, unlist)))
clusterTableMerge[is.na(clusterTableMerge)] <- 0
clusterTableMerge <- clusterTableMerge %>% rowwise() %>%
  mutate(propPerClust = cellsperclusterbc/cellsperclusterall) %>% ungroup() %>%
  mutate(propPerBarClustNorm = propPerClust/sum(propPerClust))

clusterPValues <- data.frame(cluster = c(0:(nrow(clusterMatrixRandom)-1)), pvalue = 0)
for(i in 1:nrow(clusterMatrixRandom)) {
  clusterPValues$pvalue[i] <- sum(abs(clusterMatrixRandom[i, ] - mean(clusterMatrixRandom[i, ])) > abs(clusterTableMerge$propPerBarClustNorm[i] - mean(clusterMatrixRandom[i, ])))/ncol(clusterMatrixRandom)
}

sigLabels <- c("***", "", "****", "****", "****", "**", "****", "****")
plot1 <- ggplot(clusterTableMerge, aes(x = cluster, y = propPerBarClustNorm)) +
  geom_col(fill = "#06BF06") + ylim(0, 0.2) + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  geom_hline(yintercept = mean(clusterTableMerge$cellsperclusterrandom), linetype = "dashed") +
  geom_text(aes(label = sigLabels), vjust = -0.5, size = 5, y = 0.175)
plot2 <- ggplot(clusterTableMerge, aes(x = cluster, y = cellsperclusterrandom)) +
  geom_col(fill = "#D8D8D8") + ylim(0, 0.2) + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  geom_hline(yintercept = mean(clusterTableMerge$cellsperclusterrandom), linetype = "dashed")
plot3 <- ggarrange(plot1, plot2, ncol = 1, nrow = 2)
ggsave(plot3, filename = paste0(plotDirectory, "barcodedCellClusterDistributionTables.pdf"), unit = "in", width = 4, height = 2.5)

#### expression of genes of interest in specific clusters ####
scanorama <- AddMetaData(scanorama, seuratLabelPriming$label, col.name = "priming")
DefaultAssay(scanorama) <- "RNA"
scanorama <- NormalizeData(scanorama)
Idents(scanorama) <- scanorama$seurat_clusters
scanorama.subset <- subset(scanorama, subset = seurat_clusters != 4 & seurat_clusters != 9 & seurat_clusters != 10)
DimPlot(scanorama.subset)

Idents(scanorama.subset) <- scanorama.subset$priming
markers <- FindMarkers(scanorama.subset, ident.1 = "primed", ident.2 = "nonprimed", logfc.threshold = 0)
markers <- markers %>% mutate(gene = row.names(.))
saveRDS(object = markers, file = paste0(homeDirectory, "primedMarkersAll.rds"))
markers <- readRDS(file = paste0(homeDirectory, "primedMarkersAll.rds"))

# Idents(scanorama.subset) <- scanorama.subset$seurat_clusters
# clustMarkers <- FindAllMarkers(scanorama.subset, only.pos = TRUE)
# saveRDS(object = clustMarkers, file = paste0(homeDirectory, "clusterMarkersAll.rds"))
clustMarkers <- readRDS(file = paste0(homeDirectory, "clusterMarkersAll.rds"))
pos <- position_jitter(width = 0.3, seed = 2)
ggplot(clustMarkers, aes(x = cluster, y = avg_log2FC, label = gene)) +
  geom_point(position = pos) +
  geom_text_repel(position = pos, size = 2.5)
