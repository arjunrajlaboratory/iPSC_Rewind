rm(list=ls())
gc()

library(tidyverse)
library(reshape2)
library(ggrepel)
library(ggforce)
library(ggsignif)
library(ggridges)
library(ggrastr)
library(ggpubr)
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
umapCoordinates = as_tibble(read.table(file = paste0(homeDirectory, "umapCoordinates_Scanorama_50pcs_filter.tsv"), header = TRUE, stringsAsFactors = F, sep = "\t"))
umapClusters = as_tibble(read.table(file = paste0(homeDirectory, "umapClusters_Scanorama_50pcs_filter.tsv"), header = TRUE, stringsAsFactors = F, sep = "\t"))
umapClusters <- umapClusters %>% dplyr::rename(cluster = scanorama_snn_res.0.3)
umapCoordClust <- inner_join(umapCoordinates, umapClusters, by = c("cellID", "sampleNum"))

primedCellIDList <- readRDS(file = paste0(homeDirectory, "primedCellIDList.rds"))
cutoffList <- c(10, 25, 50, 100, 150, 200, 250, 500, 1000)
i = 6
primedAllUMAP <- filter(umapCoordinates, cellID %in% c(unlist(primedCellIDList[[i]]), unlist(primedCellIDList[[i+length(cutoffList)]]), unlist(primedCellIDList[[i+2*length(cutoffList)]]))) %>%
                          dplyr::filter(sampleNum %in% c("S1", "S2", "S3"))

linCountToOverlaps <- read.table(file = paste0(homeDirectory, "filtered10XCells.txt"), sep = "\t")

#### plot UMAPs for distribution of primed cells across Seurat clusters ####
#######################################################################################################################################################
seuratLabelPriming <- left_join(umapCoordClust, primedAllUMAP, by = c("cellID", "sampleNum", "UMAP_1", "UMAP_2")) %>% mutate(label = "none")
seuratLabelPriming$label <- ifelse(seuratLabelPriming$cellID %in% linCountToOverlaps$cellID, "nonprimed", seuratLabelPriming$label)
seuratLabelPriming$label <- ifelse(seuratLabelPriming$cellID %in% primedAllUMAP$cellID, "primed", seuratLabelPriming$label)
seuratLabelPriming$cluster <- factor(seuratLabelPriming$cluster, levels = seuratLabelPriming$cluster %>% unique() %>% sort())
seuratLabelPrimingFilterPlot <- seuratLabelPriming %>% dplyr::select(cluster, UMAP_1, UMAP_2, label, sampleNum)
seuratLabelPrimingFilterPlot$text <- ""
centroids <- aggregate(cbind(UMAP_1, UMAP_2) ~ cluster, seuratLabelPrimingFilterPlot, mean) %>%
  mutate(label = "centroid") %>% mutate(text = cluster)
seuratLabelPrimingFilterPlot <- bind_rows(seuratLabelPrimingFilterPlot, centroids)

ggplot(data = seuratLabelPrimingFilterPlot %>% filter(label != "none"), aes(x = UMAP_1, y = UMAP_2, color = cluster, label = text)) +
  rasterise(geom_point(), dpi = 100) +
  geom_text_repel(color = "black", max.overlaps = Inf, size = 7.5, seed = 1234, point.size = 5, box.padding = 0.75) +
  theme(legend.position = "none") + NoAxes()
ggsave(filename = paste0(plotDirectory, "R3_UMAP_clusters.pdf"), unit = "in", height = 5, width = 5)

ggplot(data = seuratLabelPrimingFilterPlot %>% filter(label != "none"), aes(x = UMAP_1, y = UMAP_2, label = text)) +
  rasterise(geom_point(data = seuratLabelPrimingFilterPlot %>% filter(sampleNum == "S2"), aes(x = UMAP_1, y = UMAP_2), size = 2.5, alpha = 0.25, color = "black")) +
  geom_text_repel(color = "black", max.overlaps = Inf, size = 7.5, seed = 1234, point.size = 5, box.padding = 0.75) +
  theme(legend.position = "none") + NoAxes()
ggsave(filename = paste0(plotDirectory, "R3_UMAP_clusters_fast.pdf"), unit = "in", height = 5, width = 5)

ggplot(data = seuratLabelPrimingFilterPlot %>% filter(label != "none"), aes(x = UMAP_1, y = UMAP_2, label = text)) +
  rasterise(geom_point(data = seuratLabelPrimingFilterPlot %>% filter(sampleNum == "S3"), aes(x = UMAP_1, y = UMAP_2), size = 2.5, alpha = 0.25, color = "black")) +
  geom_text_repel(color = "black", max.overlaps = Inf, size = 7.5, seed = 1234, point.size = 5, box.padding = 0.75) +
  theme(legend.position = "none") + NoAxes()
ggsave(filename = paste0(plotDirectory, "R3_UMAP_clusters_slow.pdf"), unit = "in", height = 5, width = 5)

ggplot() +
  rasterise(geom_point(data = seuratLabelPriming %>% filter(label == "nonprimed"), aes(x = UMAP_1, y = UMAP_2), color = "lightgray"), dpi = 100) +
  geom_point(data = seuratLabelPriming %>% filter(label == "primed"), aes(x = UMAP_1, y = UMAP_2), color = "#ee352e", size = 2.5) +
  geom_text_repel(data = seuratLabelPrimingFilterPlot %>% filter(label != "none"), aes(x = UMAP_1, y = UMAP_2, label = text),
                  color = "black", max.overlaps = Inf, size = 7.5, seed = 1234, point.size = 5, box.padding = 0.75) + NoAxes()
ggsave(filename = paste0(plotDirectory, "R3_UMAP_primedCells.pdf"), unit = "in", height = 5, width = 5)

ggplot() +
  rasterise(geom_point(data = seuratLabelPriming %>% filter(label == "nonprimed", sampleNum == "S2"), aes(x = UMAP_1, y = UMAP_2), color = "lightgray"), dpi = 100) +
  geom_point(data = seuratLabelPriming %>% filter(label == "primed", sampleNum == "S2"), aes(x = UMAP_1, y = UMAP_2), color = "#ee352e", size = 2.5) +
  geom_text_repel(data = seuratLabelPrimingFilterPlot %>% filter(label != "none"), aes(x = UMAP_1, y = UMAP_2, label = text),
                  color = "black", max.overlaps = Inf, size = 7.5, seed = 1234, point.size = 5, box.padding = 0.75) + NoAxes()
ggsave(filename = paste0(plotDirectory, "R3_UMAP_primedFast.pdf"), unit = "in", height = 5, width = 5)

ggplot() +
  rasterise(geom_point(data = seuratLabelPriming %>% filter(label == "nonprimed", sampleNum == "S3"), aes(x = UMAP_1, y = UMAP_2), color = "lightgray"), dpi = 100) +
  geom_point(data = seuratLabelPriming %>% filter(label == "primed", sampleNum == "S3"), aes(x = UMAP_1, y = UMAP_2), color = "#ee352e", size = 2.5) +
  geom_text_repel(data = seuratLabelPrimingFilterPlot %>% filter(label != "none"), aes(x = UMAP_1, y = UMAP_2, label = text),
                  color = "black", max.overlaps = Inf, size = 7.5, seed = 1234, point.size = 5, box.padding = 0.75) + NoAxes()
ggsave(filename = paste0(plotDirectory, "R3_UMAP_primedSlow.pdf"), unit = "in", height = 5, width = 5)

FeaturePlot(scanorama_filter, features = c("SPP1", "GDF15", "CDKN1A", "FTH1", "CENPF", "MKI67", "TOP2A", "SOX21"), max.cutoff = "q80", reduction = "umap_scanorama", slot = 'scale.data', ncol = 4, raster = TRUE, raster.dpi = c(250, 250)) &
  NoAxes() & NoLegend() & scale_color_gradient2(low = "#ee352e", mid = "lightgray", high = "#0039a6", midpoint = 0.25)
ggsave(filename = paste0(plotDirectory, "R3_UMAP_primedMarkerDistribution.pdf"), unit = "in", height = 10, width = 20)

#### check distribution of all cells per condition ####
#######################################################################################################################################################
clustTableControl <- umapClusters %>% filter(sampleNum == "S1") %>% group_by(cluster) %>% summarise(cellsPerClustControl = n()) %>%
  mutate(propPerClustControl = cellsPerClustControl/sum(cellsPerClustControl)) %>% ungroup()
clustTableFast <- umapClusters %>% filter(sampleNum == "S2") %>% group_by(cluster) %>% summarise(cellsPerClustFast = n()) %>%
  mutate(propPerClustFast = cellsPerClustFast/sum(cellsPerClustFast)) %>% ungroup()
clustTableSlow <- umapClusters %>% filter(sampleNum == "S3") %>% group_by(cluster) %>% summarise(cellsPerClustSlow = n()) %>%
  mutate(propPerClustSlow = cellsPerClustSlow/sum(cellsPerClustSlow)) %>% ungroup()

clustTableMergeFast <- left_join(clustTableControl, clustTableFast, by = "cluster") %>%
  mutate(propDiff = propPerClustFast - propPerClustControl) %>% mutate(condition = "fast")
clustTableMergeSlow <- left_join(clustTableControl, clustTableSlow, by = "cluster") %>%
  mutate(propDiff = propPerClustSlow - propPerClustControl) %>% mutate(condition = "slow")
clustTableMerge <- bind_rows(clustTableMergeFast, clustTableMergeSlow)

clustTableMerge <- left_join(clustTableControl, clustTableFast, by = "cluster") %>% left_join(., clustTableSlow, by = "cluster") %>%
  mutate(propDiffFast = propPerClustFast - propPerClustControl) %>% mutate(propDiffSlow = propPerClustSlow - propPerClustControl) %>%
  mutate(propDiffDiff = propDiffFast - propDiffSlow)
clustTableMerge$cluster <- factor(clustTableMerge$cluster, levels = clustTableMerge$cluster %>% unique() %>% sort())

ggplot(clustTableMerge, aes(x = cluster, y = propDiffDiff)) +
  geom_col() +
  ylab("relative proportion in fast cells - relative proportion in slow cells") +
  geom_text(aes(x = 1, y = Inf, label = "enriched in fast cells"), vjust = 2, hjust = "left") +
  geom_text(aes(x = 1, y = -Inf, label = "enriched in slow cells"), vjust = -2, hjust = "left")
ggsave(filename = paste0(plotDirectory, "diffDiffProportionCLusters.pdf"), height = 4, width = 8)

#### check distribution of primed cells per condition ####
#######################################################################################################################################################
plotList <- list()
for(j in 1:length(cutoffList)){
  primedAllUMAP <- filter(umapCoordinates, cellID %in% c(unlist(primedCellIDList[[j]]), unlist(primedCellIDList[[j+length(cutoffList)]]), unlist(primedCellIDList[[j+2*length(cutoffList)]]))) %>%
    dplyr::filter(sampleNum %in% c("S1", "S2", "S3"))
  seuratLabelPriming <- left_join(umapClusters, primedAllUMAP, by = c("cellID", "sampleNum")) %>% mutate(label = "none")
  seuratLabelPriming$label <- ifelse(seuratLabelPriming$cellID %in% linCountToOverlaps$cellID, "nonprimed", seuratLabelPriming$label)
  seuratLabelPriming$label <- ifelse(seuratLabelPriming$cellID %in% primedAllUMAP$cellID, "primed", seuratLabelPriming$label)
  seuratLabelPriming$cluster <- factor(seuratLabelPriming$cluster, levels = seuratLabelPriming$cluster %>% unique() %>% sort())
  
  clusterTable <- seuratLabelPriming %>% dplyr::select(UMAP_1, UMAP_2, cellID, cluster, label, sampleNum)
  clusterTableAll <- clusterTable %>% group_by(cluster) %>% summarise(cellsperclusterall = n()) %>% ungroup()
  clusterTableBarcoded <- clusterTable %>% filter(label != "none") %>% group_by(cluster) %>% summarise(cellsperclusterbc = n()) %>% ungroup()
  clusterTablePrimed <- clusterTable %>% filter(label == "primed") %>% group_by(cluster) %>% summarise(cellsperclusterprimed = n()) %>% ungroup()
  
  clusterTableRandomValueList <- list()
  for(i in 1:10) {
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

plotListList <- list()
sampleList <- c("S1", "S2", "S3")
for(k in 1:length(sampleList)) {
plotList <- list()
for(j in 1:length(cutoffList)){
  primedAllUMAP <- filter(umapCoordinates, cellID %in% c(unlist(primedCellIDList[[j]]), unlist(primedCellIDList[[j+length(cutoffList)]]), unlist(primedCellIDList[[j+2*length(cutoffList)]]))) %>%
    dplyr::filter(sampleNum %in% c("S1", "S2", "S3"))
  seuratLabelPriming <- left_join(umapClusters, primedAllUMAP, by = c("cellID", "sampleNum")) %>% mutate(label = "none")
  seuratLabelPriming$label <- ifelse(seuratLabelPriming$cellID %in% linCountToOverlaps$cellID, "nonprimed", seuratLabelPriming$label)
  seuratLabelPriming$label <- ifelse(seuratLabelPriming$cellID %in% primedAllUMAP$cellID, "primed", seuratLabelPriming$label)
  seuratLabelPriming$cluster <- factor(seuratLabelPriming$cluster, levels = seuratLabelPriming$cluster %>% unique() %>% sort())
  
  clusterTable <- seuratLabelPriming %>% dplyr::select(UMAP_1, UMAP_2, cellID, cluster, label, sampleNum) %>% dplyr::filter(sampleNum == sampleList[[k]])
  clusterTableAll <- clusterTable %>% group_by(cluster) %>% summarise(cellsperclusterall = n()) %>% ungroup()
  clusterTableBarcoded <- clusterTable %>% filter(label != "none") %>% group_by(cluster) %>% summarise(cellsperclusterbc = n()) %>% ungroup()
  clusterTablePrimed <- clusterTable %>% filter(label == "primed") %>% group_by(cluster) %>% summarise(cellsperclusterprimed = n()) %>% ungroup()
  
  clusterTableRandomValueList <- list()
  for(i in 1:10) {
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
ggsave(plot = egg::ggarrange(plots = plotList, ncol = 3), filename = paste0(plotDirectory, "R3_primedDistributionWithin", sampleList[[k]], ".pdf"), units = "in", height = 10, width = 10)
}

#### analysis of number of cells initial population per condition for control versus fast versus slow ####
#######################################################################################################################################################
NCellsPerBarcodeTable <- linCountToOverlaps %>% dplyr::filter(SampleNum %in% c("S1", "S2", "S3")) %>%
  group_by(barcode, SampleNum) %>% summarize(cellNumber = n())
NCellsPerBarcodeTable$SampleNum <- factor(NCellsPerBarcodeTable$SampleNum, levels = c("S3", "S1", "S2"), labels = c("slow", "control", "fast"))

ggplot(NCellsPerBarcodeTable, aes(x = SampleNum, y = cellNumber)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "text", aes(label = round(..y.., 2)), color = "red", vjust = -1) +
  stat_compare_means(comparisons = list(c("control", "fast"), c("fast", "slow"), c("control", "slow")), method = "t.test") +
  xlab("number of cells per lineage in initial population") +
  theme(axis.title.y = element_blank()) + NoLegend()

meansProlif <- NCellsPerBarcodeTable %>% group_by(SampleNum) %>% summarise(mean = mean(cellNumber))
NCellsPerBarcodeTable <- NCellsPerBarcodeTable %>% mutate(mean = as.double(meansProlif[1, 2]))
NCellsPerBarcodeTable$mean <- ifelse(NCellsPerBarcodeTable$SampleNum == "control", as.double(meansProlif[2, 2]), NCellsPerBarcodeTable$mean)
NCellsPerBarcodeTable$mean <- ifelse(NCellsPerBarcodeTable$SampleNum == "fast", as.double(meansProlif[3, 2]), NCellsPerBarcodeTable$mean)

ggplot(NCellsPerBarcodeTable, aes(x = cellNumber)) +
  geom_histogram(aes(y = ..density..), binwidth = 1) +
  geom_vline(aes(xintercept = mean), linetype = "dashed") +
  facet_wrap(~SampleNum, nrow = 1) + xlim(0.5, 5) +
  xlab("number of cells per lineage in initial population") +
  theme(axis.title.y = element_blank()) + NoLegend()
ggsave(filename = paste0(plotDirectory, "numberCellsInitialPop_prolif.pdf"), units = "in", width = 4.5, height = 2)

primedAllBarcodes <- primedAllUMAP %>% left_join(linCountToOverlaps, by = "cellID") %>% .$barcode

NCellsPerBarcodePrimed <- NCellsPerBarcodeTable %>% mutate(label = "nonprimed")
NCellsPerBarcodePrimed$label <- ifelse(NCellsPerBarcodePrimed$barcode %in% primedAllBarcodes, "primed", NCellsPerBarcodePrimed$label)

ggplot(NCellsPerBarcodePrimed, aes(x = label, y = cellNumber)) +
  geom_boxplot(outlier.shape = NA) +
  stat_summary(fun = mean, geom = "text", aes(label = round(..y.., 2)), color = "red", vjust = -1) +
  geom_signif(comparisons = list(c("nonprimed", "primed")), test = "t.test") +
  ylab("number of cells per lineage in initial population") +
  theme(axis.title.x = element_blank()) + NoLegend()

NCellsPerBarcodePrimedPlot <- NCellsPerBarcodePrimed %>% filter(label == "primed")
NCellsPerBarcodeNonprimedPlot <- NCellsPerBarcodePrimed %>% filter(label == "nonprimed")
plot1 <- ggplot() +
  geom_histogram(data = NCellsPerBarcodePrimedPlot, aes(x = cellNumber, y = ..density..), binwidth = 1, fill = "black") +
  xlab("number of cells per barcode") + theme_classic() + xlim(0.5, 5) + ylim(0, 1) +
  geom_vline(xintercept = mean(NCellsPerBarcodePrimedPlot$cellNumber), linetype = "dashed") +
  annotate(geom = "text", label = paste0("mean = ", round(mean(NCellsPerBarcodePrimedPlot$cellNumber), 2)), y = Inf, x = mean(NCellsPerBarcodePrimedPlot$cellNumber), vjust = 1, hjust = -0.2) +
  ggtitle("distribution for primed barcoded cells")
plot2 <- ggplot() +
  geom_histogram(data = NCellsPerBarcodeNonprimedPlot, aes(x = cellNumber, y = ..density..), binwidth = 1, fill = "black") +
  xlab("number of cells per barcode") + theme_classic() + xlim(0.5, 5) + ylim(0, 1) +
  geom_vline(xintercept = mean(NCellsPerBarcodeNonprimedPlot$cellNumber), linetype = "dashed") +
  annotate(geom = "text", label = paste0("mean = ", round(mean(NCellsPerBarcodeNonprimedPlot$cellNumber), 2)), y = Inf, x = mean(NCellsPerBarcodeNonprimedPlot$cellNumber), vjust = 1, hjust = -0.2) +
  ggtitle("distribution for nonprimed barcoded cells")
plotComb <- egg::ggarrange(plot1, plot2, nrow = 2)
ggsave(plotComb, file = paste0(plotDirectory, "R3_cellNumberDistributionPlotPrimedVersusNonprimed.pdf"), height = 4, width = 2)

ggplot(NCellsPerBarcodePrimed, aes(x = label, y = cellNumber)) +
  geom_boxplot() +
  geom_jitter() +
  stat_summary(fun = mean, geom = "text", aes(label = round(..y.., 2)), color = "red", vjust = -1) +
  geom_signif(comparisons = list(c("nonprimed", "primed")), test = "t.test") +
  ylab("number of cells per lineage in initial population") +
  facet_wrap(~SampleNum) +
  theme(axis.title.x = element_blank()) + NoLegend()

for(i in c("slow", "control", "fast")) {
  NCellsPerBarcodePrimedPlot <- NCellsPerBarcodePrimed %>% filter(label == "primed", SampleNum == i)
  NCellsPerBarcodeNonprimedPlot <- NCellsPerBarcodePrimed %>% filter(label == "nonprimed", SampleNum == i)
  plot1 <- ggplot() +
    geom_histogram(data = NCellsPerBarcodePrimedPlot, aes(x = cellNumber, y = ..density..), binwidth = 1, fill = "black") +
    xlab("number of cells per barcode") + theme_classic() + xlim(0.5, 5) + ylim(0, 1) +
    geom_vline(xintercept = mean(NCellsPerBarcodePrimedPlot$cellNumber), linetype = "dashed") +
    annotate(geom = "text", label = paste0("mean = ", round(mean(NCellsPerBarcodePrimedPlot$cellNumber), 2)), y = Inf, x = mean(NCellsPerBarcodePrimedPlot$cellNumber), vjust = 1, hjust = -0.2) +
    ggtitle("distribution for primed barcoded cells")
  plot2 <- ggplot() +
    geom_histogram(data = NCellsPerBarcodeNonprimedPlot, aes(x = cellNumber, y = ..density..), binwidth = 1, fill = "black") +
    xlab("number of cells per barcode") + theme_classic() + xlim(0.5, 5) + ylim(0, 1) +
    geom_vline(xintercept = mean(NCellsPerBarcodeNonprimedPlot$cellNumber), linetype = "dashed") +
    annotate(geom = "text", label = paste0("mean = ", round(mean(NCellsPerBarcodeNonprimedPlot$cellNumber), 2)), y = Inf, x = mean(NCellsPerBarcodeNonprimedPlot$cellNumber), vjust = 1, hjust = -0.2) +
    ggtitle("distribution for nonprimed barcoded cells")
  plotComb <- egg::ggarrange(plot1, plot2, nrow = 2)
  ggsave(plotComb, file = paste0(plotDirectory, "R3_cellNumberDistributionPlotPrimedVersusNonprimed_", i, ".pdf"), height = 4, width = 2)
}

ggplot(NCellsPerBarcodePrimed, aes(x = SampleNum, y = cellNumber)) +
  geom_boxplot() +
  stat_summary(fun = mean, geom = "text", aes(label = round(..y.., 2)), color = "red", vjust = -1) +
  stat_compare_means(comparisons = list(c("control", "fast"), c("fast", "slow"), c("control", "slow")), method = "t.test") +
  xlab("number of cells per lineage in initial population") +
  theme(axis.title.y = element_blank()) + NoLegend() + facet_wrap(~label)

#### analysis of fraction of barcodes per condition for control versus fast versus slow ####
#######################################################################################################################################################
cellNumberList <- c(filter(linCountToOverlaps, SampleNum == "S1")$barcode %>% unique() %>% length(),
                    filter(linCountToOverlaps, SampleNum == "S2")$barcode %>% unique() %>% length(),
                    filter(linCountToOverlaps, SampleNum == "S3")$barcode %>% unique() %>% length())

for(i in 1:length(cutoffList)){
  primedAllUMAP <- filter(umapCoordinates, cellID %in% c(unlist(primedCellIDList[[i]]), unlist(primedCellIDList[[i+length(cutoffList)]]), unlist(primedCellIDList[[i+2*length(cutoffList)]]))) %>%
    dplyr::filter(sampleNum %in% c("S1", "S2", "S3"))
  primedAllUMAP <- inner_join(primedAllUMAP, linCountToOverlaps, by = "cellID")
  
  primedAllUMAPSummary <- primedAllUMAP %>% group_by(sampleNum) %>% summarise(count = length(unique(barcode))) %>%
    mutate(countNorm = NA)
  primedAllUMAPSummary$countNorm <- primedAllUMAPSummary$count/cellNumberList
  primedAllUMAPSummary$sampleNum <- factor(primedAllUMAPSummary$sampleNum, levels = c("S3", "S1", "S2"), labels = c("slow", "control", "fast"))
  primedAllUMAPSummary$countNormNorm <- primedAllUMAPSummary$countNorm/(primedAllUMAPSummary %>% dplyr::filter(sampleNum == "control") %>% .$countNorm)
  
  plot <- ggplot(primedAllUMAPSummary %>% dplyr::filter(sampleNum %in% c("slow", "fast")), aes(x = sampleNum, y = countNormNorm)) +
    geom_col() + geom_hline(yintercept = 1, linetype = "dashed")
  ggsave(plot = plot, file = paste0(plotDirectory, 'R2_Fast(S1)VsSlow(S3)', '_cutoff_', cutoffList[i], '.pdf'), width = 8, height = 8)
}

#### analysis of purity across samples for different lineages ####
#######################################################################################################################################################
seuratLabelPriming <- left_join(umapClusters, primedAllUMAP, by = c("cellID", "sampleNum")) %>% mutate(label = "none")
seuratLabelPriming$label <- ifelse(seuratLabelPriming$cellID %in% linCountToOverlaps$cellID, "nonprimed", seuratLabelPriming$label)
seuratLabelPriming$label <- ifelse(seuratLabelPriming$cellID %in% primedAllUMAP$cellID, "primed", seuratLabelPriming$label)
seuratLabelPriming$cluster <- factor(seuratLabelPriming$cluster, levels = seuratLabelPriming$cluster %>% unique() %>% sort())

seuratLabelPrimingBC <- inner_join(seuratLabelPriming, linCountToOverlaps, by = "cellID") %>% mutate(label = "barcoded")
seuratLabelPrimingBC$label <- ifelse(seuratLabelPrimingBC$barcode %in% primedAllBarcodes, "primed", seuratLabelPrimingBC$label)

cellsS1 <- seuratLabelPrimingBC %>% filter(sampleNum == "S1") %>% group_by(barcode, label) %>% summarize(number = n())
cellsS2 <- seuratLabelPrimingBC %>% filter(sampleNum == "S2") %>% group_by(barcode, label) %>% summarize(number = n())
cellsS3 <- seuratLabelPrimingBC %>% filter(sampleNum == "S3") %>% group_by(barcode, label) %>% summarize(number = n())

overlapS1S2 <- inner_join(cellsS1, cellsS2, by = c("barcode", "label"))
overlapS1S3 <- inner_join(cellsS1, cellsS3, by = c("barcode", "label"))
overlapS2S3 <- inner_join(cellsS2, cellsS3, by = c("barcode", "label"))
overlapAll <- inner_join(overlapS1S2, cellsS3, by = c("barcode", "label"))

library(eulerr)
venn <- euler(c(A = nrow(cellsS1),
                B = nrow(cellsS2),
                C = nrow(cellsS3),
                "A&B" = nrow(overlapS1S2),
                "A&C" = nrow(overlapS1S3),
                "B&C" = nrow(overlapS2S3),
                "A&B&C" = nrow(overlapAll)),
              shape = "ellipse")
plot(venn, labels = c("control", "fast", "slow"), quantities = list(type = "percent"), type = "percent")

cellsS1 <- seuratLabelPrimingBC %>% filter(sampleNum == "S1") %>% group_by(barcode, label) %>% summarize(number = n()) %>% filter(label == "primed")
cellsS2 <- seuratLabelPrimingBC %>% filter(sampleNum == "S2") %>% group_by(barcode, label) %>% summarize(number = n()) %>% filter(label == "primed")
cellsS3 <- seuratLabelPrimingBC %>% filter(sampleNum == "S3") %>% group_by(barcode, label) %>% summarize(number = n()) %>% filter(label == "primed")

overlapS1S2 <- inner_join(cellsS1, cellsS2, by = c("barcode", "label"))
overlapS1S3 <- inner_join(cellsS1, cellsS3, by = c("barcode", "label"))
overlapS2S3 <- inner_join(cellsS2, cellsS3, by = c("barcode", "label"))
overlapAll <- inner_join(overlapS1S2, cellsS3, by = c("barcode", "label"))

library(eulerr)
venn <- euler(c(A = nrow(cellsS1),
                B = nrow(cellsS2),
                C = nrow(cellsS3),
                "A&B" = nrow(overlapS1S2),
                "A&C" = nrow(overlapS1S3),
                "B&C" = nrow(overlapS2S3),
                "A&B&C" = nrow(overlapAll)),
              shape = "ellipse")
plot(venn, labels = c("control", "fast", "slow"), quantities = list(type = "percent"), type = "percent")

overlapPrimedCells <- linCountToOverlaps %>% filter(barcode %in% bind_rows(overlapAll)$barcode)
overlapPrimedCells <- linCountToOverlaps %>% filter(barcode %in% bind_rows(overlapAll, overlapS1S2, overlapS1S3, overlapS2S3)$barcode)
overlapPrimedCellsUMAP <- overlapPrimedCells %>% left_join(umapCoordClust, by = "cellID")

ggplot() +
  geom_point(data = umapCoordClust, aes(x = UMAP_1, y = UMAP_2), color = "lightgray") +
  geom_point(data = overlapPrimedCellsUMAP, aes(x = UMAP_1, y = UMAP_2), color = "red")

overlapPrimedCells <- linCountToOverlaps %>% filter(barcode %in% bind_rows(overlapAll, overlapS1S2, overlapS1S3, overlapS2S3)$barcode)
overlapPrimedCellsUMAP <- umapCoordClust %>% mutate(label = "none")
overlapPrimedCellsUMAP$label <- ifelse(overlapPrimedCellsUMAP$cellID %in% primedAllUMAP$cellID, "primed", overlapPrimedCellsUMAP$label)
overlapPrimedCellsUMAP$label <- ifelse(overlapPrimedCellsUMAP$cellID %in% overlapPrimedCells$cellID, "overlap", overlapPrimedCellsUMAP$label)

ggplot() +
  geom_point(data = umapCoordClust, aes(x = UMAP_1, y = UMAP_2), color = "lightgray") +
  geom_point(data = overlapPrimedCellsUMAP %>% filter(label != "none"), aes(x = UMAP_1, y = UMAP_2, color = label)) +
  facet_wrap(~sampleNum)

logNormCountsFilter <- readRDS(file = paste0(homeDirectory, "logNormalizedCountsFilter.rds"))
markerCountsOverlapPrimed <- inner_join(overlapPrimedCellsUMAP, logNormCountsFilter, by = "cellID") %>% filter(label != "none")

ggplot(markerCountsOverlapPrimed, aes(x = label, y = SPP1)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() + facet_wrap(~sampleNum.x) +
  geom_signif(comparisons = list(c("overlap", "primed")))

ggplot(markerCountsOverlapPrimed, aes(x = label, y = MKI67)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() + facet_wrap(~sampleNum.x) +
  geom_signif(comparisons = list(c("overlap", "primed")))

overlapPlotTable <- markerCountsOverlapPrimed %>% left_join(linCountToOverlaps, by = "cellID")
overlapPlot <- inner_join(overlapPlotTable %>% filter(label == "overlap", sampleNum.x == "S2"), overlapPlotTable %>% filter(label == "overlap", sampleNum.x == "S3"), by = "barcode")
ggplot(overlapPlot, aes(x = SPP1.x, y = SPP1.y)) +
  geom_point()
