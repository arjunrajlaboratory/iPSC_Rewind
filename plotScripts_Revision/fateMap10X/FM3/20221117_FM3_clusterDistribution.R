rm(list=ls())
gc()

library(Seurat)
library(tidyverse)
library(RANN)
library(tripack)
library(reshape2)
library(svglite)
library(biomaRt)
library(ggpubr)
library(ggrastr)
library(ggsignif)
library(ggrepel)

homeDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/extractedData/fateMap10X/FM3/"
plotDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/plots/fateMap10X/FM3/"

theme_set(theme_classic())

linCountToOverlapsList <- as_tibble(read.table(file = paste0(homeDirectory, "filtered10XCells.txt"), header = TRUE, stringsAsFactors = F, sep = "\t"))

#### UMAP plots for DMSO versus LSD1i ####
umapCoordinates = as_tibble(read.table(file = paste0(homeDirectory, "umapCoordinates_Scanorama_50pcs_filter.tsv"), header = TRUE, stringsAsFactors = F, sep = "\t"))
umapClusters = as_tibble(read.table(file = paste0(homeDirectory, "umapClusters_Scanorama_50pcs_filter.tsv"), header = TRUE, stringsAsFactors = F, sep = "\t"))

umapCoordClust <- inner_join(umapCoordinates, umapClusters, by = c("cellID", "sampleNum")) %>% dplyr::rename(cluster = scanorama_snn_res.0.3)
umapCoordClust$cluster <- factor(umapCoordClust$cluster, levels = c(0, 1, 2, 3, 4, 5))
rm(umapCoordinates)
rm(umapClusters)

umapCoordClustSubset <- umapCoordClust %>% dplyr::filter(cluster %in% c("0", "1", "2", "3", "4", "5"), UMAP_2 < 10)
umapCoordClustSubset$text <- ""
centroids <- aggregate(cbind(UMAP_1, UMAP_2) ~ cluster, umapCoordClustSubset, mean) %>% mutate(text = cluster)
umapCoordClustSubset <- bind_rows(umapCoordClustSubset, centroids)

umapCoordClustFilter <- umapCoordClust %>% dplyr::filter(cluster %in% c("0", "1", "2", "3", "4", "5"), UMAP_2 < 10)
umapCoordClustFilter$cluster <- factor(umapCoordClustFilter$cluster, labels = c("iPSC", "iPSC", "iPSC", "iPSC", "fibroblast", "incomplete"))
clusterTable_DMSO <- umapCoordClustFilter %>% dplyr::filter(sampleNum %in% c("S1", "S2")) %>%
  group_by(cluster) %>% summarise(nCells = n()) %>% arrange(desc(cluster)) %>%
  mutate(prop = nCells/sum(nCells)) %>% mutate(ypos = cumsum(prop) - 0.5*prop)
ggplot(clusterTable_DMSO, aes(x = "", y = prop, fill = cluster)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  coord_polar("y", start = 0) +
  theme_void() + theme(legend.position = "none") +
  geom_text(aes(y = ypos, label = paste0(cluster, "\n", round(prop, 2)*100, "%")), color = "black", size = 6)
ggsave(filename = paste0(plotDirectory, "FM3_cellTypePieChart_DMSO.pdf"), height = 3, width = 3)

clusterTable_LSD1 <- umapCoordClustFilter %>% dplyr::filter(sampleNum %in% c("S3", "S4")) %>%
  group_by(cluster) %>% summarise(nCells = n()) %>% arrange(desc(cluster)) %>%
  mutate(prop = nCells/sum(nCells)) %>% mutate(ypos = cumsum(prop) - 0.5*prop)
ggplot(clusterTable_LSD1, aes(x = "", y = prop, fill = cluster)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  coord_polar("y", start = 0) +
  theme_void() + theme(legend.position = "none") +
  geom_text(aes(y = ypos, label = paste0(cluster, "\n", round(prop, 2)*100, "%")), color = "black", size = 6)
ggsave(filename = paste0(plotDirectory, "FM3_cellTypePieChart_LSD1.pdf"), height = 3, width = 3)

clusterTable_DOT1L <- umapCoordClustFilter %>% dplyr::filter(sampleNum %in% c("S5", "S6")) %>%
  group_by(cluster) %>% summarise(nCells = n()) %>% arrange(desc(cluster)) %>%
  mutate(prop = nCells/sum(nCells)) %>% mutate(ypos = cumsum(prop) - 0.5*prop)
ggplot(clusterTable_DOT1L, aes(x = "", y = prop, fill = cluster)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  coord_polar("y", start = 0) +
  theme_void() + theme(legend.position = "none") +
  geom_text(aes(y = ypos, label = paste0(cluster, "\n", round(prop, 2)*100, "%")), color = "black", size = 6)
ggsave(filename = paste0(plotDirectory, "FM3_cellTypePieChart_DOT1L.pdf"), height = 3, width = 3)

ggplot(umapCoordClustSubset, aes(x = UMAP_1, y = UMAP_2, color = cluster, label = text)) +
  rasterize(geom_point(data = umapCoordClustSubset, aes(x = UMAP_1, y = UMAP_2), color = "lightgray"), dpi = 100) +
  rasterize(geom_point(), dpi = 100) +
  geom_text(color = "black", size = 7.5) +
  NoAxes() + NoLegend()
ggsave(filename = paste0(plotDirectory, "umapAll.pdf"), unit = "in", height = 4, width = 4)

ggplot(umapCoordClustSubset %>% dplyr::filter(sampleNum %in% c("S1", "S2", NA)), aes(x = UMAP_1, y = UMAP_2, color = cluster, label = text)) +
  rasterize(geom_point(data = umapCoordClustSubset, aes(x = UMAP_1, y = UMAP_2), color = "lightgray"), dpi = 100) +
  rasterize(geom_point(), dpi = 100) +
  geom_text(color = "black", size = 7.5) +
  NoAxes() + NoLegend()
ggsave(filename = paste0(plotDirectory, "umapDMSOAll.pdf"), unit = "in", height = 4, width = 4)

ggplot(umapCoordClustSubset %>% dplyr::filter(sampleNum %in% c("S3", "S4", NA)), aes(x = UMAP_1, y = UMAP_2, color = cluster, label = text)) +
  rasterize(geom_point(data = umapCoordClustSubset, aes(x = UMAP_1, y = UMAP_2), color = "lightgray"), dpi = 100) +
  rasterize(geom_point(), dpi = 100) +
  geom_text(color = "black", size = 7.5) +
  NoAxes() + NoLegend()
ggsave(filename = paste0(plotDirectory, "umapLSD1All.pdf"), unit = "in", height = 4, width = 4)

#### UMAP plots for twins across DMSO and LSD1i ####
allUMAP = inner_join(linCountToOverlapsList, umapCoordClustSubset, by = c("cellID")) %>% dplyr::select(-nLineages)
jointUMAP = inner_join(linCountToOverlapsList %>% dplyr::filter(SampleNum %in% c("S1", "S2", "S5", "S6")), umapCoordClustSubset, by = c("cellID")) %>% dplyr::select(-nLineages)
jointUMAP1Barcodes = jointUMAP %>% filter(SampleNum=="S1" | SampleNum=="S2") %>% dplyr::select(barcode) %>% unique()
jointUMAP2Barcodes = jointUMAP %>% filter(SampleNum=="S5" | SampleNum=="S6") %>% dplyr::select(barcode) %>% unique()

jointBarcodesOnlyBoth = inner_join(jointUMAP1Barcodes, jointUMAP2Barcodes, by = "barcode")
jointBarcodesOnlyBoth = jointBarcodesOnlyBoth %>% mutate(barcodeName = c(1:nrow(jointBarcodesOnlyBoth)))
jointBarcodesOnlyBoth$barcodeName = sub("^", "B", jointBarcodesOnlyBoth$barcodeName)

jointBarcodesOnlyBothList = jointBarcodesOnlyBoth$barcode
jointUMAPOnlyBoth = jointUMAP %>% filter(barcode %in% jointBarcodesOnlyBothList)
finalUMAPJoint = inner_join(jointUMAPOnlyBoth, jointBarcodesOnlyBoth, by = "barcode") %>% dplyr::select(-barcode,-nUMI)

finalUMAPS1 = finalUMAPJoint %>% filter(SampleNum == "S1" |
                                        SampleNum == "S2")
finalUMAPS2 = finalUMAPJoint %>% filter(SampleNum == "S5" |
                                        SampleNum == "S6")

finalUMAPS1Big = finalUMAPS1 %>% group_by(barcodeName) %>% summarise(nColony = length(barcodeName)) %>% filter(nColony > 1) %>% dplyr::select(-nColony)
finalUMAPS2Big = finalUMAPS2 %>% group_by(barcodeName) %>% summarise(nColony = length(barcodeName)) %>% filter(nColony > 1) %>% dplyr::select(-nColony)

commonS1S2Big = inner_join(finalUMAPS1Big, finalUMAPS2Big)
commonS1S2Big = commonS1S2Big$barcodeName
finalUMAPJointBig = finalUMAPJoint %>% filter(barcodeName %in% commonS1S2Big)
finalUMAPJointBig = finalUMAPJoint %>% mutate(label = "DMSO")
finalUMAPJointBig$label = ifelse(finalUMAPJointBig$SampleNum == "S5" | finalUMAPJointBig$SampleNum == "S6", "LSD1", finalUMAPJointBig$label)

finalUMAPJointPlot <- bind_rows(finalUMAPJoint, centroids)
finalUMAPJointPlot$label <- "DMSO"
finalUMAPJointPlot$label <- ifelse(finalUMAPJointPlot$SampleNum %in% c("S5", "S6"), "LSD1", finalUMAPJointPlot$label)
ggplot(finalUMAPJointPlot, aes(x = UMAP_1, y = UMAP_2, color = label, label = text)) +
  rasterize(geom_point(data = umapCoordClustSubset, aes(x = UMAP_1, y = UMAP_2), color = "lightgray"), dpi = 100) +
  rasterize(geom_point(), dpi = 100) +
  geom_text(color = "black", size = 7.5) +
  scale_color_manual(values = c("#f39cc2", "#00add0")) +
  NoAxes() + NoLegend()
ggsave(filename = paste0(plotDirectory, "umapDMSOvsLSD1TwinsOnly.pdf"), unit = "in", height = 4, width = 4)

#### fraction of cells in most dominant cluster ####
jointUMAP$label <- "DMSO"
jointUMAP$label <- ifelse(jointUMAP$SampleNum %in% c("S5", "S6"), "LSD1", jointUMAP$label)
jointBarcodesAll <- jointUMAP %>% .$barcode %>% unique()
umapClusters <- jointUMAP %>% dplyr::select(cluster) %>% unique() #### use for all barcodes in both
umapClustersTable <- tibble(cluster = rep(umapClusters$cluster, 2), label = c(rep("DMSO", nrow(umapClusters)), rep("LSD1", nrow(umapClusters))))
for(i in 1:length(jointBarcodesAll)) {
  umapBarcodeTableTemp <- jointUMAP %>% dplyr::filter(barcode == jointBarcodesAll[i])
  umapBarcodeTableDistTemp <- umapBarcodeTableTemp %>% group_by(cluster, label) %>% summarise(number = n()) %>% ungroup()
  umapBarcodeTableDistTemp <- umapBarcodeTableDistTemp %>% group_by(label) %>% mutate(fraction = number/sum(number)) %>% ungroup()
  umapBarcodeTableDistTemp <- full_join(umapClustersTable, umapBarcodeTableDistTemp, by = c("cluster", "label")) %>% replace(is.na(.), 0)
  umapBarcodeTableDistTemp <- umapBarcodeTableDistTemp %>% mutate(barcode = paste0("B", i))
  
  umapRandomTable <- sample_n(allUMAP, size = max(c(umapBarcodeTableDistTemp %>% dplyr::filter(label == "DMSO") %>% .$number %>% sum(), 
                                                umapBarcodeTableDistTemp %>% dplyr::filter(label == "LSD1") %>% .$number %>% sum())))
  umapRandomTableDistTemp <- umapRandomTable %>% group_by(cluster) %>% summarise(number = n()) %>% mutate(fraction = n()/sum(number)) %>% ungroup()
  umapRandomTableDistTemp <- full_join(umapClustersTable, umapRandomTableDistTemp, by = "cluster") %>% replace(is.na(.), 0) %>% mutate(label = "random")
  umapRandomTableDistTemp <- umapRandomTableDistTemp %>% mutate(barcode = paste0("B", i))
  if(i == 1) {
    umapBarcodeTableDist <- umapBarcodeTableDistTemp
    domClustDMSO <- dplyr::filter(umapBarcodeTableDistTemp, label == "DMSO") %>% slice(which.max(fraction))
    domClustLSD1 <- dplyr::filter(umapBarcodeTableDistTemp, label == "LSD1") %>% slice(which.max(fraction))
    domClustRandom <- umapRandomTableDistTemp %>% slice(which.max(fraction))
  } else {
    umapBarcodeTableDist <- bind_rows(umapBarcodeTableDist, umapBarcodeTableDistTemp)
    domClustDMSO <- bind_rows(domClustDMSO, dplyr::filter(umapBarcodeTableDistTemp, label == "DMSO") %>% slice(which.max(fraction)))
    domClustLSD1 <- bind_rows(domClustLSD1, dplyr::filter(umapBarcodeTableDistTemp, label == "LSD1") %>% slice(which.max(fraction)))
    domClustRandom <- bind_rows(domClustRandom, umapRandomTableDistTemp %>% slice(which.max(fraction)))
  }
}

domClustDMSO <- domClustDMSO %>% dplyr::filter(fraction > 0) %>% dplyr::filter(number > 1)
domClustLSD1 <- domClustLSD1 %>% dplyr::filter(fraction > 0) %>% dplyr::filter(number > 1)
domClustRandom <- domClustRandom %>% dplyr::filter(fraction > 0) %>% dplyr::filter(number > 1)
domClustTable <- bind_rows(domClustDMSO, domClustLSD1, domClustRandom)
domClustTable$label <- factor(domClustTable$label, levels = c("random", "DMSO", "LSD1"))

ggplot(domClustTable, aes(x = cluster, group = label)) +
  facet_wrap(~label) +
  geom_bar(aes(y = ..prop..)) +
  ylab("fraction of barcodes with dominant cluster x") + ylim(0, 1) +
  scale_x_discrete(drop = FALSE)

ggplot(domClustTable, aes(x = cluster, y = fraction)) +
  facet_wrap(~label) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.25) +
  stat_summary(fun = mean, geom = "point", size = 5) +
  ylab("fraction of barcodes with dominant cluster x") + ylim(0, 1) +
  scale_x_discrete(drop = FALSE)

ggplot(domClustTable, aes(x = label, y = fraction)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.25) +
  stat_summary(fun = mean, geom = "point", size = 5) +
  stat_compare_means(comparisons = list(c("DMSO", "random"), c("LSD1", "random"), c("DMSO", "LSD1")), method = "t.test") +
  ylab("fraction of cells in dominant cluster") +
  theme(axis.title.x = element_blank())


umapClusters <- jointUMAP %>% dplyr::select(cluster) %>% unique() #### use for barcodes only in both
umapClustersTable <- tibble(cluster = rep(umapClusters$cluster, 2), label = c(rep("DMSO", nrow(umapClusters)), rep("LSD1", nrow(umapClusters))))
for(i in 1:length(commonS1S2Big)) {
  umapBarcodeTableTemp <- finalUMAPJointBig %>% dplyr::filter(barcodeName == commonS1S2Big[i])
  umapBarcodeTableDistTemp <- umapBarcodeTableTemp %>% group_by(cluster, label) %>% summarise(number = n()) %>% ungroup()
  umapBarcodeTableDistTemp <- umapBarcodeTableDistTemp %>% group_by(label) %>% mutate(fraction = number/sum(number)) %>% ungroup()
  umapBarcodeTableDistTemp <- full_join(umapClustersTable, umapBarcodeTableDistTemp, by = c("cluster", "label")) %>% replace(is.na(.), 0)
  umapBarcodeTableDistTemp <- umapBarcodeTableDistTemp %>% mutate(barcode = paste0("B", i))
  
  umapRandomTable <- sample_n(allUMAP, size = max(c(umapBarcodeTableDistTemp %>% dplyr::filter(label == "DMSO") %>% .$number %>% sum(), 
                                                    umapBarcodeTableDistTemp %>% dplyr::filter(label == "LSD1") %>% .$number %>% sum())))
  umapRandomTableDistTemp <- umapRandomTable %>% group_by(cluster) %>% summarise(number = n()) %>% mutate(fraction = n()/sum(number)) %>% ungroup()
  umapRandomTableDistTemp <- full_join(umapClustersTable, umapRandomTableDistTemp, by = "cluster") %>% replace(is.na(.), 0) %>% mutate(label = "random")
  umapRandomTableDistTemp <- umapRandomTableDistTemp %>% mutate(barcode = paste0("B", i))
    if(i == 1) {
      umapBarcodeTableDist <- umapBarcodeTableDistTemp
      domClustDMSO <- dplyr::filter(umapBarcodeTableDistTemp, label == "DMSO") %>% slice(which.max(fraction))
      domClustLSD1 <- dplyr::filter(umapBarcodeTableDistTemp, label == "LSD1") %>% slice(which.max(fraction))
      domClustRandom <- umapRandomTableDistTemp %>% slice(which.max(fraction))
    } else {
      umapBarcodeTableDist <- bind_rows(umapBarcodeTableDist, umapBarcodeTableDistTemp)
      domClustDMSO <- bind_rows(domClustDMSO, dplyr::filter(umapBarcodeTableDistTemp, label == "DMSO") %>% slice(which.max(fraction)))
      domClustLSD1 <- bind_rows(domClustLSD1, dplyr::filter(umapBarcodeTableDistTemp, label == "LSD1") %>% slice(which.max(fraction)))
      domClustRandom <- bind_rows(domClustRandom, umapRandomTableDistTemp %>% slice(which.max(fraction)))
    }
}

domClustTable <- bind_rows(domClustDMSO, domClustLSD1, domClustRandom)
domClustTable$label <- factor(domClustTable$label, levels = c("random", "DMSO", "LSD1"))

ggplot(domClustTable %>% dplyr::filter(label != "random"), aes(x = cluster, group = label)) +
  facet_wrap(~label) +
  geom_bar(aes(y = ..prop..)) +
  ylab("fraction of barcodes with dominant cluster x") + ylim(0, 1) +
  scale_x_discrete(drop = FALSE)
ggsave(filename = paste0(plotDirectory, "domClustDistDMSOvsLSD1.pdf"), unit = "in", height = 2, width = 4)

# ggplot(domClustTable, aes(x = cluster, y = fraction)) +
#   facet_wrap(~label) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_jitter(width = 0.25) +
#   stat_summary(fun = mean, geom = "point", size = 5) +
#   ylab("fraction of barcodes with dominant cluster x") + ylim(0, 1) +
#   scale_x_discrete(drop = FALSE)

ggplot(domClustTable, aes(x = label, y = fraction)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.25) +
  stat_summary(fun = mean, geom = "point", size = 5) +
  stat_compare_means(comparisons = list(c("DMSO", "random"), c("LSD1", "random"), c("DMSO", "LSD1")), method = "t.test", paired = TRUE) +
  ylab("fraction of cells in dominant cluster") +
  theme(axis.title.x = element_blank())
ggsave(filename = paste0(plotDirectory, "fracCellsDomClustRandomvsDMSOvsLSD1.pdf"), unit = "in", height = 4, width = 4)

domClustSwitchPerBarcode <- domClustTable %>% dcast(formula = barcode ~ label, value.var = "cluster")
domClustSwitchPerBarcode <- domClustSwitchPerBarcode %>% group_by(DMSO, LSD1) %>% summarise(number = n())
domClustSwitchPerBarcode <- domClustSwitchPerBarcode %>% rowwise() %>% mutate(frac = number/sum(domClustSwitchPerBarcode$number))

domClustSwitchPerBarcodePlot <- domClustSwitchPerBarcode %>%
  inner_join(., centroids %>% dplyr::select(-text), by = c("DMSO" = "cluster")) %>% dplyr::rename(x = UMAP_1, y = UMAP_2) %>%
  inner_join(., centroids %>% dplyr::select(-text), by = c("LSD1" = "cluster")) %>% dplyr::rename(xend = UMAP_1, yend = UMAP_2)
pathSize = 1
ggplot(centroids, aes(x = UMAP_1, y = UMAP_2)) +
  rasterize(geom_point(data = umapCoordClustSubset, aes(x = UMAP_1, y = UMAP_2), color = "lightgray"), dpi = 100) +
  rasterize(geom_point(data = finalUMAPJointPlot, aes(x = UMAP_1, y = UMAP_2, color = label)), dpi = 100) +
  geom_curve(data = domClustSwitchPerBarcodePlot[1, ], aes(x = x, y = y, xend = xend, yend = yend, size = pathSize*frac)) +
  geom_curve(data = domClustSwitchPerBarcodePlot[2, ], aes(x = x, y = y, xend = xend, yend = yend, size = pathSize*frac)) +
  geom_curve(data = domClustSwitchPerBarcodePlot[3, ], aes(x = x, y = y, xend = xend, yend = yend, size = pathSize*frac)) +
  geom_curve(data = domClustSwitchPerBarcodePlot[4, ], aes(x = x, y = y, xend = xend, yend = yend, size = pathSize*frac)) +
  geom_curve(data = domClustSwitchPerBarcodePlot[6, ], aes(x = x, y = y, xend = xend, yend = yend, size = pathSize*frac)) +
  geom_curve(data = domClustSwitchPerBarcodePlot[7, ], aes(x = x, y = y, xend = xend, yend = yend, size = pathSize*frac)) +
  scale_color_manual(values = c("#f39cc2", "#00add0")) +
  geom_text(aes(label = text), color = "black", size = 7.5) +
  NoAxes() + NoLegend()
ggsave(filename = paste0(plotDirectory, "domClustSwitchPathDMSOvsLSD1.pdf"), unit = "in", height = 4, width = 4)

####normClusterTables
umapClusters <- umapClusters %>% rename(clusters = scanorama_snn_res.0.8)

barcodedCellIDsAll <- bind_rows(linCountToOverlapsListCorrected[[1]],
                                linCountToOverlapsListCorrected[[2]],
                                linCountToOverlapsListCorrected[[3]],
                                linCountToOverlapsListCorrected[[4]])$cellID
barcodedCellClustersAll <- umapClusters %>% filter(cellID %in% barcodedCellIDsAll)
clusterTableAll <- barcodedCellClustersAll %>% group_by(clusters) %>% summarise(cellsperclusterall = n()) %>% ungroup()

barcodedCellIDsDMSO <- bind_rows(linCountToOverlapsListCorrected[[1]],
                                 linCountToOverlapsListCorrected[[2]])$cellID
barcodedCellClustersDMSO <- umapClusters %>% filter(cellID %in% barcodedCellIDsDMSO)
clusterTableDMSO <- barcodedCellClustersDMSO %>% group_by(clusters) %>% summarise(cellsperclusterDMSO = n()) %>% ungroup()

barcodedCellIDsLSD1 <- bind_rows(linCountToOverlapsListCorrected[[3]],
                                 linCountToOverlapsListCorrected[[4]])$cellID
barcodedCellClustersLSD1 <- umapClusters %>% filter(cellID %in% barcodedCellIDsLSD1)
clusterTableLSD1 <- barcodedCellClustersLSD1 %>% group_by(clusters) %>% summarise(cellsperclusterLSD1 = n()) %>% ungroup()

clusterTableRandomValueDMSOList <- list()
for(i in 1:100) {
  clusterTableRandomTemp <- barcodedCellClustersAll %>% sample_n(size = nrow(barcodedCellClustersDMSO)) %>% group_by(clusters) %>% summarise(cellsperclusterrandomDMSO = n()) %>% ungroup()
  clusterTableRandomTemp <- left_join(clusterTableAll, clusterTableRandomTemp, by = "clusters")
  clusterTableRandomTemp[is.na(clusterTableRandomTemp)] <- 0
  clusterTableRandomValueDMSOList [[i]] <- clusterTableRandomTemp$cellsperclusterrandomDMSO
}
clusterMatrixRandomDMSO <- sapply(clusterTableRandomValueDMSOList, unlist)

clusterTableRandomValueLSD1List <- list()
for(i in 1:100) {
  clusterTableRandomTemp <- barcodedCellClustersAll %>% sample_n(size = nrow(barcodedCellClustersLSD1)) %>% group_by(clusters) %>% summarise(cellsperclusterrandomLSD1 = n()) %>% ungroup()
  clusterTableRandomTemp <- left_join(clusterTableAll, clusterTableRandomTemp, by = "clusters")
  clusterTableRandomTemp[is.na(clusterTableRandomTemp)] <- 0
  clusterTableRandomValueLSD1List [[i]] <- clusterTableRandomTemp$cellsperclusterrandomLSD1
}
clusterMatrixRandomLSD1 <- sapply(clusterTableRandomValueLSD1List, unlist)

clusterTableMerge <- left_join(clusterTableAll, clusterTableDMSO, by = "clusters") %>% left_join(., clusterTableLSD1, by = "clusters") %>%
  mutate(cellsperclusterrandomDMSO = rowMeans(sapply(clusterTableRandomValueDMSOList, unlist))) %>%
  mutate(cellsperclusterrandomLSD1 = rowMeans(sapply(clusterTableRandomValueLSD1List, unlist)))
clusterTableMerge[is.na(clusterTableMerge)] <- 0
clusterTableMerge <- clusterTableMerge %>% rowwise() %>%
  mutate(propPerClustDMSO = cellsperclusterDMSO/cellsperclusterall) %>% ungroup() %>%
  mutate(propPerBarClustNormDMSO = propPerClustDMSO/sum(propPerClustDMSO)) %>%
  mutate(propPerClustRandomDMSO = cellsperclusterrandomDMSO/cellsperclusterall) %>% ungroup() %>%
  mutate(propPerBarClustRandomNormDMSO = propPerClustRandomDMSO/sum(propPerClustRandomDMSO)) %>%
  mutate(propPerClustLSD1 = cellsperclusterLSD1/cellsperclusterall) %>% ungroup() %>%
  mutate(propPerBarClustNormLSD1 = propPerClustLSD1/sum(propPerClustLSD1)) %>%
  mutate(propPerClustRandomLSD1 = cellsperclusterrandomLSD1/cellsperclusterall) %>% ungroup() %>%
  mutate(propPerBarClustRandomNormLSD1 = propPerClustRandomLSD1/sum(propPerClustRandomLSD1))

clusterPValuesDMSO <- data.frame(cluster = c(0:(nrow(clusterMatrixRandomDMSO)-1)), pvalue = 0)
for(i in 1:nrow(clusterMatrixRandom)) {
  clusterPValuesDMSO$pvalue[i] <- sum(clusterMatrixRandomDMSO[i, ] >= clusterTableMerge$cellsperclusterDMSO[i])/100
}

theme_set(theme_classic())

ggplot(clusterTableMerge, aes(x = clusters, y = propPerBarClustNormDMSO)) +
  geom_col(fill = "red") + ylim(0, 0.2) + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  geom_hline(yintercept = 0.08, linetype = "dashed")
ggsave(filename = paste0(plotDirectory, "clusterPlotDMSO.pdf"), height = 1, width = 4)

ggplot(clusterTableMerge, aes(x = clusters, y = propPerBarClustNormLSD1)) +
  geom_col(fill = "red") + ylim(0, 0.2) + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  geom_hline(yintercept = 0.08, linetype = "dashed")
ggsave(filename = paste0(plotDirectory, "clusterPlotLSD1.pdf"), height = 1, width = 4)

plot2 <- ggplot(clusterTableMerge, aes(x = cluster, y = propPerBarClustRandomNorm)) +
  geom_col() + ylim(0, 0.3) + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  geom_hline(yintercept = mean(clusterTableMerge$propPerBarClustRandomNorm), linetype = "dashed")
plot3 <- ggarrange(plot1, plot2, ncol = 1, nrow = 2)
ggsave(plot3, filename = "/Users/naveenjain/Downloads/Rplot02.svg", width = 4, height = 2)

##### analyze clusters
DefaultAssay(scanorama_filter) <- "RNA"
scanorama_filter <- NormalizeData(scanorama_filter)
clusterMarkers <- FindAllMarkers(scanorama_filter, only.pos = TRUE, logfc.threshold = 0.5)
saveRDS(clusterMarkers, paste0(homeDirectory, "clusterMarkers_log2FCOver50.rds"))
clusterMarkers <- readRDS(paste0(homeDirectory, "clusterMarkers_log2FCOver50.rds"))

markersPluri <- c("NANOG", "PODXL", "ALPL", "TDGF1", "POU5F1", "SOX2")
markersFibro <- c("SNAI2", "LUM", "COL1A1", "FOSL1")

clusterMarkers <- clusterMarkers

ggplot() +
  geom_jitter(data = clusterMarkers %>% filter(!(gene %in% c(markersPluri, markersFibro))), aes(x = cluster, y = avg_log2FC, color = cluster)) +
  geom_jitter(data = clusterMarkers %>% filter(gene %in% c(markersPluri)), aes(x = cluster, y = avg_log2FC), color = "red") +
  geom_text_repel(data = clusterMarkers %>% filter(gene %in% c(markersPluri)), aes(x = cluster, y = avg_log2FC, label = gene), color = "red") +
  geom_jitter(data = clusterMarkers %>% filter(gene %in% c(markersFibro)), aes(x = cluster, y = avg_log2FC), color = "blue") +
  geom_text_repel(data = clusterMarkers %>% filter(gene %in% c(markersFibro)), aes(x = cluster, y = avg_log2FC, label = gene), color = "blue") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_classic(base_size = 20) + theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + NoLegend()
ggsave(filename = paste0(plotDirectory, "clusterLog2FCGenes.pdf"), height = 2, width = 8)

##### comparing DMSO and LSD1i conditions

twinOverlapDMSOBarcodes <- union(linCountToOverlapsListCorrected[[1]]$barcode, linCountToOverlapsListCorrected[[2]]$barcode) %>% unique()
twinOverlapDMSO <- bind_rows(linCountToOverlapsListCorrected[[1]], linCountToOverlapsListCorrected[[2]]) %>% filter(barcode %in% twinOverlapDMSOBarcodes)

twinOverlapLSD1Barcodes <- union(linCountToOverlapsListCorrected[[3]]$barcode, linCountToOverlapsListCorrected[[4]]$barcode) %>% unique()
twinOverlapLSD1 <- bind_rows(linCountToOverlapsListCorrected[[3]], linCountToOverlapsListCorrected[[4]]) %>% filter(barcode %in% twinOverlapLSD1Barcodes)

twinOverlapDMSOLSD1Barcodes <- intersect(twinOverlapDMSOBarcodes, twinOverlapLSD1Barcodes)
twinOverlapDMSOLSD1 <- bind_rows(twinOverlapDMSO, twinOverlapLSD1) %>% filter(barcode %in% twinOverlapDMSOLSD1Barcodes) %>% inner_join(umapCoordinates, by = 'cellID')

ggplot() +
  geom_point(data = umapCoordinates, aes(UMAP_1, UMAP_2), color = "gray93") +
  geom_point(data = filter(twinOverlapDMSOLSD1, SampleNum == "S1"), aes(UMAP_1, UMAP_2), color = "dodgerblue4") +
  geom_point(data = filter(twinOverlapDMSOLSD1, SampleNum == "S2"), aes(UMAP_1, UMAP_2), color = "dodgerblue1") +
  geom_point(data = filter(twinOverlapDMSOLSD1, SampleNum == "S5"), aes(UMAP_1, UMAP_2), color = "hotpink4") +
  geom_point(data = filter(twinOverlapDMSOLSD1, SampleNum == "S6"), aes(UMAP_1, UMAP_2), color = "hotpink") +
  theme_classic() + theme(legend.position = "none") + NoLegend() + NoAxes() +
  facet_wrap(~barcode)


####################20210318
###merging cellIDs and Barcodes
cellIDs = logNormalizedCountsS1$cellID
rowstokeep = which(cellIDs %in% unlist(linCountTooverlaps[,1])) #comparing cells present in both BCSeq and Seurat
logNormalizedCountsSubset = logNormalizedCountsS1[rowstokeep,]
logNormalizedCountsSubsetWBarcodes = inner_join(logNormalizedCountsSubset,linCountTooverlaps, by = c("cellID", "sampleNum"))


cutoff_ratio = 0.6
upperCutoff = 200
lowerCutoff = 50

ggplot(logNormalizedCountsSubsetWBarcodes, aes(x=SOX10)) +
  geom_density(alpha=0.5)

###For Figure 1, paper
#geneSet = c("ACTA2", "ACTG2", "MYOCD","TAGLN")
geneSet = c("ACTA2", "ACTG2", "MYOCD") 
geneSet = c("ACTA2", "MYOCD", "NGFR", "S100B", "MLANA", "SOX10") 
ColoniesSample1 = logNormalizedCountsSubsetWBarcodes %>% filter (sampleNum == "S1") %>% group_by(BC50StarcodeD8) %>% select(BC50StarcodeD8, ACTA2) %>% summarise(colonySize = length(BC50StarcodeD8), sum = sum(ACTA2 > 0.5)) %>% mutate(ratio = sum / colonySize) %>% filter(ratio > cutoff_ratio, colonySize < upperCutoff & colonySize > lowerCutoff) %>% mutate(gene = 'ACTA2')
BCsIndicesToKeep = which(logNormalizedCountsSubsetWBarcodes$BC50StarcodeD8 %in% unlist(ColoniesSample1[,1]))
BCsToKeep = logNormalizedCountsSubsetWBarcodes[BCsIndicesToKeep,]  %>% select(geneSet, cellID)
colonyS1Coordinates = inner_join(BCsToKeep, umapCoordinates, by = c("cellID"))

allCells = logNormalizedCountsSubsetWBarcodes %>% select(geneSet) %>% mutate(type = "all")
colony = colonyS1Coordinates %>% select(geneSet) %>% mutate(type = "l1")
densityPlot = bind_rows(allCells,colony)

meltData = melt(data = densityPlot, id.vars = c("type"), measure.vars = geneSet)
meltDataL1 = meltData

######UMAP PLOT
plot<- ggplot() +
  geom_point(data = umapCoordinatesS1, aes(x = UMAP_1, y = UMAP_2), color = "gray93", size = 1.5, shape = 16) +
  geom_point(data = colonyS1Coordinates, aes(x = UMAP_1, y = UMAP_2), color = "hotpink3", size = 1.5, shape = 16) +
  create_lpr_theme() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = rel(0.6)),
        legend.text = element_text(size = rel(0.6), angle = 30),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  labs(color =  "Scaled\n Log UMI counts") 
ggsave(plot, file = paste0(plotDirectory, 'FM01_ACTA2Colony.svg'), width = 6, height = 5.325)

linCountToOverlapsListCorrectedAll <- bind_rows(linCountToOverlapsListCorrected[[1]],
                                                linCountToOverlapsListCorrected[[2]],
                                                linCountToOverlapsListCorrected[[3]],
                                                linCountToOverlapsListCorrected[[4]],
                                                linCountToOverlapsListCorrected[[5]],
                                                linCountToOverlapsListCorrected[[6]],)

##### global cluster analysis for colonies
umapClusters <- as_tibble(read.table(file = paste0(homeDirectory, "umapClusters_Scanorama_50pcs_filter.tsv"), header = TRUE, stringsAsFactors = F, sep = "\t"))

conditionList <- c("DMSO A", "DMSO B", "LSD1i A", "LSD1i B", "DOT1Li A", "DOT1Li B", "all colonies")
plotList6 <- list()
plotList7 <- list()
for(i in 1:(length(linCountToOverlapsListCorrected)+1)){
  if(i != 7){
    umapClustersSample = umapClusters %>% filter(sampleNum == paste0("S", i))
    umapCoordinatesSample = umapCoordinates %>% filter(sampleNum == paste0("S", i))
    jointUMAPSample = inner_join(linCountToOverlapsListCorrected[[i]], umapCoordinatesSample, by = c("cellID")) %>% dplyr::select(-nLineages)
  } else{
    umapClustersSample = umapClusters
    umapCoordinatesSample = umapCoordinates
    jointUMAPSample = inner_join(linCountToOverlapsListCorrectedAll, umapCoordinatesSample, by = c("cellID")) %>% dplyr::select(-nLineages)
  }
  
  coloniesToAnalyze = jointUMAPSample %>% group_by(barcode) %>% summarise(nColony = length(barcode)) %>% filter(nColony >= 5)
  colonyUMAP = inner_join(coloniesToAnalyze, jointUMAPSample, by = c("barcode")) %>% select(barcode, cellID)
  colonyCluster = inner_join(colonyUMAP, umapClustersSample, by = c("cellID")) %>% select(-cellID,-sampleNum)
  
  samplingLengths = coloniesToAnalyze$nColony
  maxFraction = c()
  secondMaxFraction = c()
  
  for (j in c(1:length(samplingLengths))) {
    subSampled = sample_n(colonyCluster, samplingLengths[j])
    fraction = subSampled %>% group_by(scanorama_snn_res.0.8) %>% summarise(nFraction = length(scanorama_snn_res.0.8)/samplingLengths[j])
    maxFraction[j] = max(fraction$nFraction)
    secondMaxFraction[j] = max(fraction$nFraction[fraction$nFraction!=max(fraction$nFraction)]) 
  }
  
  secondMaxFraction[!is.finite(secondMaxFraction)] <- NA
  secondMaxFraction[which(is.na(secondMaxFraction))] = maxFraction[which(is.na(secondMaxFraction))]
  secondMaxFraction <- secondMaxFraction[!is.na(secondMaxFraction)]
  
  mean(maxFraction)
  mean(secondMaxFraction)
  
  finalFractionRandom = tibble(nfractionFinal = numeric()) %>% add_row(nfractionFinal = maxFraction) %>% mutate(type = "random")
  finalFractionRandomBoth = tibble(nfractionFinal = numeric()) %>% add_row(nfractionFinal = maxFraction + secondMaxFraction) %>% mutate(type = "random")
  
  fractionColonies = colonyCluster %>% group_by(barcode, scanorama_snn_res.0.8) %>% summarise(nCount = length(scanorama_snn_res.0.8))
  fractionColonies = fractionColonies %>% group_by(barcode) %>% summarise(nfraction = max(nCount), nfractionSecond = max(nCount[nCount!=max(nCount)]))
  fractionColonies$nfractionSecond[!is.finite(fractionColonies$nfractionSecond)] <- 0 #####correcting for second maximum which does not exist
  
  finalFractionColonies = inner_join(fractionColonies, coloniesToAnalyze, by = "barcode") %>% mutate(nfractionFinal = nfraction/nColony) %>% select(nfractionFinal) %>% mutate(type = "colonies")
  finalFractionColoniesBoth = inner_join(fractionColonies, coloniesToAnalyze, by = "barcode") %>% mutate(nfractionFinal = (nfraction+nfractionSecond)/nColony) %>% select(nfractionFinal) %>% mutate(type = "colonies")
  
  fractionFinalClusters = bind_rows(finalFractionColonies,finalFractionRandom)
  fractionFinalClustersBoth = bind_rows(finalFractionColoniesBoth,finalFractionRandomBoth)
  
  fractionFinalClusters$type <- factor(fractionFinalClusters$type, levels=c("random", "colonies"))
  fractionFinalClustersBoth$type <- factor(fractionFinalClustersBoth$type, levels=c("random", "colonies"))
  
  plotList6[[i]] <- ggplot(fractionFinalClusters, aes(x = type, y=nfractionFinal, fill = as.factor(type)), show.legend = FALSE) + 
    geom_boxplot() +
    scale_fill_manual(values=c("grey", "goldenrod3")) +
    stat_summary(fun=mean) +
    theme_classic((base_size = 20)) + theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks = element_line(colour = "black", size = 1.5), text=element_text(family="Helvetica")) + 
    ylim(0,1) + ggtitle(conditionList[[i]]) +
    if(i != 1) {theme(legend.position = "none", axis.title.x = element_blank()) + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank() )} else{theme(legend.position = "none", axis.title.x = element_blank())}
  
  plotList7[[i]] <- ggplot(fractionFinalClustersBoth, aes(x = type, y=nfractionFinal, fill = as.factor(type)), show.legend = FALSE) + 
    geom_boxplot() +
    scale_fill_manual(values=c("grey", "goldenrod3")) +
    stat_summary(fun=mean) +
    theme_classic((base_size = 20)) + theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks = element_line(colour = "black", size = 1.5), text=element_text(family="Helvetica")) + 
    ylim(0,1) + ggtitle(conditionList[[i]]) +
    if(i != 1) {theme(legend.position = "none", axis.title.x = element_blank()) + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank() )} else{theme(legend.position = "none", axis.title.x = element_blank())}
}
egg::ggarrange(plots = plotList6, nrow = 1, ncol = 7, align = 'v')
egg::ggarrange(plots = plotList7, nrow = 1, ncol = 7, align = 'v')

##### colony specific cluster analysis
umapClusters <- as_tibble(read.table(file = paste0(homeDirectory, "umapClusters_Scanorama_50pcs_filter.tsv"), header = TRUE, stringsAsFactors = F, sep = "\t"))

conditionList <- c("DMSO A", "DMSO B", "LSD1i A", "LSD1i B", "DOT1Li A", "DOT1Li B", "all colonies")
plotList8 <- list()
for(n in 1:(length(linCountToOverlapsListCorrected)+1)){
  if(n != 7){
    umapClustersSample = umapClusters %>% filter(sampleNum == paste0("S", n))
    umapCoordinatesSample = umapCoordinates %>% filter(sampleNum == paste0("S", n))
    jointUMAPSample = inner_join(linCountToOverlapsListCorrected[[n]], umapCoordinatesSample, by = c("cellID")) %>% dplyr::select(-nLineages)
  } else{
    umapClustersSample = umapClusters
    umapCoordinatesSample = umapCoordinates
    jointUMAPSample = inner_join(linCountToOverlapsListCorrectedAll, umapCoordinatesSample, by = c("cellID")) %>% dplyr::select(-nLineages)
  }
  
  coloniesToAnalyze = jointUMAPSample %>% group_by(barcode) %>% summarise(nColony = length(barcode)) %>% filter(nColony >= 5)
  colonyUMAP = inner_join(coloniesToAnalyze, jointUMAPSample, by = c("barcode")) %>% select(barcode, cellID)
  colonyCluster = inner_join(colonyUMAP, umapClustersSample, by = c("cellID")) %>% select(-cellID,-sampleNum)
  
  fractionColonies = colonyCluster %>% group_by(barcode, scanorama_snn_res.0.8) %>% summarise(nCount = length(scanorama_snn_res.0.8))
  fractionColonies <- fractionColonies %>% group_by(barcode) %>% filter(nCount == max(nCount))%>% slice(1)
  finalFractionColonies = inner_join(fractionColonies, coloniesToAnalyze, by = "barcode") %>% mutate(nfractionFinal = nCount/nColony) %>% select(nfractionFinal) %>% mutate(type = "colonies")
  
  fractionColoniesSeuratClusters = fractionColonies %>% select(barcode, scanorama_snn_res.0.8)
  fractionColoniesSeuratClustersColony = inner_join(fractionColoniesSeuratClusters, coloniesToAnalyze, by = "barcode")
  
  maxFraction = c()
  medianFraction = c()
  meanFraction = c()
  
  for (i in c(1:length(fractionColoniesSeuratClustersColony$barcode))) {
    for (j in c(1:50)) {
      subSampled = sample_n(colonyCluster, fractionColoniesSeuratClustersColony$nColony[i])
      fraction = subSampled %>% group_by(scanorama_snn_res.0.8) %>% summarise(nFraction = length(scanorama_snn_res.0.8)/fractionColoniesSeuratClustersColony$nColony[i])
      if (fractionColoniesSeuratClustersColony$scanorama_snn_res.0.8[i] %in% fraction$scanorama_snn_res.0.8) {
        fractionTemp = fraction %>% filter(scanorama_snn_res.0.8 == fractionColoniesSeuratClustersColony$scanorama_snn_res.0.8[i])
        maxFraction[j] = fractionTemp$nFraction
      } else {
        maxFraction[j] = 0
      }
    }
    meanFraction[i] = mean(maxFraction)
    medianFraction[i] = median(maxFraction) 
  }
  
  finalFractionRandom = tibble(nfractionFinal = numeric()) %>% add_row(nfractionFinal = meanFraction) %>% mutate(type = "random")
  fractionFinalClusters = bind_rows(finalFractionColonies, finalFractionRandom)
  fractionFinalClusters$type <- factor(fractionFinalClusters$type, levels=c("random", "colonies"))
  
  summaryfractionFinalClusters = fractionFinalClusters %>% group_by(type) %>% summarise(mean = mean(nfractionFinal))
  
  plotList8[[n]] = ggplot(fractionFinalClusters, aes(x = type, y=nfractionFinal, fill = as.factor(type)), show.legend = FALSE) + 
    geom_boxplot() +
    geom_jitter(width = 0.1, shape = 16) +
    scale_fill_manual(values=c("grey","goldenrod3")) +
    stat_summary(fun=mean) +
    theme_classic((base_size = 20)) + theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks = element_line(colour = "black", size = 1.5), text=element_text(family="Helvetica")) + 
    ylim(0,1) + ggtitle(conditionList[[n]]) +
    if(n != 1) {theme(legend.position = "none", axis.title.x = element_blank()) + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank() )} else{theme(legend.position = "none", axis.title.x = element_blank())}
}
egg::ggarrange(plots = plotList8, nrow = 1, ncol = 7, align = 'v')

######representative UMAPS##############
umapClustersS1UMAP = inner_join(umapClustersS1, umapCoordinatesS1, by = c("cellID", "sampleNum"))
finalFractionColoniesAll = inner_join(fractionColonies,coloniesToAnalyze, by = "BC50StarcodeD8") %>% mutate(nfractionFinal = nCount/nColony) %>% select(nCount,nfractionFinal,BC50StarcodeD8, seurat_clusters)

umapCoordinatesB_high = jointUMAP %>% filter(BC50StarcodeD8 == "ATAGATCGACAAGGAGATGTTCTTCGTGGACATCTAGATCAACATCAAGT")  
umapCoordinatesB_medium = jointUMAP %>% filter(BC50StarcodeD8 == "ATTCTAGTTGTAGTACTAGATGATCATCATGTTGTTCTTGGTCTTGTTCA")  
umapCoordinatesB_low = jointUMAP %>% filter(BC50StarcodeD8 == "ATTCCTGATGGTCATGGAGAACCTGCTGGTCTTCCACATGGACAAGGAGC")  

B_high =  finalFractionColoniesAll %>% filter(BC50StarcodeD8 == "ATAGATCGACAAGGAGATGTTCTTCGTGGACATCTAGATCAACATCAAGT") %>% arrange(-nCount)
B_medium = finalFractionColoniesAll %>% filter(BC50StarcodeD8 == "ATTCTAGTTGTAGTACTAGATGATCATCATGTTGTTCTTGGTCTTGTTCA") %>% arrange(-nCount)
B_low = finalFractionColoniesAll %>% filter(BC50StarcodeD8 == "ATTCCTGATGGTCATGGAGAACCTGCTGGTCTTCCACATGGACAAGGAGC") %>% arrange(-nCount)
toPlotBHigh = umapClustersS1UMAP %>% filter(seurat_clusters %in% B_high$seurat_clusters[1:2])
toPlotBMedium = umapClustersS1UMAP %>% filter(seurat_clusters %in% B_medium$seurat_clusters[1:2])
toPlotBLow = umapClustersS1UMAP %>% filter(seurat_clusters %in% B_low$seurat_clusters[1:2])

#toPlotBHigh = inner_join(B_high,jointUMAP, by = c("BC50StarcodeD8")) %>% select(-nUMI,sampleNum,BC50StarcodeD8) %>% filter(seurat_clusters %in% B_high$seurat_clusters[1:2])
#toPlotBMedium = inner_join(finalFractionColoniesAll,jointUMAP, by = c("BC50StarcodeD8")) %>% select(-nUMI,sampleNum,BC50StarcodeD8) %>% filter(seurat_clusters %in% B_medium$seurat_clusters[1:2])
#toPlotBLow = inner_join(B_low,jointUMAP, by = c("BC50StarcodeD8")) %>% select(-nUMI,sampleNum,BC50StarcodeD8) %>% filter(seurat_clusters %in% B_low$seurat_clusters[1:2])

plot<- ggplot() +
  geom_point(data = umapCoordinatesS1, aes(x = UMAP_1, y = UMAP_2), color = "gray93", size = 1.5, shape = 16) +
  geom_point(data = toPlotBMedium, aes(x = UMAP_1, y = UMAP_2, color = seurat_clusters), size = 1.5, shape = 16) +
  scale_color_manual(values=c("gray40", "gray60")) +
  geom_point(data = umapCoordinatesB_medium, aes(x = UMAP_1, y = UMAP_2),color = "hotpink3", size = 1.5, shape = 16) +
  create_lpr_theme() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = rel(0.6)),
        legend.text = element_text(size = rel(0.6), angle = 30),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  labs(color =  "Scaled\n Log UMI counts") 
ggsave(plot, file = paste0(plotDirectory, 'FM01_MediumClusterFraction.svg'), width = 6, height = 6.158)

plot<- ggplot() +
  geom_point(data = umapCoordinatesS1, aes(x = UMAP_1, y = UMAP_2), color = "gray93", size = 1.5, shape = 16) +
  geom_point(data = toPlotBHigh, aes(x = UMAP_1, y = UMAP_2, color = seurat_clusters), size = 1.5, shape = 16) +
  scale_color_manual(values=c("gray40", "gray60")) +
  geom_point(data = umapCoordinatesB_high, aes(x = UMAP_1, y = UMAP_2),color = "hotpink3", size = 1.5, shape = 16) +
  create_lpr_theme() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = rel(0.6)),
        legend.text = element_text(size = rel(0.6), angle = 30),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  labs(color =  "Scaled\n Log UMI counts") 
ggsave(plot, file = paste0(plotDirectory, 'FM01_HighClusterFraction.svg'), width = 6, height = 6.158)

plot<- ggplot() +
  geom_point(data = umapCoordinatesS1, aes(x = UMAP_1, y = UMAP_2), color = "gray93", size = 1.5, shape = 16) +
  geom_point(data = toPlotBLow, aes(x = UMAP_1, y = UMAP_2, color = seurat_clusters), size = 1.5, shape = 16) +
  scale_color_manual(values=c("gray40", "gray60")) +
  geom_point(data = umapCoordinatesB_low, aes(x = UMAP_1, y = UMAP_2),color = "hotpink3", size = 1.5, shape = 16) +
  create_lpr_theme() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = rel(0.6)),
        legend.text = element_text(size = rel(0.6), angle = 30),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  labs(color =  "Scaled\n Log UMI counts") 
ggsave(plot, file = paste0(plotDirectory, 'FM01_LowClusterFraction.svg'), width = 6, height = 6.158)


####################################################################################################################################
######################################## Neighbor analysis for bias #########################################################
####################################################################################################################################
library(RANN)
library(tripack)
library(reshape2)

pcaCoordinates = as_tibble(read.table(file = paste0(home1Directory, "pcaCoordinatesSCT_S1S2_50pcs_filter.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t"))
jointPCA = inner_join(linCountTooverlaps, pcaCoordinates, by = c("cellID", "sampleNum")) %>% select(-nLineages)
jointPCAS1 = jointPCA %>% filter(sampleNum=="S1") 
PCAtoAnalyze = jointPCAS1 %>% group_by(BC50StarcodeD8) %>% summarise(nColony = length(BC50StarcodeD8)) %>% filter(nColony >20) ####since for smaller colonies neighbor analysis not ideal or have to reduce the number of allowed neighbors (which is ok, maybe I can do a response curve)
jointPCAS1_onlyColonies = inner_join(jointPCAS1,PCAtoAnalyze, by = "BC50StarcodeD8")

neighborRandomFinalExtract = c();
neighborColonyFinalExtract = c();

for (i in c(1:length(PCAtoAnalyze$BC50StarcodeD8))) {
  neighborRandomFinal = c();
  neighborColonyFinal = c();
  for (j in c(1:5)) {
    subSampled1 = sample_n(jointPCAS1_onlyColonies, PCAtoAnalyze$nColony[i]) %>% mutate(name = "random1")
    subSampled2 = sample_n(jointPCAS1_onlyColonies, PCAtoAnalyze$nColony[i]) %>% mutate(name = "random2")
    random12 = bind_rows(subSampled1,subSampled2)
    random12BarcodeX = random12 %>% mutate(num = c(1:nrow(random12)))
    Random1index = random12BarcodeX %>% filter(name == "random1") %>% select(num)
    Random2index = random12BarcodeX %>% filter(name == "random2") %>% select(num)
    knnPCARandom = nn2(random12[,5:54], random12[,5:54],k = min(10, nrow(subSampled1)))
    neighborKNNRandom = as_tibble(knnPCARandom[[1]]) ####gets the neighbor
    
    barcodeColony = jointPCAS1_onlyColonies %>% filter(BC50StarcodeD8 == PCAtoAnalyze$BC50StarcodeD8[i]) %>% mutate(name = "colony")
    random1Colony = bind_rows(subSampled1,barcodeColony)
    random1ColonyBarcodeX = random1Colony %>% mutate(num = c(1:nrow(random1Colony)))
    Randomindex = random1ColonyBarcodeX %>% filter(name == "random1") %>% select(num)
    Colonyindex = random1ColonyBarcodeX %>% filter(name == "colony") %>% select(num)
    knnPCAColony = nn2(random1Colony[,5:54], random1Colony[,5:54], k = min(10, nrow(random1Colony)))
    neighborKNNColony = as_tibble(knnPCAColony[[1]]) ####gets the neighbor
    
    nRandom1 = c();
    nRandom2 = c();
    nColony1 = c();
    nColony2 = c();
    
    for (k in c(1:nrow(random12))) {
      nRandom1[k] =  sum(neighborKNNRandom[k,2:pmin(10,nrow(subSampled1))] %in% Random1index$num) 
      nRandom2[k] =  sum(neighborKNNRandom[k,2:pmin(10,nrow(subSampled1))] %in% Random2index$num)
      nColony1[k] =  sum(neighborKNNColony[k,2:pmin(10,nrow(random1Colony))] %in% Randomindex$num) 
      nColony2[k] =  sum(neighborKNNColony[k,2:pmin(10,nrow(random1Colony))] %in% Colonyindex$num) 
    }
    #neighborRandom = tibble(nRandom1,nRandom2) %>% mutate(fractionR1 = nRandom1/(nRandom1+nRandom2), fractionR2= nRandom2/(nRandom1+nRandom2)) %>% mutate(num = c(1:nrow(random12)))
    neighborRandomR1 = tibble(nRandom1,nRandom2) %>% mutate(fractionR1 = nRandom1/(nRandom1+nRandom2), fractionR2= nRandom2/(nRandom1+nRandom2)) %>% mutate(num = c(1:nrow(random12))) %>% filter(num %in% Random1index$num)
    neighborRandomR2 = tibble(nRandom1,nRandom2) %>% mutate(fractionR1 = nRandom1/(nRandom1+nRandom2), fractionR2= nRandom2/(nRandom1+nRandom2)) %>% mutate(num = c(1:nrow(random12))) %>% filter(num %in% Random2index$num)
    R1R1 = as_tibble(neighborRandomR1$fractionR1) %>% mutate(isSelf = "withSelf")
    R2R2 = as_tibble(neighborRandomR2$fractionR2) %>% mutate(isSelf = "withSelf")
    R1R2 = as_tibble(neighborRandomR1$fractionR2) %>% mutate(isSelf = "withOther")
    R2R1 = as_tibble(neighborRandomR2$fractionR1) %>% mutate(isSelf = "withOther")
    toPlotRandom = bind_rows(R1R1,R2R2,R1R2,R2R1) %>% rename(fraction = value) %>% mutate(isSelf = factor(isSelf, levels=c("withSelf","withOther")))
    withSelf = toPlotRandom %>% filter(isSelf == "withSelf") %>% select(-isSelf)
    withOther = toPlotRandom %>% filter(isSelf == "withOther") %>% select(-isSelf)
    neighborRandomFinal[j]= mean(withOther$fraction)/mean(withSelf$fraction)
    
    #neighborColony = tibble(nColony1,nColony2) %>% mutate(fractionR1 = nColony1/(nColony1+nColony2), fractionR2= nColony2/(nColony1+nColony2)) %>% mutate(num = c(1:nrow(random1Colony)))
    neighborColonyR1 = tibble(nColony1,nColony2) %>% mutate(fractionR1 = nColony1/(nColony1+nColony2), fractionR2= nColony2/(nColony1+nColony2)) %>% mutate(num = c(1:nrow(random1Colony))) %>% filter(num %in% Randomindex$num)
    neighborColonyR2 = tibble(nColony1,nColony2) %>% mutate(fractionR1 = nColony1/(nColony1+nColony2), fractionR2= nColony2/(nColony1+nColony2)) %>% mutate(num = c(1:nrow(random1Colony)))%>% filter(num %in% Colonyindex$num)
    R1R1_colony = as_tibble(neighborColonyR1$fractionR1) %>% mutate(isSelf = "withSelf")
    R2R2_colony = as_tibble(neighborColonyR2$fractionR2) %>% mutate(isSelf = "withSelf")
    R1R2_colony = as_tibble(neighborColonyR1$fractionR2) %>% mutate(isSelf = "withOther")
    R2R1_colony = as_tibble(neighborColonyR2$fractionR1) %>% mutate(isSelf = "withOther")
    toPlotColony = bind_rows(R1R1_colony,R2R2_colony,R1R2_colony,R2R1_colony) %>% rename(fraction = value) %>% mutate(isSelf = factor(isSelf, levels=c("withSelf","withOther")))
    withSelfColony = toPlotColony %>% filter(isSelf == "withSelf") %>% select(-isSelf)
    withOtherColony = toPlotColony %>% filter(isSelf == "withOther") %>% select(-isSelf)
    neighborColonyFinal[j]= mean(withOtherColony$fraction)/mean(withSelfColony$fraction)
  }
  neighborRandomFinalExtract[i] = mean(neighborRandomFinal)
  neighborColonyFinalExtract[i] = mean(neighborColonyFinal)
}

neighborRandomFinalExtract1 = tibble(mixingFinal = numeric()) %>% add_row(mixingFinal = neighborRandomFinalExtract) %>% mutate(type = "random")
neighborColonyFinalExtract1 = tibble(mixingFinal = numeric()) %>% add_row(mixingFinal =  neighborColonyFinalExtract) %>% mutate(type = "colony")

plotAll = bind_rows(neighborRandomFinalExtract1,neighborColonyFinalExtract1)
plotAll$type <- factor(plotAll$type, levels=c("random", "colony"))
plot = ggplot(plotAll, aes(x = type, y=mixingFinal, fill = as.factor(type)), show.legend = FALSE) + 
  geom_boxplot() +
  scale_fill_manual(values=c("grey","goldenrod3")) +
  stat_summary(fun=mean) +
  theme_classic((base_size = 36)) + theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks = element_line(colour = "black", size = 1.5), text=element_text(family="Helvetica")) + 
  theme(legend.position = "none") +
  ylim(0,1.2)
ggsave(plot, file = paste0(plotDirectory, 'FM01_neighborAnalysisToRandom.svg'), width = 4, height = 8)

plot = ggplot(plotAll, aes(x = type, y=mixingFinal, fill = as.factor(type)), show.legend = FALSE) + 
  geom_boxplot() +
  geom_jitter(width = 0.1, shape = 16) +
  scale_fill_manual(values=c("grey","goldenrod3")) +
  stat_summary(fun=mean) +
  theme_classic((base_size = 36)) + theme(axis.line = element_line(colour = 'black', size = 1.5), axis.ticks = element_line(colour = "black", size = 1.5), text=element_text(family="Helvetica")) + 
  theme(legend.position = "none") +
  ylim(0,1.2)
ggsave(plot, file = paste0(plotDirectory, 'FM01_neighborAnalysisToRandomWDots.svg'), width = 4, height = 8)

####################################################################################################################################
######################################## ColonySizeAnalysis #########################################################
####################################################################################################################################
###MLANA: 0; NGFR = 2, 12; ACTA2: 7; VCAM1 : 5,13; 10: IFIT2 10 (snn 0.4)

################################################
############this block only once################
################################################
sample1_2 <- FindClusters(object=sample1_2, resolution = 0.6, verbose = FALSE)
jointUMAP = inner_join(linCountTooverlaps, umapCoordinatesS1, by = c("cellID", "sampleNum")) %>% select(-nLineages)
umapClusters = (sample1_2[['seurat_clusters']])
cells_Clusters = rownames(umapClusters) #CellIds with Sample number as prefix
cells_Clusters_cellID = sub("S\\d_", "", cells_Clusters)
cells_Clusters_sampleNum = gsub("[^S12]", "", cells_Clusters)
umapClusters = as_tibble(umapClusters)
umapClusters = umapClusters %>% mutate(cellID = cells_Clusters_cellID,
                                       sampleNum = cells_Clusters_sampleNum)
write.table(umapClusters, file=paste0(home1Directory,'umapClusters_s1s2Scanorama_50pcs_filter_snn0_6.tsv'), col.names = TRUE, sep='\t')
###this block ^^ only once###########

umapClusters = as_tibble(read.table(file = paste0(home1Directory, "umapClusters_s1s2Scanorama_50pcs_filter_snn0_6.tsv"), header = TRUE, stringsAsFactors=F, sep = "\t"))

clusterMLANACellIDs = umapClusters %>% filter(seurat_clusters %in% c(0,3))
clusterMLANAUMAP = inner_join(clusterMLANACellIDs,umapCoordinates, by = c("cellID", "sampleNum"))
clusterMLANAUMAPBarcodes = inner_join(linCountTooverlaps, clusterMLANAUMAP, by = c("cellID", "sampleNum")) %>% select(-nLineages,seurat_clusters)
colonySizeClusterMLANA = clusterMLANAUMAPBarcodes %>% group_by(BC50StarcodeD8,sampleNum) %>% summarise(nColony = length(BC50StarcodeD8))
colonySizeAll = jointUMAP %>% group_by(BC50StarcodeD8,sampleNum) %>% summarise(nColony = length(BC50StarcodeD8))
fractionColony = inner_join(colonySizeClusterMLANA,colonySizeAll, by = c("BC50StarcodeD8", "sampleNum")) %>% mutate(fractionClusterMLANA = nColony.x/nColony.y) %>% filter(fractionClusterMLANA >0.5) %>% select(-nColony.x)
colonySizeClusterMLANASummary = fractionColony %>% group_by(nColony.y,sampleNum) %>% summarise(total = length(nColony.y))

colonySizeClusterMLANASummary2 = colonySizeClusterMLANASummary %>% group_by(sampleNum) %>% summarise(singlet = sum(total[nColony.y == 1]),
                                                                                                     small = sum(total[nColony.y >1 & nColony.y <4]),
                                                                                                     large = sum(total[nColony.y >3]))
colonySizeClusterMLANASummary2 = colonySizeClusterMLANASummary2 %>% gather("type", "value", -sampleNum)
MLANASummary = colonySizeClusterMLANASummary2 %>% mutate(cluster = "mlana") %>% mutate(percent = 100*value/(sum(value)))

####NGFR
clusterMLANACellIDs = umapClusters %>% filter(seurat_clusters %in% c(7))
clusterMLANAUMAP = inner_join(clusterMLANACellIDs,umapCoordinates, by = c("cellID", "sampleNum"))
clusterMLANAUMAPBarcodes = inner_join(linCountTooverlaps, clusterMLANAUMAP, by = c("cellID", "sampleNum")) %>% select(-nLineages,seurat_clusters)
colonySizeClusterMLANA = clusterMLANAUMAPBarcodes %>% group_by(BC50StarcodeD8,sampleNum) %>% summarise(nColony = length(BC50StarcodeD8))
colonySizeAll = jointUMAP %>% group_by(BC50StarcodeD8,sampleNum) %>% summarise(nColony = length(BC50StarcodeD8))
fractionColony = inner_join(colonySizeClusterMLANA,colonySizeAll, by = c("BC50StarcodeD8", "sampleNum")) %>% mutate(fractionClusterMLANA = nColony.x/nColony.y) %>% filter(fractionClusterMLANA >0.5) %>% select(-nColony.x)
colonySizeClusterMLANASummary = fractionColony %>% group_by(nColony.y,sampleNum) %>% summarise(total = length(nColony.y))

colonySizeClusterMLANASummary2 = colonySizeClusterMLANASummary %>% group_by(sampleNum) %>% summarise(singlet = sum(total[nColony.y == 1]),
                                                                                                     small = sum(total[nColony.y >1 & nColony.y <4]),
                                                                                                     large = sum(total[nColony.y >3]))
colonySizeClusterMLANASummary2 = colonySizeClusterMLANASummary2 %>% gather("type", "value", -sampleNum) %>% mutate(percent = 100*value/(sum(value)))
NGFRSummary = colonySizeClusterMLANASummary2 %>% mutate(cluster = "ngfr")

####IFIT2
clusterMLANACellIDs = umapClusters %>% filter(seurat_clusters %in% c(12))
clusterMLANAUMAP = inner_join(clusterMLANACellIDs,umapCoordinates, by = c("cellID", "sampleNum"))
clusterMLANAUMAPBarcodes = inner_join(linCountTooverlaps, clusterMLANAUMAP, by = c("cellID", "sampleNum")) %>% select(-nLineages,seurat_clusters)
colonySizeClusterMLANA = clusterMLANAUMAPBarcodes %>% group_by(BC50StarcodeD8,sampleNum) %>% summarise(nColony = length(BC50StarcodeD8))
colonySizeAll = jointUMAP %>% group_by(BC50StarcodeD8,sampleNum) %>% summarise(nColony = length(BC50StarcodeD8))
fractionColony = inner_join(colonySizeClusterMLANA,colonySizeAll, by = c("BC50StarcodeD8", "sampleNum")) %>% mutate(fractionClusterMLANA = nColony.x/nColony.y) %>% filter(fractionClusterMLANA >0.5) %>% select(-nColony.x)
colonySizeClusterMLANASummary = fractionColony %>% group_by(nColony.y,sampleNum) %>% summarise(total = length(nColony.y))

colonySizeClusterMLANASummary2 = colonySizeClusterMLANASummary %>% group_by(sampleNum) %>% summarise(singlet = sum(total[nColony.y == 1]),
                                                                                                     small = sum(total[nColony.y >1 & nColony.y <4]),
                                                                                                     large = sum(total[nColony.y >3]))
colonySizeClusterMLANASummary2 = colonySizeClusterMLANASummary2 %>% gather("type", "value", -sampleNum) %>% mutate(percent = 100*value/(sum(value)))
IFITSummary = colonySizeClusterMLANASummary2 %>% mutate(cluster = "IFIT")

####VCAM1

clusterMLANACellIDs = umapClusters %>% filter(seurat_clusters %in% c(15,4,6))
clusterMLANAUMAP = inner_join(clusterMLANACellIDs,umapCoordinates, by = c("cellID", "sampleNum"))
clusterMLANAUMAPBarcodes = inner_join(linCountTooverlaps, clusterMLANAUMAP, by = c("cellID", "sampleNum")) %>% select(-nLineages,seurat_clusters)
colonySizeClusterMLANA = clusterMLANAUMAPBarcodes %>% group_by(BC50StarcodeD8,sampleNum) %>% summarise(nColony = length(BC50StarcodeD8))
colonySizeAll = jointUMAP %>% group_by(BC50StarcodeD8,sampleNum) %>% summarise(nColony = length(BC50StarcodeD8))
fractionColony = inner_join(colonySizeClusterMLANA,colonySizeAll, by = c("BC50StarcodeD8", "sampleNum")) %>% mutate(fractionClusterMLANA = nColony.x/nColony.y) %>% filter(fractionClusterMLANA >0.5) %>% select(-nColony.x)
colonySizeClusterMLANASummary = fractionColony %>% group_by(nColony.y,sampleNum) %>% summarise(total = length(nColony.y))

colonySizeClusterMLANASummary2 = colonySizeClusterMLANASummary %>% group_by(sampleNum) %>% summarise(singlet = sum(total[nColony.y == 1]),
                                                                                                     small = sum(total[nColony.y >1 & nColony.y <4]),
                                                                                                     large = sum(total[nColony.y >3]))
colonySizeClusterMLANASummary2 = colonySizeClusterMLANASummary2 %>% gather("type", "value", -sampleNum) %>% mutate(percent = 100*value/(sum(value)))
VCAM1Summary = colonySizeClusterMLANASummary2 %>% mutate(cluster = "VCAM1")

####ACTA2
clusterMLANACellIDs = umapClusters %>% filter(seurat_clusters %in% c(8))
clusterMLANAUMAP = inner_join(clusterMLANACellIDs,umapCoordinates, by = c("cellID", "sampleNum"))
clusterMLANAUMAPBarcodes = inner_join(linCountTooverlaps, clusterMLANAUMAP, by = c("cellID", "sampleNum")) %>% select(-nLineages,seurat_clusters)
colonySizeClusterMLANA = clusterMLANAUMAPBarcodes %>% group_by(BC50StarcodeD8,sampleNum) %>% summarise(nColony = length(BC50StarcodeD8))
colonySizeAll = jointUMAP %>% group_by(BC50StarcodeD8,sampleNum) %>% summarise(nColony = length(BC50StarcodeD8))
fractionColony = inner_join(colonySizeClusterMLANA,colonySizeAll, by = c("BC50StarcodeD8", "sampleNum")) %>% mutate(fractionClusterMLANA = nColony.x/nColony.y) %>% filter(fractionClusterMLANA >0.5) %>% select(-nColony.x)
colonySizeClusterMLANASummary = fractionColony %>% group_by(nColony.y,sampleNum) %>% summarise(total = length(nColony.y))

colonySizeClusterMLANASummary2 = colonySizeClusterMLANASummary %>% group_by(sampleNum) %>% summarise(singlet = sum(total[nColony.y == 1]),
                                                                                                     small = sum(total[nColony.y >1 & nColony.y <4]),
                                                                                                     large = sum(total[nColony.y >3]))
colonySizeClusterMLANASummary2 = colonySizeClusterMLANASummary2 %>% gather("type", "value", -sampleNum)
ACTA2Summary = colonySizeClusterMLANASummary2 %>% mutate(cluster = "ACTA2") %>% mutate(percent = 100*value/(sum(value)))


####AXL, SERPINE1
clusterMLANACellIDs = umapClusters %>% filter(seurat_clusters %in% c(11,1,2))
clusterMLANAUMAP = inner_join(clusterMLANACellIDs,umapCoordinates, by = c("cellID", "sampleNum"))
clusterMLANAUMAPBarcodes = inner_join(linCountTooverlaps, clusterMLANAUMAP, by = c("cellID", "sampleNum")) %>% select(-nLineages,seurat_clusters)
colonySizeClusterMLANA = clusterMLANAUMAPBarcodes %>% group_by(BC50StarcodeD8,sampleNum) %>% summarise(nColony = length(BC50StarcodeD8))
colonySizeAll = jointUMAP %>% group_by(BC50StarcodeD8,sampleNum) %>% summarise(nColony = length(BC50StarcodeD8))
fractionColony = inner_join(colonySizeClusterMLANA,colonySizeAll, by = c("BC50StarcodeD8", "sampleNum")) %>% mutate(fractionClusterMLANA = nColony.x/nColony.y) %>% filter(fractionClusterMLANA >0.5) %>% select(-nColony.x)
colonySizeClusterMLANASummary = fractionColony %>% group_by(nColony.y,sampleNum) %>% summarise(total = length(nColony.y))

colonySizeClusterMLANASummary2 = colonySizeClusterMLANASummary %>% group_by(sampleNum) %>% summarise(singlet = sum(total[nColony.y == 1]),
                                                                                                     small = sum(total[nColony.y >1 & nColony.y <4]),
                                                                                                     large = sum(total[nColony.y >3]))
colonySizeClusterMLANASummary2 = colonySizeClusterMLANASummary2 %>% gather("type", "value", -sampleNum)
AXLSummary = colonySizeClusterMLANASummary2 %>% mutate(cluster = "AXL") %>% mutate(percent = 100*value/(sum(value)))

#########
All = bind_rows(AXLSummary, ACTA2Summary,MLANASummary, NGFRSummary,IFITSummary, VCAM1Summary)
All$type <- factor(All$type, levels=c("large", "small", "singlet"))
All$cluster <- factor(All$cluster, levels=c("mlana", "ngfr", "IFIT", "ACTA2", "AXL", "VCAM1"))

write.table(All, file=paste0(plotDirectory,'FM01_AllColoniesSummary_PercentValue.tsv'), col.names = TRUE, sep='\t')

plot = ggplot(data=All, aes(x=cluster, y=percent, fill=type)) +
  geom_bar(stat="identity") +
  create_lpr_theme() + theme(legend.position = "none")
ggsave(plot, file = paste0(plotDirectory, 'FM01_AllColoniesSummary_Percent.svg'), width = 8, height = 5)

plot = ggplot(data=All, aes(x=cluster, y=value, fill=type)) +
  geom_bar(stat="identity") +
  create_lpr_theme() + theme(legend.position = "none")
ggsave(plot, file = paste0(plotDirectory, 'FM01_AllColoniesSummary_Value.svg'), width = 8, height = 5)
