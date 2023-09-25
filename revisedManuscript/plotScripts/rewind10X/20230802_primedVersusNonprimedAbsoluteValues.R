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

homeDirectory <- "/Users/naveenj/Dropbox (RajLab)/Shared_Naveen/Paper/extractedData/rewind10X/R1/"
plotDirectory <- "/Users/naveenj/Dropbox (RajLab)/Shared_Naveen/Revisions/plots/rewind10X/"

scanorama <- readRDS(paste0(homeDirectory, "scTransform.integrated"))
scanorama <- FindClusters(object = scanorama, resolution = 0.45)
Idents(scanorama) <- scanorama$seurat_clusters
DimPlot(scanorama)

barcodesR1 <- read.table("/Users/naveenj/Dropbox (RajLab)/Shared_Naveen/Paper/extractedData/rewind10X/R1/stepThreeStarcodeShavedReads_BC_10XAndGDNA.txt", header = TRUE, sep = '\t')
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

logNormCounts <- read.table(paste0(homeDirectory, "logNormalizedCounts_Scanorama_50pcs_filterRound.tsv"), header = T, sep = "\t")

logNormPlot <- inner_join(logNormCounts, seuratLabelPriming, by = "cellID") %>% dplyr::filter(label %in% c("nonprimed", "primed"))

logNormPlotSelect <- logNormPlot %>% dplyr::select(cellID, label, NANOG, POU5F1, SOX2, COL1A1, COL1A2, S100A4,
                                   SERPINE2, THBS1, ACTA2, TPM2,
                                   SPP1, FTH1, GDF15, CDKN1A, CENPF, TOP2A, ASPM, MKI67, UBC, GAPDH)

logNormPlotSelectMelt <- melt(logNormPlotSelect, id.vars = c("cellID", "label"))

ggplot(logNormPlotSelectMelt, aes(x = label, y = value)) + 
  stat_summary(fun = mean, geom = "col", aes(fill = label)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.25) +
  facet_wrap(~variable, nrow = 2) + NoLegend() + ylab("log(normalized RNA counts)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), axis.title.x = element_blank())
ggsave(filename = paste0(plotDirectory, "primingMarkersAbsoluteValues.pdf"), height = 4, width = 8)
