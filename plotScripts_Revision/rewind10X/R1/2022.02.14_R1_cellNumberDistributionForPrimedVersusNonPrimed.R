library(tidyverse)
library(fuzzyjoin)
library(reshape2)
library(egg)
library(ggsignif)
library(ggrepel)
library(ggpubr)
library(scales)
library(spgs)
library(Seurat)
library(clustree)
options(future.globals.maxSize = 8000 * 1024 ^ 2)

dataDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/extractedData/rewind10X/R1/"
plotDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/plots/rewind10X/R1/"

filteredBarcodes <- read.table(paste0(dataDirectory, 'filtered10XCells.txt'), header = TRUE, sep = '\t')
filteredBarcodes <- filteredBarcodes %>% filter(BC50StarcodeD8 != 'ATTCTAGTTGTAGTACGAGTAGCACATGTTCTACGTGGAGGACGAGAACG')
nCellsPerBarcodeOverlapAll <- filteredBarcodes %>% group_by(BC50StarcodeD8) %>% summarise(nCells = length(cellID))

barcodesR1 <- read.table("/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/extractedData/rewind10X/R1/stepThreeStarcodeShavedReads_BC_10XAndGDNA.txt", header = TRUE, sep = '\t')
probedBarcodesR1 <- barcodesR1 %>% filter(cellID == "dummy") %>% dplyr::select(UMI, BC50StarcodeD8, SampleNum)
probedBarcodesR1$UMI <- as.numeric(probedBarcodesR1$UMI)
probedBarcodesR1 <- probedBarcodesR1 %>% group_by(BC50StarcodeD8, SampleNum) %>% summarise(nUMI = sum(UMI))
probedBarcodesTopR1 <- probedBarcodesR1 %>% ungroup() %>% slice_max(., order_by = nUMI, n = 100)

primedCells <- inner_join(filteredBarcodes, probedBarcodesTopR1, by = "BC50StarcodeD8") %>% mutate(label = "primed")
nCellsPerBarcodeOverlapPrimed <- primedCells %>% group_by(BC50StarcodeD8) %>% summarise(nCells = length(cellID)) %>% mutate(label = "primed")
nonprimedCells <- anti_join(filteredBarcodes, probedBarcodesTopR1, by = "BC50StarcodeD8") %>% mutate(label = "nonprimed")
nCellsPerBarcodeOverlapNonprimed <- nonprimedCells %>% group_by(BC50StarcodeD8) %>% summarise(nCells = length(cellID)) %>% mutate(label = "nonprimed")

plot1 <- ggplot() +
  geom_histogram(data = nCellsPerBarcodeOverlapAll, aes(x = nCells, y = ..density..), binwidth = 1, fill = "black") +
  geom_histogram(data = nCellsPerBarcodeOverlapPrimed, aes(x = nCells, y = ..density..), binwidth = 1, fill = "red", alpha = 0.5) +
  xlab("number of cells per barcode") + theme_classic() + xlim(0, max(nCellsPerBarcodeOverlapAll$nCells)) + ylim(0, 1) +
  ggtitle("distribution for all barcoded cells")

nCellsPerBarcodeOverlapUnprimedTemp <- nCellsPerBarcodeOverlapNonprimed #%>% sample_n(nrow(nCellsPerBarcodeOverlapPrimed))

plot5 <- ggplot() +
  geom_histogram(data = nCellsPerBarcodeOverlapPrimed, aes(x = nCells, y = ..density..), binwidth = 1, fill = "black") +
  xlab("number of cells per barcode") + theme_classic() + xlim(0.5, 5) + ylim(0, 1) +
  geom_vline(xintercept = mean(nCellsPerBarcodeOverlapPrimed$nCells), linetype = "dashed") +
  annotate(geom = "text", label = paste0("mean = ", round(mean(nCellsPerBarcodeOverlapPrimed$nCells), 2)), y = Inf, x = mean(nCellsPerBarcodeOverlapPrimed$nCells), vjust = 1, hjust = -0.2) +
  ggtitle("distribution for primed barcoded cells")

plot6 <- ggplot() +
  geom_histogram(data = nCellsPerBarcodeOverlapUnprimedTemp, aes(x = nCells, y = ..density..), binwidth = 1, fill = "black") +
  xlab("number of cells per barcode") + theme_classic() + ylim(0, 1) +
  geom_vline(xintercept = mean(nCellsPerBarcodeOverlapUnprimedTemp$nCells), linetype = "dashed") +
  annotate(geom = "text", label = paste0("mean = ", round(mean(nCellsPerBarcodeOverlapUnprimedTemp$nCells), 2)), y = Inf, x = mean(nCellsPerBarcodeOverlapUnprimedTemp$nCells), vjust = 1, hjust = -0.2) +
  ggtitle("randomly sampled distribution of same size for nonprimed barcoded cells") +
  scale_x_continuous(breaks = seq(1, 5, by = 1), limits = c(0.5, 5))

plot7 <- egg::ggarrange(plot5, plot6, nrow = 2)
ggsave(plot7, file = paste0(plotDirectory, "cellNumberDistributionPlotPrimedVersusNonprimed.pdf"), height = 4, width = 2)

nCellsPerBarcodeOverlapComb <- bind_rows(nCellsPerBarcodeOverlapPrimed, nCellsPerBarcodeOverlapNonprimed)
ggplot(nCellsPerBarcodeOverlapComb, aes(x = label, y = nCells)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  geom_signif(comparisons = list(c("nonprimed", "primed")), test = "t.test")

meanList <- c()
for(i in 1:1000){
  nCellsPerBarcodeOverlapUnprimedTemp <- nCellsPerBarcodeOverlapNonprimed %>% sample_n(nrow(nCellsPerBarcodeOverlapPrimed))
  meanList <- c(meanList, mean(nCellsPerBarcodeOverlapUnprimedTemp$nCells))
}

meanList <- as.data.frame(meanList)
pValue <- round(sum(meanList$meanList > mean(nCellsPerBarcodeOverlapPrimed$nCells))/1000, 3)

plot2 <- ggplot() +
  geom_histogram(data = nCellsPerBarcodeOverlapPrimed, aes(x = nCells, y = ..density..), binwidth = 1, fill = "black") +
  xlab("number of cells per barcode") + theme_classic() + xlim(0, max(nCellsPerBarcodeOverlapAll$nCells)) +
  ggtitle("distribution for primed barcoded cells") + ylim(0, 1)

plot3 <- ggplot() +
  geom_histogram(data = meanList, aes(x = meanList), fill = "black") +
  geom_vline(xintercept = mean(nCellsPerBarcodeOverlapPrimed$nCells), color = "blue") +
  annotate(geom = "text", label = paste0("p-value = ", pValue), x = Inf, y = Inf, hjust = 1, vjust = 1) +
  xlab("mean number of cells per barcode per set") + ylab("count") + theme_classic() +
  ggtitle("distribution of 1000 randomly sampled\nsets of nonprimed barcoded cells")

plot4 <- egg::ggarrange(plot2, plot3, ncol = 2)
ggsave(plot4, file = paste0(plotDirectory, "cellNumberDistributionPlot_Cutoff", j, ".pdf"), height = 4, width = 8)