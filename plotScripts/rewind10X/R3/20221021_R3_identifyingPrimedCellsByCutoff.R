rm(list=ls())
gc()

library(tidyverse)
library(reshape2)
library(ggrepel)
library(ggforce)
library(concaveman)
library(egg)
library(RColorBrewer)
library(Seurat)
library(biomaRt)
library(spgs)
library(ggsignif)
library(ggridges)
library(viridis)

theme_set(theme_classic())

homeDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/extractedData/rewind10X/R3/"
plotDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/plots/rewind10X/R3/"

scanorama_filter <- readRDS(paste0(homeDirectory, "scanorma_filter.rds"))
barcodes = as_tibble(read.table(paste0(homeDirectory, 'stepThreeStarcodeShavedReads_BC_10XAndGDNA.txt'), stringsAsFactors = F, header = T))
linCountToOverlaps <- read.table(file = paste0(homeDirectory, "filtered10XCells.txt"), sep = "\t")

barcodesMergeEnd <- barcodes %>% filter(cellID == "dummy")
barcodesMergeEnd$UMI <- as.numeric(barcodesMergeEnd$UMI)

#### identifying gDNA barcodes for primed lineages w/ 3 independent transductions ####
#######################################################################################################################################################
sampleNames <- c("FS_1A", "FS_1B", "FS_2A", "FS_2B", "FS_3A", "FS_3B")

barcodesMergeEndCount <- barcodesMergeEnd %>% filter(SampleNum %in% sampleNames) %>%
  group_by(BC50StarcodeD8, SampleNum) %>% summarise(nUMI = sum(UMI))

spikeBC1 <- 'TCCAGGTCCTCCTACTTGTACAACACCTTGTACAGCTGCTAGTGGTAGAAGAGGTACAACAACAACACGAGCATCATGAGGATCTACAGCATCAAGAACA' %>% reverseComplement(., content = "dna", case = "as is") %>% substr(0,50)
spikeBC2 <- 'ACGTTGTGCATGACCTTGATCACCAGCTCGATGTCGAACATCACGAGCTCGTTCTGCATCTGCAAGAACACCTCGTCCTTGAACTGCTCGACGTCCATGA' %>% reverseComplement(., content = "dna", case = "as is") %>% substr(0,50)
nBC1 <- 20000
nBC2 <- 5000

standardTableAll = list()
lmr = list()
for (i in 1:length(sampleNames)) {
  standardTable <- filter(barcodesMergeEndCount, SampleNum == sampleNames[i]) %>% dplyr::filter(., BC50StarcodeD8 == spikeBC1 | BC50StarcodeD8 == spikeBC2) %>% arrange(-nUMI) %>% mutate(ratio = .$nUMI[1]/.$nUMI[2])
  lmr[i] = coef(lm(c(nBC1, nBC2) ~ 0 + standardTable$nUMI))
  if(is.null(dim(standardTableAll))){
    standardTableAll = standardTable
  } else {
    standardTableAll = bind_rows(standardTableAll, standardTable)
  }
}
barcodesMergeEndCountFilter <- barcodesMergeEndCount %>% filter(BC50StarcodeD8 != spikeBC1) %>% filter(BC50StarcodeD8 != spikeBC2)

#### overlap scatterplots for different cutoffs and different replicates ####
#######################################################################################################################################################
overlapCellIDList <- list()
overlapTableList <- list()
sampleList <- c(1, 3, 5)
cutoffList <- c(10, 25, 50, 100, 150, 200, 250, 500, 1000)
k <- 1
for(i in 1:length(sampleList)){
  sampleA <- barcodesMergeEndCountFilter %>% filter(SampleNum == sampleNames[sampleList[i]]) %>% 
    rowwise() %>% mutate(nUMINorm = nUMI*lmr[[i]]) %>% ungroup() %>% dplyr::select(-SampleNum)
  sampleB <- barcodesMergeEndCountFilter %>% filter(SampleNum == sampleNames[sampleList[i] + 1]) %>% 
    rowwise() %>% mutate(nUMINorm = nUMI*lmr[[i + 1]]) %>% ungroup() %>% dplyr::select(-SampleNum)
  sampleOverlap <- full_join(sampleA, sampleB, by = "BC50StarcodeD8") %>% replace(is.na(.), 0)
  overlapTableList[[i]] <- sampleOverlap
  overlapValue <- nrow(inner_join(sampleA, sampleB, by = "BC50StarcodeD8")) / max(nrow(sampleA), nrow(sampleB))
  plot <- ggplot(sampleOverlap, aes(x = log2(nUMINorm.x + 1), y = log2(nUMINorm.y + 1))) +
    geom_point() + geom_abline(slope = 1, intercept = 1) +
    annotate("text", x = Inf, y = Inf, label = paste0(round(overlapValue, 2)), hjust = 1, vjust = 1) +
    theme_classic()
  ggsave(plot = plot, file = paste0(plotDirectory, 'R3_heritability_', sampleNames[sampleList[i]], "vs", sampleNames[sampleList[i] + 1], '.pdf'), width = 8, height = 8)
  for(j in 1:length(cutoffList)){
    sampleOverlapTop <- sampleOverlap %>% rowwise() %>% mutate(nUMIMax = max(nUMINorm.x, nUMINorm.y)) %>%
      ungroup() %>% slice_max(nUMIMax, n = cutoffList[j], with_ties = FALSE)
    overlapCellIDList[[k]] <- filter(linCountToOverlaps, barcode %in% sampleOverlapTop$BC50StarcodeD8)$cellID %>% unique() %>% list(.)
    k <- k + 1
  }
}

for(i in 1:length(overlapCellIDList)){
  print(length(unlist(overlapCellIDList[[i]])))
}

saveRDS(overlapCellIDList, file = paste0(homeDirectory, "primedCellIDList.rds"))
saveRDS(overlapTableList, file = paste0(homeDirectory, "overalapTableList.rds"))
