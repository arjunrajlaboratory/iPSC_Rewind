library(tidyverse)
library(reshape2)
library(fuzzyjoin)
library(spgs)
library(stringdist)
library(stringr)

theme_set(theme_classic())

dataDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/extractedData/prolifOverlap/"
plotDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/plots/prolifOverlap/"

#### load in groupd shaved reads table from 10X barcode matching pipeline ####
barcodeTable <- as_tibble(read.table(file = paste0(dataDirectory, "exp3/stepThreeStarcodeShavedReads.txt"), header = TRUE, stringsAsFactors = FALSE, sep = "\t"))

#### normalize read counts based on spike-in barcodes ####
barcodeTableGroup <- barcodeTable %>% dplyr::select(UMI, BC30StarcodeD8, SampleName)
length(unique(barcodeTableGroup$BC30StarcodeD8))
barcodeTableGroup <- barcodeTableGroup %>% group_by(BC30StarcodeD8, SampleName) %>% mutate(nUMI = sum(UMI)) %>% dplyr::select(-UMI) %>% unique() %>% dplyr::rename(BC = BC30StarcodeD8)

sampleNames <- unique(barcodeTable$SampleName)

spikeBC1 <- 'TCCAGGTCCTCCTACTTGTACAACACCTTGTACAGCTGCTAGTGGTAGAAGAGGTACAACAACAACACGAGCATCATGAGGATCTACAGCATCAAGAACA' %>% substr(0,30)
spikeBC2 <- 'ACGTTGTGCATGACCTTGATCACCAGCTCGATGTCGAACATCACGAGCTCGTTCTGCATCTGCAAGAACACCTCGTCCTTGAACTGCTCGACGTCCATGA' %>% substr(0,30)
nBC1 <- 20000
nBC2 <- 5000

sampleNamesNorm <- c(sampleNames[[6]], sampleNames[[7]])
standardTableAll = list()
lmr = list()
for (i in 1:length(sampleNamesNorm)) { #only need to do heritability samples because did not add to overlap sorted samples
  standardTable <- filter(barcodeTableGroup, SampleName == sampleNamesNorm[i]) %>% filter(BC == spikeBC1 | BC == spikeBC2) %>% arrange(-nUMI) %>% mutate(ratio = .$nUMI[1]/.$nUMI[2])
  lmr[i] = coef(lm(c(nBC1, nBC2) ~ 0 + standardTable$nUMI))
  if(is.null(dim(standardTableAll))){
    standardTableAll = standardTable
  } else {
    standardTableAll = bind_rows(standardTableAll, standardTable)
  }
}

barcodeTableGroupFilter <- barcodeTableGroup %>% filter(BC != spikeBC1) %>% filter(BC != spikeBC2)

#### compare overlap for heritability splits in gDNA barcodes #####
heritA <- filter(barcodeTableGroupFilter, SampleName == sampleNames[6]) %>% mutate(nUMI = nUMI*lmr[[1]]*1000) %>% ungroup() %>% dplyr::select(-SampleName)
heritB <- filter(barcodeTableGroupFilter, SampleName == sampleNames[7]) %>% mutate(nUMI = nUMI*lmr[[2]]*1000) %>% ungroup() %>% dplyr::select(-SampleName)

overlap <- full_join(heritA, heritB, by = "BC") %>% replace(is.na(.), 0)

ggplot() +
  geom_point(data = overlap, aes(x = log10(nUMI.x), y = log10(nUMI.y))) +
  geom_abline(slope = 1, intercept = 0) +
  theme_classic() + xlab("normalized counts in split A") + ylab("normalized counts in split B")

overlap <- overlap %>% rowwise() %>% mutate(nUMI = mean(nUMI.x, nUMI.y)) %>% ungroup()

#### compare overlap for individual overlap sorts with iPSCs ####
sampleNamesOverlap <- c(sampleNames[2], sampleNames[1], sampleNames[3])
for(i in 1:length(sampleNamesOverlap)) {
  normalizeTableTemp <- filter(barcodeTableGroupFilter, SampleName == sampleNamesOverlap[i]) %>% rowwise() %>% mutate(nUMI = round(nUMI/sum(.$nUMI)*10^6)) %>% ungroup()
  barcodeTableGroupFilter[barcodeTableGroupFilter$SampleName == sampleNamesOverlap[i], ] <- normalizeTableTemp
}

colonyList <- c(10, 25, 50, 100, 250, 500)
cutoffList <- c(0.1, 0.25, 0.5, 0.9)
for(n in 1:length(colonyList)) {
  colonyTable <- overlap %>% rowwise() %>% mutate(nUMI = mean(nUMI.x, nUMI.y)) %>% ungroup() %>% slice_max(nUMI, n = colonyList[n], with_ties = FALSE)
  overlapMatrix <- matrix(nrow = (length(sampleNamesOverlap)), ncol = length(cutoffList))
  for(i in 1:length(sampleNamesOverlap)) {
    for(j in 1:length(cutoffList)) {
      overlapTableTemp <- filter(barcodeTableGroupFilter, SampleName == sampleNamesOverlap[i]) %>% ungroup() %>% slice_max(nUMI, prop = cutoffList[j], with_ties = FALSE)
      overlapMatrix[i, j] <- nrow(inner_join(overlapTableTemp, colonyTable, by = 'BC'))/nrow(overlapTableTemp)
    }
  }
  rownames(overlapMatrix) <- sampleNamesOverlap
  colnames(overlapMatrix) <- cutoffList
  for(k in 1:ncol(overlapMatrix)) {
    overlapMatrix[, k] <- overlapMatrix[, k] / overlapMatrix[sampleNamesOverlap[2], k]
  }
  overlapMatrix <- as.data.frame(overlapMatrix) %>% mutate(name = rownames(overlapMatrix)) %>% melt(., id.vars = 'name')
  overlapMatrix$name <- factor(overlapMatrix$name, levels = sampleNamesOverlap, labels = c("slow", "control", "fast"))
  overlapMatrix <- overlapMatrix %>% mutate(colonyThresh = paste0(colonyList[n]))
  if(n == 1) {
    overlapMatrixPlot <- overlapMatrix
  } else{
    overlapMatrixPlot <- bind_rows(overlapMatrixPlot, overlapMatrix)
  }
}
overlapMatrixPlot$colonyThresh <- factor(overlapMatrixPlot$colonyThresh, levels = as.character(colonyList))
plot1 <- ggplot(overlapMatrixPlot, aes(y = value, x = name)) +
  geom_col(aes(fill = name)) + facet_grid(variable ~ colonyThresh, scales = "fixed") + ggtitle("exp 3") +
  ylab('normalized overlap') + theme_classic() + theme(axis.title.x = element_blank(), legend.position = 'none')

saveRDS(overlapMatrixPlot, file = paste0(dataDirectory, "overlapMatrixTable_exp3.rds"))
ggsave(filename = paste0(plotDirectory, "overlapMatrixPlot_exp3.pdf"), plot = plot1, height = 5, width = 5, unit = "in")

#### compare overlap across individual overlap sorts ####
sampleNamesOverlap <- c(sampleNames[2], sampleNames[1], sampleNames[3])
cutoffList <- c(0.025, 0.05, 0.1, 0.25, 0.5, 0.9)
plotList <- list()

for(m in 1:length(sampleNamesOverlap)) {
  for(n in 1:length(sampleNamesOverlap)) {
    sampleOverlap <- data_frame(overlap = NA, cutoff = cutoffList)
    for(j in 1:length(cutoffList)) {
      sampleA <- barcodeTableGroupFilter %>% filter(SampleName == sampleNamesOverlap[m]) %>% ungroup %>% slice_max(order_by = nUMI, prop = cutoffList[j], with_ties = FALSE)
      sampleB <- barcodeTableGroupFilter %>% filter(SampleName == sampleNamesOverlap[n]) %>% ungroup %>% slice_max(order_by = nUMI, prop = cutoffList[j], with_ties = FALSE)
      overlap <- nrow(inner_join(sampleA, sampleB, by = "BC"))/min(nrow(sampleA), nrow(sampleB))
      sampleOverlap$overlap[j] <- overlap
    }
    sampleOverlap <- sampleOverlap %>% mutate(sample = sampleNamesOverlap[n])
    if(n == 1) {
      sampleOverlapPlot <- sampleOverlap
    } else{
      sampleOverlapPlot <- bind_rows(sampleOverlapPlot, sampleOverlap)
    }
  }
  sampleOverlapPlot$cutoff <- factor(sampleOverlapPlot$cutoff, levels = (as.character(cutoffList)))
  sampleOverlapPlot$sample <- factor(sampleOverlapPlot$sample, levels = sampleNamesOverlap, labels = c("H", "C", "L"))
  plotList[[m]] <- ggplot(sampleOverlapPlot, aes(x = cutoff, y = overlap, group = sample, color = sample)) +
    geom_point() + geom_line() + ylim(0, 1) + theme_classic() + ggtitle(paste0("overlap with ", sampleNamesOverlap[m], "\nsample with changing cutoff"))
}
plot2 <- egg::ggarrange(plots = plotList, ncol = 3)

ggsave(filename = paste0(plotDirectory, "sampleOverlapPlot_exp3.pdf"), plot = plot2, height = 3, width = 9, unit = "in")

for(n in 1:length(sampleNamesOverlap)) {
  sampleOverlap <- data_frame(overlap = NA, cutoff = cutoffList)
  for(j in 1:length(cutoffList)) {
    sampleA <- barcodeTableGroupFilter %>% filter(SampleName == sampleNamesOverlap[2]) %>% ungroup %>% slice_max(order_by = nUMI, prop = cutoffList[j], with_ties = FALSE)
    sampleB <- barcodeTableGroupFilter %>% filter(SampleName == sampleNamesOverlap[n]) %>% ungroup %>% slice_max(order_by = nUMI, prop = cutoffList[j], with_ties = FALSE)
    overlap <- nrow(inner_join(sampleA, sampleB, by = "BC"))/min(nrow(sampleA), nrow(sampleB))
    sampleOverlap$overlap[j] <- overlap
  }
  sampleOverlap <- sampleOverlap %>% mutate(sample = sampleNamesOverlap[n])
  if(n == 1) {
    sampleOverlapPlot <- sampleOverlap
  } else{
    sampleOverlapPlot <- bind_rows(sampleOverlapPlot, sampleOverlap)
  }
}
sampleOverlapPlot$cutoff <- factor(sampleOverlapPlot$cutoff, levels = (as.character(cutoffList)))
sampleOverlapPlot$sample <- factor(sampleOverlapPlot$sample, levels = sampleNamesOverlap, labels = c("H", "C", "L"))
saveRDS(sampleOverlapPlot, file = paste0(dataDirectory, "sampleOverlapPlot_exp3.rds"))
