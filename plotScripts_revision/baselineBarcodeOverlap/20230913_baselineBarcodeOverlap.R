rm(list=ls())
gc()

library(tidyverse)
library(spgs)
library(concaveman)
library(egg)
library(reshape2)
library(ggrepel)
library(ggforce)
library(ggsignif)
library(ggridges)
library(ggpmisc)
library(viridis)
library(VennDiagram)

theme_set(theme_classic())

homeDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Revisions/extractedData/baselineBarcodeOverlap/"
plotDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Revisions/plots/baselineBarcodeOverlap/"

barcodes = as_tibble(read.table(paste0(homeDirectory, 'stepThreeStarcodeShavedReads.txt'), stringsAsFactors = F, header = T))
barcodes$UMI <- as.numeric(barcodes$UMI)
barcodesCondense <- barcodes %>% group_by(BC50StarcodeD8, SampleNum) %>% summarise(nUMI = sum(UMI))

barcodesList <- barcodesCondense$BC50StarcodeD8 %>% unique()
count_wsn <- function(sequence) {
  return(length(gregexpr("[AT]([GC])[ATGC]", sequence)[[1]]))
}
counts <- sapply(barcodesList, count_wsn)
barcodeFilterList <- names(counts[counts >= 14])

barcodesFilter <- barcodesCondense %>% ungroup(BC50StarcodeD8, SampleNum) %>% dplyr::filter(BC50StarcodeD8 %in% barcodeFilterList)
barcodesFilter <- barcodesFilter %>% group_by(SampleNum) %>% mutate(nUMIFrac = nUMI / sum(nUMI) * 10^6) %>% ungroup()

sampleTable <- barcodesFilter %>% dplyr::filter(SampleNum == "1A")
ggplot(sampleTable, aes(x = nUMIFrac)) +
  geom_histogram()

#### PLOT OVERLAP SCATTERPLOTS FOR 2 WAY HERITABILITY SPLITS ###############################################################################
sampleNames <- c("1A", "1B", "2A", "2B", "3A", "3B", "4A", "4B")
sampleList <- c(1, 3, 5, 7)
cutoffList <- c(0, 10, 50, 100, 500, 1000)

for(i in 1:length(sampleList)){
  sampleA <- barcodesFilter %>% filter(SampleNum == sampleNames[sampleList[i]]) %>% dplyr::select(BC50StarcodeD8, nUMIFrac)
  sampleB <- barcodesFilter %>% filter(SampleNum == sampleNames[sampleList[i] + 1]) %>% dplyr::select(BC50StarcodeD8, nUMIFrac)
  sampleOverlap <- full_join(sampleA, sampleB, by = "BC50StarcodeD8") %>% replace(is.na(.), 0)
  plot <- ggplot(sampleOverlap, aes(x = nUMIFrac.x, y = nUMIFrac.y)) +
    geom_point() + geom_abline(slope = 1, intercept = 1) +
    theme_classic()
  ggsave(plot = plot, file = paste0(plotDirectory, 'bcOverlapPlot_', sampleNames[sampleList[i]], "vs", sampleNames[sampleList[i] + 1], '.pdf'), width = 8, height = 8)
  
  overlapFrac <- inner_join(sampleA, sampleB, by = "BC50StarcodeD8") %>% nrow()
  grid.newpage()
  venndiagram <- draw.pairwise.venn(area1 = nrow(sampleA), area2 = nrow(sampleB), cross.area = overlapFrac,
                                    euler.d = TRUE, fill = c(alpha("gray5", 0.5), alpha("gray5", 0.5)),
                                    fontfamily = "sans", cex = 2.5, print.mode = "percent", sigdigs = 3)
  ggsave(plot = venndiagram, file = paste0(plotDirectory, 'bcOverlapVenn_', sampleNames[sampleList[i]], "vs", sampleNames[sampleList[i] + 1], '.pdf'), width = 4, height = 4)
  
  for(j in 1:length(cutoffList)) {
    sampleA <- barcodesFilter %>% filter(SampleNum == sampleNames[sampleList[i]]) %>% dplyr::select(BC50StarcodeD8, nUMI) %>% dplyr::filter(nUMI > cutoffList[[j]])
    sampleB <- barcodesFilter %>% filter(SampleNum == sampleNames[sampleList[i] + 1]) %>% dplyr::select(BC50StarcodeD8, nUMI) %>% dplyr::filter(nUMI > cutoffList[[j]])
    overlapValue <- nrow(inner_join(sampleA, sampleB, by = "BC50StarcodeD8")) / (nrow(sampleA) + nrow(sampleB) - nrow(inner_join(sampleA, sampleB, by = "BC50StarcodeD8")))
    
    overlapTableEntry <- data_frame(A = sampleNames[sampleList[i]], B = sampleNames[sampleList[i] + 1], value = overlapValue, cutoff = cutoffList[j], number = mean(nrow(sampleA), nrow(sampleB)))
    
    if(i == 1 & j == 1){
      overlapTable <- overlapTableEntry
    } else{
      overlapTable <- bind_rows(overlapTable, overlapTableEntry)
    }
  }
}

overlapTable$cond <- ifelse(overlapTable$A %in% c("1A", "2A"), "2-3", "4-5")

ggplot(overlapTable, aes(x = cond, y = value)) +
  stat_summary(fun = mean, geom = "col") +
  geom_jitter(height = 0, width = 0.1) +
  facet_grid(~cutoff) +
  ylim(0, 1) + ylab("clone barcode overlap between 2 splits") + xlab("number of divisions before splitting")
ggsave(file = paste0(plotDirectory, 'bcOverlapAcrossThresholds_2Way.pdf'), width = 4, height = 4)

#### PLOT OVERLAP SCATTERPLOTS FOR 3 WAY HERITABILITY SPLITS ###############################################################################
sampleNames <- c("5A", "5B", "5C", "6A", "6B", "6C")
sampleList <- c(1, 4)
cutoffList <- c(0, 10, 50, 100, 500, 1000)

for(i in 1:length(sampleList)){
  sampleA <- barcodesFilter %>% filter(SampleNum == sampleNames[sampleList[i]]) %>% dplyr::select(BC50StarcodeD8, nUMIFrac)
  sampleB <- barcodesFilter %>% filter(SampleNum == sampleNames[sampleList[i] + 1]) %>% dplyr::select(BC50StarcodeD8, nUMIFrac)
  sampleC <- barcodesFilter %>% filter(SampleNum == sampleNames[sampleList[i] + 2]) %>% dplyr::select(BC50StarcodeD8, nUMIFrac)
  sampleOverlap <- full_join(sampleA, sampleB, by = "BC50StarcodeD8") %>% replace(is.na(.), 0)
  plot <- ggplot(sampleOverlap, aes(x = nUMIFrac.x, y = nUMIFrac.y)) +
    geom_point() + geom_abline(slope = 1, intercept = 1) +
    theme_classic()
  ggsave(plot = plot, file = paste0(plotDirectory, 'bcOverlapPlot_', sampleNames[sampleList[i]], "vs", sampleNames[sampleList[i] + 1], '.pdf'), width = 8, height = 8)
  
  overlapFrac <- inner_join(sampleA, sampleB, by = "BC50StarcodeD8") %>% nrow()
  grid.newpage()
  venndiagram <- draw.pairwise.venn(area1 = nrow(sampleA), area2 = nrow(sampleB), cross.area = overlapFrac,
                                    euler.d = TRUE, fill = c(alpha("gray5", 0.5), alpha("gray5", 0.5)),
                                    fontfamily = "sans", cex = 2.5, print.mode = "percent", sigdigs = 3)
  ggsave(plot = venndiagram, file = paste0(plotDirectory, 'bcOverlapVenn_', sampleNames[sampleList[i]], "vs", sampleNames[sampleList[i] + 1], '.pdf'), width = 4, height = 4)
  
  for(j in 1:length(cutoffList)) {
    sampleA <- barcodesFilter %>% filter(SampleNum == sampleNames[sampleList[i]]) %>% dplyr::select(BC50StarcodeD8, nUMI) %>% dplyr::filter(nUMI > cutoffList[[j]])
    sampleB <- barcodesFilter %>% filter(SampleNum == sampleNames[sampleList[i] + 1]) %>% dplyr::select(BC50StarcodeD8, nUMI) %>% dplyr::filter(nUMI > cutoffList[[j]])
    sampleC <- barcodesFilter %>% filter(SampleNum == sampleNames[sampleList[i] + 2]) %>% dplyr::select(BC50StarcodeD8, nUMI) %>% dplyr::filter(nUMI > cutoffList[[j]])
    overlapValue1 <- nrow(inner_join(sampleA, sampleB, by = "BC50StarcodeD8")) / (nrow(sampleA) + nrow(sampleB) - nrow(inner_join(sampleA, sampleB, by = "BC50StarcodeD8")))
    overlapValue2 <- nrow(inner_join(sampleA, sampleC, by = "BC50StarcodeD8")) / (nrow(sampleA) + nrow(sampleC) - nrow(inner_join(sampleA, sampleC, by = "BC50StarcodeD8")))
    overlapValue3 <- nrow(inner_join(sampleB, sampleC, by = "BC50StarcodeD8")) / (nrow(sampleB) + nrow(sampleC) - nrow(inner_join(sampleB, sampleC, by = "BC50StarcodeD8")))
    
    overlapTableEntry <- data_frame(A = rep(sampleNames[sampleList[i]], 3),
                                    B = rep(sampleNames[sampleList[i] + 1], 3),
                                    C = rep(sampleNames[sampleList[i] + 2], 3),
                                    value = c(overlapValue1, overlapValue2, overlapValue3),
                                    cutoff = rep(cutoffList[j], 3),
                                    number = rep(mean(nrow(sampleA), nrow(sampleB)), 3))
    
    if(i == 1 & j == 1){
      overlapTable3W <- overlapTableEntry
    } else{
      overlapTable3W <- bind_rows(overlapTable3W, overlapTableEntry)
    }
  }
}

overlapTable3W$cond <- "4-5,"

ggplot(overlapTable3W, aes(x = cond, y = value)) +
  stat_summary(fun = mean, geom = "col") +
  geom_jitter(height = 0, width = 0.1) +
  facet_grid(~cutoff) +
  ylim(0, 1) + ylab("clone barcode overlap between 3 splits") + xlab("number of divisions before splitting")
ggsave(file = paste0(plotDirectory, 'bcOverlapAcrossThresholds_3Way.pdf'), width = 4, height = 4)

overlapComb <- bind_rows(overlapTable, overlapTable3W)
ggplot(overlapComb, aes(x = cond, y = value)) +
  stat_summary(fun = mean, geom = "col") +
  geom_jitter(height = 0, width = 0.1) +
  facet_grid(~cutoff) +
  ylim(0, 1) + ylab("clone barcode overlap between splits") + xlab("number of divisions before splitting") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave(file = paste0(plotDirectory, 'bcOverlapAcrossThresholds_Comb.pdf'), width = 6, height = 3)

ggplot(overlapComb, aes(x = cond, y = number)) +
  stat_summary(fun = mean, geom = "col") +
  geom_jitter(height = 0, width = 0.1) +
  facet_grid(~cutoff) +
  ylab("clone barcode overlap between splits") + xlab("number of divisions before splitting")
