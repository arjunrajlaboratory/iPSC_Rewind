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
library(ggrastr)
library(viridis)
library(VennDiagram)
library(ggalluvial)

theme_set(theme_classic())

homeDirectory <- "/Users/naveenj/Dropbox (RajLab)/Shared_Naveen/Revisions/extractedData/continuousCulture/exp1/"
plotDirectory <- "/Users/naveenj/Dropbox (RajLab)/Shared_Naveen/Revisions/plots/continuousCulture/exp1/"

barcodes = as_tibble(read.table(paste0(homeDirectory, 'stepThreeStarcodeShavedReads.txt'), stringsAsFactors = F, header = T))
barcodes$UMI <- as.numeric(barcodes$UMI)
sampleNames <- c("5", "6", "7", "8", "9", "10", "11", "12")

barcodesFibroblasts <- barcodes %>% filter(!(SampleNum %in% sampleNames)) %>%
  group_by(BC50StarcodeD8, SampleNum) %>% summarise(nUMI = sum(UMI))
barcodesFibroblasts <- barcodesFibroblasts %>% group_by(SampleNum) %>% mutate(nUMIFrac = nUMI / sum(nUMI) * n())

#### NORMALIZE COLONY READ COUNTS BY SPIKE-INS #######################################################################################
barcodesColonies <- barcodes %>% filter(SampleNum %in% sampleNames) %>%
  group_by(BC50StarcodeD8, SampleNum) %>% summarise(nUMI = sum(UMI))

spikeBC1 <- 'TCCAGGTCCTCCTACTTGTACAACACCTTGTACAGCTGCTAGTGGTAGAAGAGGTACAACAACAACACGAGCATCATGAGGATCTACAGCATCAAGAACA' %>% reverseComplement(., content = "dna", case = "as is") %>% substr(0, 50)
spikeBC2 <- 'ACGTTGTGCATGACCTTGATCACCAGCTCGATGTCGAACATCACGAGCTCGTTCTGCATCTGCAAGAACACCTCGTCCTTGAACTGCTCGACGTCCATGA' %>% reverseComplement(., content = "dna", case = "as is") %>% substr(0, 50)
nBC1 <- 20000
nBC2 <- 5000

standardTableAll = list()
lmr = list()
for (i in 1:length(sampleNames)) {
  standardTable <- filter(barcodesColonies, SampleNum == sampleNames[i]) %>% dplyr::filter(., BC50StarcodeD8 == spikeBC1 | BC50StarcodeD8 == spikeBC2) %>% arrange(-nUMI) %>% mutate(ratio = .$nUMI[1]/.$nUMI[2])
  lmr[i] = coef(lm(c(nBC1, nBC2) ~ 0 + standardTable$nUMI))
  if(is.null(dim(standardTableAll))){
    standardTableAll = standardTable
  } else {
    standardTableAll = bind_rows(standardTableAll, standardTable)
  }
}
barcodesColoniesFilter <- barcodesColonies %>% filter(BC50StarcodeD8 != spikeBC1) %>% filter(BC50StarcodeD8 != spikeBC2)

#### PLOT OVERLAP SCATTERPLOTS FOR HERITABILITY SPLITS ###############################################################################
overlapTableList <- list()
sampleList <- c(1, 3, 5, 7)
reprogEffic <- c(170, 323, 169, 83)

for(i in 1:length(sampleList)){
  sampleA <- barcodesColoniesFilter %>% filter(SampleNum == sampleNames[sampleList[i]]) %>% 
    rowwise() %>% mutate(nUMINorm = nUMI*lmr[[i]]) %>% ungroup() %>% dplyr::select(-SampleNum)
  sampleB <- barcodesColoniesFilter %>% filter(SampleNum == sampleNames[sampleList[i] + 1]) %>% 
    rowwise() %>% mutate(nUMINorm = nUMI*lmr[[i + 1]]) %>% ungroup() %>% dplyr::select(-SampleNum)
  sampleOverlap <- full_join(sampleA, sampleB, by = "BC50StarcodeD8") %>% replace(is.na(.), 0)
  overlapTableList[[i]] <- sampleOverlap
  overlapValue <- nrow(inner_join(sampleA, sampleB, by = "BC50StarcodeD8")) / (nrow(sampleA) + nrow(sampleB) - nrow(inner_join(sampleA, sampleB, by = "BC50StarcodeD8")))
  plot <- ggplot(sampleOverlap, aes(x = nUMINorm.x, y = nUMINorm.y + 1)) +
    geom_point() + geom_abline(slope = 1, intercept = 1) +
    annotate("text", x = Inf, y = Inf, label = paste0(round(overlapValue, 2)), hjust = 1, vjust = 1) +
    theme_classic()
  ggsave(plot = plot, file = paste0(plotDirectory, 'exp1_heritability_', sampleNames[sampleList[i]], "vs", sampleNames[sampleList[i] + 1], '.pdf'), width = 8, height = 8)
  
  overlapFrac <- inner_join(sampleA %>% slice_max(nUMINorm, n = reprogEffic[i]),
                            sampleB %>% slice_max(nUMINorm, n = reprogEffic[i]),
                            by = "BC50StarcodeD8") %>% nrow()
  grid.newpage()
  venndiagram <- draw.pairwise.venn(area1 = reprogEffic[i], area2 = reprogEffic[i], cross.area = overlapFrac,
                     euler.d = TRUE, fill = c(alpha("gray5", 0.5), alpha("gray5", 0.5)),
                     fontfamily = "sans", cex = 2.5, print.mode = "percent", sigdigs = 3)
  ggsave(plot = venndiagram, file = paste0(plotDirectory, 'exp1_heritabilityVenn_', sampleNames[sampleList[i]], "vs", sampleNames[sampleList[i] + 1], '.pdf'), width = 4, height = 4)
}

#### PLOT OVERLAP BAR GRAPHS FOR EACH TIME POINT #####################################################################################
cutoffLists <- list(c(170, 323, 169, 83), c(200, 200, 200, 200))
cutoffType <- c("observed", "same")

for(i in 1:length(cutoffLists)) {
cutoffList <- cutoffLists[[i]]
  
primed_p1 <- overlapTableList[[1]] %>% rowwise() %>% mutate(nUMINorm = max(nUMINorm.x, nUMINorm.y)) %>% ungroup() %>% slice_max(nUMINorm, n = cutoffList[1]) %>% dplyr::select(BC50StarcodeD8) %>% mutate(p1 = 1)
primed_p2 <- overlapTableList[[2]] %>% rowwise() %>% mutate(nUMINorm = max(nUMINorm.x, nUMINorm.y)) %>% ungroup() %>% slice_max(nUMINorm, n = cutoffList[2]) %>% dplyr::select(BC50StarcodeD8) %>% mutate(p2 = 1)
primed_p3 <- overlapTableList[[3]] %>% rowwise() %>% mutate(nUMINorm = max(nUMINorm.x, nUMINorm.y)) %>% ungroup() %>% slice_max(nUMINorm, n = cutoffList[3]) %>% dplyr::select(BC50StarcodeD8) %>% mutate(p3 = 1)
primed_p4 <- overlapTableList[[4]] %>% rowwise() %>% mutate(nUMINorm = max(nUMINorm.x, nUMINorm.y)) %>% ungroup() %>% slice_max(nUMINorm, n = cutoffList[4]) %>% dplyr::select(BC50StarcodeD8) %>% mutate(p4 = 1)

primedP1P2 <- inner_join(primed_p1, primed_p2, by = "BC50StarcodeD8")
primedP1P3 <- inner_join(primed_p1, primed_p3, by = "BC50StarcodeD8")
primedP1P4 <- inner_join(primed_p1, primed_p4, by = "BC50StarcodeD8")
primedP2P4 <- inner_join(primed_p2, primed_p4, by = "BC50StarcodeD8")
primedP3P4 <- inner_join(primed_p3, primed_p4, by = "BC50StarcodeD8")
primedP2P3 <- inner_join(primed_p2, primed_p3, by = "BC50StarcodeD8")

overlapTableReprog <- data.frame(cond = c("P1vsP2", "P1vsP3", "P1vsP4", "P2vsP4", "P3vsP4"),
                           overlap = c(nrow(primedP1P2)/(nrow(primed_p1) + nrow(primed_p2) - nrow(primedP1P2)),
                                       nrow(primedP1P3)/(nrow(primed_p1) + nrow(primed_p3) - nrow(primedP1P3)),
                                       nrow(primedP1P4)/(nrow(primed_p1) + nrow(primed_p4) - nrow(primedP1P4)),
                                       nrow(primedP2P4)/(nrow(primed_p2) + nrow(primed_p4) - nrow(primedP2P4)),
                                       nrow(primedP3P4)/(nrow(primed_p3) + nrow(primed_p4) - nrow(primedP3P4))))

ggplot(overlapTableReprog, aes(x = cond, y = overlap)) +
  geom_col() +
  theme(axis.title.x = element_blank())
ggsave(filename = paste0(plotDirectory, "primedOverlap_", cutoffType[[i]], ".pdf"), height = 2, width = 4)

saveRDS(object = overlapTableReprog, file = paste0(homeDirectory, "overlapTablePrimed_", cutoffType[[i]], ".rds"))

bcFibro_p1 <- barcodesFibroblasts %>% dplyr::filter(SampleNum == 1, nUMI > 0)
bcFibro_p2 <- barcodesFibroblasts %>% dplyr::filter(SampleNum == 2, nUMI > 0)
bcFibro_p3 <- barcodesFibroblasts %>% dplyr::filter(SampleNum == 3, nUMI > 0)
bcFibro_p4 <- barcodesFibroblasts %>% dplyr::filter(SampleNum == 4, nUMI > 0)

bcFibroP1P2 <- inner_join(bcFibro_p1, bcFibro_p2, by = "BC50StarcodeD8")
bcFibroP1P3 <- inner_join(bcFibro_p1, bcFibro_p3, by = "BC50StarcodeD8")
bcFibroP1P4 <- inner_join(bcFibro_p1, bcFibro_p4, by = "BC50StarcodeD8")
bcFibroP2P4 <- inner_join(bcFibro_p2, bcFibro_p4, by = "BC50StarcodeD8")
bcFibroP3P4 <- inner_join(bcFibro_p3, bcFibro_p4, by = "BC50StarcodeD8")

overlapTableInit <- data.frame(cond = c("P1vsP2", "P1vsP3", "P1vsP4", "P2vsP4", "P3vsP4"),
                           overlap = c(nrow(bcFibroP1P2)/(nrow(bcFibro_p1) + nrow(bcFibro_p2) - nrow(bcFibroP1P2)),
                                       nrow(bcFibroP1P3)/(nrow(bcFibro_p1) + nrow(bcFibro_p3) - nrow(bcFibroP1P3)),
                                       nrow(bcFibroP1P4)/(nrow(bcFibro_p1) + nrow(bcFibro_p4) - nrow(bcFibroP1P4)),
                                       nrow(bcFibroP2P4)/(nrow(bcFibro_p2) + nrow(bcFibro_p4) - nrow(bcFibroP2P4)),
                                       nrow(bcFibroP3P4)/(nrow(bcFibro_p3) + nrow(bcFibro_p4) - nrow(bcFibroP3P4))))
ggplot(overlapTableInit, aes(x = cond, y = overlap)) +
  geom_col() +
  theme(axis.title.x = element_blank())
ggsave(filename = paste0(plotDirectory, "baseOverlap_", cutoffType[[i]], ".pdf"), height = 2, width = 4)

saveRDS(object = overlapTableInit, file = paste0(homeDirectory, "overlapTableBase_", cutoffType[[i]], ".rds"))

overlapTableNorm <- tibble(cond = overlapTableInit$cond, overlap = overlapTableReprog$overlap / overlapTableInit$overlap)
overlapTableNorm$overlap <- overlapTableNorm$overlap / sum(overlapTableNorm$overlap)
ggplot(overlapTableNorm, aes(x = cond, y = overlap)) +
  geom_col() +
  theme(axis.title.x = element_blank())
ggsave(filename = paste0(plotDirectory, "normalizedOverlap_", cutoffType[[i]], ".pdf"), height = 2, width = 4)

saveRDS(object = overlapTableNorm, file = paste0(homeDirectory, "overlapTableNorm_", cutoffType[[i]], ".rds"))

primed_combined <- full_join(full_join(full_join(primed_p1, primed_p2, by = "BC50StarcodeD8"), primed_p3, by = "BC50StarcodeD8"), primed_p4, by = "BC50StarcodeD8")
primed_combined <- replace_na(primed_combined, replace = list(p1 = 0, p2 = 0, p3 = 0, p4 = 0))

#### PLOT CORRELATION BETWEEN POPULATION FRACTION AND COLONY SIZE ####################################################################
barcodesFibroblastsCast <- dcast(barcodesFibroblasts, id.vars = "BC50StarcodeD8", formula = BC50StarcodeD8 ~ SampleNum, value.var = "nUMIFrac")
barcodesFibroblastsCast[is.na(barcodesFibroblastsCast)] <- 0
barcodesFibroblastsMelt <- melt(barcodesFibroblastsCast, id.vars = "BC50StarcodeD8") %>% dplyr::rename(SampleNum = variable, nUMIFrac = value)

w1_primed <- barcodesFibroblastsMelt %>% dplyr::filter(BC50StarcodeD8 %in% dplyr::filter(primed_combined, p1 == 1, p2 == 0, p3 == 0, p4 == 0)$BC50StarcodeD8) %>% group_by(SampleNum) %>% summarize(nUMIFracMean = median(nUMIFrac)) %>% mutate(label = "week 1 primed")
w2_primed <- barcodesFibroblastsMelt %>% dplyr::filter(BC50StarcodeD8 %in% dplyr::filter(primed_combined, p1 == 0, p2 == 1, p3 == 0, p4 == 0)$BC50StarcodeD8) %>% group_by(SampleNum) %>% summarize(nUMIFracMean = median(nUMIFrac)) %>% mutate(label = "week 2 primed")
w3_primed <- barcodesFibroblastsMelt %>% dplyr::filter(BC50StarcodeD8 %in% dplyr::filter(primed_combined, p1 == 0, p2 == 0, p3 == 1, p4 == 0)$BC50StarcodeD8) %>% group_by(SampleNum) %>% summarize(nUMIFracMean = median(nUMIFrac)) %>% mutate(label = "week 3 primed")
w4_primed <- barcodesFibroblastsMelt %>% dplyr::filter(BC50StarcodeD8 %in% dplyr::filter(primed_combined, p1 == 0, p2 == 0, p3 == 0, p4 == 1)$BC50StarcodeD8) %>% group_by(SampleNum) %>% summarize(nUMIFracMean = median(nUMIFrac)) %>% mutate(label = "week 4 primed")
all_primed <- bind_rows(w1_primed, w2_primed, w3_primed, w4_primed)

saveRDS(object = all_primed, file = paste0(homeDirectory, "primedAbundanceCorrelation_", cutoffType[[i]], ".rds"))

ggplot(all_primed, aes(x = SampleNum, y = nUMIFracMean, group = label)) +
  geom_line() + facet_grid(~label)
ggsave(filename = paste0(plotDirectory, "cloneAbundancePlot_", cutoffType[[i]], ".pdf"), height = 2, width = 5)

#### PLOT SANKEY DIAGRAM TO DEMONSTRATE FLOW #########################################################################################
nAll = primed_combined$BC50StarcodeD8 %>% unique() %>% length()

nP1 = nrow(dplyr::filter(primed_combined, p1 == 1 & p2 == 0 & p3 == 0 & p4 == 0)) / nAll
nP2 = nrow(dplyr::filter(primed_combined, p1 == 0 & p2 == 1 & p3 == 0 & p4 == 0)) / nAll
nP3 = nrow(dplyr::filter(primed_combined, p1 == 0 & p2 == 0 & p3 == 1 & p4 == 0)) / nAll
nP4 = nrow(dplyr::filter(primed_combined, p1 == 0 & p2 == 0 & p3 == 0 & p4 == 1)) / nAll
nP1P2 = nrow(dplyr::filter(primed_combined, p1 == 1 & p2 == 1 & p3 == 0 & p4 == 0)) / nAll
nP1P3 = nrow(dplyr::filter(primed_combined, p1 == 1 & p2 == 0 & p3 == 1 & p4 == 0)) / nAll
nP1P4 = nrow(dplyr::filter(primed_combined, p1 == 1 & p2 == 0 & p3 == 0 & p4 == 1)) / nAll
nP2P3 = nrow(dplyr::filter(primed_combined, p1 == 0 & p2 == 1 & p3 == 1 & p4 == 0)) / nAll
nP2P4 = nrow(dplyr::filter(primed_combined, p1 == 0 & p2 == 1 & p3 == 0 & p4 == 1)) / nAll
nP3P4 = nrow(dplyr::filter(primed_combined, p1 == 0 & p2 == 0 & p3 == 1 & p4 == 1)) / nAll
nP1P2P3 = nrow(dplyr::filter(primed_combined, p1 == 1 & p2 == 1 & p3 == 1 & p4 == 0)) / nAll
nP1P2P4 = nrow(dplyr::filter(primed_combined, p1 == 1 & p2 == 1 & p3 == 0 & p4 == 1)) / nAll
nP1P3P4 = nrow(dplyr::filter(primed_combined, p1 == 1 & p2 == 0 & p3 == 1 & p4 == 1)) / nAll
nP2P3P4 = nrow(dplyr::filter(primed_combined, p1 == 0 & p2 == 1 & p3 == 1 & p4 == 1)) / nAll
nP1P2P3P4 = nrow(dplyr::filter(primed_combined, p1 == 1 & p2 == 1 & p3 == 1 & p4 == 1)) / nAll

onP1P2P3P4 = nP1P2P3P4 + nP1P3P4 + nP1P2P4 + nP1P4
onP1 = nP1
onP2 = nP2
onP3 = nP3
onP4 = nP4
onP1P2 = nP1P2
onP2P3 = nP2P3
onP3P4 = nP3P4
onP1P2P3 = nP1P2P3 + nP1P3
onP2P3P4 = nP2P3P4 + nP2P4

freq = c(onP1P2P3P4, onP1, onP2, onP3, onP4, onP1P2, onP2P3, onP3P4, onP1P2P3, onP2P3P4) * 100
w1 = c("primed", "primed", "nonprimed", "nonprimed", "nonprimed", "primed", "nonprimed", "nonprimed", "primed", "nonprimed")
w2 = c("primed", "nonprimed", "primed", "nonprimed", "nonprimed", "primed", "primed", "nonprimed", "primed", "primed")
w3 = c("primed", "nonprimed", "nonprimed", "primed", "nonprimed", "nonprimed", "primed", "primed", "primed", "primed")
w4 = c("primed", "nonprimed", "nonprimed", "nonprimed", "primed", "nonprimed", "nonprimed", "primed", "nonprimed", "primed")
alluvial <- data.frame(freq = freq, w1 = w1, w2 = w2, w3 = w3, w4 = w4)

ggplot(alluvial, aes(y = freq, axis1 = w1, axis2 = w2, axis3 = w3, axis4 = w4, fill = w1)) +
  geom_flow() +
  geom_stratum() +
  scale_fill_manual(values = c(nonprimed = "gray", primed = "red")) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), angle = 90) +
  ylab("percentage of clone barcodes primed at some point") +
  theme(legend.position = "none")
ggsave(filename = paste0(plotDirectory, "sankey_W1_", cutoffType[[i]], ".pdf"), height = 4, width = 4)

ggplot(alluvial, aes(y = freq, axis1 = w1, axis2 = w2, axis3 = w3, axis4 = w4, fill = w2)) +
  geom_flow() +
  geom_stratum() +
  scale_fill_manual(values = c(nonprimed = "gray", primed = "red")) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), angle = 90) +
  ylab("percentage of clone barcodes primed at some point") +
  theme(legend.position = "none")
ggsave(filename = paste0(plotDirectory, "sankey_W2_", cutoffType[[i]], ".pdf"), height = 4, width = 4)

ggplot(alluvial, aes(y = freq, axis1 = w1, axis2 = w2, axis3 = w3, axis4 = w4, fill = w3)) +
  geom_flow() +
  geom_stratum() +
  scale_fill_manual(values = c(nonprimed = "gray", primed = "red")) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), angle = 90) +
  ylab("percentage of clone barcodes primed at some point") +
  theme(legend.position = "none")
ggsave(filename = paste0(plotDirectory, "sankey_W3_", cutoffType[[i]], ".pdf"), height = 4, width = 4)

ggplot(alluvial, aes(y = freq, axis1 = w1, axis2 = w2, axis3 = w3, axis4 = w4, fill = w4)) +
  geom_flow() +
  geom_stratum() +
  scale_fill_manual(values = c(nonprimed = "gray", primed = "red")) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), angle = 90) +
  ylab("percentage of clone barcodes primed at some point") +
  theme(legend.position = "none")
ggsave(filename = paste0(plotDirectory, "sankey_W4_", cutoffType[[i]], ".pdf"), height = 4, width = 4)

freq = rep(freq, 4)
week = c(rep("week 1", 10), rep("week 2", 10), rep("week 3", 10), rep("week 4", 10))
priming = c(w1, w2, w3, w4)
alluvium = c(rep(1:10, 4))
alluvial <- data.frame(freq = freq, week = week, priming = priming, alluvium = alluvium)

ggplot(alluvial, aes(y = freq, x = week, stratum = priming, alluvium = alluvium, fill = week)) +
  geom_flow() +
  geom_stratum() +
  scale_fill_manual(values = c(nonprimed = "gray", primed = "red")) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), angle = 90) +
  ylab("percentage of clone barcodes primed at some point")
ggsave(filename = paste0(plotDirectory, "sankey_comb_", cutoffType[[i]], ".pdf"), height = 4, width = 4)
}

#### PLOT CORRELATION BETWEEN POPULATION FRACTION AND COLONY SIZE ####################################################################
cutoffList <- c(170, 333, 169, 83)
cutoffList <- c(200, 200, 200, 200)

primed_p1 <- overlapTableList[[1]] %>% rowwise() %>% mutate(nUMINorm = max(nUMINorm.x, nUMINorm.y)) %>% ungroup() %>% slice_max(nUMINorm, n = cutoffList[1]) %>% dplyr::select(BC50StarcodeD8, nUMINorm) %>% mutate(split = 1)
primed_p2 <- overlapTableList[[2]] %>% rowwise() %>% mutate(nUMINorm = max(nUMINorm.x, nUMINorm.y)) %>% ungroup() %>% slice_max(nUMINorm, n = cutoffList[2]) %>% dplyr::select(BC50StarcodeD8, nUMINorm) %>% mutate(split = 2)
primed_p3 <- overlapTableList[[3]] %>% rowwise() %>% mutate(nUMINorm = max(nUMINorm.x, nUMINorm.y)) %>% ungroup() %>% slice_max(nUMINorm, n = cutoffList[3]) %>% dplyr::select(BC50StarcodeD8, nUMINorm) %>% mutate(split = 3)
primed_p4 <- overlapTableList[[4]] %>% rowwise() %>% mutate(nUMINorm = max(nUMINorm.x, nUMINorm.y)) %>% ungroup() %>% slice_max(nUMINorm, n = cutoffList[4]) %>% dplyr::select(BC50StarcodeD8, nUMINorm) %>% mutate(split = 4)

barcodes_p1 <- barcodesFibroblasts %>% dplyr::filter(SampleNum == 1) %>% mutate(priming = "nonprimed")
barcodes_p1$priming <- ifelse(barcodes_p1$BC50StarcodeD8 %in% primed_p1$BC50StarcodeD8, "primed", barcodes_p1$priming)
ggplot(barcodes_p1, aes(x = priming, y = nUMIFrac)) +
  rasterize(geom_jitter(), dpi = 300) +
  geom_boxplot(outlier.shape = NA, width = 0.25) +
  theme(axis.title.x = element_blank()) + ylab("normalized relative abundance per barcode")
ggsave(filename = paste0(plotDirectory, "W1_abundanceVersusColonySizeBoxplot.pdf"), height = 4, width = 4)
overlap_p1 <- left_join(primed_p1, barcodes_p1, by = "BC50StarcodeD8")
ggplot(overlap_p1, aes(x = nUMIFrac, y = nUMINorm)) +
  rasterize(geom_point(), dpi = 300) +
  geom_smooth(method = "lm", se = FALSE) + stat_poly_eq(use_label(c("R2"))) +
  xlab("normalized relative abundance per barcode in initial population") + ylab("normalized cells per barcode in reprogrammed population")
ggsave(filename = paste0(plotDirectory, "W1_abundanceVersusColonySize.pdf"), height = 4, width = 4)

barcodes_p2 <- barcodesFibroblasts %>% dplyr::filter(SampleNum == 2) %>% mutate(priming = "nonprimed")
barcodes_p2$priming <- ifelse(barcodes_p2$BC50StarcodeD8 %in% primed_p2$BC50StarcodeD8, "primed", barcodes_p2$priming)
ggplot(barcodes_p2, aes(x = priming, y = nUMIFrac)) +
  rasterize(geom_jitter(), dpi = 300) +
  geom_boxplot(outlier.shape = NA, width = 0.25) +
  theme(axis.title.x = element_blank()) + ylab("normalized relative abundance per barcode")
ggsave(filename = paste0(plotDirectory, "W2_abundanceVersusColonySizeBoxplot.pdf"), height = 4, width = 4)
overlap_p2 <- left_join(primed_p2, barcodes_p2, by = "BC50StarcodeD8")
ggplot(overlap_p2, aes(x = nUMIFrac, y = nUMINorm)) +
  rasterize(geom_point(), dpi = 300) +
  geom_smooth(method = "lm", se = FALSE) + stat_poly_eq(use_label(c("R2"))) +
  xlab("normalized relative abundance per barcode in initial population") + ylab("normalized cells per barcode in reprogrammed population")
ggsave(filename = paste0(plotDirectory, "W2_abundanceVersusColonySize.pdf"), height = 4, width = 4)

barcodes_p3 <- barcodesFibroblasts %>% dplyr::filter(SampleNum == 3) %>% mutate(priming = "nonprimed")
barcodes_p3$priming <- ifelse(barcodes_p3$BC50StarcodeD8 %in% primed_p3$BC50StarcodeD8, "primed", barcodes_p3$priming)
ggplot(barcodes_p3, aes(x = priming, y = nUMIFrac)) +
  rasterize(geom_jitter(), dpi = 300) +
  geom_boxplot(outlier.shape = NA, width = 0.25) +
  theme(axis.title.x = element_blank()) + ylab("normalized relative abundance per barcode")
ggsave(filename = paste0(plotDirectory, "W3_abundanceVersusColonySizeBoxplot.pdf"), height = 4, width = 4)
overlap_p3 <- left_join(primed_p3, barcodes_p3, by = "BC50StarcodeD8")
ggplot(overlap_p3, aes(x = nUMIFrac, y = nUMINorm)) +
  rasterize(geom_point(), dpi = 300) +
  geom_smooth(method = "lm", se = FALSE) + stat_poly_eq(use_label(c("R2"))) +
  xlab("normalized relative abundance per barcode in initial population") + ylab("normalized cells per barcode in reprogrammed population")
ggsave(filename = paste0(plotDirectory, "W3_abundanceVersusColonySize.pdf"), height = 4, width = 4)

barcodes_p4 <- barcodesFibroblasts %>% dplyr::filter(SampleNum == 4) %>% mutate(priming = "nonprimed")
barcodes_p4$priming <- ifelse(barcodes_p4$BC50StarcodeD8 %in% primed_p4$BC50StarcodeD8, "primed", barcodes_p4$priming)
ggplot(barcodes_p4, aes(x = priming, y = nUMIFrac)) +
  rasterize(geom_jitter(), dpi = 300) +
  geom_boxplot(outlier.shape = NA, width = 0.25) +
  theme(axis.title.x = element_blank()) + ylab("normalized relative abundance per barcode")
ggsave(filename = paste0(plotDirectory, "W4_abundanceVersusColonySizeBoxplot.pdf"), height = 4, width = 4)
overlap_p4 <- left_join(primed_p4, barcodes_p4, by = "BC50StarcodeD8")
ggplot(overlap_p4, aes(x = nUMIFrac, y = nUMINorm)) +
  rasterize(geom_point(), dpi = 300) +
  geom_smooth(method = "lm", se = FALSE) + stat_poly_eq(use_label(c("R2"))) +
  xlab("normalized relative abundance per barcode in initial population") + ylab("normalized cells per barcode in reprogrammed population")
ggsave(filename = paste0(plotDirectory, "W4_abundanceVersusColonySize.pdf"), height = 4, width = 4)

#### PLOT HISTOGRAM OF COLONY SIZES FOR EACH TIME POINT ##############################################################################
primed_p1 <- overlapTableList[[1]] %>% rowwise() %>% mutate(nUMINorm = max(nUMINorm.x, nUMINorm.y)) %>% mutate(week = "1")
primed_p2 <- overlapTableList[[2]] %>% rowwise() %>% mutate(nUMINorm = max(nUMINorm.x, nUMINorm.y)) %>% mutate(week = "2")
primed_p3 <- overlapTableList[[3]] %>% rowwise() %>% mutate(nUMINorm = max(nUMINorm.x, nUMINorm.y)) %>% mutate(week = "3")
primed_p4 <- overlapTableList[[4]] %>% rowwise() %>% mutate(nUMINorm = max(nUMINorm.x, nUMINorm.y)) %>% mutate(week = "4")
primed_all <- bind_rows(primed_p1, primed_p2, primed_p3, primed_p4)

ggplot(primed_all, aes(x = nUMINorm, fill = week)) +
  geom_histogram(alpha = 0.5, position = "identity", binwidth = 100) +
  xlim(200, 5000) + xlab("normalized reads per barcode")