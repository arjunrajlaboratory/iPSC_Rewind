rm(list=ls())
gc()

library(tidyverse)
library(spgs)
library(Seurat)
library(ggrepel)
library(reshape2)
library(ggrastr)

theme_set(theme_classic())

homeDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Original Manuscript/extractedData/fateMap10X/FM3/"
plotDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Original Manuscript/plots/fateMap10X/FM3/"

linCountToOverlaps <- as_tibble(read.table(file = paste0(homeDirectory, "filtered10XCells.txt"), header = TRUE, stringsAsFactors = F, sep = "\t"))
barcodes <- as_tibble(read.table(file = "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Original Manuscript/extractedData/fateMap10X/FM3/stepThreeStarcodeShavedReads_BC_10XAndGDNA.txt", header = TRUE, stringsAsFactors = F, sep = "\t"))

#### normalize gDNA barcode reads by spike-ins ####
#######################################################################################################################################################
probedBarcodes <- barcodes %>% filter(cellID == "dummy") %>% dplyr::select(UMI, BC50StarcodeD8, SampleNum) %>% rename(barcode = BC50StarcodeD8)
probedBarcodes$UMI <- as.numeric(probedBarcodes$UMI)
probedBarcodes <- probedBarcodes %>% group_by(barcode, SampleNum) %>% summarise(nUMI = sum(UMI))

standardTableAll = list()
lmr = list()

spikeBC1 <- 'TCCAGGTCCTCCTACTTGTACAACACCTTGTACAGCTGCTAGTGGTAGAAGAGGTACAACAACAACACGAGCATCATGAGGATCTACAGCATCAAGAACA' %>% reverseComplement(., content = "dna", case = "as is") %>% substr(0,48) %>% paste0("AT", .)
spikeBC2 <- 'ACGTTGTGCATGACCTTGATCACCAGCTCGATGTCGAACATCACGAGCTCGTTCTGCATCTGCAAGAACACCTCGTCCTTGAACTGCTCGACGTCCATGA' %>% reverseComplement(., content = "dna", case = "as is") %>% substr(0,48) %>% paste0("AT", .)
nBC1 <- 20000
nBC2 <- 5000

sampleNames <- probedBarcodes$SampleNum %>% unique()
standardTableAll = list()
lmr = list()
for (i in 1:length(sampleNames)) {
  standardTable <- filter(probedBarcodes, SampleNum == sampleNames[i]) %>% dplyr::filter(., barcode == spikeBC1 | barcode == spikeBC2) %>% arrange(-nUMI) %>% mutate(ratio = .$nUMI[1]/.$nUMI[2])
  lmr[i] = coef(lm(c(nBC1, nBC2) ~ 0 + standardTable$nUMI))
  if(is.null(dim(standardTableAll))){
    standardTableAll = standardTable
  } else {
    standardTableAll = bind_rows(standardTableAll, standardTable)
  }
}

probedBarcodesFilter <- probedBarcodes %>% filter(barcode != spikeBC1) %>% filter(barcode != spikeBC2)

#### compare within condition scatterplots and across condition scatterplots ####
#######################################################################################################################################################
sample_A <- probedBarcodesFilter %>% dplyr::filter(SampleNum == sampleNames[1]) %>%
  rowwise() %>% mutate(nUMINorm = nUMI*lmr[[1]]) %>% ungroup() %>% dplyr::select(-SampleNum)
sample_B <- probedBarcodesFilter %>% dplyr::filter(SampleNum == sampleNames[2]) %>%
  rowwise() %>% mutate(nUMINorm = nUMI*lmr[[2]]) %>% ungroup() %>% dplyr::select(-SampleNum)
barcodesDMSO_1 <- full_join(sample_A, sample_B, by = "barcode") %>% replace(is.na(.), 0) %>% mutate(condition = "DMSO_1")
ggplot(barcodesDMSO_1, aes(x = nUMINorm.x, y = nUMINorm.y)) +
  geom_point() + geom_abline(slope = 1, intercept = 0)
rm(sample_A)
rm(sample_B)

sample_A <- probedBarcodesFilter %>% dplyr::filter(SampleNum == sampleNames[9]) %>%
  rowwise() %>% mutate(nUMINorm = nUMI*lmr[[9]]) %>% ungroup() %>% dplyr::select(-SampleNum)
sample_B <- probedBarcodesFilter %>% dplyr::filter(SampleNum == sampleNames[4]) %>%
  rowwise() %>% mutate(nUMINorm = nUMI*lmr[[4]]) %>% ungroup() %>% dplyr::select(-SampleNum)
barcodesDMSO_2 <- full_join(sample_A, sample_B, by = "barcode") %>% replace(is.na(.), 0) %>% mutate(condition = "DMSO_2")
ggplot(barcodesDMSO_2, aes(x = nUMINorm.x, y = nUMINorm.y)) +
  geom_point() + geom_abline(slope = 1, intercept = 0)
rm(sample_A)
rm(sample_B)

barcodesDMSO <- bind_rows(barcodesDMSO_1 %>% rowwise() %>% mutate(nUMIMax = max(nUMINorm.x, nUMINorm.y)) %>% dplyr::select(-nUMI.x, -nUMI.y, -nUMINorm.x, -nUMINorm.y),
                          barcodesDMSO_2 %>% rowwise() %>% mutate(nUMIMax = max(nUMINorm.x, nUMINorm.y)) %>% dplyr::select(-nUMI.x, -nUMI.y, -nUMINorm.x, -nUMINorm.y))

sample_A <- probedBarcodesFilter %>% dplyr::filter(SampleNum == sampleNames[6]) %>%
  rowwise() %>% mutate(nUMINorm = nUMI*lmr[[6]]) %>% ungroup() %>% dplyr::select(-SampleNum)
sample_B <- probedBarcodesFilter %>% dplyr::filter(SampleNum == sampleNames[10]) %>%
  rowwise() %>% mutate(nUMINorm = nUMI*lmr[[10]]) %>% ungroup() %>% dplyr::select(-SampleNum)
barcodesLSD1_1 <- full_join(sample_A, sample_B, by = "barcode") %>% replace(is.na(.), 0) %>% mutate(condition = "LSD1_1")
ggplot(barcodesLSD1_1, aes(x = nUMINorm.x, y = nUMINorm.y)) +
  geom_point() + geom_abline(slope = 1, intercept = 0)
rm(sample_A)
rm(sample_B)

sample_A <- probedBarcodesFilter %>% dplyr::filter(SampleNum == sampleNames[5]) %>%
  rowwise() %>% mutate(nUMINorm = nUMI*lmr[[5]]) %>% ungroup() %>% dplyr::select(-SampleNum)
sample_B <- probedBarcodesFilter %>% dplyr::filter(SampleNum == sampleNames[7]) %>%
  rowwise() %>% mutate(nUMINorm = nUMI*lmr[[7]]) %>% ungroup() %>% dplyr::select(-SampleNum)
barcodesLSD1_2 <- full_join(sample_A, sample_B, by = "barcode") %>% replace(is.na(.), 0) %>% mutate(condition = "LSD1_2")
ggplot(barcodesLSD1_2, aes(x = nUMINorm.x, y = nUMINorm.y)) +
  geom_point() + geom_abline(slope = 1, intercept = 0)
rm(sample_A)
rm(sample_B)

barcodesLSD1 <- bind_rows(barcodesLSD1_1 %>% rowwise() %>% mutate(nUMIMax = max(nUMINorm.x, nUMINorm.y)) %>% dplyr::select(-nUMI.x, -nUMI.y, -nUMINorm.x, -nUMINorm.y),
                          barcodesLSD1_2 %>% rowwise() %>% mutate(nUMIMax = max(nUMINorm.x, nUMINorm.y)) %>% dplyr::select(-nUMI.x, -nUMI.y, -nUMINorm.x, -nUMINorm.y))

sample_A <- probedBarcodesFilter %>% dplyr::filter(SampleNum == sampleNames[11]) %>%
  rowwise() %>% mutate(nUMINorm = nUMI*lmr[[11]]) %>% ungroup() %>% dplyr::select(-SampleNum)
sample_B <- probedBarcodesFilter %>% dplyr::filter(SampleNum == sampleNames[12]) %>%
  rowwise() %>% mutate(nUMINorm = nUMI*lmr[[12]]) %>% ungroup() %>% dplyr::select(-SampleNum)
barcodesDOT1L_1 <- full_join(sample_A, sample_B, by = "barcode") %>% replace(is.na(.), 0) %>% mutate(condition = "DOT1L_1")
ggplot(barcodesDMSO_1, aes(x = nUMINorm.x, y = nUMINorm.y)) +
  geom_point() + geom_abline(slope = 1, intercept = 0)
rm(sample_A)
rm(sample_B)

sample_A <- probedBarcodesFilter %>% dplyr::filter(SampleNum == sampleNames[3]) %>%
  rowwise() %>% mutate(nUMINorm = nUMI*lmr[[3]]) %>% ungroup() %>% dplyr::select(-SampleNum)
sample_B <- probedBarcodesFilter %>% dplyr::filter(SampleNum == sampleNames[8]) %>%
  rowwise() %>% mutate(nUMINorm = nUMI*lmr[[8]]) %>% ungroup() %>% dplyr::select(-SampleNum)
barcodesDOT1L_2 <- full_join(sample_A, sample_B, by = "barcode") %>% replace(is.na(.), 0) %>% mutate(condition = "DOT1L_2")
ggplot(barcodesDMSO_1, aes(x = nUMINorm.x, y = nUMINorm.y)) +
  geom_point() + geom_abline(slope = 1, intercept = 0)
rm(sample_A)
rm(sample_B)

barcodesDOT1L <- bind_rows(barcodesDOT1L_1 %>% rowwise() %>% mutate(nUMIMax = max(nUMINorm.x, nUMINorm.y)) %>% dplyr::select(-nUMI.x, -nUMI.y, -nUMINorm.x, -nUMINorm.y),
                           barcodesDOT1L_2 %>% rowwise() %>% mutate(nUMIMax = max(nUMINorm.x, nUMINorm.y)) %>% dplyr::select(-nUMI.x, -nUMI.y, -nUMINorm.x, -nUMINorm.y))

overlapLSD1 <- full_join(barcodesDMSO %>% dplyr::select(-condition), barcodesLSD1 %>% dplyr::select(-condition), by = "barcode") %>% replace(is.na(.), 0)
overlapLSD1 <- overlapLSD1 %>% rowwise() %>% mutate(foldchange = log((nUMIMax.y + 1)/(nUMIMax.x + 1)))
ggplot() +
  geom_point(data = overlapLSD1, aes(x = nUMIMax.x, y = nUMIMax.y)) +
  geom_abline(slope = 1, intercept = 0) +
  xlab("normalized cell number in DMSO") + ylab("normalized cell number in LSD1i")

overlapDOT1L <- full_join(barcodesDMSO %>% dplyr::select(-condition), barcodesDOT1L %>% dplyr::select(-condition), by = "barcode") %>% replace(is.na(.), 0)
overlapDOT1L <- overlapDOT1L %>% rowwise() %>% mutate(foldchange = log((nUMIMax.y + 1)/(nUMIMax.x + 1)))
ggplot() +
  geom_point(data = overlapDOT1L, aes(x = nUMIMax.x, y = nUMIMax.y)) +
  geom_abline(slope = 1, intercept = 0) +
  xlab("normalized cell number in DMSO") + ylab("normalized cell number in DOT1Li")

#######################################################################################################################################################
#### DON'T NEED TO RUN AGAIN ####
barcodesOverlap <- overlapLSD1
cutoff = c(100, 250, 500, 1000, 2500, 5000, 10000, 25000)
foldChangecutoffDep = 2
foldChangecutoffInd = 2

ratioDepXDepY = c()
ratioDepXInd = c()
ratioDepYInd = c()
ratioIndDepYIndDepX = c()
depX = c()
depY = c()
ind = c()
for (i in 1:length(cutoff)) {
  drugDepX <- barcodesOverlap %>% filter(foldchange < -foldChangecutoffDep | foldchange > foldChangecutoffDep) %>% filter(nUMIMax.x > cutoff[i] & nUMIMax.y < cutoff[i])
  drugDepY <- barcodesOverlap %>% filter(foldchange < -foldChangecutoffDep | foldchange > foldChangecutoffDep) %>% filter(nUMIMax.y > cutoff[i] & nUMIMax.x < cutoff[i])
  drugInd <- barcodesOverlap %>% filter(foldchange > -foldChangecutoffInd & foldchange < foldChangecutoffInd) %>% filter(nUMIMax.x > cutoff[i] & nUMIMax.y > cutoff[i])
  
  depX[i] = nrow(drugDepX)
  depY[i] = nrow(drugDepY)
  ind[i] = nrow(drugInd)
  
  ratioDepXDepY[i] = nrow(drugDepY)/nrow(drugDepX)
  ratioDepXInd[i] = nrow(drugInd)/nrow(drugDepX)
  ratioDepYInd[i] = nrow(drugInd)/nrow(drugDepY)
  ratioIndDepYIndDepX[i] = (nrow(drugInd) + nrow(drugDepY))/(nrow(drugInd) + nrow(drugDepX))
  
  plot <- ggplot() +
    geom_point(aes(x=barcodesOverlap$nUMIMax.x,y=barcodesOverlap$nUMIMax.y), alpha = 0.4, size = 4, shape = 16) +
    geom_point(aes(x=drugInd$nUMIMax.x,y=drugInd$nUMIMax.y), col="orange", size = 4, shape = 16) +
    geom_point(aes(x=drugDepY$nUMIMax.x,y=drugDepY$nUMIMax.y), col="springgreen3", size = 4, shape = 16) +
    geom_point(aes(x=drugDepX$nUMIMax.x,y=drugDepX$nUMIMax.y), col="turquoise3", size = 4, shape = 16) +
    geom_abline(slope=1,intercept=0) +
    geom_vline(xintercept = cutoff[i]) +
    geom_hline(yintercept = cutoff[i]) +
    xlab(paste0("DMSO")) + ylab(paste0("LSD1i"))
  ggsave(plot = plot, file = paste0(plotDirectory, 'DMSOvsLSD1i/FM3_DMSOvsLSD1i_cutoff_', cutoff[i] ,'.pdf'), width = 4, height = 4)
}

barcodesOverlap <- overlapDOT1L
cutoff = c(100, 250, 500, 1000, 2500, 5000, 10000, 25000)
foldChangecutoffDep = 2
foldChangecutoffInd = 2

ratioDepXDepY = c()
ratioDepXInd = c()
ratioDepYInd = c()
ratioIndDepYIndDepX = c()
depX = c()
depY = c()
ind = c()
for (i in 1:length(cutoff)) {
  drugDepX <- barcodesOverlap %>% filter(foldchange < -foldChangecutoffDep | foldchange > foldChangecutoffDep) %>% filter(nUMIMax.x > cutoff[i] & nUMIMax.y < cutoff[i])
  drugDepY <- barcodesOverlap %>% filter(foldchange < -foldChangecutoffDep | foldchange > foldChangecutoffDep) %>% filter(nUMIMax.y > cutoff[i] & nUMIMax.x < cutoff[i])
  drugInd <- barcodesOverlap %>% filter(foldchange > -foldChangecutoffInd & foldchange < foldChangecutoffInd) %>% filter(nUMIMax.x > cutoff[i] & nUMIMax.y > cutoff[i])
  
  depX[i] = nrow(drugDepX)
  depY[i] = nrow(drugDepY)
  ind[i] = nrow(drugInd)
  
  ratioDepXDepY[i] = nrow(drugDepY)/nrow(drugDepX)
  ratioDepXInd[i] = nrow(drugInd)/nrow(drugDepX)
  ratioDepYInd[i] = nrow(drugInd)/nrow(drugDepY)
  ratioIndDepYIndDepX[i] = (nrow(drugInd) + nrow(drugDepY))/(nrow(drugInd) + nrow(drugDepX))
  
  plot <- ggplot() +
    geom_point(aes(x=barcodesOverlap$nUMIMax.x,y=barcodesOverlap$nUMIMax.y), alpha = 0.4, size = 4, shape = 16) +
    geom_point(aes(x=drugInd$nUMIMax.x,y=drugInd$nUMIMax.y), col="orange", size = 4, shape = 16) +
    geom_point(aes(x=drugDepY$nUMIMax.x,y=drugDepY$nUMIMax.y), col="springgreen3", size = 4, shape = 16) +
    geom_point(aes(x=drugDepX$nUMIMax.x,y=drugDepX$nUMIMax.y), col="turquoise3", size = 4, shape = 16) +
    geom_abline(slope=1,intercept=0) +
    geom_vline(xintercept = cutoff[i]) +
    geom_hline(yintercept = cutoff[i]) +
    xlab(paste0("DMSO")) + ylab(paste0("LSD1i"))
  ggsave(plot = plot, file = paste0(plotDirectory, 'DMSOvsLSD1i/FM3_DMSOvsDOT1Li_cutoff_', cutoff[i] ,'.pdf'), width = 4, height = 4)
}
#### DON'T NEED TO RUN AGAIN ####
#######################################################################################################################################################
#### plotting scatterplot and Venn diagram for LSD1i independent versus dependent ####
#######################################################################################################################################################
umapCoordinates = as_tibble(read.table(file = paste0(homeDirectory, "umapCoordinates_Scanorama_50pcs_filter.tsv"), header = TRUE, stringsAsFactors = F, sep = "\t"))

barcodesOverlap <- overlapLSD1
cutoff = 5000
foldChangecutoffDep = 2
foldChangecutoffInd = 2

drugDepX <- barcodesOverlap %>% filter(foldchange < -foldChangecutoffDep | foldchange > foldChangecutoffDep) %>% filter(nUMIMax.x > cutoff & nUMIMax.y < cutoff)
drugDepY <- barcodesOverlap %>% filter(foldchange < -foldChangecutoffDep | foldchange > foldChangecutoffDep) %>% filter(nUMIMax.y > cutoff & nUMIMax.x < cutoff)
drugInd <- barcodesOverlap %>% filter(nUMIMax.x > cutoff & nUMIMax.y > cutoff)

ggplot() +
  rasterise(geom_point(aes(x=barcodesOverlap$nUMIMax.x,y=barcodesOverlap$nUMIMax.y), color = "lightgray", shape = 16), dpi = 150) +
  geom_point(aes(x=drugInd$nUMIMax.x,y=drugInd$nUMIMax.y), col="orange", shape = 16) +
  geom_point(aes(x=drugDepY$nUMIMax.x,y=drugDepY$nUMIMax.y), col="springgreen3", shape = 16) +
  geom_point(aes(x=drugDepX$nUMIMax.x,y=drugDepX$nUMIMax.y), col="turquoise3", shape = 16) +
  geom_abline(slope=1, intercept=0, linetype = "dashed") +
  # geom_vline(xintercept = cutoff, linetype = "dashed") +
  # geom_hline(yintercept = cutoff, linetype = "dashed") +
  xlab(paste0("DMSO")) + ylab(paste0("LSD1i")) + coord_fixed(ratio = 1, xlim = c(0, 70000), ylim = c(0, 70000))
ggsave(filename = paste0(plotDirectory, 'barcodeScatterPlot_DMSOvsLSD1_10X.pdf'), width = 4, height = 4, useDingbats = F)

library(VennDiagram)
grid.newpage()
venndiagram <- draw.pairwise.venn(area1 = nrow(drugInd) + nrow(drugDepX), area2 = nrow(drugInd) + nrow(drugDepY), cross.area = nrow(drugInd),
                                  euler.d = TRUE, fill = c(alpha("gray5", 0.5), alpha("gray5", 0.5)),
                                  fontfamily = "sans", cex = 2.5, print.mode = "percent", sigdigs = 3)
ggsave(plot = venndiagram, filename = paste0(plotDirectory, 'overlapVennDiagram_LSD1_10X.pdf'), width = 5, height = 5, useDingbats = F)

#### analysis comparing DNA-sequencing and single-cell RNA-sequencing ####
#######################################################################################################################################################
umapCoordinates <- umapCoordinates %>% mutate(label = "DMSO")
umapCoordinates$label <- ifelse(umapCoordinates$sampleNum %in% c("S3", "S4"), "LSD1", umapCoordinates$label)
umapCoordinates$label <- ifelse(umapCoordinates$sampleNum %in% c("S5", "S6"), "DOT1L", umapCoordinates$label)

labelProp <- umapCoordinates %>% group_by(label) %>% summarize(nCells = n())

barcodeAllProp <- labelProp <- inner_join(umapCoordinates, linCountToOverlaps, by = "cellID") %>% dplyr::filter(!is.na(barcode)) %>%
  group_by(label) %>% mutate(nCellsCond = n()) %>% ungroup() %>%
  group_by(barcode, label, nCellsCond) %>% summarize(nCells = n()) %>% ungroup() %>%
  mutate(prop = nCells/nCellsCond)
barcodeAllPropCast <- dcast(barcodeAllProp, formula = barcode ~ label) %>% replace(is.na(.), 0)
ggplot(barcodeAllPropCast, aes(x = DMSO + 0.000001, y = LSD1 + 0.000001)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10")

drugDepXCellID <- inner_join(drugDepX, linCountToOverlaps, by = "barcode") %>% dplyr::select(cellID) %>% unique() %>% filter(cellID %in% linCountToOverlaps$cellID)
drugDepYCellID <- inner_join(drugDepY, linCountToOverlaps, by = "barcode") %>% dplyr::select(cellID) %>% unique() %>% filter(cellID %in% linCountToOverlaps$cellID)
drugIndCellID <- inner_join(drugInd, linCountToOverlaps, by = "barcode") %>% dplyr::select(cellID) %>% unique() %>% filter(cellID %in% linCountToOverlaps$cellID)

drugDepXUMAP <- filter(umapCoordinates, cellID %in% drugDepXCellID$cellID) %>% dplyr::filter(label != "DOT1L")
drugDepYUMAP <- filter(umapCoordinates, cellID %in% drugDepYCellID$cellID) %>% dplyr::filter(label != "DOT1L")
drugIndUMAP <- filter(umapCoordinates, cellID %in% drugIndCellID$cellID) %>% dplyr::filter(label != "DOT1L")

drugDepXUMAP %>% group_by(label) %>% summarize(nCells = n())
drugDepYUMAP %>% group_by(label) %>% summarize(nCells = n())
drugIndUMAP %>% group_by(label) %>% summarize(nCells = n())

drugDepXUMAP <- left_join(drugDepXUMAP, linCountToOverlaps %>% dplyr::select(cellID, barcode), by = "cellID")
drugDepYUMAP <- left_join(drugDepYUMAP, linCountToOverlaps %>% dplyr::select(cellID, barcode), by = "cellID")
drugIndUMAP <- left_join(drugIndUMAP, linCountToOverlaps %>% dplyr::select(cellID, barcode), by = "cellID")

drugDepX_barcodeProp <- drugDepXUMAP %>% left_join(., labelProp, by = c("label", "barcode")) %>% rename(nCellsCondJoin = nCellsCond) %>%
  group_by(barcode, label, nCellsCondJoin) %>% summarize(nCells = n()) %>% ungroup() %>%
  mutate(prop = nCells/nCellsCondJoin)
drugDepX_barcodePropCast <- dcast(drugDepX_barcodeProp, formula = barcode ~ label) %>% replace(is.na(.), 0) %>% mutate(cond = "depX")

drugDepY_barcodeProp <- drugDepYUMAP %>% left_join(., labelProp, by = c("label", "barcode")) %>% rename(nCellsCondJoin = nCellsCond) %>%
  group_by(barcode, label, nCellsCondJoin) %>% summarize(nCells = n()) %>% ungroup() %>%
  mutate(prop = nCells/nCellsCondJoin)
drugDepY_barcodePropCast <- dcast(drugDepY_barcodeProp, formula = barcode ~ label) %>% replace(is.na(.), 0) %>% mutate(cond = "depY")

drugInd_barcodeProp <- drugIndUMAP %>% left_join(., labelProp, by = c("label", "barcode")) %>% rename(nCellsCondJoin = nCellsCond) %>%
  group_by(barcode, label, nCellsCondJoin) %>% summarize(nCells = n()) %>% ungroup() %>%
  mutate(prop = nCells/nCellsCondJoin)
drugInd_barcodePropCast <- dcast(drugInd_barcodeProp, formula = barcode ~ label) %>% replace(is.na(.), 0) %>% mutate(cond = "ind")

barcodePropDrugAll <- bind_rows(drugDepX_barcodePropCast, drugDepY_barcodePropCast, drugInd_barcodePropCast)
barcodePropDrugComb <- full_join(barcodePropDrugAll, dplyr::select(barcodeAllPropCast, -DOT1L), by = c("barcode", "DMSO", "LSD1"))
barcodePropDrugComb[is.na(barcodePropDrugComb)] <- "none"
barcodePropDrugComb$cond <- factor(barcodePropDrugComb$cond, levels = c("none", "depX", "depY", "ind"))

ggplot(barcodePropDrugAll, aes(x = DMSO + 0.000001, y = LSD1 + 0.000001, color = cond)) +
  geom_jitter(width = 0.1, height = 0.1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  coord_fixed(ratio = 1) +
  scale_x_continuous(trans = "log10", limits = c(0.000001, 1)) +
  scale_y_continuous(trans = "log10", limits = c(0.000001, 1)) +
  scale_color_manual(values = c("turquoise3", "springgreen3", "orange")) + NoLegend()
ggsave(filename = paste0(plotDirectory, 'barcodeScatterPlot_DMSOvsLSD1_10X_Overlap.pdf'), width = 4, height = 4, useDingbats = F)

drugDepYCorrBarcode <- drugDepY_barcodePropCast %>% dplyr::filter(DMSO == 0)
drugIndCorrBarcode <- drugInd_barcodePropCast %>% dplyr::filter(DMSO != 0)

drugDepYCorr <- drugDepY %>% dplyr::filter(barcode %in% drugDepYCorrBarcode$barcode)
drugIndCorr <- drugInd %>% dplyr::filter(barcode %in% drugIndCorrBarcode$barcode)

#### look at differentially expressed genes in bulk across independent and dependent lineages ####
#######################################################################################################################################################
drugDepY <- drugDepYCorr
drugInd <- drugIndCorr

drugDepXCellID <- inner_join(drugDepX, linCountToOverlaps, by = "barcode") %>% dplyr::select(cellID) %>% unique() %>% filter(cellID %in% linCountToOverlaps$cellID)
drugDepYCellID <- inner_join(drugDepY, linCountToOverlaps, by = "barcode") %>% dplyr::select(cellID) %>% unique() %>% filter(cellID %in% linCountToOverlaps$cellID)
drugIndCellID <- inner_join(drugInd, linCountToOverlaps, by = "barcode") %>% dplyr::select(cellID) %>% unique() %>% filter(cellID %in% linCountToOverlaps$cellID)

# drugDepXUMAP <- filter(umapCoordinates, cellID %in% drugDepXCellID$cellID) %>% dplyr::filter(label != "DOT1L")
# drugDepYUMAP <- filter(umapCoordinates, cellID %in% drugDepYCellID$cellID) %>% dplyr::filter(label != "DOT1L")
# drugIndUMAP <- filter(umapCoordinates, cellID %in% drugIndCellID$cellID) %>% dplyr::filter(label != "DOT1L")
# 
# plot1 <- ggplot() +
#   geom_point(data = umapCoordinates, aes(x = UMAP_1, y = UMAP_2), shape = 16) +
#   geom_point(data = drugDepXUMAP, aes(x = UMAP_1, y = UMAP_2), shape = 16, size = 2.5, color = "turquoise3") +
#   facet_wrap(~sampleNum)
# plot2 <- ggplot() +
#   geom_point(data = umapCoordinates, aes(x = UMAP_1, y = UMAP_2), shape = 16) +
#   geom_point(data = drugDepYUMAP, aes(x = UMAP_1, y = UMAP_2), shape = 16, size = 2.5, color = "springgreen3") +
#   facet_wrap(~sampleNum)
# plot3 <- ggplot() +
#   geom_point(data = umapCoordinates, aes(x = UMAP_1, y = UMAP_2), shape = 16) +
#   geom_point(data = drugIndUMAP, aes(x = UMAP_1, y = UMAP_2), shape = 16, size = 2.5, color = "orange") +
#   facet_wrap(~sampleNum)
# egg::ggarrange(plot1, plot2, plot3, ncol = 3)

umapCoordinates = as_tibble(read.table(file = paste0(homeDirectory, "umapCoordinates_Scanorama_50pcs_filter.tsv"), header = TRUE, stringsAsFactors = F, sep = "\t"))
umapClusters <- as_tibble(read.table(file = paste0(homeDirectory, "umapClusters_Scanorama_50pcs_filter.tsv"), header = TRUE, stringsAsFactors = F, sep = "\t"))
umapClusters <- umapClusters %>% dplyr::rename(cluster = 1)
umapCoordinates <- inner_join(umapCoordinates, umapClusters, by = c("cellID", "sampleNum"))

labelsToAdd <- umapCoordinates %>% mutate(label = "none")
labelsToAdd$label <- ifelse(labelsToAdd$cellID %in% (linCountToOverlaps %>% .$cellID) &
                              labelsToAdd$sampleNum %in% c("S1", "S2") &
                              labelsToAdd$cluster %in% c("0", "1", "2", "3"), "barcoded iPSC", labelsToAdd$label)
labelsToAdd$label <- ifelse(labelsToAdd$cellID %in% drugDepYCellID$cellID &
                              labelsToAdd$sampleNum %in% c("S3", "S4") &
                              labelsToAdd$cluster %in% c("0", "1", "2", "3"), "LSD1i dep", labelsToAdd$label)
labelsToAdd$label <- ifelse(labelsToAdd$cellID %in% drugIndCellID$cellID &
                              labelsToAdd$sampleNum %in% c("S3", "S4") &
                              labelsToAdd$cluster %in% c("0", "1", "2", "3"), "LSD1i ind", labelsToAdd$label)
labelsToAddFinal <- labelsToAdd %>% dplyr::select(label)

scanorama_filter <- readRDS(paste0(homeDirectory, "scanorama_filter.rds"))
scanorama_filter_labels <- AddMetaData(scanorama_filter, metadata = labelsToAddFinal)
scanorama_filter_labels$label <- labelsToAddFinal$label
Idents(scanorama_filter_labels) <- scanorama_filter_labels$label
DimPlot(scanorama_filter_labels, split.by = "label")

#### detour for PCA plot for aggregate expression ####
scanorama_filter_labels_PCA <- subset(scanorama_filter_labels, subset = label %in% c("LSD1i ind", "LSD1i dep"))
DimPlot(scanorama_filter_labels_PCA, reduction = "pca_scanorama")

labelsToAdd <- umapCoordinates %>% mutate(label = "nonbarcoded")
labelsToAdd$label <- ifelse(labelsToAdd$cellID %in% (linCountToOverlaps %>% .$cellID) &
                              labelsToAdd$cluster %in% c("4"), "fibroblast", labelsToAdd$label)
labelsToAdd$label <- ifelse(labelsToAdd$cellID %in% (linCountToOverlaps %>% .$cellID) &
                              labelsToAdd$cluster %in% c("5"), "incomplete", labelsToAdd$label)
labelsToAdd$label <- ifelse(labelsToAdd$cellID %in% (linCountToOverlaps %>% .$cellID) &
                              labelsToAdd$sampleNum %in% c("S1", "S2") &
                              labelsToAdd$cluster %in% c("0", "1", "2", "3"), "DMSO iPSC", labelsToAdd$label)
labelsToAdd$label <- ifelse(labelsToAdd$cellID %in% (linCountToOverlaps %>% .$cellID) &
                              labelsToAdd$sampleNum %in% c("S3", "S4") &
                              labelsToAdd$cluster %in% c("0", "1", "2", "3"), "LSD1i iPSC", labelsToAdd$label)
labelsToAdd$label <- ifelse(labelsToAdd$cellID %in% drugDepYCellID$cellID &
                              labelsToAdd$sampleNum %in% c("S3", "S4") &
                              labelsToAdd$cluster %in% c("0", "1", "2", "3"), "LSD1i dep", labelsToAdd$label)
labelsToAdd$label <- ifelse(labelsToAdd$cellID %in% drugIndCellID$cellID &
                              labelsToAdd$sampleNum %in% c("S3", "S4") &
                              labelsToAdd$cluster %in% c("0", "1", "2", "3"), "LSD1i ind", labelsToAdd$label)
labelsToAddFinal <- labelsToAdd %>% dplyr::select(label)

scanorama_filter_labels <- AddMetaData(scanorama_filter, metadata = labelsToAddFinal)
scanorama_filter_labels$label <- labelsToAddFinal$label
Idents(scanorama_filter_labels) <- scanorama_filter_labels$label
DimPlot(scanorama_filter_labels, split.by = "label")

scanorama_aggregate <- scanorama_filter_labels %>% AggregateExpression(return.seurat = TRUE)
scanorama_aggregate <- NormalizeData(scanorama_aggregate)
scanorama_aggregate <- FindVariableFeatures(scanorama_aggregate, nfeatures = 8000)
all.genes <- rownames(scanorama_aggregate)
scanorama_aggregate <- ScaleData(scanorama_aggregate, features = all.genes)
scanorama_aggregate <- RunPCA(scanorama_aggregate, npcs = 4, approx = FALSE)
DimPlot(scanorama_aggregate, reduction = "pca", dims = c(1, 2), pt.size = 5)

# scanorama_aggregate_subset <- subset(scanorama_aggregate, idents = c("fibroblast", "incomplete", "DMSO iPSC", "LSD1i iPSC"))
scanorama_aggregate_subset <- subset(scanorama_aggregate, idents = c("fibroblast", "incomplete", "DMSO iPSC", "LSD1i dep", "LSD1i ind"))
scanorama_aggregate_subset <- NormalizeData(scanorama_aggregate_subset)
scanorama_aggregate_subset <- FindVariableFeatures(scanorama_aggregate_subset, nfeatures = 8000)
all.genes <- rownames(scanorama_aggregate_subset)
scanorama_aggregate_subset <- ScaleData(scanorama_aggregate_subset, features = all.genes)
scanorama_aggregate_subset <- RunPCA(scanorama_aggregate_subset, npcs = 3, approx = FALSE)
DimPlot(scanorama_aggregate_subset, reduction = "pca", pt.size = 5) & NoLegend()
ggsave(filename = paste0(plotDirectory, "FM3_PCAFinalState_PC1vsPC2.pdf"), height = 3, width = 3)
pcLoadings <- as_tibble(scanorama_aggregate_subset[["pca"]][, 1:2])
pcLoadings$gene <- rownames(scanorama_aggregate_subset[["pca"]][, 1:2])
pcLoadingsMelt <- melt(pcLoadings, id.vars = "gene")

pcLoadingsMelt <- pcLoadingsMelt %>% dplyr::filter(variable == "PC_1")

ggplot() +
  rasterise(geom_point(data = pcLoadingsMelt, aes(x = variable, y = value), color = "lightgray", position = position_jitter(seed = 1234)), dpi = 300) +
  geom_point(data = pcLoadingsMelt %>% dplyr::filter(gene %in% genesPluri),
             aes(x = variable, y = value), color = "red", position = position_jitter(seed = 1234)) +
  geom_text_repel(data = pcLoadingsMelt %>% dplyr::filter(gene %in% genesPluri),
                  aes(x = variable, y = value, label = gene), color = "red", position = position_jitter(seed = 1234)) +
  # geom_point(data = pcLoadingsMelt %>% dplyr::filter(gene %in% genesFibro),
  #            aes(x = variable, y = value), color = "blue", position = position_jitter(seed = 1234)) +
  # geom_text_repel(data = pcLoadingsMelt %>% dplyr::filter(gene %in% genesFibro),
  #                 aes(x = variable, y = value, label = gene), color = "blue", position = position_jitter(seed = 1234)) +
  geom_hline(yintercept = 0, linetype = "dashed")
ggsave(filename = paste0(plotDirectory, "FM3_PCALoadings_PC1.pdf"), height = 4, width = 4)

mat <- Seurat::GetAssayData(scanorama_aggregate_subset, assay = "scanorama", slot = "scale.data")
pca <- scanorama_aggregate_subset[["pca"]]

total_variance <- sum(matrixStats::rowVars(mat))

eigValues = (pca@stdev)^2  ## EigenValues
varExplained = eigValues / total_variance

ggplot(pcLoadingsMelt, aes(x = variable, y = value)) +
  geom_jitter()
######################################################

DefaultAssay(scanorama_filter_labels) <- "RNA"
scanorama_filter_labels <- NormalizeData(scanorama_filter_labels)

markersLSD1IndVsDep <- FindMarkers(object = scanorama_filter_labels, ident.1 = "LSD1i ind", ident.2 = "LSD1i dep", logfc.threshold = 0)
markersLSD1IndVsDep <- markersLSD1IndVsDep %>% mutate(gene = rownames(.))

# markersLSD1Ind <- FindMarkers(object = scanorama_filter_labels, ident.1 = "LSD1i ind", ident.2 = "barcoded iPSC", logfc.threshold = 0.2)
# markersLSD1Ind <- markersLSD1Ind %>% mutate(gene = rownames(.))
# markersLSD1Dep <- FindMarkers(object = scanorama_filter_labels, ident.1 = "LSD1i dep", ident.2 = "barcoded iPSC", logfc.threshold = 0.2)
# markersLSD1Dep <- markersLSD1Dep %>% mutate(gene = rownames(.))

# library(gprofiler2)
# gostres <- gost(query = markersDMSOvsLSD1 %>% dplyr::filter(avg_log2FC > 0.5) %>% arrange(., desc(abs(avg_log2FC))) %>% .$gene,
#                 organism = "hsapiens", ordered_query = TRUE, sources = "GO:BP", evcodes = TRUE)
# gostres <- gost(query = markersDMSOvsLSD1 %>% dplyr::filter(avg_log2FC < -0.5) %>% arrange(., desc(abs(avg_log2FC))) %>% .$gene,
#                 organism = "hsapiens", ordered_query = TRUE, sources = "GO:BP", evcodes = TRUE)
# res <- gostres$result
# 
naiveSCMarkers <- c("KLF2", "ZPF42", "ESRRB", "DPPA3", "TFCP2L1", "FGF4", "TBX3", "CDH1")
primedSCMarkers <- c("CD24", "OTX2", "GDF3", "DNMT3B", "FGF5", "MEIS1", "POU3F1", "CER1")

# markersScatterPlot <- full_join(markersLSD1Ind %>% dplyr::select(gene, avg_log2FC), markersLSD1Dep %>% dplyr::select(gene, avg_log2FC), by = "gene") %>% replace(is.na(.), 0)
# # saveRDS(markersScatterPlot, file = paste0(plotDirectory, "markersScatterPlotLSD1iIndVsDep.pdf"))
# markersScatterPlot <- readRDS(file = paste0(plotDirectory, "markersScatterPlotLSD1iIndVsDep.pdf"))
# 
# ggplot(markersScatterPlot, aes(x = avg_log2FC.x, y = avg_log2FC.y)) +
#   rasterise(geom_point(), dpi = 150) + geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
#   geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
#   geom_rect(aes(xmin = -0.5, xmax = 0.5, ymin = -0.5, ymax = 0.5), linetype = "dashed", fill = NA, color = "black") +
#   geom_point(data = markersScatterPlot %>% dplyr::filter(gene %in% genesPluri), color = "red") +
#   geom_text_repel(data = markersScatterPlot %>% dplyr::filter(gene %in% genesPluri), aes(label = gene), color = "red") 
# 
# ggplot(markersScatterPlot, aes(x = avg_log2FC.x, y = avg_log2FC.y)) +
#   rasterise(geom_point(), dpi = 150) + geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
#   geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
#   geom_rect(aes(xmin = -0.5, xmax = 0.5, ymin = -0.5, ymax = 0.5), linetype = "dashed", fill = NA, color = "black") +
#   geom_point(data = markersScatterPlot %>% dplyr::filter(gene %in% primedSCMarkers), color = "red") +
#   geom_text_repel(data = markersScatterPlot %>% dplyr::filter(gene %in% primedSCMarkers), aes(label = gene), color = "red") +
#   geom_point(data = markersScatterPlot %>% dplyr::filter(gene %in% naiveSCMarkers), color = "blue") +
#   geom_text_repel(data = markersScatterPlot %>% dplyr::filter(gene %in% naiveSCMarkers), aes(label = gene), color = "blue") 
# 
# ggplot(markersScatterPlot, aes(x = avg_log2FC.x, y = avg_log2FC.y)) +
#   rasterise(geom_point(), dpi = 150) + geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
#   geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
#   geom_rect(aes(xmin = -0.5, xmax = 0.5, ymin = -0.5, ymax = 0.5), linetype = "dashed", fill = NA, color = "black") +
#   geom_point(data = markersScatterPlot %>% dplyr::filter(gene %in% genesNeg), color = "blue") +
#   geom_text_repel(data = markersScatterPlot %>% dplyr::filter(gene %in% genesNeg), aes(label = gene), color = "blue")
# 
# ggplot(markersScatterPlot, aes(x = avg_log2FC.x, y = avg_log2FC.y)) +
#   rasterise(geom_point(), dpi = 150) + geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
#   geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
#   geom_rect(aes(xmin = -0.5, xmax = 0.5, ymin = -0.5, ymax = 0.5), linetype = "dashed", fill = NA, color = "black") +
#   geom_text_repel(data = markersScatterPlot %>% dplyr::filter(abs(avg_log2FC.x) > 0.5 | abs(avg_log2FC.y > 0.5)), aes(label = gene))

#### plots incoporating previous marker categories
# genesPluri <- c("NANOG", "ALPL", "POU5F1", "SOX2", "PODXL", "LIN28A", "UTF1", "TERT", "ZFP42", "DNMT3B")
# genesFibro <- c("LUM", "S100A4", "THY1", "PDGFRA", "COL1A1", "COL5A1", "LOXL1", "FBLN1", "FBLN2", "VTN")
# genesEpi <- c("CDH1", "CLDN3", "KRT3", "OCLN", "EPCAM", "ANPEP", "MUC1", "CD24")
# genesMes <- c("CDH2", "VIM", "FN1", "ZEB1", "SNAI2", "TWIST1", "TWIST2", "TGFB1")
# genesMyo <- c("ACTA2", "TPX2", "TAGLN", "MYL9", "CDKN1A", "CDKN2A", "CCN2", "GLI2")

markersPrev <- readRDS("/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/extractedData/rewind10X/markersComb.rds")
genesPluri <- c("NANOG", "ALPL", "POU5F1", "SOX2", "PODXL", "LIN28A", "UTF1", "TERT", "ZFP42", "DNMT3B")
genesFibro <- c("LUM", "S100A4", "THY1", "PDGFRA", "COL1A1", "COL1A2", "LOXL1", "FBLN1", "FBLN2", "VTN")
genesEpi <- c("CDH1", "CLDN3", "KRT3", "OCLN", "EPCAM", "ANPEP", "MUC1", "CD24")
genesMes <- c("CDH2", "VIM", "ENTPD1", "SNAI1", "SNAI2", "TWIST1", "TWIST2", "TGFB1")
genesActFibro <- c("FAP", "TAGLN", "TNC", "SPARC", "LGALS1", "SERPINE2", "FBN1", "THBS1", "SHOX2", "TBX2")
genesMyoFibro <- c("ACTA2", "TPM2", "POSTN", "IGB1", "MMP1", "MMP14", "FN1", "AOC3", "NKX2-3", "LRRC17")

genesPos <- markersPrev %>% slice_max(., order_by = mean, n = 25) %>% .$gene
genesNeg <- markersPrev %>% slice_min(., order_by = mean, n = 25) %>% .$gene

# genesMarkerList <- list(genesPluri, genesFibro, genesEpi, genesMes, genesActFibro, genesMyoFibro)
# genesNamesList <- c("pluri", "fibro", "epi", "mes", "act", "myo")

genesMarkerList <- list(genesPluri, genesFibro)
genesNamesList <- c("pluri", "fibro")

# markersPrevious <- genesPluri
# 
# VlnPlot(scanorama_filter_labels, features = markersPrevious, ncol = 4) &
#   stat_summary(fun = mean, geom = "bar", color = "black") &
#   stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.25)

# markersLSD1IndVsDep <- FindMarkers(object = scanorama_filter_labels, ident.1 = "LSD1i ind", ident.2 = "LSD1i dep", logfc.threshold = 0.2)
# markersLSD1IndVsDep <- markersLSD1IndVsDep %>% mutate(gene = rownames(.))
markersLSD1IndVsDep <- markersLSD1IndVsDep %>% arrange(., avg_log2FC) %>% mutate(number = seq.int(nrow(markersLSD1IndVsDep)))
markersLSD1IndVsDep <- markersLSD1IndVsDep %>% mutate(markers = "none")
for(i in 1:length(genesMarkerList)) {
  markersLSD1IndVsDep$markers <- ifelse(markersLSD1IndVsDep$gene %in% genesMarkerList[[i]], genesNamesList[i], markersLSD1IndVsDep$markers)
}
# markersLSD1IndVsDep$markers <- factor(markersLSD1IndVsDep$markers, levels = c("none", "mes", "myo", "act", "fibro", "epi", "pluri"))
markersLSD1IndVsDep$markers <- factor(markersLSD1IndVsDep$markers, levels = c("none", "pluri", "fibro"))
markersLSD1IndVsDep$label <- ifelse(markersLSD1IndVsDep$markers == "none", "", markersLSD1IndVsDep$gene)
markersLSD1IndVsDep_bulk <- markersLSD1IndVsDep

# markersLSD1Ind <- FindMarkers(object = scanorama_filter_labels, ident.1 = "LSD1i ind", ident.2 = "barcoded iPSC", logfc.threshold = 0.2)
# markersLSD1Ind <- markersLSD1Ind %>% mutate(gene = rownames(.))
# markersLSD1Ind <- markersLSD1Ind %>% arrange(., avg_log2FC) %>% mutate(number = seq.int(nrow(markersLSD1Ind)))
# markersLSD1Ind <- markersLSD1Ind %>% mutate(markers = "none")
# for(i in 1:length(genesMarkerList)) {
#   markersLSD1Ind$markers <- ifelse(markersLSD1Ind$gene %in% genesMarkerList[[i]], genesNamesList[i], markersLSD1Ind$markers)
# }
# markersLSD1Ind$markers <- factor(markersLSD1Ind$markers, levels = c("none", "mes", "myo", "fibro", "epi", "pluri"))
# markersLSD1Ind$label <- ifelse(markersLSD1Ind$markers == "none", "", markersLSD1Ind$gene)
# 
# markersLSD1Dep <- FindMarkers(object = scanorama_filter_labels, ident.1 = "LSD1i dep", ident.2 = "barcoded iPSC", logfc.threshold = 0.2)
# markersLSD1Dep <- markersLSD1Dep %>% mutate(gene = rownames(.))
# markersLSD1Dep <- markersLSD1Dep %>% arrange(., avg_log2FC) %>% mutate(number = seq.int(nrow(markersLSD1Dep)))
# markersLSD1Dep <- markersLSD1Dep %>% mutate(markers = "none")
# for(i in 1:length(genesMarkerList)) {
#   markersLSD1Dep$markers <- ifelse(markersLSD1Dep$gene %in% genesMarkerList[[i]], genesNamesList[i], markersLSD1Dep$markers)
# }
# markersLSD1Dep$markers <- factor(markersLSD1Dep$markers, levels = c("none", "mes", "myo", "fibro", "epi", "pluri"))
# markersLSD1Dep$label <- ifelse(markersLSD1Dep$markers == "none", "", markersLSD1Dep$gene)

ggplot(markersLSD1IndVsDep, aes(x = markers, y = avg_log2FC, fill = markers)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(seed = 1234)) +
  geom_text_repel(aes(label = label), position = position_jitter(seed = 1234)) +
  geom_hline(yintercept = 0, linetype = "dashed") + NoLegend()

# ggplot(markersLSD1Ind, aes(x = markers, y = avg_log2FC, label = label, fill = markers)) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_point(position = position_jitter(seed = 1234)) +
#   geom_text_repel(position = position_jitter(seed = 1234)) +
#   geom_hline(yintercept = 0, linetype = "dashed") + NoLegend()
# 
# ggplot(markersLSD1Dep, aes(x = markers, y = avg_log2FC, label = label, fill = markers)) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_point(position = position_jitter(seed = 1234)) +
#   geom_text_repel(position = position_jitter(seed = 1234)) +
#   geom_hline(yintercept = 0, linetype = "dashed") + NoLegend()

#### calculate neighbor matrix between LSD1 independent and LSD1 dependent ####
#######################################################################################################################################################
library(RANN)
library(tripack)
library(svglite)

umapCoordinates = as_tibble(read.table(file = paste0(homeDirectory, "umapCoordinates_Scanorama_50pcs_filter.tsv"), header = TRUE, stringsAsFactors = F, sep = "\t"))
umapClusters <- as_tibble(read.table(file = paste0(homeDirectory, "umapClusters_Scanorama_50pcs_filter.tsv"), header = TRUE, stringsAsFactors = F, sep = "\t"))
umapClusters <- umapClusters %>% dplyr::rename(clusters = 1)
umapCoordClust <- inner_join(umapCoordinates, umapClusters, by = c("cellID", "sampleNum"))

pcaCoordinates = as_tibble(read.table(file = paste0(homeDirectory, "pcaCoordinates_Scanorama_50pcs_filter.tsv"), header = TRUE, stringsAsFactors = F, sep = "\t"))

drugDepYCellID <- inner_join(drugDepY, linCountToOverlaps %>% dplyr::filter(SampleNum %in% c("S5", "S6")), by = "barcode") %>% dplyr::select(cellID) %>% unique() %>% filter(cellID %in% linCountToOverlaps$cellID)
drugIndCellID <- inner_join(drugInd, linCountToOverlaps %>% dplyr::filter(SampleNum %in% c("S5", "S6")), by = "barcode") %>% dplyr::select(cellID) %>% unique() %>% filter(cellID %in% linCountToOverlaps$cellID)

umapCoordinatesIndVsDep <- umapCoordClust %>% dplyr::filter(sampleNum %in% c("S3", "S4"), clusters %in% c("0", "1", "2", "3", "4", "5"), UMAP_2 < 10)
drugDepYUMAP <- dplyr::filter(umapCoordinatesIndVsDep, cellID %in% drugDepYCellID$cellID)
drugIndUMAP <- dplyr::filter(umapCoordinatesIndVsDep, cellID %in% drugIndCellID$cellID)

ggplot() +
  rasterise(geom_point(data = umapCoordinatesIndVsDep, aes(x = UMAP_1, y = UMAP_2), shape = 16, color = "lightgray"), dpi = 150) +
  geom_point(data = drugDepYUMAP, aes(x = UMAP_1, y = UMAP_2), shape = 16, size = 2.5, color = "springgreen3") + NoAxes()
ggsave(filename = paste0(plotDirectory, "UMAP_LSD1InhibInd.pdf"), height = 4, width = 4)

ggplot() +
  rasterise(geom_point(data = umapCoordinatesIndVsDep, aes(x = UMAP_1, y = UMAP_2), shape = 16, color = "lightgray"), dpi = 150) +
  geom_point(data = drugIndUMAP, aes(x = UMAP_1, y = UMAP_2), shape = 16, size = 2.5, color = "orange") + NoAxes()
ggsave(filename = paste0(plotDirectory, "UMAP_LSD1InhibDep.pdf"), height = 4, width = 4)

mixingCoeffList <- c()
for(n in 1:100) {
  iBarcodeBx = pcaCoordinates %>% dplyr::filter(cellID %in% drugDepYCellID$cellID) %>% mutate(label = "dep")
  jBarcodeBx = pcaCoordinates %>% dplyr::filter(cellID %in% drugIndCellID$cellID) %>% sample_n(., size = nrow(iBarcodeBx), seed = n) %>% mutate(label = "ind")
  ijBarcodeBx = bind_rows(iBarcodeBx, jBarcodeBx)
  BarcodeRef = ijBarcodeBx %>% mutate(num = c(1:nrow(ijBarcodeBx)))
  S1index = BarcodeRef %>% filter(label == "dep") %>% dplyr::select(num)
  S2index = BarcodeRef %>% filter(label == "ind") %>% dplyr::select(num)
  
  knnPCA = nn2(ijBarcodeBx[, 1:50], ijBarcodeBx[, 1:50])
  neighborKNN = as_tibble(knnPCA[[1]])
  nS1 = c();
  nS2 = c();
  fraction2 = (nrow(S2index)/nrow(S1index))
  for (k in c(1:nrow(ijBarcodeBx))) {
    nS1[k] =  (sum(neighborKNN[k,2:10] %in% S1index$num))*(nrow(S2index)/nrow(S1index))
    nS2[k] =  sum(neighborKNN[k,2:10] %in% S2index$num)
  }
  neighbor = tibble(nS1,nS2) %>% mutate(fractionS1 = nS1/(nS1+nS2), fractionS2= nS2/(nS1+nS2)) %>% mutate(num = c(1:nrow(ijBarcodeBx)))
  neighbor = inner_join(neighbor, BarcodeRef, by = "num")
  neighborS1 = tibble(nS1,nS2) %>% mutate(fractionS1 = nS1/(nS1+nS2), fractionS2= nS2/(nS1+nS2)) %>% mutate(num = c(1:nrow(ijBarcodeBx))) %>% filter(num %in% S1index$num)
  neighborS2 = tibble(nS1,nS2) %>% mutate(fractionS1 = nS1/(nS1+nS2), fractionS2= nS2/(nS1+nS2)) %>% mutate(num = c(1:nrow(ijBarcodeBx))) %>% filter(num %in% S2index$num)
  S1S1 = as_tibble(neighborS1$fractionS1) %>% mutate(isSelf = "withSelf") 
  S2S2 = as_tibble(neighborS2$fractionS2) %>% mutate(isSelf = "withSelf")
  S1S2 = as_tibble((neighborS1$fractionS2)) %>% mutate(isSelf ="withOther")
  S2S1 = as_tibble((neighborS2$fractionS1)) %>% mutate(isSelf = "withOther")
  toPlot = bind_rows(S1S1,S2S2,S1S2,S2S1) %>% rename(fraction = value) %>% mutate(isSelf = factor(isSelf, levels=c("withSelf","withOther")))
  withSelf = toPlot %>% filter(isSelf == "withSelf") %>% dplyr::select(-isSelf)
  withOther = toPlot %>% filter(isSelf == "withOther") %>% dplyr::select(-isSelf)
  mixingCoeffList[n] = mean(withOther$fraction)/mean(withSelf$fraction)
}

ggplot(data = NULL, mapping = aes(x = "", y = mixingCoeffList)) +
  stat_summary(fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.25)

ggplot(data = NULL, mapping = aes(x = "", y = mixingCoeffList)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.25) + ylim(0, 1.5)
ggsave(filename = paste0(plotDirectory, 'mixingCoeff_LSD1IndVsDep.pdf'), width = 3, height = 3, useDingbats = F)

library(DescTools)
KLD <- function(x,y) sum(x * log2(x/y))
JSD <- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
mixedrank = function(x) order(gtools::mixedorder(x))

umapTable <- as_tibble(read.table(file = paste0(homeDirectory, "umapCoordinates_Scanorama_50pcs_filter.tsv"), header = TRUE, stringsAsFactors = F, sep = "\t"))
clusterTable <- as_tibble(read.table(file = paste0(homeDirectory, "umapClusters_Scanorama_50pcs_filter.tsv"), header = TRUE, stringsAsFactors = F, sep = "\t"))
filteredBarcodeTable <- as_tibble(read.table(file = paste0(homeDirectory, "filtered10XCells.txt"), header = TRUE, stringsAsFactors = F, sep = "\t"))

jointUMAP = umapCoordClust %>% inner_join(linCountToOverlaps, ., by = c("cellID")) %>% dplyr::select(-SampleNum, -nUMI, -fracUMI, -nLineages)
jointUMAPFilter <- jointUMAP %>% dplyr::filter(!(clusters %in% c("4", "5", "7", "8"))) %>% dplyr::filter(UMAP_2 < 10) %>% dplyr::filter(sampleNum %in% c("S3", "S4"))

clusterList <- as.tibble(clusters = c("0", "1", "2", "3"))
nCellsPerBarcodeFilter <- drugDepY %>% .$barcode

jointUMAPFilterDepY <- jointUMAPFilter %>% dplyr::filter(barcode %in% drugDepY$barcode)
jointUMAPFilterInd <- jointUMAPFilter %>% dplyr::filter(barcode %in% drugInd$barcode)

ggplot() +
  geom_point(data = jointUMAPFilter, aes(x = UMAP_1, y = UMAP_2), shape = 16, color = "lightgray") +
  geom_point(data = jointUMAPFilterDepY, aes(x = UMAP_1, y = UMAP_2), shape = 16, size = 2.5, color = "springgreen3")

ggplot() +
  geom_point(data = jointUMAPFilter, aes(x = UMAP_1, y = UMAP_2), shape = 16, color = "lightgray") +
  geom_point(data = jointUMAPFilterInd, aes(x = UMAP_1, y = UMAP_2), shape = 16, size = 2.5, color = "orange")

sibACells <- jointUMAPFilter %>% dplyr::filter(barcode %in% nCellsPerBarcodeFilter)
barPerClustA <- jointUMAPFilterInd %>% dplyr::select(clusters) %>% group_by(clusters) %>% summarise(allBarPerClust = n())
sibACellsDistTable <- sibACells %>% group_by(clusters) %>% summarise(cellsPerCluster = n()) %>% ungroup() %>% right_join(., barPerClustA, by = "clusters") %>%
  replace(., is.na(.), 0) %>%
  mutate(propBarPerClust = cellsPerCluster / allBarPerClust) %>% ungroup() %>%
  mutate(prop = propBarPerClust / sum(propBarPerClust))

n = 100
seeds = 1:n
JSTest = rep(0, n)
JS_A = rep(0, nrow(jointUMAPFilterDepY))
meanJS_A = rep(0, nrow(jointUMAPFilterDepY))
sdJS_A = rep(0, nrow(jointUMAPFilterDepY))
setwd(paste0(plotDirectory, "constraint_LSD1iIndVsDep/"))

bartitle = "dep vs. ind"
barcodeclusters = sibACellsDistTable
allrandomprop = barPerClustA
#generate random samples
for (j in seeds) {
  j = 1
  set.seed(j)
  randomsibA = suppressMessages(
    jointUMAPFilterInd %>% sample_n(size = nrow(jointUMAPFilterDepY)) %>%
      dplyr::select(clusters) %>%
      group_by(clusters) %>%
      summarise(randomcellspercluster = n()) %>%
      full_join(barPerClustA, by = c('clusters')) %>%
      mutate_if(is.numeric, ~ replace(., is.na(.), 0)) %>%
      mutate(randomprop = randomcellspercluster /
               allBarPerClust) %>%
      mutate(randomprop = randomprop / sum(randomprop)) %>%
      arrange(clusters)
  )
  allrandomprop[, paste0("randomprop", j)] = randomsibA$randomprop
}
#for entropy analysis
allrandomprop = allrandomprop %>% rowwise() %>%
  mutate(avgrandomprop = mean(c_across(starts_with("randomprop"))))
summaryrandomprop = allrandomprop %>%
  inner_join(barcodeclusters, by = 'clusters') %>%
  dplyr::select(clusters, prop, avgrandomprop)
#calculate entropy values
p_observed = summaryrandomprop %>% ungroup() %>%
  arrange(clusters) %>%
  summarise(prob_observed = prop + 0.000001) #add pseudocount to avoid zero in the numerator and/or denominator of KL
p_random = summaryrandomprop %>% ungroup() %>%
  arrange(clusters) %>%
  summarise(prob_random = avgrandomprop + 0.000001)

observedentropy = Entropy(p_observed, base = 2)
randomentropy = Entropy(p_random, base = 2)
JS_A[i] <- JSD(p_random, p_observed)

plotsummaryrandomprop = summaryrandomprop %>% transmute(clusters, observed_prop = prop, randomavg_prop = avgrandomprop) %>%
  pivot_longer(cols = -clusters, names_to = c("sample", ".value"), names_pattern = "(.+)_(.+)")
sample.labs <- c(paste0("observed ", bartitle, ": entropy = ", round(observedentropy, 3)),
                 paste0("random sample: entropy = ", round(randomentropy, 3)))
names(sample.labs) <- c("observed", "randomavg")
entropyplot_prop = ggplot(plotsummaryrandomprop, aes(x = as.factor(clusters), y = prop, fill = sample)) +
  geom_bar(stat = 'identity') +
  facet_wrap( ~ sample, ncol = 1, labeller = labeller(sample = sample.labs)) +
  scale_fill_manual(breaks = c("observed", "randomavg"), values = c("#009292", "gray80")) +
  ggtitle(paste(bartitle, ': ', nrow(jointUMAPFilterDepY), ' cells\n', 'Jensen-Shannon distance = ', round(JS_A[i], 3), sep = "")) +
  xlab("Seurat clusters, res=0.5") +
  ylab(paste0("normalized cell proportion")) +
  theme_classic() +
  theme(legend.position = "none")
pdf(paste(bartitle, "entropyPlot.pdf", sep="_"), width = 4, height = 4)
print(entropyplot_prop)
dev.off()

#evaluate JSD significance
for (k in seeds + n) {
  set.seed(k)
  JSDsibA = jointUMAPFilter %>% sample_n(size = nrow(jointUMAPFilterDepY)) %>%
    dplyr::select(clusters) %>%
    full_join(barPerClustA, by = 'clusters') %>%
    group_by(clusters, allBarPerClust) %>%
    summarise(cellspercluster = n(), .groups = 'drop_last') %>%
    mutate(propperbarclust = cellspercluster / allBarPerClust) %>%
    ungroup() %>%
    mutate(prop = propperbarclust / sum(propperbarclust))
  p_test = JSDsibA %>% ungroup() %>%
    arrange(clusters) %>%
    summarise(prob_test = prop + 0.000001)
  JSTest[k - n] <- JSD(p_random, p_test)
}
meanJS_A[i] = mean(JSTest)
sdJS_A[i] = sd(JSTest)
normp = pnorm(JS_A[i], mean = meanJS_A[i], sd = sdJS_A[i], lower.tail = FALSE)
if (round(normp, 3) == 0) {
  normp3 <- "< 0.001"
} else{
  normp3 <- paste("=", round(normp, 3))
}
JSplot = ggplot() +
  geom_histogram(aes(x = JSTest), binwidth = 0.01, fill = 'gray80') +
  geom_vline(xintercept = JS_A[i], color = '#009292') +
  ggtitle(paste0(bartitle, ': ', nrow(jointUMAPFilterDepY), ' cells', '\np ', normp3)) +
  xlab("Jensen-Shannon distance") +
  coord_cartesian(xlim = c(0, 1)) + #without coord_cartesian throws an error about rows being removed
  theme_classic()
pdf(paste(bartitle, "JSDPlot.pdf", sep="_"), width = 4, height = 4)
print(JSplot)
dev.off()

#### calculate DE of pluripotency markers by pairwise comparisons of twins ####
#######################################################################################################################################################
# logNormalizedCounts = scanorama_filter_labels@assays$RNA@data
# logNormalizedCountsRound = round(logNormalizedCounts, 4)
# cells_count = sub("-1", "", colnames(logNormalizedCountsRound))
# cells_count_cellID = sub("S\\d_", "", cells_count)
# cells_count_sampleNum = gsub("[^S123456]", "", cells_count)
# logNormalizedCountsRound = as_tibble(as.data.frame((t(as.matrix(logNormalizedCountsRound)))))
# logNormalizedCountsRound = logNormalizedCountsRound %>% mutate(cellID = cells_count_cellID, sampleNum = cells_count_sampleNum)
# saveRDS(object = logNormalizedCountsRound, file = paste0(homeDirectory, "logNormRoundRNA.rds"))
# 
# logNormPluri <- logNormalizedCountsRound %>% dplyr::select(cellID, sampleNum, NANOG, ALPL, POU5F1, SOX2, PODXL, LIN28A, UTF1, TERT, ZFP42, DNMT3B, GAS5, TERF1, TRIM28, ZNF483, NLRP7, FGF2, DPPA5, EPCAM, CDH1, NODAL, LEFTY1)
# rm(logNormalizedCounts)
# rm(logNormalizedCountsRound)

# cellsInDMSOSC <- linCountToOverlaps %>% dplyr::filter(SampleNum %in% c("S1", "S2"))
# cellsInLSD1SC <- linCountToOverlaps %>% dplyr::filter(SampleNum %in% c("S5", "S6"))
# 
# barcodesinDMSO <- cellsInDMSOSC$barcode %>% unique()
# barcodesinLSD1 <- cellsInLSD1SC$barcode %>% unique()
# 
# bcLSD1independent <- intersect(barcodesinDMSO, barcodesinLSD1)
# bcLSD1dependent <- setdiff(barcodesinLSD1, barcodesinDMSO)
# 
# umapCoordinates = as_tibble(read.table(file = paste0(homeDirectory, "umapCoordinates_Scanorama_50pcs_filter.tsv"), header = TRUE, stringsAsFactors = F, sep = "\t"))
# umapClusters = as_tibble(read.table(file = paste0(homeDirectory, "umapClusters_Scanorama_50pcs_filter.tsv"), header = TRUE, stringsAsFactors = F, sep = "\t")) %>%
#   rename(clusters = 1)
# umapCoordClust <- inner_join(umapCoordinates, umapClusters, by = c("cellID", "sampleNum"))
# 
# jointUMAP <- left_join(umapCoordClust, linCountToOverlaps %>% dplyr::select(cellID, barcode, SampleNum), by = "cellID")
# jointUMAP <- jointUMAP %>% mutate(barcoded = "none", condition = "none") %>% dplyr::filter(SampleNum %in% c("S5", "S6"))
# jointUMAP$barcoded <- ifelse(!is.na(jointUMAP$barcoded), "barcoded", jointUMAP$barcoded)
# jointUMAP$condition <- ifelse(jointUMAP$barcode %in% bcLSD1independent, "LSD1 ind", jointUMAP$condition)
# jointUMAP$condition <- ifelse(jointUMAP$barcode %in% bcLSD1dependent, "LSD1 dep", jointUMAP$condition)
# 
# ggplot(jointUMAP, aes(x = UMAP_1, y = UMAP_2)) +
#   geom_point() + facet_wrap(~condition)
# 
# cellIDLSD1ind <- jointUMAP %>% dplyr::filter(condition == "LSD1 ind" & clusters %in% c(0, 1, 2, 3, 5)) %>% .$cellID
# cellIDLSD1dep <- jointUMAP %>% dplyr::filter(condition == "LSD1 dep" & clusters %in% c(0, 1, 2, 3, 5)) %>% .$cellID

# umapCoordinates = as_tibble(read.table(file = paste0(homeDirectory, "umapCoordinates_Scanorama_50pcs_filter.tsv"), header = TRUE, stringsAsFactors = F, sep = "\t"))
# umapClusters <- as_tibble(read.table(file = paste0(homeDirectory, "umapClusters_Scanorama_50pcs_filter.tsv"), header = TRUE, stringsAsFactors = F, sep = "\t"))
# umapClusters <- umapClusters %>% dplyr::rename(clusters = 1)
# umapCoordClust <- inner_join(umapCoordinates, umapClusters, by = c("cellID", "sampleNum"))
# 
# jointUMAP <- left_join(umapCoordClust, linCountToOverlaps %>% dplyr::select(cellID, barcode, SampleNum), by = "cellID")
# jointUMAP <- jointUMAP %>% mutate(plot = "nonbarcoded iPSC")
# jointUMAP$plot <- ifelse(!is.na(jointUMAP$barcode), "barcoded iPSC", jointUMAP$plot)
# jointUMAP$plot <- ifelse(jointUMAP$clusters == "4", "fibroblast", jointUMAP$plot)
# jointUMAP$plot <- ifelse(jointUMAP$clusters == "5", "incomplete", jointUMAP$plot)
# jointUMAP$plot <- ifelse(jointUMAP$cellID %in% drugIndCellID$cellID&
#                            jointUMAP$sampleNum %in% c("S3", "S4") &
#                            jointUMAP$cluster %in% c("0", "1", "2", "3"), "LSD1 ind", jointUMAP$plot)
# jointUMAP$plot <- ifelse(jointUMAP$cellID %in% drugDepYCellID$cellID &
#                            jointUMAP$sampleNum %in% c("S3", "S4") &
#                            jointUMAP$cluster %in% c("0", "1", "2", "3"), "LSD1 dep", jointUMAP$plot)
# 
# logNormPluri <- readRDS(file = "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/extractedData/fateMap10X/FM3/logNormSubsetPluri.rds")
# 
# pluriMarkerPlotTable <- left_join(jointUMAP %>% dplyr::select(UMAP_1, UMAP_2, cellID, barcode, clusters, sampleNum, plot), logNormPluri, by = c("cellID", "sampleNum")) %>% dplyr::filter(sampleNum %in% c("S1", "S2", "S3", "S4"))
# 
# pluriMarkerPlotTable$plot <- factor(pluriMarkerPlotTable$plot, levels = c("nonbarcoded iPSC", "barcoded iPSC", "LSD1 ind", "LSD1 dep", "incomplete", "fibroblast"))
# ggplot(pluriMarkerPlotTable, aes(x = plot, y = NANOG)) +
#   geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), scale = "area") +
#   stat_summary(fun = mean, geom = "point", size = 2)
# 
# pluriMarkerPlotTableMelt <- melt(pluriMarkerPlotTable, id.vars = c("UMAP_1", "UMAP_2", "cellID", "barcode", "clusters", "sampleNum", "plot"))
# ggplot(pluriMarkerPlotTableMelt, aes(x = plot, y = value, fill = plot)) +
#   geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), scale = "area") +
#   stat_summary(fun = mean, geom = "point", size = 2) + facet_wrap(~variable) & NoLegend()

# jointUMAP <- left_join(umapCoordClust, linCountToOverlaps %>% dplyr::select(cellID, barcode, SampleNum), by = "cellID")
# logNormPluri <- readRDS(file = paste0(homeDirectory, "logNormRoundRNA.rds"))
# pluriMarkerPlotTable <- left_join(jointUMAP %>% dplyr::select(UMAP_1, UMAP_2, cellID, barcode, clusters, sampleNum), logNormPluri, by = c("cellID", "sampleNum")) %>% dplyr::filter(sampleNum %in% c("S1", "S2", "S3", "S4"))
# 
# geneAllList <- colnames(logNormPluri)
# 
# pluriMarkerCorrPlot_DMSO <- pluriMarkerPlotTable %>% dplyr::filter(sampleNum %in% c("S1", "S2") & clusters %in% c("0", "1", "2", "3")) %>% group_by(barcode) %>%
#   summarise(across(-c(UMAP_1, UMAP_2, cellID, clusters, sampleNum), ~ mean(.x, na.rm = TRUE)))
# pluriMarkerCorrPlot_LSD1 <- pluriMarkerPlotTable %>% dplyr::filter(sampleNum %in% c("S3", "S4") & clusters %in% c("0", "1", "2", "3")) %>% group_by(barcode) %>%
#   summarise(across(-c(UMAP_1, UMAP_2, cellID, clusters, sampleNum), ~ mean(.x, na.rm = TRUE)))
# pluriMarkerCorrPlot <- inner_join(pluriMarkerCorrPlot_DMSO, pluriMarkerCorrPlot_LSD1, by = "barcode") %>% dplyr::filter(!is.na(barcode))
# 
# saveRDS(object = pluriMarkerCorrPlot, file = paste0(homeDirectory, "corrPlotDMSOvsLSD1AllGenes.rds"))

pluriMarkerCorrPlot <- readRDS(file = paste0(homeDirectory, "corrPlotDMSOvsLSD1AllGenes.rds"))

geneList <- c("NANOG", "ALPL", "POU5F1", "SOX2", "PODXL", "LIN28A", "UTF1", "TERT", "ZFP42", "DNMT3B", "GAS5", "TERF1", "TRIM28", "ZNF483", "NLRP7", "FGF2", "DPPA5", "EPCAM", "CDH1", "NODAL", "LEFTY1", "GDF3", "FOXD3", "CYP26A1", "OTX2")
for(i in 1:length(geneList)) {
  ggplot(pluriMarkerCorrPlot, aes_string(x = paste0(geneList[i], ".x"), y = paste0(geneList[i], ".y"))) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed")
  # ggsave(filename = paste0(plotDirectory, "FM3_corrPlot_", geneList[i], ".pdf"), height = 4, width = 4)
  
  pluriMarkerCorrTable <- pluriMarkerCorrPlot %>% dplyr::select(barcode, paste0(geneList[i], ".x"), paste0(geneList[i], ".y"))
  pluriMarkerCorrTable <- pluriMarkerCorrTable %>% dplyr::rename(x = 2, y = 3)
  pluriMarkerCorrTable$foldchange <- pluriMarkerCorrTable$x + 0.001 - pluriMarkerCorrTable$y
  
  foldchangeTableTemp <- pluriMarkerCorrTable %>% mutate(gene = geneList[i])
  if(i == 1) {
    foldchangeTable <- foldchangeTableTemp
  } else{
    foldchangeTable <- bind_rows(foldchangeTable, foldchangeTableTemp)
  }
}

foldchangeTable$gene <- factor(foldchangeTable$gene, levels = c("NANOG", "POU5F1", "SOX2", "ALPL", "PODXL",
                                                          "DNMT3B", "LIN28A", "UTF1", "TERT", "ZFP42",
                                                          "GDF3", "FOXD3", "CYP26A1", "OTX2",
                                                          "GAS5", "TERF1",
                                                          "TRIM28", "ZNF483", "NLRP7", "FGF2", "DPPA5",
                                                          "EPCAM", "CDH1", "NODAL", "LEFTY1"))

ggplot(foldchangeTable, aes(x = gene, y = foldchange)) +
  geom_boxplot(fill = "gray") +
  #stat_summary(fun = mean, geom = "bar", fill = "lightgray") +
  stat_summary(fun = mean, geom = "point", size = 2.5, color = "blue") +
  #stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.25, color = "blue") +
  #geom_jitter(width = 0.25) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = c(-0.5, 0.5), linetype = "dashed")
ggsave(filename = paste0(plotDirectory, "pairwisePluriMarkerDEDMSOvsLSD1i.pdf"), height = 3, width = 9)

foldchangeTable$category <- "none"
foldchangeTable$category <- ifelse(foldchangeTable$gene %in% c("NANOG", "POU5F1", "SOX2", "ALPL", "PODXL", "DNMT3B", "LIN28A", "UTF1", "TERT", "ZFP42"), "common", foldchangeTable$category)
foldchangeTable$category <- ifelse(foldchangeTable$gene %in% c("GAS5", "TERF1"), "cluster 0", foldchangeTable$category)
foldchangeTable$category <- ifelse(foldchangeTable$gene %in% c("NLRP7", "DPPA5"), "cluster 1", foldchangeTable$category)
foldchangeTable$category <- ifelse(foldchangeTable$gene %in% c("EPCAM", "CDH1"), "cluster 3", foldchangeTable$category)

foldchangeTablePlot <- foldchangeTable %>% group_by(gene, category) %>% summarise(mean = mean(foldchange), sd = sd(foldchange))
foldchangeTablePlot$category <- factor(foldchangeTablePlot$category, levels = c("core", "common", "cluster 0", "cluster 1", "cluster 3"))

ggplot(foldchangeTablePlot %>% dplyr::filter(category != "none"), aes(x = "", y = mean)) +
  stat_summary(fun = median, geom = "bar", fill = "lightgray", width = 0.25) +
  geom_point(position = position_jitter(seed = 1234)) +
  geom_text_repel(position = position_jitter(seed = 1234), aes(label = gene)) +
  #geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.25, position = position_jitter(seed = 1234)) +
  facet_wrap(~category, ncol = 5) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = c(-0.5, 0.5), linetype = "dashed")

ggplot(foldchangeTable %>% dplyr::filter(category != "none"), aes(x = gene, y = foldchange)) +
  geom_boxplot(fill = "gray") +
  #stat_summary(fun = mean, geom = "bar", fill = "lightgray") +
  stat_summary(fun = mean, geom = "point", size = 2.5, color = "blue") +
  #stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.25, color = "blue") +
  #geom_jitter(width = 0.25) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = c(-0.5, 0.5), linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave(filename = paste0(plotDirectory, "pairwisePluriMarkerDEDMSOvsLSD1i.pdf"), height = 4, width = 4)
