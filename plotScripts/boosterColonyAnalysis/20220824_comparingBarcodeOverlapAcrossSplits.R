rm(list=ls())
gc()

#### load required libraries and functions----
library(tidyverse)
library(reshape2)
library(ggpubr)
library(fuzzyjoin)
library(VennDiagram)
library(ggforce)
library(spgs)

theme_set(theme_classic())

homeDirectory <- '/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Original Manuscript/extractedData/barcodeOverlap/'
plotDirectory <- '/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Original Manuscript/plots/barcodeOverlap/'

#### load dataset into data frames----
combDir14 <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Original Manuscript/rawData/boosters/barcodes/combined/hiFTTM14/"
combDir16 <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Original Manuscript/rawData/boosters/barcodes/combined/hiFTTM16/"

DOT1LDir13 <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Original Manuscript/rawData/boosters/barcodes/DOT1Li/hiFTTM13/"

LSD1Dir6 <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Original Manuscript/rawData/boosters/barcodes/LSD1i/hiFTTM6/"
LSD1Dir28 <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Original Manuscript/rawData/boosters/barcodes/LSD1i/hiFTTM28/"

sampleDirList <- c(combDir14, combDir16, DOT1LDir13, LSD1Dir6)

sampleTablesList <- list(rep('NA', length(sampleDirList)))
sampleFoldersList <- list(rep('NA', length(sampleDirList)))

for(j in 1:length(sampleDirList)) {
  sampleFolders <- list.files(path = sampleDirList[[j]], pattern = '*d*.txt', recursive = TRUE)
  sampleFoldersList[[j]] <- sampleFolders
  sampleTables <- list(rep("NA", length(sampleFolders)))
  for (i in 1:length(sampleFolders)) {
    sampleTable = read.table(paste0(sampleDirList[[j]], sampleFolders[i]), header = FALSE, sep = "\t")
    sampleTables[[i]] <- sampleTable
  }
  sampleTablesList[[j]] <- sampleTables
}

sampleDirListSC <- c(LSD1Dir28)

sampleTablesListSC <- list(rep('NA', length(sampleDirListSC)))
sampleFoldersListSC <- list(rep('NA', length(sampleDirListSC)))

for(j in 1:length(sampleDirListSC)) {
  sampleTable = read.table(paste0(sampleDirListSC[[j]], "stepThreeStarcodeShavedReads.txt"), header = TRUE, sep = "\t") %>% dplyr::select(2, 4, 8) %>%
    dplyr::rename(BC50 = 2, sampleNum = 3) %>% group_by(BC50, sampleNum) %>% summarize(nUMI = sum(UMI)) %>% dplyr::filter(nUMI > 10)
  sampleNames <- sampleTable$sampleNum %>% unique()
  sampleFoldersListSC[[j]] <- sampleNames
  sampleTables <- list(rep("NA", length(sampleNames)))
  for (i in 1:length(sampleNames)) {
    sampleTables[[i]] <- sampleTable %>% dplyr::filter(sampleNum == sampleNames[[i]])
  }
  sampleTablesListSC[[j]] <- sampleTables
}

#### work with original barcode outputs ####
#### LSD1 ##########################################################################################################################
rep1 <- list(c(1, 1), c(1, 2), c(1, 5), c(1, 6))
rep2 <- list(c(2, 1), c(2, 3))
rep3 <- list(c(4, 1), c(4, 2))

repList <- list(rep1, rep2, rep3)

cutoff = c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)
foldChangecutoffDep = 2

numIndList <- list()
numDepXList <- list()
numDepYList <- list()
ratioDepYIndList <- list()
scatterPlotsList <- list()
for(j in 1:length(repList)) {
  if(length(repList[[j]]) == 4) {
    DMSO_A_Temp <- sampleTablesList[[repList[[j]][[1]][1]]][[repList[[j]][[1]][2]]] %>% mutate(V2 = V2/(sum(V2)) * 10^6)
    DMSO_B_Temp <- sampleTablesList[[repList[[j]][[2]][1]]][[repList[[j]][[2]][2]]] %>% mutate(V2 = V2/(sum(V2)) * 10^6)
    DMSO_Temp <- full_join(DMSO_A_Temp, DMSO_B_Temp, by = "V1") %>% replace(is.na(.), 0) %>% rowwise() %>% mutate(max = max(V2.x, V2.y))
    
    if(length(repList[[j]]) == 3) {
      LSD1_Temp <- sampleTablesList[[repList[[j]][[3]][1]]][[repList[[j]][[3]][2]]] %>% mutate(V2 = V2/(sum(V2)) * 10^6) %>% dplyr::rename(max = V2)
    } else {
      LSD1_A_Temp <- sampleTablesList[[repList[[j]][[3]][1]]][[repList[[j]][[3]][2]]] %>% mutate(V2 = V2/(sum(V2)) * 10^6)
      LSD1_B_Temp <- sampleTablesList[[repList[[j]][[4]][1]]][[repList[[j]][[4]][2]]] %>% mutate(V2 = V2/(sum(V2)) * 10^6)
      LSD1_Temp <- full_join(LSD1_A_Temp, LSD1_B_Temp, by = "V1") %>% replace(is.na(.), 0) %>% rowwise() %>% mutate(max = max(V2.x, V2.y))
      
      # ggplot(DMSO_Temp, aes(x = V2.x, y = V2.y)) +
      #   geom_point()
      # 
      # ggplot(LSD1_Temp, aes(x = V2.x, y = V2.y)) +
      #   geom_point()
    }
    
  } else {
    DMSO_Temp <- sampleTablesList[[repList[[j]][[1]][1]]][[repList[[j]][[1]][2]]] %>% mutate(V2 = V2/(sum(V2)) * 10^6) %>% dplyr::rename(max = V2)
    LSD1_Temp <- sampleTablesList[[repList[[j]][[2]][1]]][[repList[[j]][[2]][2]]] %>% mutate(V2 = V2/(sum(V2)) * 10^6) %>% dplyr::rename(max = V2)
  }
  
  barcodesOverlap <- full_join(DMSO_Temp %>% dplyr::select(V1, max),
                               LSD1_Temp %>% dplyr::select(V1, max), by = "V1") %>% replace(is.na(.), 0)
  
  correctionFactor <- barcodesOverlap %>% ungroup() %>% slice_max(order_by = max.x, n = 5) %>%
    rowwise() %>% mutate(foldchange = max.x/max.y) %>% .$foldchange %>% mean()
  
  barcodesOverlap <- barcodesOverlap %>% rowwise() %>% mutate(max.y_corr = max.y * correctionFactor) %>%
    mutate(foldchange = log((max.y_corr + 1) / (max.x + 1)))
  
  # ggplot(barcodesOverlap, aes(x = max.x, y = max.y_corr)) +
  #   geom_point() +
  #   geom_abline(slope = 1, intercept = 0)
  
  numInd <- list()
  numDepX <- list()
  numDepY <- list()
  ratioDepYInd <- list()
  scatterPlots <- list()
  for(i in 1:length(cutoff)) {
    drugDepX <- barcodesOverlap %>% filter(foldchange < -foldChangecutoffDep | foldchange > foldChangecutoffDep) %>% filter(max.x > cutoff[i] & max.y_corr < cutoff[i])
    drugDepY <- barcodesOverlap %>% filter(foldchange < -foldChangecutoffDep | foldchange > foldChangecutoffDep) %>% filter(max.y_corr > cutoff[i] & max.x < cutoff[i])
    drugInd <- barcodesOverlap %>% filter(foldchange > -foldChangecutoffDep & foldchange < foldChangecutoffDep) %>% filter(max.x > cutoff[i] & max.y_corr > cutoff[i])
    
    scatterPlots[[i]] <- ggplot() +
      geom_point(data = barcodesOverlap, aes(x = max.x, y = max.y_corr), color = "gray") +
      geom_point(data = drugDepX, aes(x = max.x, y = max.y_corr), color = "red") +
      geom_point(data = drugDepY, aes(x = max.x, y = max.y_corr), color = "green") +
      geom_point(data = drugInd, aes(x = max.x, y = max.y_corr), color = "black") +
      geom_abline(intercept = 0, slope = 1)
    
    numInd[[i]] <- nrow(drugInd)
    numDepX[[i]] <- nrow(drugDepX)
    numDepY[[i]] <- nrow(drugDepY)
    ratioDepYInd[[i]] <- nrow(drugDepY)/(nrow(drugInd) + nrow(drugDepX))
  }
  numIndList[[j]] <- numInd
  numDepXList[[j]] <- numDepX
  numDepYList[[j]] <- numDepY
  ratioDepYIndList[[j]] <- ratioDepYInd
  scatterPlotsList[[j]] <- scatterPlots
}

####plot example scatterplot for figure########################################################################################################################################
###############################################################################################################################################################################
j = 1

if(length(repList[[j]]) == 4) {
  DMSO_A_Temp <- sampleTablesList[[repList[[j]][[1]][1]]][[repList[[j]][[1]][2]]] %>% mutate(V2 = V2/(sum(V2)) * 10^6)
  DMSO_B_Temp <- sampleTablesList[[repList[[j]][[2]][1]]][[repList[[j]][[2]][2]]] %>% mutate(V2 = V2/(sum(V2)) * 10^6)
  DMSO_Temp <- full_join(DMSO_A_Temp, DMSO_B_Temp, by = "V1") %>% replace(is.na(.), 0) %>% rowwise() %>% mutate(max = max(V2.x, V2.y))
  
  if(length(repList[[j]]) == 3) {
    LSD1_Temp <- sampleTablesList[[repList[[j]][[3]][1]]][[repList[[j]][[3]][2]]] %>% mutate(V2 = V2/(sum(V2)) * 10^6) %>% dplyr::rename(max = V2)
  } else {
    LSD1_A_Temp <- sampleTablesList[[repList[[j]][[3]][1]]][[repList[[j]][[3]][2]]] %>% mutate(V2 = V2/(sum(V2)) * 10^6)
    LSD1_B_Temp <- sampleTablesList[[repList[[j]][[4]][1]]][[repList[[j]][[4]][2]]] %>% mutate(V2 = V2/(sum(V2)) * 10^6)
    LSD1_Temp <- full_join(LSD1_A_Temp, LSD1_B_Temp, by = "V1") %>% replace(is.na(.), 0) %>% rowwise() %>% mutate(max = max(V2.x, V2.y))
    
    # ggplot(DMSO_Temp, aes(x = V2.x, y = V2.y)) +
    #   geom_point()
    # 
    # ggplot(LSD1_Temp, aes(x = V2.x, y = V2.y)) +
    #   geom_point()
  }
  
} else {
  DMSO_Temp <- sampleTablesList[[repList[[j]][[1]][1]]][[repList[[j]][[1]][2]]] %>% mutate(V2 = V2/(sum(V2)) * 10^6) %>% dplyr::rename(max = V2)
  LSD1_Temp <- sampleTablesList[[repList[[j]][[2]][1]]][[repList[[j]][[2]][2]]] %>% mutate(V2 = V2/(sum(V2)) * 10^6) %>% dplyr::rename(max = V2)
}

barcodesOverlap <- full_join(DMSO_Temp %>% dplyr::select(V1, max),
                             LSD1_Temp %>% dplyr::select(V1, max), by = "V1") %>% replace(is.na(.), 0)

correctionFactor <- barcodesOverlap %>% ungroup() %>% slice_max(order_by = max.x, n = 5) %>%
  rowwise() %>% mutate(foldchange = max.x/max.y) %>% .$foldchange %>% mean()

barcodesOverlap <- barcodesOverlap %>% rowwise() %>% mutate(max.y_corr = max.y * correctionFactor) %>%
  mutate(foldchange = log((max.y_corr + 1) / (max.x + 1)))

i = 5

drugDepX <- barcodesOverlap %>% filter(foldchange < -foldChangecutoffDep | foldchange > foldChangecutoffDep) %>% filter(max.x > cutoff[i] & max.y_corr < cutoff[i])
drugDepY <- barcodesOverlap %>% filter(foldchange < -foldChangecutoffDep | foldchange > foldChangecutoffDep) %>% filter(max.y_corr > cutoff[i] & max.x < cutoff[i])
drugInd <- barcodesOverlap %>% filter(foldchange > -foldChangecutoffDep & foldchange < foldChangecutoffDep) %>% filter(max.x > cutoff[i] & max.y_corr > cutoff[i])

library(ggrastr)
ggplot() +
  rasterise(geom_point(data = barcodesOverlap, aes(x = max.x, y = max.y_corr), color = "lightgray", shape = 16), dpi = 300) +
  geom_point(data = drugDepX, aes(x = max.x, y = max.y_corr), color = "red", shape = 16) +
  geom_point(data = drugDepY, aes(x = max.x, y = max.y_corr), color = "green", shape = 16) +
  geom_point(data = drugInd, aes(x = max.x, y = max.y_corr), color = "black", shape = 16) +
  geom_abline(intercept = 0, slope = 1) + xlim(0, 5000) + ylim(0, 5000)
ggsave(filename = paste0(plotDirectory, 'DMSOvsLSD1i_scatterplotExampleForFigure_corrected.pdf'), width = 4, height = 4, useDingbats = F)

i = 5

drugDepX <- barcodesOverlap %>% filter(foldchange < -foldChangecutoffDep | foldchange > foldChangecutoffDep) %>% filter(max.x > cutoff[i] & max.y < cutoff[i])
drugDepY <- barcodesOverlap %>% filter(foldchange < -foldChangecutoffDep | foldchange > foldChangecutoffDep) %>% filter(max.y > cutoff[i] & max.x < cutoff[i])
drugInd <- barcodesOverlap %>% filter(foldchange > -foldChangecutoffDep & foldchange < foldChangecutoffDep) %>% filter(max.x > cutoff[i] & max.y > cutoff[i])

ggplot() +
  rasterise(geom_point(data = barcodesOverlap, aes(x = max.x, y = max.y), color = "lightgray", shape = 16), dpi = 300) +
  geom_point(data = drugDepX, aes(x = max.x, y = max.y), color = "red", shape = 16) +
  geom_point(data = drugDepY, aes(x = max.x, y = max.y), color = "green", shape = 16) +
  geom_point(data = drugInd, aes(x = max.x, y = max.y), color = "black", shape = 16) +
  geom_abline(intercept = 0, slope = 1) + xlim(0, 5000) + ylim(0, 5000)
ggsave(filename = paste0(plotDirectory, 'DMSOvsLSD1i_scatterplotExampleForFigure_uncorrected.pdf'), width = 4, height = 4, useDingbats = F)

# grid.newpage()
# venndiagram <- draw.pairwise.venn(area1 = nrow(drugInd) + nrow(drugDepY), area2 = nrow(drugInd) + nrow(drugDepY), cross.area = nrow(drugInd),
#                                   euler.d = TRUE, fill = c(alpha("gray5", 0.5), alpha("gray5", 0.5)),
#                                   fontfamily = "sans", cex = 2.5, print.mode = "percent", sigdigs = 3)
# ggsave(plot = venndiagram, filename = paste0(plotDirectory, 'overlapVennDiagram_LSD1_1.pdf'), width = 5, height = 5, useDingbats = F)
###############################################################################################################################################################################

aggregateSampleTable <- tibble(rep = c(rep("1", length(cutoff)), rep("2", length(cutoff)), rep("3", length(cutoff))),
                               cutoff = rep(cutoff, 3),
                               ratio = c(ratioDepYIndList[[1]], ratioDepYIndList[[2]], ratioDepYIndList[[3]]),
                               numInd = c(numIndList[[1]], numIndList[[2]], numIndList[[3]]),
                               numDepX = c(numDepXList[[1]], numDepXList[[2]], numDepXList[[3]]),
                               numDepY = c(numDepYList[[1]], numDepYList[[2]], numDepYList[[3]]))
aggregateSampleTable <- aggregateSampleTable %>% rowwise() %>% mutate(nCells = numInd + numDepX + numDepY)

ggplot(aggregateSampleTable, aes(x = cutoff, y = as.double(ratio))) +
  geom_point() +
  stat_summary(fun = mean, geom = "point", size = 2, color = "blue") +
  stat_summary(fun.data = mean_se, geom = "errorbar", color = "blue") + ylim(0, 9) +
  geom_hline(yintercept = 1, linetype = "dashed")
ggsave(filename = paste0(plotDirectory, 'DMSOvsLSD1i_foldchangeNumLineagesVaryingCutoff.pdf'), width = 5, height = 2, useDingbats = F)

aggregateSampleTableFilter <- aggregateSampleTable %>% dplyr::filter(cutoff == 500)
ggplot(aggregateSampleTableFilter, aes(x = as.factor(cutoff), y = as.double(ratio))) +
  geom_jitter(width = 0.1, height = 0.1) +
  stat_summary(fun = mean, geom = "point", size = 4) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1) +
  ylim(0, 8.1)
ggsave(filename = paste0(plotDirectory, 'DMSOvsLSD1i_foldchangeNumLineages.pdf'), width = 1.5, height = 3, useDingbats = F)

grid.newpage()
venndiagram <- draw.pairwise.venn(area1 = aggregateSampleTableFilter$numInd[[1]] + aggregateSampleTableFilter$numDepY[[1]], area2 = aggregateSampleTableFilter$numInd[[1]] + aggregateSampleTableFilter$numDepX[[1]], cross.area = aggregateSampleTableFilter$numInd[[1]],
                                  euler.d = TRUE, fill = c(alpha("gray5", 0.5), alpha("gray5", 0.5)),
                                  fontfamily = "sans", cex = 2.5, print.mode = "percent", sigdigs = 3)
ggsave(plot = venndiagram, filename = paste0(plotDirectory, 'overlapVennDiagram_LSD1_1.pdf'), width = 5, height = 5, useDingbats = F)

grid.newpage()
venndiagram <- draw.pairwise.venn(area1 = aggregateSampleTableFilter$numInd[[3]] + aggregateSampleTableFilter$numDepY[[3]], area2 = aggregateSampleTableFilter$numInd[[3]] + aggregateSampleTableFilter$numDepX[[3]], cross.area = aggregateSampleTableFilter$numInd[[3]],
                                  euler.d = TRUE, fill = c(alpha("gray5", 0.5), alpha("gray5", 0.5)),
                                  fontfamily = "sans", cex = 2.5, print.mode = "percent", sigdigs = 3)
ggsave(plot = venndiagram, filename = paste0(plotDirectory, 'overlapVennDiagram_LSD1_2.pdf'), width = 5, height = 5, useDingbats = F)

grid.newpage()
venndiagram <- draw.pairwise.venn(area1 = aggregateSampleTableFilter$numInd[[2]] + aggregateSampleTableFilter$numDepY[[2]], area2 = aggregateSampleTableFilter$numInd[[2]] + aggregateSampleTableFilter$numDepX[[2]], cross.area = aggregateSampleTableFilter$numInd[[2]],
                                  euler.d = TRUE, fill = c(alpha("gray5", 0.5), alpha("gray5", 0.5)),
                                  fontfamily = "sans", cex = 2.5, print.mode = "percent", sigdigs = 3)
ggsave(plot = venndiagram, filename = paste0(plotDirectory, 'overlapVennDiagram_LSD1_3.pdf'), width = 5, height = 5, useDingbats = F)

#### make pie chart plot of number of reads per LSD1-dependent and LSD1-independent ####
barcodesOverlapList <- list()
for(j in 1:length(repList)) {
  if(length(repList[[j]]) == 4) {
    DMSO_A_Temp <- sampleTablesList[[repList[[j]][[1]][1]]][[repList[[j]][[1]][2]]] %>% dplyr::mutate(V2 = V2/(sum(V2)) * 10^6)
    DMSO_B_Temp <- sampleTablesList[[repList[[j]][[2]][1]]][[repList[[j]][[2]][2]]] %>% dplyr::mutate(V2 = V2/(sum(V2)) * 10^6)
    DMSO_Temp <- full_join(DMSO_A_Temp, DMSO_B_Temp, by = "V1") %>% replace(is.na(.), 0) %>% rowwise() %>% dplyr::mutate(max = max(V2.x, V2.y))
    
    if(length(repList[[j]]) == 3) {
      LSD1_Temp <- sampleTablesList[[repList[[j]][[3]][1]]][[repList[[j]][[3]][2]]] %>% dplyr::mutate(V2 = V2/(sum(V2)) * 10^6) %>% dplyr::rename(max = V2)
    } else {
      LSD1_A_Temp <- sampleTablesList[[repList[[j]][[3]][1]]][[repList[[j]][[3]][2]]] %>% dplyr::mutate(V2 = V2/(sum(V2)) * 10^6)
      LSD1_B_Temp <- sampleTablesList[[repList[[j]][[4]][1]]][[repList[[j]][[4]][2]]] %>% dplyr::mutate(V2 = V2/(sum(V2)) * 10^6)
      LSD1_Temp <- full_join(LSD1_A_Temp, LSD1_B_Temp, by = "V1") %>% replace(is.na(.), 0) %>% rowwise() %>% dplyr::mutate(max = max(V2.x, V2.y))
      
      # ggplot(DMSO_Temp, aes(x = V2.x, y = V2.y)) +
      #   geom_point()
      # 
      # ggplot(LSD1_Temp, aes(x = V2.x, y = V2.y)) +
      #   geom_point()
    }
    
  } else {
    DMSO_Temp <- sampleTablesList[[repList[[j]][[1]][1]]][[repList[[j]][[1]][2]]] %>% dplyr::mutate(V2 = V2/(sum(V2)) * 10^6) %>% dplyr::rename(max = V2)
    LSD1_Temp <- sampleTablesList[[repList[[j]][[2]][1]]][[repList[[j]][[2]][2]]] %>% dplyr::mutate(V2 = V2/(sum(V2)) * 10^6) %>% dplyr::rename(max = V2)
  }
  
  barcodesOverlap <- full_join(DMSO_Temp %>% dplyr::select(V1, max),
                               LSD1_Temp %>% dplyr::select(V1, max), by = "V1") %>% replace(is.na(.), 0)
  
  correctionFactor <- barcodesOverlap %>% ungroup() %>% slice_max(order_by = max.x, n = 5) %>%
    rowwise() %>% dplyr::mutate(foldchange = max.x/max.y) %>% .$foldchange %>% mean()
  
  barcodesOverlap <- barcodesOverlap %>% rowwise() %>% dplyr::mutate(max.y_corr = max.y * correctionFactor) %>%
    dplyr::mutate(foldchange = log((max.y_corr + 1) / (max.x + 1)))
  
  barcodesOverlapList[[j]] <- barcodesOverlap
}

for(j in 1:length(barcodesOverlapList)) {
  dataTableTemp <- barcodesOverlapList[[j]]
  
  drugDepX <- dataTableTemp %>% filter(foldchange < -foldChangecutoffDep | foldchange > foldChangecutoffDep) %>% filter(max.x > 500 & max.y_corr < 500)
  drugDepY <- dataTableTemp %>% filter(foldchange < -foldChangecutoffDep | foldchange > foldChangecutoffDep) %>% filter(max.y_corr > 500 & max.x < 500)
  drugInd <- dataTableTemp %>% filter(foldchange > -foldChangecutoffDep & foldchange < foldChangecutoffDep) %>% filter(max.x > 500 & max.y_corr > 500)
  
  drugDepXReads1 <- drugDepX %>% .$max.x %>% sum()
  drugDepXReads2 <- drugDepX %>% .$max.y_corr %>% sum()
  drugDepYReads1 <- drugDepY %>% .$max.x %>% sum()
  drugDepYReads2 <- drugDepY %>% .$max.y_corr %>% sum()
  drugIndReads1 <- drugInd %>% .$max.x %>% sum()
  drugIndReads2 <- drugInd %>% .$max.y_corr %>% sum()
  
  readFrac1 <- data.frame(group = c("X", "Y", "IND"),
                          value = c(drugDepXReads1, drugDepYReads1, drugIndReads1))
  readFrac1 <- readFrac1 %>%
    arrange(desc(group)) %>%
    mutate(prop = value / sum(readFrac1$value) * 100) %>%
    mutate(ypos = cumsum(prop) - 0.5 * prop)
  
  ggplot(readFrac1, aes(x = "", y = prop, fill = group)) +
    geom_bar(stat = "identity", width = 1, color = "black") +
    theme_void() + theme(legend.position = "none") +
    geom_text(aes(y = ypos, label = round(prop, 0)), color = "black", size = 6) +
    coord_polar("y", start = 0)
  ggsave(filename = paste0(plotDirectory, "fracReadsPieChart_DMSO_", j, ".pdf"), height = 4, width = 4)
  rm(readFrac1)
  
  readFrac2 <- data.frame(group = c("X", "Y", "IND"),
                          value = c(drugDepXReads2, drugDepYReads2, drugIndReads2))
  readFrac2 <- readFrac2 %>%
    arrange(desc(group)) %>%
    mutate(prop = value / sum(readFrac2$value) * 100) %>%
    mutate(ypos = cumsum(prop) - 0.5 * prop)
  
  ggplot(readFrac2, aes(x = "", y = prop, fill = group)) +
    geom_bar(stat = "identity", width = 1, color = "black") +
    theme_void() + theme(legend.position = "none") +
    geom_text(aes(y = ypos, label = round(prop, 0)), color = "black", size = 6) +
    coord_polar("y", start = 0)
  ggsave(filename = paste0(plotDirectory, "fracReadsPieChart_LSD1_", j, ".pdf"), height = 4, width = 4)
  rm(readFrac2)
}

#### DOT1L #########################################################################################################################
rep1 <- list(c(1, 1), c(1, 2), c(1, 3))
rep2 <- list(c(2, 1), c(2, 2))
rep3 <- list(c(3, 1), c(3, 2), c(3, 3), c(3, 4))

repList <- list(rep1, rep2, rep3)

cutoff = c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)
foldChangecutoffDep = 2

numIndList <- list()
numDepXList <- list()
numDepYList <- list()
ratioDepYIndList <- list()
scatterPlotsList <- list()
for(j in 1:length(repList)) {
  if(length(repList[[j]]) == 4) {
    DMSO_A_Temp <- sampleTablesList[[repList[[j]][[1]][1]]][[repList[[j]][[1]][2]]] %>% mutate(V2 = V2/(sum(V2)) * 10^6)
    DMSO_B_Temp <- sampleTablesList[[repList[[j]][[2]][1]]][[repList[[j]][[2]][2]]] %>% mutate(V2 = V2/(sum(V2)) * 10^6)
    DMSO_Temp <- full_join(DMSO_A_Temp, DMSO_B_Temp, by = "V1") %>% replace(is.na(.), 0) %>% rowwise() %>% mutate(max = max(V2.x, V2.y))
    
    if(length(repList[[j]]) == 3) {
      LSD1_Temp <- sampleTablesList[[repList[[j]][[3]][1]]][[repList[[j]][[3]][2]]] %>% mutate(V2 = V2/(sum(V2)) * 10^6) %>% dplyr::rename(max = V2)
    } else {
      LSD1_A_Temp <- sampleTablesList[[repList[[j]][[3]][1]]][[repList[[j]][[3]][2]]] %>% mutate(V2 = V2/(sum(V2)) * 10^6)
      LSD1_B_Temp <- sampleTablesList[[repList[[j]][[4]][1]]][[repList[[j]][[4]][2]]] %>% mutate(V2 = V2/(sum(V2)) * 10^6)
      LSD1_Temp <- full_join(LSD1_A_Temp, LSD1_B_Temp, by = "V1") %>% replace(is.na(.), 0) %>% rowwise() %>% mutate(max = max(V2.x, V2.y))
      
      # ggplot(DMSO_Temp, aes(x = V2.x, y = V2.y)) +
      #   geom_point()
      # 
      # ggplot(LSD1_Temp, aes(x = V2.x, y = V2.y)) +
      #   geom_point()
    }
    
  } else {
    DMSO_Temp <- sampleTablesList[[repList[[j]][[1]][1]]][[repList[[j]][[1]][2]]] %>% mutate(V2 = V2/(sum(V2)) * 10^6) %>% dplyr::rename(max = V2)
    LSD1_Temp <- sampleTablesList[[repList[[j]][[2]][1]]][[repList[[j]][[2]][2]]] %>% mutate(V2 = V2/(sum(V2)) * 10^6) %>% dplyr::rename(max = V2)
  }
  
  barcodesOverlap <- full_join(DMSO_Temp %>% dplyr::select(V1, max),
                               LSD1_Temp %>% dplyr::select(V1, max), by = "V1") %>% replace(is.na(.), 0)
  
  correctionFactor <- barcodesOverlap %>% ungroup() %>% slice_max(order_by = max.x, n = 5) %>%
    rowwise() %>% mutate(foldchange = max.x/max.y) %>% .$foldchange %>% mean()
  
  barcodesOverlap <- barcodesOverlap %>% rowwise() %>% mutate(max.y_corr = max.y * correctionFactor) %>%
    mutate(foldchange = log((max.y_corr + 1) / (max.x + 1)))
  
  # ggplot(barcodesOverlap, aes(x = max.x, y = max.y_corr)) +
  #   geom_point() +
  #   geom_abline(slope = 1, intercept = 0)
  
  numInd <- list()
  numDepX <- list()
  numDepY <- list()
  ratioDepYInd <- list()
  scatterPlots <- list()
  for(i in 1:length(cutoff)) {
    drugDepX <- barcodesOverlap %>% filter(foldchange < -foldChangecutoffDep | foldchange > foldChangecutoffDep) %>% filter(max.x > cutoff[i] & max.y_corr < cutoff[i])
    drugDepY <- barcodesOverlap %>% filter(foldchange < -foldChangecutoffDep | foldchange > foldChangecutoffDep) %>% filter(max.y_corr > cutoff[i] & max.x < cutoff[i])
    drugInd <- barcodesOverlap %>% filter(foldchange > -foldChangecutoffDep & foldchange < foldChangecutoffDep) %>% filter(max.x > cutoff[i] & max.y_corr > cutoff[i])
    
    scatterPlots[[i]] <- ggplot() +
      geom_point(data = barcodesOverlap, aes(x = max.x, y = max.y_corr), color = "gray") +
      geom_point(data = drugDepX, aes(x = max.x, y = max.y_corr), color = "red") +
      geom_point(data = drugDepY, aes(x = max.x, y = max.y_corr), color = "green") +
      geom_point(data = drugInd, aes(x = max.x, y = max.y_corr), color = "black") +
      geom_abline(intercept = 0, slope = 1)
    
    numInd[[i]] <- nrow(drugInd)
    numDepX[[i]] <- nrow(drugDepX)
    numDepY[[i]] <- nrow(drugDepY)
    ratioDepYInd[[i]] <- nrow(drugDepY)/(nrow(drugInd) + nrow(drugDepX))
  }
  numIndList[[j]] <- numInd
  numDepXList[[j]] <- numDepX
  numDepYList[[j]] <- numDepY
  ratioDepYIndList[[j]] <- ratioDepYInd
  scatterPlotsList[[j]] <- scatterPlots
}
  
ggpubr::ggarrange(plotlist = scatterPlotsList[[3]])

aggregateSampleTable <- tibble(rep = c(rep("1", length(cutoff)), rep("2", length(cutoff)), rep("3", length(cutoff))),
                               cutoff = rep(cutoff, 3),
                               ratio = c(ratioDepYIndList[[1]], ratioDepYIndList[[2]], ratioDepYIndList[[3]]),
                               numInd = c(numIndList[[1]], numIndList[[2]], numIndList[[3]]),
                               numDepX = c(numDepXList[[1]], numDepXList[[2]], numDepXList[[3]]),
                               numDepY = c(numDepYList[[1]], numDepYList[[2]], numDepYList[[3]]))
aggregateSampleTable <- aggregateSampleTable %>% rowwise() %>% mutate(nCells = numInd + numDepX + numDepY)

ggplot(aggregateSampleTable, aes(x = cutoff, y = as.double(ratio))) +
  geom_point() +
  stat_summary(fun = mean, geom = "point", size = 2, color = "blue") +
  stat_summary(fun.data = mean_se, geom = "errorbar", color = "blue") + ylim(0, 9) +
  geom_hline(yintercept = 1, linetype = "dashed")
ggsave(filename = paste0(plotDirectory, 'DMSOvsDOT1Li_foldchangeNumLineagesVaryingCutoff.pdf'), width = 5, height = 2, useDingbats = F)

aggregateSampleTableFilter <- aggregateSampleTable %>% dplyr::filter(cutoff == 500)
ggplot(aggregateSampleTableFilter, aes(x = as.factor(cutoff), y = as.double(ratio))) +
  geom_jitter(width = 0.1, height = 0.1) +
  stat_summary(fun = mean, geom = "point", size = 4) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1) +
  ylim(0, 8.1)
ggsave(filename = paste0(plotDirectory, 'DMSOvsDOT1Li_foldchangeNumLineages.pdf'), width = 1.5, height = 3, useDingbats = F)

grid.newpage()
venndiagram <- draw.pairwise.venn(area1 = aggregateSampleTableFilter$numInd[[1]] + aggregateSampleTableFilter$numDepY[[1]], area2 = aggregateSampleTableFilter$numInd[[1]] + aggregateSampleTableFilter$numDepX[[1]], cross.area = aggregateSampleTableFilter$numInd[[1]],
                                  euler.d = TRUE, fill = c(alpha("gray5", 0.5), alpha("gray5", 0.5)),
                                  fontfamily = "sans", cex = 2.5, print.mode = "percent", sigdigs = 3)
ggsave(plot = venndiagram, filename = paste0(plotDirectory, 'overlapVennDiagram_DOT1L_1.pdf'), width = 5, height = 5, useDingbats = F)

grid.newpage()
venndiagram <- draw.pairwise.venn(area1 = aggregateSampleTableFilter$numInd[[3]] + aggregateSampleTableFilter$numDepY[[3]], area2 = aggregateSampleTableFilter$numInd[[3]] + aggregateSampleTableFilter$numDepX[[3]], cross.area = aggregateSampleTableFilter$numInd[[3]],
                                  euler.d = TRUE, fill = c(alpha("gray5", 0.5), alpha("gray5", 0.5)),
                                  fontfamily = "sans", cex = 2.5, print.mode = "percent", sigdigs = 3)
ggsave(plot = venndiagram, filename = paste0(plotDirectory, 'overlapVennDiagram_DOT1L_2.pdf'), width = 5, height = 5, useDingbats = F)

grid.newpage()
venndiagram <- draw.pairwise.venn(area1 = aggregateSampleTableFilter$numInd[[2]] + aggregateSampleTableFilter$numDepY[[2]], area2 = aggregateSampleTableFilter$numInd[[2]] + aggregateSampleTableFilter$numDepX[[2]], cross.area = aggregateSampleTableFilter$numInd[[2]],
                                  euler.d = TRUE, fill = c(alpha("gray5", 0.5), alpha("gray5", 0.5)),
                                  fontfamily = "sans", cex = 2.5, print.mode = "percent", sigdigs = 3)
ggsave(plot = venndiagram, filename = paste0(plotDirectory, 'overlapVennDiagram_DOT1L_3.pdf'), width = 5, height = 5, useDingbats = F)
  
#### input samples into condensed format ####
# spikeBC1 <- 'TCCAGGTCCTCCTACTTGTACAACACCTTGTACAGCTGCTAGTGGTAGAAGAGGTACAACAACAACACGAGCATCATGAGGATCTACAGCATCAAGAACA' %>% substr(0, 50)
# spikeBC2 <- 'ACGTTGTGCATGACCTTGATCACCAGCTCGATGTCGAACATCACGAGCTCGTTCTGCATCTGCAAGAACACCTCGTCCTTGAACTGCTCGACGTCCATGA' %>% substr(0, 50)
# nBC1 <- 20000
# nBC2 <- 5000
# 
# standardTableAll = list()
# lmr = list()
# for(j in 1:length(sampleTablesListSC)) {
#   for(i in 1:length(sampleTablesListSC[[j]])) {
#     standardTable <- sampleTablesListSC[[j]][[i]] %>% dplyr::filter(., BC50 == spikeBC1 | BC50 == spikeBC2) %>% arrange(-nUMI) %>% mutate(ratio = .$nUMI[1]/.$nUMI[2])
#       if(is.null(dim(standardTableAll))) {
#         standardTableAll = standardTable
#         lmr = tibble(name = standardTable$sampleNum, coef = coef(lm(c(nBC1, nBC2) ~ 0 + standardTable$nUMI)))
#       } else {
#         standardTableAll = bind_rows(standardTableAll, standardTable)
#         lmr = bind_rows(lmr, tibble(name = standardTable$sampleNum, coef = coef(lm(c(nBC1, nBC2) ~ 0 + standardTable$nUMI))))
#       }
#   }
# }
# 
# lmr <- unique(lmr)
# 
# combinedSampleTable <- bind_rows(sampleTablesListSC) %>% filter(BC50 != spikeBC1) %>% filter(BC50 != spikeBC2)
# 
# rep1 <- c("TM28_DMSO_2A", "TM28_DMSO_2B", "TM28_LSD1i_2A", "TM28_LSD1i_2B", "TM28_DOT1Li_2A", "TM28_DOT1Li_2B")
# rep2 <- c("TM28_DMSO_3A", "TM28_DMSO_3B", "TM28_LSD1i_3A", "TM28_LSD1i_3B", "TM28_DOT1Li_3A", "TM28_DOT1Li_3B")
# 
# repList <- list(rep1, rep2)
# DMSOTablesList <- list()
# LSD1TablesList <- list()
# DOT1LTablesList <- list()
# 
# for(j in 1:length(repList)) {
#   print(j)
#   DMSO_A_Temp <- combinedSampleTable %>% dplyr::filter(sampleNum == repList[[j]][[1]]) %>% rowwise() %>% mutate(nCells = nUMI * (dplyr::filter(lmr, name == repList[[j]][[1]]) %>% .$coef))
#   DMSO_B_Temp <- combinedSampleTable %>% dplyr::filter(sampleNum == repList[[j]][[2]]) %>% rowwise() %>% mutate(nCells = nUMI * (dplyr::filter(lmr, name == repList[[j]][[2]]) %>% .$coef))
# 
#   LSD1_A_Temp <- combinedSampleTable %>% dplyr::filter(sampleNum == repList[[j]][[3]]) %>% rowwise() %>% mutate(nCells = nUMI * (dplyr::filter(lmr, name == repList[[j]][[3]]) %>% .$coef))
#   LSD1_B_Temp <- combinedSampleTable %>% dplyr::filter(sampleNum == repList[[j]][[4]]) %>% rowwise() %>% mutate(nCells = nUMI * (dplyr::filter(lmr, name == repList[[j]][[4]]) %>% .$coef))
# 
#   DOT1L_A_Temp <- combinedSampleTable %>% dplyr::filter(sampleNum == repList[[j]][[5]]) %>% rowwise() %>% mutate(nCells = nUMI * (dplyr::filter(lmr, name == repList[[j]][[5]]) %>% .$coef))
#   DOT1L_B_Temp <- combinedSampleTable %>% dplyr::filter(sampleNum == repList[[j]][[6]]) %>% rowwise() %>% mutate(nCells = nUMI * (dplyr::filter(lmr, name == repList[[j]][[6]]) %>% .$coef))
# 
#   DMSOTablesList[[j]] <- full_join(DMSO_A_Temp, DMSO_B_Temp, by = "BC50") %>% mutate(nCells.x = ifelse(is.na(nCells.x), 0, nCells.x), nCells.y = ifelse(is.na(nCells.y), 0, nCells.y)) %>% rowwise() %>% mutate(nCellsMax = max(nCells.x, nCells.y))
#   LSD1TablesList[[j]] <- full_join(LSD1_A_Temp, LSD1_B_Temp, by = "BC50") %>% mutate(nCells.x = ifelse(is.na(nCells.x), 0, nCells.x), nCells.y = ifelse(is.na(nCells.y), 0, nCells.y)) %>% rowwise() %>% mutate(nCellsMax = max(nCells.x, nCells.y))
#   DOT1LTablesList[[j]] <- full_join(DOT1L_A_Temp, DOT1L_B_Temp, by = "BC50") %>% mutate(nCells.x = ifelse(is.na(nCells.x), 0, nCells.x), nCells.y = ifelse(is.na(nCells.y), 0, nCells.y)) %>% rowwise() %>% mutate(nCellsMax = max(nCells.x, nCells.y))
# }
# 
# saveRDS(DMSOTablesList, file = paste0(homeDirectory, "DMSOTablesList.rds"))
# saveRDS(LSD1TablesList, file = paste0(homeDirectory, "LSD1TablesList.rds"))
# saveRDS(DOT1LTablesList, file = paste0(homeDirectory, "DOT1LTablesList.rds"))

rep1 <- c("TM28_DMSO_2A", "TM28_DMSO_2B", "TM28_LSD1i_2A", "TM28_LSD1i_2B", "TM28_DOT1Li_2A", "TM28_DOT1Li_2B")
rep2 <- c("TM28_DMSO_3A", "TM28_DMSO_3B", "TM28_LSD1i_3A", "TM28_LSD1i_3B", "TM28_DOT1Li_3A", "TM28_DOT1Li_3B")
repList <- list(rep1, rep2)

DMSOTablesList <- readRDS(file = paste0(homeDirectory, "DMSOTablesList.rds"))
LSD1TablesList <- readRDS(file = paste0(homeDirectory, "LSD1TablesList.rds"))
DOT1LTablesList <- readRDS(file = paste0(homeDirectory, "DOT1LTablesList.rds"))

#### compare samples for each replicate ####
#### LSD1 ##########################################################################################################################
cutoff = c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)
foldChangecutoffDep = 2
foldChangecutoffInd = 2

numIndList <- list()
numDepXList <- list()
numDepYList <- list()
ratioDepYIndList <- list()

for(j in 1:length(repList)) {
  DMSO_Temp <- DMSOTablesList[[j]]
  LSD1_Temp <- LSD1TablesList[[j]]
  DOT1L_Temp <- DOT1LTablesList[[j]]
  
  barcodesOverlap <- full_join(DMSO_Temp %>% dplyr::select(BC50, nCellsMax),
                               LSD1_Temp %>% dplyr::select(BC50, nCellsMax), by = "BC50") %>% replace(is.na(.), 0) %>%
    mutate(foldchange = log((nCellsMax.y + 1) / (nCellsMax.x + 1)))
  
  numInd <- list()
  numDepX <- list()
  numDepY <- list()
  ratioDepYInd <- list()
  for(i in 1:length(cutoff)) {
    drugDepX <- barcodesOverlap %>% filter(foldchange < -foldChangecutoffDep | foldchange > foldChangecutoffDep) %>% filter(nCellsMax.x > cutoff[i] & nCellsMax.y < cutoff[i])
    drugDepY <- barcodesOverlap %>% filter(foldchange < -foldChangecutoffDep | foldchange > foldChangecutoffDep) %>% filter(nCellsMax.y > cutoff[i] & nCellsMax.x < cutoff[i])
    drugInd <- barcodesOverlap %>% filter(foldchange > -foldChangecutoffInd & foldchange < foldChangecutoffInd) %>% filter(nCellsMax.x > cutoff[i] & nCellsMax.y > cutoff[i])
    
    numInd[[i]] <- nrow(drugInd)
    numDepX[[i]] <- nrow(drugDepX)
    numDepY[[i]] <- nrow(drugDepY)
    ratioDepYInd[[i]] <- nrow(drugDepY)/(nrow(drugInd) + nrow(drugDepX))
  }
  
  numIndList[[j]] <- numInd
  numDepXList[[j]] <- numDepX
  numDepYList[[j]] <- numDepY
  ratioDepYIndList[[j]] <- ratioDepYInd
}

####plot example scatterplot for figure########################################################################################################################################
###############################################################################################################################################################################
j = 1
i = 2

DMSO_Temp <- DMSOTablesList[[j]]
LSD1_Temp <- LSD1TablesList[[j]]

barcodesOverlap <- full_join(DMSO_Temp %>% dplyr::select(BC50, nCellsMax),
                             LSD1_Temp %>% dplyr::select(BC50, nCellsMax), by = "BC50") %>% replace(is.na(.), 0) %>%
  mutate(foldchange = log((nCellsMax.y + 1) / (nCellsMax.x + 1)))

drugDepX <- barcodesOverlap %>% filter(foldchange < -foldChangecutoffDep | foldchange > foldChangecutoffDep) %>% filter(nCellsMax.x > cutoff[i] & nCellsMax.y < cutoff[i])
drugDepY <- barcodesOverlap %>% filter(foldchange < -foldChangecutoffDep | foldchange > foldChangecutoffDep) %>% filter(nCellsMax.y > cutoff[i] & nCellsMax.x < cutoff[i])
drugInd <- barcodesOverlap %>% filter(foldchange > -foldChangecutoffInd & foldchange < foldChangecutoffInd) %>% filter(nCellsMax.x > cutoff[i] & nCellsMax.y > cutoff[i])

ggplot() +
  rasterise(geom_point(data = barcodesOverlap, aes(x = nCellsMax.x, y = nCellsMax.y), color = "lightgray", shape = 16), dpi = 300) +
  geom_point(data = drugDepX, aes(x = nCellsMax.x, y = nCellsMax.y), color = "red", shape = 16) +
  geom_point(data = drugDepY, aes(x = nCellsMax.x, y = nCellsMax.y), color = "green", shape = 16) +
  geom_point(data = drugInd, aes(x = nCellsMax.x, y = nCellsMax.y), color = "black", shape = 16) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") + coord_fixed(ratio = 1, xlim = c(0, 175000), ylim = c(0, 175000))
  # scale_x_continuous(trans = "log10", limits = c(0.1, 175000)) + scale_y_continuous(trans = "log10", limits = c(0.1, 175000))
ggsave(filename = paste0(plotDirectory, 'DMSOvsLSD1i_spikeIn_scatterplotExampleForFigure_corrected.pdf'), width = 4, height = 4, useDingbats = F)

ggplot() +
  rasterise(geom_point(data = barcodesOverlap, aes(x = nCellsMax.x, y = nCellsMax.y), color = "lightgray", shape = 16), dpi = 300) +
  geom_point(data = drugDepX, aes(x = nCellsMax.x, y = nCellsMax.y), color = "red", shape = 16) +
  geom_point(data = drugDepY, aes(x = nCellsMax.x, y = nCellsMax.y), color = "green", shape = 16) +
  geom_point(data = drugInd, aes(x = nCellsMax.x, y = nCellsMax.y), color = "black", shape = 16) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") + coord_fixed(ratio = 1, xlim = c(0, 30000), ylim = c(0, 30000))
# scale_x_continuous(trans = "log10", limits = c(0.1, 175000)) + scale_y_continuous(trans = "log10", limits = c(0.1, 175000))
ggsave(filename = paste0(plotDirectory, 'DMSOvsLSD1i_spikeIn_scatterplotExampleForFigure_corrected_zoomIn.pdf'), width = 4, height = 4, useDingbats = F)

DMSO_Temp$nUMI.x[is.na(DMSO_Temp$nUMI.x)] <- 0
DMSO_Temp$nUMI.y[is.na(DMSO_Temp$nUMI.y)] <- 0
DMSO_Temp <- DMSO_Temp %>% mutate(nUMIRPM.x = nUMI.x/sum(DMSO_Temp$nUMI.x)*10^6,
                                  nUMIRPM.y = nUMI.y/sum(DMSO_Temp$nUMI.y)*10^6) %>%
  mutate(nUMIRPMMax = max(nUMIRPM.x, nUMIRPM.y))

LSD1_Temp$nUMI.x[is.na(LSD1_Temp$nUMI.x)] <- 0
LSD1_Temp$nUMI.y[is.na(LSD1_Temp$nUMI.y)] <- 0
LSD1_Temp <- LSD1_Temp %>% mutate(nUMIRPM.x = nUMI.x/sum(LSD1_Temp$nUMI.x)*10^6,
                                  nUMIRPM.y = nUMI.y/sum(LSD1_Temp$nUMI.y)*10^6) %>%
  mutate(nUMIRPMMax = max(nUMIRPM.x, nUMIRPM.y))

barcodesOverlap <- full_join(DMSO_Temp %>% dplyr::select(BC50, nUMIRPMMax),
                             LSD1_Temp %>% dplyr::select(BC50, nUMIRPMMax), by = "BC50") %>% replace(is.na(.), 0) %>%
  mutate(foldchange = log((nUMIRPMMax.y + 1) / (nUMIRPMMax.x + 1)))

drugDepX <- barcodesOverlap %>% filter(foldchange < -foldChangecutoffDep | foldchange > foldChangecutoffDep) %>% filter(nUMIRPMMax.x > cutoff[i] & nUMIRPMMax.y < cutoff[i])
drugDepY <- barcodesOverlap %>% filter(foldchange < -foldChangecutoffDep | foldchange > foldChangecutoffDep) %>% filter(nUMIRPMMax.y > cutoff[i] & nUMIRPMMax.x < cutoff[i])
drugInd <- barcodesOverlap %>% filter(foldchange > -foldChangecutoffDep & foldchange < foldChangecutoffDep) %>% filter(nUMIRPMMax.x > cutoff[i] & nUMIRPMMax.y > cutoff[i])

ggplot() +
  rasterise(geom_point(data = barcodesOverlap, aes(x = nUMIRPMMax.x, y = nUMIRPMMax.y), color = "lightgray", shape = 16), dpi = 300) +
  geom_point(data = drugDepX, aes(x = nUMIRPMMax.x, y = nUMIRPMMax.y), color = "red", shape = 16) +
  geom_point(data = drugDepY, aes(x = nUMIRPMMax.x, y = nUMIRPMMax.y), color = "green", shape = 16) +
  geom_point(data = drugInd, aes(x = nUMIRPMMax.x, y = nUMIRPMMax.y), color = "black", shape = 16) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") + coord_fixed(ratio = 1, xlim = c(0, 120000), ylim = c(0, 120000))
  # scale_x_continuous(trans = "log10", limits = c(0.1, 175000)) + scale_y_continuous(trans = "log10", limits = c(0.1, 175000))
ggsave(filename = paste0(plotDirectory, 'DMSOvsLSD1i_spikeIn_scatterplotExampleForFigure_uncorrected.pdf'), width = 4, height = 4, useDingbats = F)

ggplot() +
  rasterise(geom_point(data = barcodesOverlap, aes(x = nUMIRPMMax.x, y = nUMIRPMMax.y), color = "lightgray", shape = 16), dpi = 300) +
  geom_point(data = drugDepX, aes(x = nUMIRPMMax.x, y = nUMIRPMMax.y), color = "red", shape = 16) +
  geom_point(data = drugDepY, aes(x = nUMIRPMMax.x, y = nUMIRPMMax.y), color = "green", shape = 16) +
  geom_point(data = drugInd, aes(x = nUMIRPMMax.x, y = nUMIRPMMax.y), color = "black", shape = 16) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") + coord_fixed(ratio = 1, xlim = c(0, 30000), ylim = c(0, 30000))
# scale_x_continuous(trans = "log10", limits = c(0.1, 175000)) + scale_y_continuous(trans = "log10", limits = c(0.1, 175000))
ggsave(filename = paste0(plotDirectory, 'DMSOvsLSD1i_spikeIn_scatterplotExampleForFigure_uncorrected_zoomIn.pdf'), width = 4, height = 4, useDingbats = F)
###############################################################################################################################################################################

aggregateSampleTable <- tibble(rep = c(rep("1", length(cutoff)), rep("2", length(cutoff))),
                               cutoff = rep(cutoff, 2),
                               ratio = c(ratioDepYIndList[[1]], ratioDepYIndList[[2]]),
                               numInd = c(numIndList[[1]], numIndList[[2]]),
                               numDepX = c(numDepXList[[1]], numDepXList[[2]]),
                               numDepY = c(numDepYList[[1]], numDepYList[[2]]))
aggregateSampleTable <- aggregateSampleTable %>% rowwise() %>% mutate(nCells = numInd + numDepX + numDepY, area1 = numInd + numDepY, area2 = numInd + numDepX)

ggplot(aggregateSampleTable, aes(x = cutoff, y = as.double(ratio))) +
  geom_point() +
  stat_summary(fun = mean, geom = "point", size = 2, color = "blue") +
  stat_summary(fun.data = mean_se, geom = "errorbar", color = "blue") + ylim(0, 9) +
  geom_hline(yintercept = 1, linetype = "dashed")
ggsave(filename = paste0(plotDirectory, 'DMSOvsLSD1i_spikeIn_foldchangeNumLineagesVaryingCutoff.pdf'), width = 5, height = 2, useDingbats = F)

aggregateSampleTableFilter <- aggregateSampleTable %>% dplyr::filter(cutoff == 200)
ggplot(aggregateSampleTableFilter, aes(x = as.factor(1), y = as.double(ratio))) +
  geom_jitter(width = 0.1, height = 0.1) +
  stat_summary(fun = mean, geom = "point", size = 4) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.1) +
  ylim(0, 8.1) + geom_hline(yintercept = 1, linetype = "dashed")
ggsave(filename = paste0(plotDirectory, 'DMSOvsLSD1i_spikeIn_foldchangeNumLineages.pdf'), width = 1.5, height = 3, useDingbats = F)

grid.newpage()
venndiagram <- draw.pairwise.venn(area1 = aggregateSampleTableFilter$numInd[[1]] + aggregateSampleTableFilter$numDepY[[1]], area2 = aggregateSampleTableFilter$numInd[[1]] + aggregateSampleTableFilter$numDepX[[1]], cross.area = aggregateSampleTableFilter$numInd[[1]],
                                  euler.d = TRUE, fill = c(alpha("gray5", 0.5), alpha("gray5", 0.5)),
                                  fontfamily = "sans", cex = 2.5, print.mode = "percent", sigdigs = 3)
ggsave(plot = venndiagram, filename = paste0(plotDirectory, 'overlapVennDiagram_LSD1_spikeIn_1.pdf'), width = 5, height = 5, useDingbats = F)

grid.newpage()
venndiagram <- draw.pairwise.venn(area1 = aggregateSampleTableFilter$numInd[[2]] + aggregateSampleTableFilter$numDepY[[2]], area2 = aggregateSampleTableFilter$numInd[[2]] + aggregateSampleTableFilter$numDepX[[2]], cross.area = aggregateSampleTableFilter$numInd[[2]],
                                  euler.d = TRUE, fill = c(alpha("gray5", 0.5), alpha("gray5", 0.5)),
                                  fontfamily = "sans", cex = 2.5, print.mode = "percent", sigdigs = 3)
ggsave(plot = venndiagram, filename = paste0(plotDirectory, 'overlapVennDiagram_LSD1_spikeIn_2.pdf'), width = 5, height = 5, useDingbats = F)

#### demonstrate being LSD1i dependent is heritable ################################################################################
####################################################################################################################################
rep1 <- 1
rep2 <- 2

repList <- list(rep1, rep2)
cutoff <- c(200, 200)
foldChangecutoffDep <- c(2, 2)

for(j in 1:length(repList)) {
  if(length(repList[[j]]) == 4) {
    DMSO_A_Temp <- sampleTablesList[[repList[[j]][[1]][1]]][[repList[[j]][[1]][2]]] %>% mutate(V2 = V2/(sum(V2)) * 10^6)
    DMSO_B_Temp <- sampleTablesList[[repList[[j]][[2]][1]]][[repList[[j]][[2]][2]]] %>% mutate(V2 = V2/(sum(V2)) * 10^6)
    DMSO_Temp <- full_join(DMSO_A_Temp, DMSO_B_Temp, by = "V1") %>% replace(is.na(.), 0) %>% rowwise() %>% mutate(max = max(V2.x, V2.y))
    
    LSD1_A_Temp <- sampleTablesList[[repList[[j]][[3]][1]]][[repList[[j]][[3]][2]]] %>% mutate(V2 = V2/(sum(V2)) * 10^6) %>% dplyr::rename(max = V2)
    LSD1_B_Temp <- sampleTablesList[[repList[[j]][[4]][1]]][[repList[[j]][[4]][2]]] %>% mutate(V2 = V2/(sum(V2)) * 10^6) %>% dplyr::rename(max = V2)
  } else{
    DMSO_Temp <- DMSOTablesList[[repList[[j]]]] %>% dplyr::rename(V1 = BC50, max = nCellsMax)
    DMSO_sampleNames <- c(DMSO_Temp$sampleNum.x %>% unique(), DMSO_Temp$sampleNum.y %>% unique()) %>% .[!is.na(.)]
    DMSO_A_Temp <- DMSO_Temp %>% dplyr::filter(sampleNum.x == DMSO_sampleNames[1])
    DMSO_B_Temp <- DMSO_Temp %>% dplyr::filter(sampleNum.y == DMSO_sampleNames[2])
    
    LSD1_Temp <- LSD1TablesList[[repList[[j]]]] %>% dplyr::rename(V1 = BC50)
    LSD1_sampleNames <- c(LSD1_Temp$sampleNum.x %>% unique(), LSD1_Temp$sampleNum.y %>% unique()) %>% .[!is.na(.)]
    LSD1_A_Temp <- LSD1_Temp %>% dplyr::filter(sampleNum.x == LSD1_sampleNames[1]) %>% dplyr::rename(max = nCells.x)
    LSD1_B_Temp <- LSD1_Temp %>% dplyr::filter(sampleNum.y == LSD1_sampleNames[2]) %>% dplyr::rename(max = nCells.y)
  }
  DMSO_memoryOverlap <- tibble(ind = nrow(inner_join(DMSO_A_Temp, DMSO_B_Temp, by = "V1")),
                               depX = nrow(left_join(DMSO_A_Temp, DMSO_B_Temp, by = "V1")),
                               depY = nrow(left_join(DMSO_B_Temp, DMSO_A_Temp, by = "V1"))) %>%
    rowwise() %>% mutate(nCells = sum(ind, depX, depY), area1 = ind + depX, area2 = ind + depY, overlap = ind/nCells, cond = "DMSO")
  
  barcodesOverlap_A <- full_join(DMSO_Temp %>% dplyr::select(V1, max),
                                 LSD1_A_Temp %>% dplyr::select(V1, max), by = "V1") %>% replace(is.na(.), 0)
  if(length(repList[[j]]) == 4) {
    correctionFactor_A <- barcodesOverlap_A %>% ungroup() %>% slice_max(order_by = max.x, n = 5) %>%
      rowwise() %>% mutate(foldchange = max.x/max.y) %>% .$foldchange %>% mean()
    barcodesOverlap_A <- barcodesOverlap_A %>% rowwise() %>% mutate(max.y_corr = max.y * correctionFactor_A) %>%
      mutate(foldchange = log((max.y_corr + 1) / (max.x + 1)))
    drugDep_A <- barcodesOverlap_A %>% filter(foldchange < -foldChangecutoffDep[j] | foldchange > foldChangecutoffDep[j]) %>% filter(max.y_corr > cutoff[j] & max.x < cutoff[j])
  } else{
    barcodesOverlap_A <- barcodesOverlap_A %>% rowwise() %>%mutate(foldchange = log((max.y + 1) / (max.x + 1)))
    drugDep_A <- barcodesOverlap_A %>% filter(foldchange < -foldChangecutoffDep[j] | foldchange > foldChangecutoffDep[j]) %>% filter(max.y > cutoff[j] & max.x < cutoff[j])
  }
  barcodesOverlap_B <- full_join(DMSO_Temp %>% dplyr::select(V1, max),
                                 LSD1_B_Temp %>% dplyr::select(V1, max), by = "V1") %>% replace(is.na(.), 0)
  if(length(repList[[j]]) == 4) {
    correctionFactor_B <- barcodesOverlap_B %>% ungroup() %>% slice_max(order_by = max.x, n = 5) %>%
      rowwise() %>% mutate(foldchange = max.x/max.y) %>% .$foldchange %>% mean()
    barcodesOverlap_B <- barcodesOverlap_B %>% rowwise() %>% mutate(max.y_corr = max.y * correctionFactor_B) %>%
      mutate(foldchange = log((max.y_corr + 1) / (max.x + 1)))
    drugDep_B <- barcodesOverlap_B %>% filter(foldchange < -foldChangecutoffDep[j] | foldchange > foldChangecutoffDep[j]) %>% filter(max.y_corr > cutoff[j] & max.x < cutoff[j])
  } else{
    barcodesOverlap_B <- barcodesOverlap_B %>% rowwise() %>%mutate(foldchange = log((max.y + 1) / (max.x + 1)))
    drugDep_B <- barcodesOverlap_B %>% filter(foldchange < -foldChangecutoffDep[j] | foldchange > foldChangecutoffDep[j]) %>% filter(max.y > cutoff[j] & max.x < cutoff[j])
  }
    
  LSD1_memoryOverlap <- tibble(ind = nrow(inner_join(drugDep_A, drugDep_B, by = "V1")),
                               depX = nrow(left_join(drugDep_A, drugDep_B, by = "V1")),
                               depY = nrow(left_join(drugDep_B, drugDep_A, by = "V1"))) %>%
    rowwise() %>% mutate(nCells = sum(ind, depX, depY), area1 = ind + depX, area2 = ind + depY, overlap = ind/nCells, cond = "LSD1")
  
  memoryOverlapTableTemp <- bind_rows(DMSO_memoryOverlap, LSD1_memoryOverlap) %>% mutate(rep = paste0(j))
  if(j == 1) {
    memoryOverlapTable1 <- memoryOverlapTableTemp
  } else{
    memoryOverlapTable1 <- bind_rows(memoryOverlapTable1, memoryOverlapTableTemp)
  }
}

ggplot(memoryOverlapTable1, aes(x = cond, y = overlap)) +
  stat_summary(fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.25) +
  geom_signif(comparisons = list(c("DMSO", "LSD1"))) +
  geom_jitter(width = 0.1)

#### demonstrate being LSD1i dependent and DOT1Li dependent are distinct ###########################################################
####################################################################################################################################
rep1 <- list(c(1, 1), c(1, 2), c(1, 5), c(1, 6), c(1, 3))
rep2 <- list(c(2, 1), c(2, 3), c(2, 2))
rep3 <- 1

repList <- list(rep1, rep2, rep3)
cutoff <- c(500, 500, 200)
foldChangecutoffDep <- c(2, 2, 2)

for(j in 1:length(repList)) {
  if(j %in% c(1, 2)) {
    if(length(repList[[j]]) == 5) {
      DMSO_A_Temp <- sampleTablesList[[repList[[j]][[1]][1]]][[repList[[j]][[1]][2]]] %>% mutate(V2 = V2/(sum(V2)) * 10^6)
      DMSO_B_Temp <- sampleTablesList[[repList[[j]][[2]][1]]][[repList[[j]][[2]][2]]] %>% mutate(V2 = V2/(sum(V2)) * 10^6)
      DMSO_Temp <- full_join(DMSO_A_Temp, DMSO_B_Temp, by = "V1") %>% replace(is.na(.), 0) %>% rowwise() %>% mutate(max = max(V2.x, V2.y))
      
      LSD1_A_Temp <- sampleTablesList[[repList[[j]][[3]][1]]][[repList[[j]][[3]][2]]] %>% mutate(V2 = V2/(sum(V2)) * 10^6)
      LSD1_B_Temp <- sampleTablesList[[repList[[j]][[4]][1]]][[repList[[j]][[4]][2]]] %>% mutate(V2 = V2/(sum(V2)) * 10^6)
      LSD1_Temp <- full_join(LSD1_A_Temp, LSD1_B_Temp, by = "V1") %>% replace(is.na(.), 0) %>% rowwise() %>% mutate(max = max(V2.x, V2.y))
      
      DOT1L_Temp <- sampleTablesList[[repList[[j]][[5]][1]]][[repList[[j]][[5]][2]]] %>% mutate(V2 = V2/(sum(V2)) * 10^6) %>% dplyr::rename(max = V2)
    } else{
      DMSO_Temp <- sampleTablesList[[repList[[j]][[1]][1]]][[repList[[j]][[1]][2]]] %>% mutate(V2 = V2/(sum(V2)) * 10^6) %>% dplyr::rename(max = V2)
      LSD1_Temp <- sampleTablesList[[repList[[j]][[2]][1]]][[repList[[j]][[2]][2]]] %>% mutate(V2 = V2/(sum(V2)) * 10^6) %>% dplyr::rename(max = V2)
      DOT1L_Temp <- sampleTablesList[[repList[[j]][[3]][1]]][[repList[[j]][[3]][2]]] %>% mutate(V2 = V2/(sum(V2)) * 10^6) %>% dplyr::rename(max = V2)
    }
  } else{
    DMSO_Temp <- DMSOTablesList[[repList[[j]]]] %>% dplyr::rename(V1 = BC50, max = nCellsMax)
    LSD1_Temp <- LSD1TablesList[[repList[[j]]]] %>% dplyr::rename(V1 = BC50, max = nCellsMax)
    DOT1L_Temp <- DOT1LTablesList[[repList[[j]]]] %>% dplyr::rename(V1 = BC50, max = nCellsMax)
  }
  barcodesOverlap_A <- full_join(DMSO_Temp %>% dplyr::select(V1, max),
                                 LSD1_Temp %>% dplyr::select(V1, max), by = "V1") %>% replace(is.na(.), 0)
  if(j %in% c(1, 2)) {
    correctionFactor_A <- barcodesOverlap_A %>% ungroup() %>% slice_max(order_by = max.x, n = 5) %>%
      rowwise() %>% mutate(foldchange = max.x/max.y) %>% .$foldchange %>% mean()
    barcodesOverlap_A <- barcodesOverlap_A %>% rowwise() %>% mutate(max.y_corr = max.y * correctionFactor_A) %>%
      mutate(foldchange = log((max.y_corr + 1) / (max.x + 1)))
    drugDep_A <- barcodesOverlap_A %>% filter(foldchange < -foldChangecutoffDep[j] | foldchange > foldChangecutoffDep[j]) %>% filter(max.y_corr > cutoff[j] & max.x < cutoff[j])
  } else{
    barcodesOverlap_A <- barcodesOverlap_A %>% rowwise() %>%mutate(foldchange = log((max.y + 1) / (max.x + 1)))
    drugDep_A <- barcodesOverlap_A %>% filter(foldchange < -foldChangecutoffDep[j] | foldchange > foldChangecutoffDep[j]) %>% filter(max.y > cutoff[j] & max.x < cutoff[j])
  }
  barcodesOverlap_B <- full_join(DMSO_Temp %>% dplyr::select(V1, max),
                                 DOT1L_Temp %>% dplyr::select(V1, max), by = "V1") %>% replace(is.na(.), 0)
  if(j %in% c(1, 2)) {
    correctionFactor_B <- barcodesOverlap_B %>% ungroup() %>% slice_max(order_by = max.x, n = 5) %>%
      rowwise() %>% mutate(foldchange = max.x/max.y) %>% .$foldchange %>% mean()
    barcodesOverlap_B <- barcodesOverlap_B %>% rowwise() %>% mutate(max.y_corr = max.y * correctionFactor_B) %>%
      mutate(foldchange = log((max.y_corr + 1) / (max.x + 1)))
    drugDep_B <- barcodesOverlap_B %>% filter(foldchange < -foldChangecutoffDep[j] | foldchange > foldChangecutoffDep[j]) %>% filter(max.y_corr > cutoff[j] & max.x < cutoff[j])
  } else{
    barcodesOverlap_B <- barcodesOverlap_B %>% rowwise() %>%mutate(foldchange = log((max.y + 1) / (max.x + 1)))
    drugDep_B <- barcodesOverlap_B %>% filter(foldchange < -foldChangecutoffDep[j] | foldchange > foldChangecutoffDep[j]) %>% filter(max.y > cutoff[j] & max.x < cutoff[j])
  }
  
  drug_memoryOverlap <- tibble(ind = nrow(inner_join(drugDep_A, drugDep_B, by = "V1")),
                               depX = nrow(left_join(drugDep_A, drugDep_B, by = "V1")),
                               depY = nrow(left_join(drugDep_B, drugDep_A, by = "V1"))) %>%
    rowwise() %>% mutate(nCells = sum(ind, depX, depY), area1 = ind + depX, area2 = ind + depY, overlap = ind/nCells, cond = "DOT1L")
  
  memoryOverlapTableTemp <- drug_memoryOverlap %>% mutate(rep = paste0(j))
  if(j == 1) {
    memoryOverlapTable2 <- memoryOverlapTableTemp
  } else{
    memoryOverlapTable2 <- bind_rows(memoryOverlapTable2, memoryOverlapTableTemp)
  }
}

memoryOverlapTable <- bind_rows(memoryOverlapTable1, memoryOverlapTable2)
memoryOverlapTable$cond <- factor(memoryOverlapTable$cond, levels = c("DMSO", "LSD1", "DOT1L"))

ggplot(memoryOverlapTable, aes(x = cond, y = overlap)) +
  stat_summary(fun = mean, geom = "bar") +
  geom_jitter(height = 0, width = 0.25, size = 2) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.25) +
  #stat_summary(fun = mean, geom = "point", size = 5) +
  geom_signif(comparisons = list(c("DMSO", "LSD1"), c("LSD1", "DOT1L"), c("DMSO", "DOT1L")), test = "t.test")
ggsave(file = paste0(plotDirectory, "barcodeOverlapForLSD1iConditions.pdf"), height = 3, width = 3)
