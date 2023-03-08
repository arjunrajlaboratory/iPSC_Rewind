#### load required libraries and functions----
library(tidyverse)
library(reshape2)
library(ggpubr)
library(fuzzyjoin)
library(VennDiagram)
library(ggforce)

theme_set(theme_classic())

#### load dataset into data frames----
sampleDir1 <- '/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/extractedData/barcodeOverlap/hiFTTM1/'
sampleDir2 <- '/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/extractedData/barcodeOverlap/hiFTTM4/'
sampleDir3 <- '/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/extractedData/barcodeOverlap/hiFTTM21/'
sampleDirList <- c(sampleDir1, sampleDir2, sampleDir3)

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

outputDir <- '/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/plots/barcodeOverlap/'

#### compare heritability overlap----
colonyNumberList <- c(300, 400, 300)
overlapScatterList <- list()
numOverlap <- rep(0, length(sampleDirList))
numA <- rep(0, length(sampleDirList))
numB <- rep(0, length(sampleDirList))

for (i in 1:length(sampleDirList)) {
    colonyNumber <- colonyNumberList[i]
    filterThreshold = 0
    a <- sampleTablesList[[i]][[1]] %>% filter(V2>filterThreshold) %>% mutate(V2=V2/sum(V2)*10^6) %>% mutate(V2=V2/sum(V2)*10^6)
    b <- sampleTablesList[[i]][[2]] %>% filter(V2>filterThreshold) %>% mutate(V2=V2/sum(V2)*10^6) %>% mutate(V2=V2/sum(V2)*10^6)
    
    aTopTable <- slice_max(a, n = colonyNumber, V2, with_ties = FALSE)
    bTopTable <- slice_max(b, n = colonyNumber, V2, with_ties = FALSE)
    
    overlapScatter <- full_join(aTopTable, bTopTable, by = "V1") %>% replace(is.na(.), 0) %>%
      mutate(foldchange = log2(V2.y/V2.x))
    
    foldChangecutoffInd <- 2.5
    
    drugDepX <- overlapScatter %>% filter(foldchange < -foldChangecutoffInd)
    drugDepY <- overlapScatter %>% filter(foldchange > foldChangecutoffInd)
    drugInd <- overlapScatter %>% filter(foldchange > -foldChangecutoffInd & foldchange < foldChangecutoffInd)
    
    numOverlap[i] <- nrow(drugInd)
    numA[i] <- nrow(drugDepX) + nrow(drugInd)
    numB[i] <- nrow(drugDepY) + nrow(drugInd)
    
    overlapScatterList[[i]] <- ggplot() +
      geom_point(data = overlapScatter, aes(x = V2.x, y = V2.y), color = "#CACACA") +
      geom_point(data = drugInd, aes(x = V2.x, y = V2.y), color = "#9F9F9F") +
      coord_fixed(ratio = 1) +
      geom_abline(slope = 1, intercept = 0, linetype = 'dotted', color = 'black') +
      xlab('barcode reads per million in split A') + ylab('barcode reads per million in split B')
}

#### compare heritability overlap for different values----
colonyCutoff <- c(25, 50, 100, 250, 500, 1000)
overlapValMatrix <- matrix(ncol = length(sampleDirList), nrow = length(colonyCutoff))

for (i in 1:length(sampleDirList)) {
  filterThreshold = 0
  a <- sampleTablesList[[i]][[1]] %>% filter(V2>filterThreshold) %>% mutate(V2=V2/sum(V2)*10^6) %>% mutate(V2=V2/sum(V2)*10^6)
  b <- sampleTablesList[[i]][[2]] %>% filter(V2>filterThreshold) %>% mutate(V2=V2/sum(V2)*10^6) %>% mutate(V2=V2/sum(V2)*10^6)
  
  for(j in 1:length(colonyCutoff)){
    colonyNumber <- colonyCutoff[j]
    
    aTopTable <- slice_max(a, n = colonyNumber, V2, with_ties = FALSE)
    bTopTable <- slice_max(b, n = colonyNumber, V2, with_ties = FALSE)
    
    overlapScatter <- full_join(aTopTable, bTopTable, by = "V1") %>% replace(is.na(.), 0) %>%
      mutate(foldchange = log2(V2.y/V2.x))
    
    foldChangecutoffInd <- 2.5
    
    drugDepX <- overlapScatter %>% filter(foldchange < -foldChangecutoffInd)
    drugDepY <- overlapScatter %>% filter(foldchange > foldChangecutoffInd)
    drugInd <- overlapScatter %>% filter(foldchange >= -foldChangecutoffInd & foldchange <= foldChangecutoffInd)
    
    overlapValMatrix[j, i] <- nrow(drugInd) / (nrow(drugDepX) + nrow(drugDepY) + nrow(drugInd))
  }
}

overlapValTable <- overlapValMatrix %>% as.data.frame() %>% mutate(cutoff = as.character(colonyCutoff)) %>% melt(id = "cutoff")
overlapValTable$cutoff <- factor(overlapValTable$cutoff, levels = c("25", "50", "100", "250", "500", "1000"))
overlapPlot <- ggplot(overlapValTable, aes(x = cutoff, y = value)) +
  geom_jitter(height = 0, width = 0) + geom_line(aes(group = variable)) + ylim(0, 1) +
  stat_summary(fun = mean, geom = "point", color = "blue", size = 2.5) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", color = "blue", size = 1, width = 0.25) +
  ylab("fraction of overlapping barcodes between split A and split B") +
  xlab("cutoff for number of barcoded cells to be called an iPSC colony")
ggsave(plot = overlapPlot, filename = paste0(outputDir, 'overlapAcrossSamplesWithChangingCutoff.pdf'), width = 5, height = 5, useDingbats = F)

#### plot VennDiagram of average barcode overlap ----
grid.newpage()
venndiagram <- draw.pairwise.venn(area1 = numA[[1]], area2 = numB[[1]], cross.area = numOverlap[[1]],
                                  euler.d = TRUE, fill = c(alpha("gray5", 0.5), alpha("gray5", 0.5)),
                                  fontfamily = "sans", cex = 2.5, print.mode = "percent", sigdigs = 3)
ggsave(plot = venndiagram, filename = paste0(outputDir, 'overlapVennDiagram_2.pdf'), width = 5, height = 5, useDingbats = F)

grid.newpage()
venndiagram <- draw.pairwise.venn(area1 = numA[[2]], area2 = numB[[2]], cross.area = numOverlap[[2]],
                                  euler.d = TRUE, fill = c(alpha("gray5", 0.5), alpha("gray5", 0.5)),
                                  fontfamily = "sans", cex = 2.5, print.mode = "percent", sigdigs = 3)
ggsave(plot = venndiagram, filename = paste0(outputDir, 'overlapVennDiagram_3.pdf'), width = 5, height = 5, useDingbats = F)

grid.newpage()
venndiagram <- draw.pairwise.venn(area1 = numA[[3]], area2 = numB[[3]], cross.area = numOverlap[[3]],
                                  euler.d = TRUE, fill = c(alpha("gray5", 0.5), alpha("gray5", 0.5)),
                                  fontfamily = "sans", cex = 2.5, print.mode = "percent", sigdigs = 3)
ggsave(plot = venndiagram, filename = paste0(outputDir, 'overlapVennDiagram_1.pdf'), width = 5, height = 5, useDingbats = F)

#### graph scatterplots of heritability controls A vs. B----
ggarrange(overlapScatterList[[1]], overlapScatterList[[2]], overlapScatterList[[3]], ncol = 3)
#ggsave(filename = paste0(outputDir, 'overlapScatterPlot.pdf'), width = 25, height = 5, useDingbats = F)

scatterPlot1 <- overlapScatterList[[1]] + ylim(0, 15000) + xlim(0, 15000)
scatterPlot3 <- overlapScatterList[[2]] + 
  facet_zoom(x = V2.x < 17500, y = V2.y < 17500)
scatterPlot5 <- overlapScatterList[[3]] + 
  facet_zoom(x = V2.x < 25000, y = V2.y < 25000)

ggsave(scatterPlot1, filename = paste0(outputDir, 'overlapScatterPlot_2.pdf'), width = 5, height = 5, useDingbats = F)
ggsave(scatterPlot3, filename = paste0(outputDir, 'overlapScatterPlot_3.pdf'), width = 25/3, height = 5, useDingbats = F)
ggsave(scatterPlot5, filename = paste0(outputDir, 'overlapScatterPlot_1.pdf'), width = 25/3, height = 5, useDingbats = F)

#### compare size correlation----
colonyCutoff <- c(25, 50, 100, 250, 500, 1000)
coefficientList <- matrix(ncol = length(sampleDirList), nrow = length(colonyCutoff))
rsquaredList <- matrix(ncol = length(sampleDirList), nrow = length(colonyCutoff))

for (i in 1:length(sampleDirList)) {
  for(j in 1:length(colonyCutoff)) {
    colonyNumber <- colonyCutoff[[j]]
    filterThreshold = 0
    a <- sampleTablesList[[i]][[1]] %>% filter(V2>filterThreshold) %>% mutate(V2=V2/sum(V2)*10^6) %>% mutate(V2=V2/sum(V2)*10^6)
    b <- sampleTablesList[[i]][[2]] %>% filter(V2>filterThreshold) %>% mutate(V2=V2/sum(V2)*10^6) %>% mutate(V2=V2/sum(V2)*10^6)
    
    aTopTable <- slice_max(a, n = colonyNumber, V2, with_ties = FALSE)
    bTopTable <- slice_max(b, n = colonyNumber, V2, with_ties = FALSE)
    
    overlapScatter <- inner_join(aTopTable, bTopTable, by = "V1") %>% replace(is.na(.), 0)
    
    model <- summary(lm(V2.y ~ 0 + V2.x, data = overlapScatter))
    
    coefficientList[j, i] <- model$coefficients[[1]]
    rsquaredList[j, i] <- model$r.squared
  }
}

coefficientTable <- coefficientList %>% as.data.frame() %>% mutate(cutoff = as.character(colonyCutoff)) %>% melt(id = "cutoff")
coefficientTable$cutoff <- factor(coefficientTable$cutoff, levels = c("25", "50", "100", "250", "500", "1000"))
cofficientPlot <- ggplot(coefficientTable, aes(x = cutoff, y = value)) +
  geom_jitter(height = 0, width = 0) + geom_line(aes(group = variable)) + ylim(0, 2) +
  stat_summary(fun = mean, geom = "point", color = "blue", size = 2.5) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", color = "blue", size = 1, width = 0.25) +
  ylab("fraction of overlapping barcodes between split A and split B") +
  xlab("cutoff for number of barcoded cells to be called an iPSC colony")
ggsave(plot = cofficientPlot, filename = paste0(outputDir, 'coefficientAcrossSamplesWithChangingCutoff.pdf'), width = 5, height = 5, useDingbats = F)

rsquaredTable <- rsquaredList %>% as.data.frame() %>% mutate(cutoff = as.character(colonyCutoff)) %>% melt(id = "cutoff")
rsquaredTable$cutoff <- factor(rsquaredTable$cutoff, levels = c("25", "50", "100", "250", "500", "1000"))
rsquaredPlot <- ggplot(rsquaredTable, aes(x = cutoff, y = value^(1/2))) +
  geom_jitter(height = 0, width = 0) + geom_line(aes(group = variable)) + ylim(0, 1) +
  stat_summary(fun = mean, geom = "point", color = "blue", size = 2.5) + 
  stat_summary(fun.data = mean_se, geom = "errorbar", color = "blue", size = 1, width = 0.25) +
  ylab(expression("R value for overlapping barcodes between split A and split B")) +
  xlab("cutoff for number of barcoded cells to be called an iPSC colony")
ggsave(plot = rsquaredPlot, filename = paste0(outputDir, 'rSquaredAcrossSamplesWithChangingCutoff.pdf'), width = 5, height = 5, useDingbats = F)
