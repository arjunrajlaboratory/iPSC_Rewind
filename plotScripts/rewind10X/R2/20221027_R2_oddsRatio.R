rm(list=ls())
gc()

library(tidyverse)
library(reshape2)
library(ggrepel)
library(ggforce)
library(ggsignif)
library(ggridges)
library(RColorBrewer)
library(viridis)
library(egg)
library(Seurat)
library(biomaRt)
library(spgs)

theme_set(theme_classic())

homeDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/extractedData/rewind10X/R2/"
plotDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/plots/rewind10X/R2/"

umapCoordinates = as_tibble(read.table(file = paste0(homeDirectory, "umapCoordinates_Scanorama_50pcs_filter.tsv"), header = TRUE, stringsAsFactors = F, sep = "\t"))
umapClusters = as_tibble(read.table(file = paste0(homeDirectory, "umapClusters_Scanorama_50pcs_filter.tsv"), header = TRUE, stringsAsFactors = F, sep = "\t"))
umapClusters <- umapClusters %>% dplyr::rename(cluster = integrated_snn_res.0.45)
primedCells <- readRDS(file = paste0(homeDirectory, "primedCellsInd.rds")) %>% mutate(label = "primed")
primedAllUMAP <- filter(umapCoordinates, cellID %in% primedCells$cellID)
linCountToOverlaps <- as_tibble(read.table(file = paste0(homeDirectory, "filtered10XCells.txt"), header = TRUE, stringsAsFactors = F, sep = "\t"))

#### check odds ratios of being primed for each gene compared with proliferation speed ####
#######################################################################################################################################################
# logNormCounts <- as_tibble(read.table(file = paste0(homeDirectory, "logNormalizedCounts_Scanorama_50pcs_filterRound.tsv"), header = TRUE, stringsAsFactors = F, sep = "\t"))
# logNormCountsFilter <- logNormCounts %>% dplyr::select(cellID, sampleNum, SPP1, GDF15, CDKN1A, FTH1, TOP2A, MKI67, CENPF, SOX21, GAPDH, UBC)
# logNormCountsFilter$cellID <- gsub(pattern = "*-1", x = logNormCountsFilter$cellID, replacement = "")
# saveRDS(logNormCountsFilter, file = paste0(homeDirectory, "logNormalizedCountsFilter.rds"))
logNormCountsFilter <- readRDS(file = paste0(homeDirectory, "logNormalizedCountsFilter.rds"))

labelsToAdd <- umapCoordinates %>% mutate(label = "none")
labelsToAdd$label <- ifelse(labelsToAdd$cellID %in% linCountToOverlaps$cellID, "barcoded", labelsToAdd$label)
labelsToAdd$label <- ifelse(labelsToAdd$cellID %in% primedAllUMAP$cellID, "primed", labelsToAdd$label)
oddsTable <- inner_join(labelsToAdd, logNormCountsFilter, by = "cellID") %>% filter(label != "none")

priming <- c("primed", "nonprimed")
geneLevel <- c("high", "low")
data <- matrix(nrow = 2, ncol = 2)
dimnames(data) <- list("gene level" = geneLevel, "priming" = priming)

oddsRatioList <- data.frame(estimate = NA, log_estimate = NA, se = NA, gene = NA, prop = NA)
geneList <- c("CENPF", "MKI67", "TOP2A", "SOX21", "SPP1", "GDF15", "CDKN1A", "FTH1", "GAPDH", "UBC")
propList <- c(0.01, 0.05, 0.10, 0.15, 0.25)

for(i in 1:length(geneList)) {
    for(j in 1:length(propList)) {
      oddsTable$geneLevel <- "no"
      oddsTable$geneLevel <- ifelse(oddsTable[[geneList[i]]] > min(slice_max(oddsTable, oddsTable[[geneList[i]]], prop = propList[j], with_ties = FALSE) %>% .[[geneList[i]]]), "yes", oddsTable$geneLevel)
      data[1, 1] <- nrow(oddsTable %>% dplyr::filter(label == "primed" & geneLevel == "yes")) + 0.5
      data[2, 1] <- nrow(oddsTable %>% dplyr::filter(label == "primed" & geneLevel == "no")) + 0.5
      data[1, 2] <- nrow(oddsTable %>% dplyr::filter(label == "barcoded" & geneLevel == "yes")) + 0.5
      data[2, 2] <- nrow(oddsTable %>% dplyr::filter(label == "barcoded" & geneLevel == "no")) + 0.5
      or_table <- data.frame(estimate = (data[1, 1] * data[2, 2])/(data[1, 2] * data[2, 1]),
                             log_estimate = log((data[1, 1] * data[2, 2])/(data[1, 2] * data[2, 1])),
                             se = sqrt(1/data[1, 1] + 1/data[1, 2] + 1/data[2, 1] + 1/data[2, 2]),
                             gene = geneList[i],
                             prop = propList[j])
      oddsRatioList <- bind_rows(oddsRatioList, or_table)
    }
  }

oddsRatioListFilter <- oddsRatioList[-1, ]
oddsRatioListFilter$prop <- factor(oddsRatioListFilter$prop, levels = propList %>% sort())
oddsRatioListFilter$gene <- factor(oddsRatioListFilter$gene, levels = geneList)

saveRDS(object = oddsRatioListFilter, file = paste0(homeDirectory, "R2_oddsRatioListFilter.rds"))

ggplot(oddsRatioListFilter, aes(x = prop, y = log_estimate)) +
  geom_errorbar(aes(ymin = log_estimate - se, ymax = log_estimate + se), width = 0.25) +
  geom_point(size = 5) + facet_wrap(~gene) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("log(odds ratio) of being primed with high expression of each gene") + theme(axis.title.x = element_blank())

ggplot(oddsRatioListFilter %>% filter(gene %in% geneList[1:4], prop == 0.1), aes(x = gene, y = log_estimate)) +
  geom_errorbar(aes(ymin = log_estimate - se, ymax = log_estimate + se), width = 0.25) +
  geom_point(size = 5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("log(odds ratio) of being primed with high expression of each gene") + theme(axis.title.x = element_blank())
ggsave(filename = paste0(plotDirectory, "R2_oddsRatioMarkers_positive.pdf"), unit = "in", height = 2, width = 5)

ggplot(oddsRatioListFilter %>% filter(gene %in% geneList[5:8], prop == 0.1), aes(x = gene, y = log_estimate)) +
  geom_errorbar(aes(ymin = log_estimate - se, ymax = log_estimate + se), width = 0.25) +
  geom_point(size = 5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("log(odds ratio) of being primed with high expression of each gene") + theme(axis.title.x = element_blank())
ggsave(filename = paste0(plotDirectory, "R2_oddsRatioMarkers_negative.pdf"), unit = "in", height = 2, width = 5)