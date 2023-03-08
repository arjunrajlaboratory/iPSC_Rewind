rm(list=ls())
gc()

library(tidyverse)
library(reshape2)
library(ggrepel)
library(ggforce)
library(ggsignif)
library(ggridges)
library(ggrastr)
library(RColorBrewer)
library(viridis)
library(egg)
library(Seurat)
library(biomaRt)
library(spgs)

theme_set(theme_classic())

homeDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/extractedData/rewind10X/R3/"
plotDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/plots/rewind10X/R3/"

# scanorama_filter <- readRDS(paste0(homeDirectory, "scanorma_filter.rds"))
umapCoordinates = as_tibble(read.table(file = paste0(homeDirectory, "umapCoordinates_Scanorama_50pcs_filter.tsv"), header = TRUE, stringsAsFactors = F, sep = "\t"))
umapClusters = as_tibble(read.table(file = paste0(homeDirectory, "umapClusters_Scanorama_50pcs_filter.tsv"), header = TRUE, stringsAsFactors = F, sep = "\t"))
umapClusters <- umapClusters %>% dplyr::rename(cluster = scanorama_snn_res.0.3)

primedCellIDList <- readRDS(file = paste0(homeDirectory, "primedCellIDList.rds"))
cutoffList <- c(10, 25, 50, 100, 150, 200, 250, 500, 1000)
i = 6
primedAllUMAP <- filter(umapCoordinates, cellID %in% c(unlist(primedCellIDList[[i]]), unlist(primedCellIDList[[i+length(cutoffList)]]), unlist(primedCellIDList[[i+2*length(cutoffList)]]))) %>%
  dplyr::filter(sampleNum %in% c("S1", "S2", "S3"))

linCountToOverlaps = as_tibble(read.table(file = paste0(homeDirectory, "filtered10XCells.txt"), header = TRUE, stringsAsFactors = F, sep = "\t"))

#### check odds ratios of being primed for each gene compared with proliferation speed ####
#######################################################################################################################################################
# logNormCounts <- as_tibble(read.table(file = paste0(homeDirectory, "logNormalizedCounts_Scanorama_50pcs_filterRound.tsv"), header = TRUE, stringsAsFactors = F, sep = "\t"))
# logNormCountsFilter <- logNormCounts %>% dplyr::select(cellID, sampleNum, SPP1, GDF15, CDKN1A, FTH1, TOP2A, MKI67, CENPF, SOX21, GAPDH, UBC)
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
geneList <- c("S2", "CENPF", "MKI67", "TOP2A", "SOX21", "S3", "SPP1", "GDF15", "CDKN1A", "FTH1", "GAPDH", "UBC")
propList <- c(0.01, 0.05, 0.10, 0.15, 0.25)

for(i in 1:length(geneList)) {
  i = 7
  if(geneList[i] %in% c("S2", "S3")) {
    oddsTable$geneLevel <- "no"
    oddsTable$geneLevel <- ifelse(oddsTable$sampleNum.x == geneList[i], "yes", oddsTable$geneLevel)
    data[1, 1] <- nrow(oddsTable %>% filter(label == "primed" & geneLevel == "yes")) + 0.5
    data[2, 1] <- nrow(oddsTable %>% filter(label == "primed" & geneLevel == "no")) + 0.5
    data[1, 2] <- nrow(oddsTable %>% filter(label == "barcoded" & geneLevel == "yes")) + 0.5
    data[2, 2] <- nrow(oddsTable %>% filter(label == "barcoded" & geneLevel == "no")) + 0.5
    or_table <- data.frame(estimate = (data[1, 1] * data[2, 2])/(data[1, 2] * data[2, 1]),
                           log_estimate = log((data[1, 1] * data[2, 2])/(data[1, 2] * data[2, 1])),
                           se = sqrt(1/data[1, 1] + 1/data[1, 2] + 1/data[2, 1] + 1/data[2, 2]),
                           gene = geneList[i],
                           prop = NA)
    oddsRatioList <- bind_rows(oddsRatioList, or_table)
  } else {
    for(j in 1:length(propList)) {
      oddsTable$geneLevel <- "no"
      oddsTable$geneLevel <- ifelse(oddsTable[[geneList[i]]] > min(slice_max(oddsTable, oddsTable[[geneList[i]]], prop = propList[j], with_ties = FALSE) %>% .[[geneList[i]]]), "yes", oddsTable$geneLevel)
      data[1, 1] <- nrow(oddsTable %>% filter(label == "primed" & geneLevel == "yes")) + 0.5
      data[2, 1] <- nrow(oddsTable %>% filter(label == "primed" & geneLevel == "no")) + 0.5
      data[1, 2] <- nrow(oddsTable %>% filter(label == "barcoded" & geneLevel == "yes")) + 0.5
      data[2, 2] <- nrow(oddsTable %>% filter(label == "barcoded" & geneLevel == "no")) + 0.5
      or_table <- data.frame(estimate = (data[1, 1] * data[2, 2])/(data[1, 2] * data[2, 1]),
                             log_estimate = log((data[1, 1] * data[2, 2])/(data[1, 2] * data[2, 1])),
                             se = sqrt(1/data[1, 1] + 1/data[1, 2] + 1/data[2, 1] + 1/data[2, 2]),
                             gene = geneList[i],
                             prop = propList[j])
      oddsRatioList <- bind_rows(oddsRatioList, or_table)
    }
  }
}

oddsRatioListFilter <- oddsRatioList[-1, ]
oddsRatioListFilter$prop <- factor(oddsRatioListFilter$prop, levels = propList %>% sort())
oddsRatioListFilter$gene <- factor(oddsRatioListFilter$gene, levels = geneList)

saveRDS(object = oddsRatioListFilter, file = paste0(homeDirectory, "R3_oddsRatioListFilter.rds"))

oddsRatioListFilterP <- oddsRatioListFilter %>% filter(!(gene %in% c("S2", "S3"))) %>% filter(prop == 0.1) %>%
  rowwise() %>% mutate(pval1 = pnorm(abs(abs(log_estimate) - abs(0.47449380))/sqrt(se^2 + 0.1465549^2),
                                     mean = 0, sd = 1, lower.tail = FALSE)* 2) %>% # comparing with fast
  rowwise() %>% mutate(pval2 = pnorm(abs(abs(log_estimate) - abs(-0.35874367))/sqrt(se^2 + 0.1676660^2),
                                     mean = 0, sd = 1, lower.tail = FALSE)* 2) # comparing with slow

oddsRatioListFilterPlot <- oddsRatioListFilter %>% filter(prop %in% c(NA, 0.1))
oddsRatioListFilterPlot$pval <- c(NA, oddsRatioListFilterP$pval1[1:4], NA, oddsRatioListFilterP$pval2[5:10])

ggplot(oddsRatioListFilterPlot %>% filter(gene %in% geneList[1:5]), aes(x = gene, y = log_estimate)) +
  geom_errorbar(aes(ymin = log_estimate - se, ymax = log_estimate + se), width = 0.25) +
  geom_text(aes(label = round(pval, 2)), hjust = -1) +
  geom_point(size = 5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("log(odds ratio) of being primed with high expression of each gene") + theme(axis.title.x = element_blank())
ggsave(filename = paste0(plotDirectory, "R3_oddsRatioMarkers_positive.pdf"), unit = "in", height = 2, width = 5)

ggplot(oddsRatioListFilterPlot %>% filter(gene %in% geneList[6:10]), aes(x = gene, y = log_estimate)) +
  geom_errorbar(aes(ymin = log_estimate - se, ymax = log_estimate + se), width = 0.25) +
  geom_text(aes(label = round(pval, 2)), hjust = -1) +
  geom_point(size = 5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("log(odds ratio) of being primed with high expression of each gene") + theme(axis.title.x = element_blank())
ggsave(filename = paste0(plotDirectory, "R3_oddsRatioMarkers_negative.pdf"), unit = "in", height = 2, width = 5)

#### check odds ratio without dichotimizing gene expression ####
#######################################################################################################################################################
oddsTableNB <- oddsTable %>% mutate(outcome = 0) %>% mutate(prolif = 0)
oddsTableNB$outcome <- ifelse(oddsTableNB$label == "primed", 1, oddsTableNB$outcome)
oddsTableNB$prolif <- ifelse(oddsTable$sampleNum.x == "S2", 1, oddsTableNB$prolif)
oddsTableNB$prolif <- ifelse(oddsTable$sampleNum.x == "S1", 0.5, oddsTableNB$prolif)

ggplot(oddsTableNB, aes(x = SPP1, y = outcome)) +
  geom_jitter()

geneList <- c("CENPF", "MKI67", "TOP2A", "SOX21", "SPP1", "GDF15", "CDKN1A", "FTH1", "GAPDH", "UBC")
for(i in 1:length(geneList)) {
  model <- glm(formula = formula(paste("outcome ~ ", geneList[i], "* prolif")), family = "binomial", data = oddsTableNB)
  modelSum <- summary(model, test = "LRT")
  coeffTableTemp <- modelSum$coefficients %>% as_tibble(rownames = "variable") %>% mutate(gene = geneList[i])
  if(i == 1) {
    coeffTable <- coeffTableTemp
  } else {
    coeffTable <- bind_rows(coeffTable, coeffTableTemp)
  }
}

coeffTable <- coeffTable %>% dplyr::filter(variable != "(Intercept)") %>% mutate(coefficient = "interaction")
coeffTable$coefficient <- ifelse(coeffTable$variable %in% geneList, "gene", coeffTable$coefficient)
coeffTable$coefficient <- ifelse(coeffTable$variable == "prolif", "prolif", coeffTable$coefficient)
coeffTable$coefficient <- factor(coeffTable$coefficient, levels = c("gene", "prolif", "interaction"))
coeffTable$gene <- factor(coeffTable$gene, levels = geneList)

coeffTable <- coeffTable %>% dplyr::rename(., pval = 'Pr(>|z|)')
ggplot(coeffTable, aes(x = coefficient, y = -log(pval))) +
  geom_col() + facet_wrap(~gene, ncol = 4) +
  geom_hline(yintercept = -log(0.05), linetype = "dashed")

library(scales)
ggplot(coeffTable %>% dplyr::filter(gene %in% geneList[1:4]), aes(x = coefficient, y = pval)) +
  geom_col() + facet_wrap(~gene, ncol = 4) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  scale_y_continuous(trans = c("log", "reverse"), labels = label_number(accuracy = 0.01), limits = c(1, 0.0075), breaks = c(exp(0), exp(-1), exp(-2), exp(-3), exp(-4), exp(-5)))
ggsave(filename = paste0(plotDirectory, "logisticModelSum_positiveMarkers.pdf"), units = "in", height = 2.5, width = 7.5)

ggplot(coeffTable %>% dplyr::filter(gene %in% geneList[5:8]), aes(x = coefficient, y = pval)) +
  geom_col() + facet_wrap(~gene, ncol = 4) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  scale_y_continuous(trans = c("log", "reverse"), labels = label_number(accuracy = 0.01), limits = c(1, 0.0075), breaks = c(exp(0), exp(-1), exp(-2), exp(-3), exp(-4), exp(-5)))
ggsave(filename = paste0(plotDirectory, "logisticModelSum_negativeMarkers.pdf"), units = "in", height = 2.5, width = 7.5)

ggplot(coeffTable %>% dplyr::filter(gene %in% geneList[9:10]), aes(x = coefficient, y = pval)) +
  geom_col() + facet_wrap(~gene, ncol = 4) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  scale_y_continuous(trans = c("log", "reverse"), labels = label_number(accuracy = 0.01), limits = c(1, 0.0075), breaks = c(exp(0), exp(-1), exp(-2), exp(-3), exp(-4), exp(-5)))
ggsave(filename = paste0(plotDirectory, "logisticModelSum_controlMarkers.pdf"), units = "in", height = 2.5, width = 3.75)

geneList <- c("prolif", "CENPF", "MKI67", "TOP2A", "SOX21", "SPP1", "GDF15", "CDKN1A", "FTH1", "GAPDH", "UBC")
for(i in 1:length(geneList)) {
  model <- glm(formula = formula(paste("outcome ~ ", geneList[i])), family = "binomial", data = oddsTableNB)
  modelSum <- summary(model, test = "LRT")
  rSquareTableTemp <- tibble(gene = geneList[i], rSquare = with(modelSum, 1-deviance/null.deviance))
  if(i == 1) {
    rSquareTable <- rSquareTableTemp
  } else {
    rSquareTable <- bind_rows(rSquareTable, rSquareTableTemp)
  }
}

rSquareTable$gene <- factor(rSquareTable$gene, levels = geneList)
ggplot(rSquareTable %>% dplyr::filter(gene != "prolif"), aes(x = gene, y = rSquare * 100)) +
  geom_col() +
  geom_hline(yintercept = rSquareTable %>% dplyr::filter(gene == "prolif") %>% .$rSquare * 100, linetype = "dashed") +
  ylab("percent variation in priming explained by each gene") +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
ggsave(filename = paste0(plotDirectory, "logisticModelRSquare.pdf"), units = "in", height = 3, width = 3)

nullmodel <- glm(outcome ~ 1, data = oddsTableNB)
fullmodel <- glm(outcome ~ prolif + CENPF + MKI67 + TOP2A + SOX21 + SPP1 + GDF15 + CDKN1A + FTH1, data = oddsTableNB)
summary(fullmodel)

library(stats)
forward <- step(nullmodel, direction = "both", scope = formula(fullmodel))

#### check odds ratios of being primed for each gene given fast versus slow versus control ####
#######################################################################################################################################################
i = 7
primedAllUMAP <- filter(umapCoordinates, cellID %in% c(unlist(primedCellIDList[[i]]), unlist(primedCellIDList[[i+length(cutoffList)]]), unlist(primedCellIDList[[i+2*length(cutoffList)]]))) %>%
  dplyr::filter(sampleNum %in% c("S1", "S2", "S3"))
labelsToAdd <- umapCoordinates %>% mutate(label = "none")
labelsToAdd$label <- ifelse(labelsToAdd$cellID %in% linCountToOverlaps$cellID, "barcoded", labelsToAdd$label)
labelsToAdd$label <- ifelse(labelsToAdd$cellID %in% primedAllUMAP$cellID, "primed", labelsToAdd$label)
oddsTable <- inner_join(labelsToAdd, logNormCountsFilter, by = "cellID") %>% filter(label != "none")

oddsRatioList <- data.frame(estimate = NA, log_estimate = NA, se = NA, gene = NA, prop = NA)
# geneList <- c("TOP2A", "MKI67", "CENPF", "SOX21", "SPP1", "GDF15", "CDKN1A", "FTH1")
geneList <- c("SPP1", "GDF15", "CDKN1A", "FTH1")
propList <- c(0.01, 0.05, 0.10, 0.15, 0.25)
condList <- c("S1", "S2", "S3")

for(i in 1:length(geneList)) {
  for(j in 1:length(propList)) {
    for(k in 1:length(condList)) {
      oddsTableSub <- oddsTable %>% dplyr::filter(sampleNum.x == condList[k])
      oddsTableSub$geneLevel <- "no"
      oddsTableSub$geneLevel <- ifelse(oddsTableSub[[geneList[i]]] > min(slice_max(oddsTableSub, oddsTableSub[[geneList[i]]], prop = propList[j], with_ties = FALSE) %>% .[[geneList[i]]]), "yes", oddsTableSub$geneLevel)
      data[1, 1] <- nrow(oddsTableSub %>% filter(label == "primed" & geneLevel == "yes")) + 0.5
      data[2, 1] <- nrow(oddsTableSub %>% filter(label == "primed" & geneLevel == "no")) + 0.5
      data[1, 2] <- nrow(oddsTableSub %>% filter(label == "barcoded" & geneLevel == "yes")) + 0.5
      data[2, 2] <- nrow(oddsTableSub %>% filter(label == "barcoded" & geneLevel == "no")) + 0.5
      or_table <- data.frame(estimate = (data[1, 1] * data[2, 2])/(data[1, 2] * data[2, 1]),
                             log_estimate = log((data[1, 1] * data[2, 2])/(data[1, 2] * data[2, 1])),
                             se = sqrt(1/data[1, 1] + 1/data[1, 2] + 1/data[2, 1] + 1/data[2, 2]),
                             cond = condList[k],
                             gene = geneList[i],
                             prop = propList[j])
      oddsRatioList <- bind_rows(oddsRatioList, or_table)
    }
  }
}

oddsRatioListFilter <- oddsRatioList[-1, ]
oddsRatioListFilter$prop <- factor(oddsRatioListFilter$prop, levels = propList %>% sort())
oddsRatioListFilter$gene <- factor(oddsRatioListFilter$gene, levels = geneList)
oddsRatioListFilter$cond <- factor(oddsRatioListFilter$cond, levels = c("S3", "S1", "S2"), labels = c("slow", "control", "fast"))

ggplot(oddsRatioListFilter %>% filter(prop == 0.1), aes(x = cond, y = log_estimate)) +
  geom_errorbar(aes(ymin = log_estimate - se*1.645, ymax = log_estimate + se*1.645)) +
  geom_point(size = 5) +
  geom_hline(yintercept = 0, linetype = "dashed") + facet_wrap(~gene)

ggplot(oddsRatioListFilter %>% filter(prop == 0.15) %>% filter(cond == "slow"), aes(x = gene, y = log_estimate)) +
  geom_errorbar(aes(ymin = log_estimate - se*1.96, ymax = log_estimate + se*1.96), width = 0.25) +
  geom_point(size = 5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.title.x = element_blank()) + ylab("log(odds ratio) of being primed with high expression of each gene if slow")
ggsave(filename = paste0(plotDirectory, "oddsRatioMarkersVsCycling.pdf"), unit = "in", height = 2, width = 4)

#### check expression of certain primed state markers in different primed populations ####
#######################################################################################################################################################
i = 6
linCountToOverlaps = as_tibble(read.table(file = paste0(homeDirectory, "filtered10XCells.txt"), header = TRUE, stringsAsFactors = F, sep = "\t"))
primedAllUMAP <- filter(umapCoordinates, cellID %in% c(unlist(primedCellIDList[[i]]), unlist(primedCellIDList[[i+length(cutoffList)]]), unlist(primedCellIDList[[i+2*length(cutoffList)]]))) %>%
  dplyr::filter(sampleNum %in% c("S1", "S2", "S3"))
labelsToAdd <- umapCoordinates %>% mutate(label = "none")
labelsToAdd$label <- ifelse(labelsToAdd$cellID %in% linCountToOverlaps$cellID, "barcoded", labelsToAdd$label)
labelsToAdd$label <- ifelse(labelsToAdd$cellID %in% primedAllUMAP$cellID, "primed", labelsToAdd$label)
logNormCountsFilter <- readRDS(file = paste0(homeDirectory, "logNormalizedCountsFilter.rds"))
oddsTable <- inner_join(labelsToAdd, logNormCountsFilter, by = "cellID") %>% filter(label != "none")
oddsTable$sampleNum.x <- factor(oddsTable$sampleNum.x, levels = c("S3", "S1", "S2"), labels = c("slow", "control", "fast"))

ggplot(oddsTable, aes(x = sampleNum.x, y = SPP1)) +
  rasterise(geom_jitter(size = 0.25), dpi = 100) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_signif(comparisons = list(c("slow", "control"), c("control", "fast"))) +
  geom_signif(comparisons = list(c("slow", "fast")), position = position_nudge(x = 0, y = 0.5)) +
  theme(axis.title.x = element_blank())
ggsave(filename = paste0(plotDirectory, "R3_SPP1ExpressionByCategory.pdf"), unit = "in", height = 3, width = 3)

ggplot(oddsTable %>% dplyr::filter(label == "primed"), aes(x = sampleNum.x, y = SPP1)) +
  rasterise(geom_jitter(size = 0.25), dpi = 100) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_signif(comparisons = list(c("slow", "control"), c("fast", "control"))) +
  geom_signif(comparisons = list(c("slow", "fast")), position = position_nudge(x = 0, y = 0.5)) +
  theme(axis.title.x = element_blank())
ggsave(filename = paste0(plotDirectory, "R3_SPP1ExpressionByCategoryPrimed.pdf"), unit = "in", height = 3, width = 3)

ggplot(oddsTable, aes(x = label, y = exp(1)^SPP1)) +
  geom_boxplot(outlier.shape = NA) +
  rasterise(geom_jitter(size = 0.25, width = 0.25), dpi = 300) +
  facet_wrap(~sampleNum.x) +
  scale_y_continuous(trans = "log", breaks = c(0, 1, 10, 100, 1000, 10000)) +
  geom_signif(comparisons = list(c("barcoded", "primed")))
ggsave(filename = paste0(plotDirectory, "R3_SPP1ExpressionByAllCategories.pdf"), unit = "in", height = 4, width = 7)
