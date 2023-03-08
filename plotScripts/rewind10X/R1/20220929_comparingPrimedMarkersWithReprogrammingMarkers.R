library(tidyverse)
library(biomaRt)
library(ggrastr)
library(ggsignif)
library(ggrepel)
library(egg)
library(reshape2)

theme_set(theme_classic())

plotDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/plots/rewind10X/R1/"

SCResults <- readRDS('/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/extractedData/rewind10X/R1/primedMarkersAll.rds') %>%
  dplyr::select(gene, avg_log2FC, p_val_adj)

#### pull out foldchange values for protein-coding genes during reprogramming
repTable <- read.csv(file = "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/extractedData/rewind10X/previousDatasets/bulkRNASeqDuringReprogramming.csv", header = TRUE)

repTablePlot <- repTable %>% dplyr::select(gene_short_name, 4, 5, 6, 7, 8, 9, 10, 12, 13) %>%
  dplyr::rename(day_0 = 2, day_2 = 3, day_5 = 4, day_8 = 5, day_10 = 6, day_14 = 7, day_20 = 8, day_24 = 9, ipsc = 10)
repTablePlotMelt <- melt(repTablePlot, id.vars = c("gene_short_name"))
repTablePlotMeltNorm <- repTablePlotMelt

mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = 'https://www.ensembl.org')
genes <- biomaRt::getBM(attributes = c("external_gene_name", "transcript_biotype"), filters = c("transcript_biotype"), values = list("protein_coding"), mart = mart)

repTablePC <- repTable %>% filter(gene_short_name %in% genes$external_gene_name)

repTablePC <- repTablePC %>% rowwise() %>% 
  mutate(earlyFC = log2((X2dd_DOX_plus + 0.01)/(hiF.T_P14 + 0.01))) %>%
  mutate(lateFC = log2((hIPSC.T_P10 + 0.01)/(hiF.T_P14 + 0.01))) %>%
  dplyr::select(gene_short_name, earlyFC, lateFC)

#### compare foldchange values for primed markers versus reprogramming
overlapTable <- full_join(SCResults, repTablePC, by = c("gene" = "gene_short_name"))
overlapTable[is.na(overlapTable)] <- 0

theme_set(theme_classic())

genesPluri <- c("NANOG", "ALPL", "POU5F1", "SOX2", "PODXL", "LIN28A", "UTF1", "TERT", "ZFP42", "DNMT3B")
genesFibro <- c("LUM", "S100A4", "THY1", "PDGFRA", "COL1A1", "COL5A1", "LOXL1", "FBLN1", "FBLN2", "VTN")
genesEpi <- c("CDH1", "CLDN3", "KRT3", "OCLN", "EPCAM", "ANPEP", "MUC1", "CD24")
genesMes <- c("CDH2", "VIM", "FN1", "ZEB1", "SNAI2", "TWIST1", "TWIST2", "TGFB1")

geneList1 <- c(genesPluri)
geneList2 <- c(genesFibro)

# geneList1 <- c(genesEpi)
# geneList2 <- c(genesMes)

genesPos <- c("MKI67", "TOP2A", "CENPF", "SMC4", "CDC20", "ASPM")
genesNeg <- c("CDKN1A", "CDKN2A", "SPP1", "SERPINE2", "THBS1", "TAGLN", "MYL9", "ACTA2", "TPM2", "POSTN", "GDF15", "IGFBP7", "GAS5")

plot1 <- ggplot() + 
  rasterize(geom_point(data = overlapTable, aes(x = avg_log2FC, y = lateFC), alpha = 0.25), dpi = 300) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(data = filter(overlapTable, gene %in% genesPos), aes(x = avg_log2FC, y = lateFC), color = "#EE352E") +
  geom_point(data = filter(overlapTable, gene %in% genesNeg), aes(x = avg_log2FC, y = lateFC), color = "blue") +
  geom_text_repel(data = filter(overlapTable, gene %in% genesPos), aes(x = avg_log2FC, y = lateFC, label = gene), max.overlaps = 30, xlim = c(0.25, NA), color = "#EE352E") +
  geom_text_repel(data = filter(overlapTable, gene %in% genesNeg), aes(x = avg_log2FC, y = lateFC, label = gene), max.overlaps = 15, xlim = c(NA, -0.25), color = "blue") +
  ylab("log2(foldchange) in iPSCs / fibroblasts") + xlab("log2(foldchange) in primed cells / nonprimed cells") +
  xlim(-1.4, 1.4)
ggsave(plot = plot1, filename = paste0(plotDirectory, "markersInPrimedCellsVersusiPSCs_expMarkers.pdf"), unit = "in", height = 5, width = 5, useDingbats = F)

give.n <- function(x){
  return(c(y = mean(x), label = round(mean(x), 2)))
}

plot1sub <- ggplot() +
  geom_violin(data = filter(overlapTable, avg_log2FC > 0.25), aes(x = 2, y = lateFC), fill = "#EE352E", draw_quantiles = c(0.25, 0.75), linetype = "dashed") +
  geom_violin(data = filter(overlapTable, avg_log2FC > 0.25), aes(x = 2, y = lateFC), fill = "transparent", draw_quantiles = c(0.5)) +
  stat_summary(data = filter(overlapTable, avg_log2FC > 0.25), fun.data = give.n, geom = "text", hjust = 0,
               aes(x = 2.05, y = lateFC)) +
  stat_summary(data = filter(overlapTable, avg_log2FC > 0.25), fun = mean, geom = "point", aes(x = 2, y = lateFC), size = 5) +
  geom_violin(data = filter(overlapTable, avg_log2FC < -0.25), aes(x = 1, y = lateFC), fill = "#A7A9AC", draw_quantiles = c(0.25, 0.75), linetype = "dashed") +
  geom_violin(data = filter(overlapTable, avg_log2FC < -0.25), aes(x = 1, y = lateFC), fill = "transparent", draw_quantiles = c(0.5)) +
  stat_summary(data = filter(overlapTable, avg_log2FC < -0.25), fun.data = give.n, geom = "text", hjust = 0,
               aes(x = 1.05, y = lateFC)) +
  stat_summary(data = filter(overlapTable, avg_log2FC < -0.25), fun = mean, geom = "point", aes(x = 1, y = lateFC), size = 5) + ylim(-12, 12)
ggsave(plot = plot1sub, filename = paste0(plotDirectory, "markersInPrimedCellsVersusiPSCs_expMarkersviolin.pdf"), unit = "in", height = 3, width = 6, useDingbats = F)

plot2 <- ggplot() + 
  rasterize(geom_point(data = overlapTable, aes(x = avg_log2FC, y = lateFC), alpha = 0.05), dpi = 100) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(data = dplyr::filter(overlapTable, gene %in% geneList1), aes(x = avg_log2FC, y = lateFC), color = "#0039a6") +
  geom_text_repel(data = dplyr::filter(overlapTable, gene %in% geneList1), aes(x = avg_log2FC, y = lateFC, label = gene),
                  force = 0.5,
                  nudge_x = 0.45-dplyr::filter(overlapTable, gene %in% geneList1)$avg_log2FC,
                  direction = "y",
                  hjust = 0,
                  segment.size = 0.2,
                  max.overlaps = Inf,
                  color = "#0039a6",
                  size = 2.5) +
  geom_point(data = dplyr::filter(overlapTable, gene %in% geneList2), aes(x = avg_log2FC, y = lateFC), color = "#fccc0a") +
  geom_text_repel(data = dplyr::filter(overlapTable, gene %in% geneList2), aes(x = avg_log2FC, y = lateFC, label = gene),
                  force = 0.5,
                  nudge_x = -0.45-dplyr::filter(overlapTable, gene %in% geneList2)$avg_log2FC,
                  direction = "y",
                  hjust = 1,
                  segment.size = 0.2,
                  max.overlaps = Inf,
                  color = "#fccc0a",
                  size = 2.5) +
  ylab("log2(foldchange) in iPSCs / fibroblasts") + xlab("log2(foldchange) in primed cells / nonprimed cells") +
  xlim(-0.5, 0.5)
ggsave(plot = plot2, filename = paste0(plotDirectory, "markersInPrimedCellsVersusiPSCs_prevMarkers.pdf"), unit = "in", height = 5, width = 2, useDingbats = F)

plot2sub <- ggplot() +
  geom_violin(data = filter(overlapTable, gene %in% geneList1), aes(x = 2, y = avg_log2FC), fill = "#0039a6", draw_quantiles = c(0.25, 0.75), linetype = "dashed", color = "black") +
  geom_violin(data = filter(overlapTable, gene %in% geneList1), aes(x = 2, y = avg_log2FC), fill = "transparent", draw_quantiles = c(0.5), color = "black") +
  stat_summary(data = filter(overlapTable, gene %in% geneList1), fun.data = give.n, geom = "text", hjust = 0,
               aes(x = 2.05, y = avg_log2FC)) +
  stat_summary(data = filter(overlapTable, gene %in% geneList1), fun = mean, geom = "point", aes(x = 2, y = avg_log2FC), size = 5) +
  geom_violin(data = filter(overlapTable, gene %in% geneList2), aes(x = 1, y = avg_log2FC), fill = "#fccc0a", draw_quantiles = c(0.25, 0.75), linetype = "dashed", color = "black") +
  geom_violin(data = filter(overlapTable, gene %in% geneList2), aes(x = 1, y = avg_log2FC), fill = "transparent", draw_quantiles = c(0.5), color = "black") +
  stat_summary(data = filter(overlapTable, gene %in% geneList2), fun.data = give.n, geom = "text", hjust = 0,
               aes(x = 1.05, y = avg_log2FC)) +
  stat_summary(data = filter(overlapTable, gene %in% geneList2), fun = mean, geom = "point", aes(x = 1, y = avg_log2FC), size = 5) + ylim(-0.25, 0.25)
ggsave(plot = plot2sub, filename = paste0(plotDirectory, "markersInPrimedCellsVersusiPSCs_prevMarkersviolin.pdf"), unit = "in", height = 3, width = 6, useDingbats = F)

#### make boxplots of different categories of genes based on markers ###############################################################################################################
genesPluri <- c("NANOG", "ALPL", "POU5F1", "SOX2", "PODXL", "LIN28A", "UTF1", "TERT", "ZFP42", "DNMT3B")
genesFibro <- c("LUM", "S100A4", "THY1", "PDGFRA", "COL1A1", "COL1A2", "LOXL1", "FBLN1", "FBLN2", "VTN")
genesEpi <- c("CDH1", "CLDN3", "KRT3", "OCLN", "EPCAM", "ANPEP", "MUC1", "CD24")
genesMes <- c("CDH2", "VIM", "ENTPD1", "SNAI1", "SNAI2", "TWIST1", "TWIST2", "TGFB1")
genesActFibro <- c("FAP", "TAGLN", "TNC", "SPARC", "LGALS1", "SERPINE2", "FBN1", "THBS1", "SHOX2", "TBX2")
genesMyoFibro <- c("ACTA2", "TPM2", "POSTN", "IGB1", "MMP1", "MMP14", "FN1", "AOC3", "NKX2-3", "LRRC17")

markersList <- c(genesPluri, genesFibro, genesEpi, genesMes, genesActFibro, genesMyoFibro)
markersList <- data.frame(gene = markersList)

markersPlot <- full_join(markersList, SCResults, by = "gene") %>% replace(is.na(.), 0)
markersPlot$label <- "none"
markersPlot$label <- ifelse(markersPlot$gene %in% genesPluri, "pluripotency", markersPlot$label)
markersPlot$label <- ifelse(markersPlot$gene %in% genesFibro, "fibroblast", markersPlot$label)
markersPlot$label <- ifelse(markersPlot$gene %in% genesEpi, "epithelial", markersPlot$label)
markersPlot$label <- ifelse(markersPlot$gene %in% genesMes, "mesenchymal", markersPlot$label)
markersPlot$label <- ifelse(markersPlot$gene %in% genesActFibro, "activated fibroblast", markersPlot$label)
markersPlot$label <- ifelse(markersPlot$gene %in% genesMyoFibro, "myofibroblast", markersPlot$label)

meanLog2FCTable <- c()
medianLog2FCTable <- c()
for(i in 1:1000) {
  seed <- i
  sampleResultsTable <- sample_n(SCResults, size = 10, seed = 1)
  meanLog2FCTable[i] <- mean(sampleResultsTable$avg_log2FC)
  medianLog2FCTable[i] <- median(sampleResultsTable$avg_log2FC)
}
meanLog2FCTable <- data.frame(value = meanLog2FCTable)
medianLog2FCTable <- data.frame(value = medianLog2FCTable)

ggplot(meanLog2FCTable, aes(y = value)) +
  geom_histogram(bins = 100) + ylim(-1, 1)

ggplot(medianLog2FCTable, aes(y = value)) +
  geom_histogram(bins = 100) + ylim(-1, 1)

markersPlot$label <- factor(markersPlot$label, levels = c("none", "pluripotency", "epithelial", "mesenchymal", "fibroblast", "activated fibroblast", "myofibroblast"))
markersPlotSummary <- markersPlot %>% group_by(label) %>% summarise(mean = mean(avg_log2FC), median = median(avg_log2FC)) %>%
  rowwise() %>% mutate(p_mean = sum(meanLog2FCTable$value > abs(mean)) * 2 / length(meanLog2FCTable$value)) %>%
  mutate(p_median = sum(medianLog2FCTable$value > abs(mean)) * 2 / length(medianLog2FCTable$value))

ggplot(markersPlot %>% dplyr::filter(label %in% c("pluripotency", "fibroblast", "activated fibroblast", "myofibroblast")), aes(x = "", y = avg_log2FC, label = gene)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  stat_summary(fun = median, geom = "bar", width = 0.25) +
  geom_point(position = position_jitter(width = 0.25, seed = 1234)) +
  geom_text_repel(max.overlaps = 15, position = position_jitter(width = 0.25, seed = 1234), size = 3) +
  facet_wrap(~label, ncol = 6) + theme_bw() + ylab("average log2(foldchange) in primed / nonprimed cells") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.title.x = element_blank())
ggsave(filename = paste0(plotDirectory, "markersCategoryBoxplots.pdf"), unit = "in", height = 3.5, width = 5, useDingbats = F)

ggplot(markersPlot %>% dplyr::filter(label %in% c("pluripotency", "fibroblast")), aes(x = "", y = avg_log2FC, label = gene)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  stat_summary(fun = median, geom = "bar", width = 0.25) +
  geom_point(position = position_jitter(width = 0.25, seed = 1234)) +
  geom_text_repel(max.overlaps = 15, position = position_jitter(width = 0.25, seed = 1234), size = 3) +
  facet_wrap(~label, ncol = 6) + theme_bw() + ylab("average log2(foldchange) in primed / nonprimed cells") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.title.x = element_blank()) +
  ylim(-1, 1)
ggsave(filename = paste0(plotDirectory, "markersCategoryBoxplots_PluriFibro.pdf"), unit = "in", height = 4, width = 4, useDingbats = F)

#### check expression of primed factors during different points of reprogramming ###################################################################################################
repTablePC <- repTablePlotMeltNorm %>% filter(gene_short_name %in% genes$external_gene_name)
repTablePCSample <- repTablePC %>% dplyr::filter(gene_short_name %in% SCResults$gene)
repTablePCScaled <- repTablePCSample %>% group_by(gene_short_name) %>% mutate(scaled = (value - mean(value)) / sd(value)) %>% group_by(gene_short_name, variable) %>% summarise(mean_scaled = mean(scaled), nrep = n()) %>% ungroup()

repTableMatrix <- repTablePlot %>% dplyr::select(-1) %>% as.matrix
rownames(repTableMatrix) <- repTablePlot$gene_short_name
hclust_matrix <- repTableMatrix %>% t() %>% scale() %>% t()
hclust_matrix <- hclust_matrix[is.finite(hclust_matrix[ , "day_0"]), ]

gene_dist <- dist(hclust_matrix)
gene_hclust <- hclust(gene_dist, method = "complete")
gene_cluster <- cutree(gene_hclust, k = 6) %>% enframe() %>% rename(gene = name, cluster = value)

repTableClusterPlot <- inner_join(repTablePCScaled, gene_cluster, by = c("gene_short_name" = "gene"))
repTableClusterPlot$cluster <- ifelse(repTableClusterPlot$cluster %in% c("2", "3"), "2", repTableClusterPlot$cluster)
ggplot(repTableClusterPlot, aes(x = variable, y = mean_scaled, group = gene_short_name)) +
  rasterise(geom_line(alpha = 0.05), dpi = 150) + facet_wrap(~cluster, ncol = 5) +
  geom_line(stat = "summary", fun = "median", size = 1.5, aes(group = 1), color = "red")
ggsave(filename = paste0(plotDirectory, "reprogExpressionClustersAll.pdf"), height = 2, width = 8, useDingbats = FALSE)
repTableAllProp <- repTableClusterPlot %>% group_by(cluster) %>% summarise(count = n()) %>% ungroup() %>% mutate(prop = count/sum(count))

markersSelect <- c("CENPF", "MKI67", "TOP2A", "SOX21", "SPP1", "GDF15", "CDKN1A", "FTH1")

repTablePCUp <- repTableClusterPlot %>% dplyr::filter(gene_short_name %in% slice_max(SCResults, avg_log2FC, n = 100)$gene)
ggplot(repTablePCUp %>% dplyr::filter(cluster == "1"), aes(x = variable, y = mean_scaled, group = gene_short_name)) +
  geom_line() + facet_wrap(~cluster, ncol = 5) +
  geom_line(stat = "summary", fun = "median", size = 1.5, aes(group = 1), color = "red")
ggsave(filename = paste0(plotDirectory, "reprogExpressionClustersExample_Up.pdf"), height = 3, width = 3, useDingbats = FALSE)

repTablePCUpProp <- repTablePCUp %>% group_by(cluster) %>% summarise(count = n()) %>% ungroup() %>% mutate(prop = count/sum(count))
repTablePCUpProp <- repTablePCUpProp %>% inner_join(repTableAllProp, by = "cluster")
repTablePCUpProp <- repTablePCUpProp %>% rowwise() %>%
  mutate(propPerClust = count.x/count.y) %>% ungroup() %>%
  mutate(propPerClustNorm = propPerClust/sum(propPerClust))
ggplot(repTablePCUpProp, aes(x = cluster, y = propPerClustNorm)) +
  geom_col() + geom_hline(yintercept = 1/nrow(repTablePCUpProp), linetype = "dashed") + ylim(0, 1)
ggsave(filename = paste0(plotDirectory, "reprogExpressionClustersProp_Up.pdf"), height = 3, width = 3, useDingbats = FALSE)
repTablePCUpRank <- repTablePCUp %>% dplyr::select(gene_short_name, cluster) %>% unique() %>% inner_join(., SCResults, by = c("gene_short_name" = "gene"))
ggplot(repTablePCUpRank, aes(x = cluster, y = avg_log2FC, label = gene_short_name)) +
  #geom_boxplot(outlier.shape = NA) +
  geom_point() +
  geom_text_repel() + ylim(0.25, 1.5) +
  geom_point(data = repTablePCUpRank %>% dplyr::filter(gene_short_name %in% markersSelect), color = "blue") +
  geom_text_repel(data = repTablePCUpRank %>% dplyr::filter(gene_short_name %in% markersSelect), color = "blue")
ggsave(filename = paste0(plotDirectory, "reprogExpressionClustersFC_Up.pdf"), height = 4, width = 4, useDingbats = FALSE)

repTablePCDown <- repTableClusterPlot %>% dplyr::filter(gene_short_name %in% slice_min(SCResults, avg_log2FC, n = 100)$gene)
ggplot(repTablePCDown %>% dplyr::filter(cluster == "2"), aes(x = variable, y = mean_scaled, group = gene_short_name)) +
  geom_line() + facet_wrap(~cluster, ncol = 5) +
  geom_line(stat = "summary", fun = "median", size = 1.5, aes(group = 1), color = "red")
ggsave(filename = paste0(plotDirectory, "reprogExpressionClustersExample_Down.pdf"), height = 3, width = 3, useDingbats = FALSE)

repTablePCDownProp <- repTablePCDown %>% group_by(cluster) %>% summarise(count = n()) %>% ungroup() %>% mutate(prop = count/sum(count))
repTablePCDownProp <- repTablePCDownProp %>% inner_join(repTableAllProp, by = "cluster")
repTablePCDownProp <- repTablePCDownProp %>% rowwise() %>%
  mutate(propPerClust = count.x/count.y) %>% ungroup() %>%
  mutate(propPerClustNorm = propPerClust/sum(propPerClust))
ggplot(repTablePCDownProp, aes(x = cluster, y = propPerClustNorm)) +
  geom_col() + geom_hline(yintercept = 1/nrow(repTablePCUpProp), linetype = "dashed") + ylim(0, 1)
ggsave(filename = paste0(plotDirectory, "reprogExpressionClustersProp_Down.pdf"), height = 3, width = 3, useDingbats = FALSE)
repTablePCDownRank <- repTablePCDown %>% dplyr::select(gene_short_name, cluster) %>% unique() %>% inner_join(., SCResults, by = c("gene_short_name" = "gene"))
ggplot(repTablePCDownRank, aes(x = cluster, y = avg_log2FC, label = gene_short_name)) +
  #geom_boxplot(outlier.shape = NA) +
  geom_point() +
  geom_text_repel() + ylim(-1.5, -0.25) +
  geom_point(data = repTablePCDownRank %>% dplyr::filter(gene_short_name %in% markersSelect), color = "blue") +
  geom_text_repel(data = repTablePCDownRank %>% dplyr::filter(gene_short_name %in% markersSelect), color = "blue")
ggsave(filename = paste0(plotDirectory, "reprogExpressionClustersFC_Down.pdf"), height = 4, width = 4, useDingbats = FALSE)
