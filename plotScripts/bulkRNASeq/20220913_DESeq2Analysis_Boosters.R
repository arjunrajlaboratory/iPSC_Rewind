library(tidyverse)
library(DESeq2)
library(reshape2)
library(biomaRt)
library(ggrepel)
library(stringr)
library(RColorBrewer)
library(pheatmap)
library(ggsignif)

theme_set(theme_classic())

dataDirectory <- '/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/extractedData/bulkRNASeq/'
plotDirectory <- '/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/plots/bulkRNASeq/'

countMatrix <- read.table(file = paste0(dataDirectory, '20220914_RNASeq_ProlifSort_Boosters_finalData.tsv'), sep = '\t', header = TRUE)
countMatrixFormat <- countMatrix %>% dplyr::select(sampleID, counts, gene_id) %>%
  dplyr::filter(sampleID %in% c("1_DMSO", "1_LSD1i", "1_DOT1Li", "2_noDrug", "2_LSD1i")) %>% 
  dcast(., gene_id ~ sampleID, value.var = "counts")
rownames(countMatrixFormat) <- countMatrixFormat$gene_id
countMatrixFormat <- countMatrixFormat %>% dplyr::select(-gene_id)

geneDic <- countMatrix %>% dplyr::select(gene_id, gene_name) %>% unique()
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
proteinCodingGenes <- getBM(attributes = c("ensembl_gene_id"), filters = 'biotype', values = c('protein_coding'), mart = ensembl)
colnames(proteinCodingGenes)[1] <- "gene_id"

condition <- c("control", "DOT1L", "LSD1", "LSD1", "control")
replicate <- c(rep("1", 3), rep("2", 2))
colData <- data_frame(cond = condition,
                      rep = replicate)
rownames(colData) <- colnames(countMatrixFormat)

dds <- DESeqDataSetFromMatrix(countData = countMatrixFormat,
                              colData = colData,
                              design = ~ cond + rep)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds <- DESeq(dds)

#### check differentially expressed genes for LSD1 and DOT1L versus control ####
#######################################################################################################################################################
res <- results(dds, contrast = c("cond", "LSD1", "control")) %>% as.data.frame() %>%
  replace(is.na(.), 0) %>%
  mutate(gene_id = rownames(.)) %>%
  filter(gene_id %in% proteinCodingGenes$gene_id) %>%
  filter(baseMean > 10) %>%
  inner_join(., geneDic, by = "gene_id")

res <- results(dds, contrast = c("cond", "DOT1L", "control")) %>% as.data.frame() %>%
  replace(is.na(.), 0) %>%
  mutate(gene_id = rownames(.)) %>%
  filter(gene_id %in% proteinCodingGenes$gene_id) %>%
  filter(baseMean > 10) %>%
  inner_join(., geneDic, by = "gene_id")

ggplot(res, aes(x = log2FoldChange, y = -log10(padj), label = gene_name)) +
  geom_point() +
  geom_text_repel(data = res %>% filter(abs(log2FoldChange) > 1),
                  aes(x = log2FoldChange, y = -log10(padj), label = gene_name), size = 2)

#### check differentially expressed genes for boosters for priming ####
#######################################################################################################################################################
markersPrevious <- readRDS("/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/extractedData/rewind10X/R1/primedMarkersAll.rds")

markersOverlap <- inner_join(markersPrevious, res, by = c("gene" = "gene_name"))
ggplot(markersOverlap, aes(x = avg_log2FC, y = log2FoldChange, label = gene)) +
  geom_point() +
  geom_text_repel()

markersPrevHigh <- markersPrevious %>% ungroup %>% slice_max(order_by = avg_log2FC, n = 50)
markersPrevLow <- markersPrevious %>% ungroup %>% slice_min(order_by = avg_log2FC, n = 50)

markersTreatHigh <- res %>% ungroup %>% slice_max(order_by = log2FoldChange, n = 50) %>% mutate(cond = "drug")
markersOverlapHigh <- res %>% dplyr::filter(gene_name %in% markersPrevHigh$gene) %>% mutate(cond = "high")
markersOverlapLow <- res %>% dplyr::filter(gene_name %in% markersPrevLow$gene) %>% mutate(cond = "low")
markersOverlapCont <- res %>% dplyr::filter(gene_name %in% c("UBC", "GAPDH", "PGH1", "ACTB")) %>% mutate(cond = "cont")

markersOverlapPlot <- bind_rows(markersTreatHigh, markersOverlapHigh, markersOverlapLow, markersOverlapCont)
genesPrime <- c("EPCAM", "KRT19", "ITGB4", "BMP4", "SPP1", "GDF15", "SQSTM1", "FTH1", "TOP2A", "MKI67", "CENPF", "ASPM", "GAPDH", "UBC")
markersOverlapPlot$cond <- factor(markersOverlapPlot$cond, levels = c("drug", "low", "cont", "high"))
markersOverlapPlot <- markersOverlapPlot %>% mutate(text = '')
markersOverlapPlot$text <- ifelse(markersOverlapPlot$gene_name %in% genesPrime, markersOverlapPlot$gene_name, markersOverlapPlot$text)

ggplot(markersOverlapPlot, aes(x = cond, y = log2FoldChange, label = text, group = cond)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(seed = 1234, width = 0.25)) +
  geom_text_repel(position = position_jitter(seed = 1234, width = 0.25), max.overlaps = 50) +
  geom_signif(comparisons = list(c("drug", "cont"), c("low", "cont"), c("cont", "high")), step_increase = 0.05) +
  geom_hline(yintercept = 0, linetype = "dashed")
ggsave(filename = paste0(plotDirectory, "previousMarkersBulkLSD1i.pdf"), height = 4, width = 5, useDingbats = FALSE)

#### check variance stabilizing distribution ####
#######################################################################################################################################################
vsd <- vst(dds, blind = FALSE)

plotPCA(vsd, intgroup = c("cond"))
plotPCA(vsd, intgroup = c("rep"))
assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$rep, design = model.matrix(~ vsd$cond))
plotPCA(vsd, intgroup = c("cond"))

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$cond, vsd$rep, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

#### check TPM values for select genes ####
#######################################################################################################################################################
genesPluri <- c("NANOG", "ALPL", "POU5F1", "SOX2", "PODXL", "LIN28A", "UTF1", "TERT", "ZFP42", "DNMT3B", "SALL4")
genesFibro <- c("LUM", "S100A4", "THY1", "PDGFRA", "COL1A1", "COL5A1", "LOXL1", "FBLN1", "FBLN2", "VTN")
genesEpi <- c("CDH1", "CLDN3", "KRT7", "OCLN", "EPCAM", "ANPEP", "MUC1", "CD24")
genesMes <- c("CDH2", "VIM", "FN1", "ZEB1", "SNAI2", "TWIST1", "TWIST2", "TGFB1")
genesPrime <- c("SPP1", "GDF15", "CDKN1A", "FTH1", "TOP2A", "MKI67", "CENPF", "SOX21", "GAPDH", "UBC")

genesList <- genesEpi

countMatrixPlot <- countMatrix %>% dplyr::select(sampleID, tpm, gene_name) %>%
  dplyr::filter(sampleID %in% c("1_DMSO", "1_LSD1i", "1_DOT1Li", "2_noDrug", "2_LSD1i")) %>%
  dplyr::filter(gene_name %in% c(genesList)) %>%
  mutate(rep = word(.$sampleID, 1, sep = "_")) %>%
  mutate(cond = word(.$sampleID, 2, sep = "_"))

countMatrixPlot$cond <- ifelse(countMatrixPlot$cond == "DMSO", "noDrug", countMatrixPlot$cond)

countMatrixPlot$cond <- factor(countMatrixPlot$cond, levels = c("noDrug", "LSD1i", "DOT1Li"))
ggplot(countMatrixPlot, aes(x = cond, y = tpm)) +
  geom_point() +
  stat_summary(fun.data = mean_se, geom = "errorbar") +
  facet_wrap(~ gene_name, scales = "free")

countMatrixPlotFormat <- data.frame(sampleID = character(),
                                    tpm = integer(),
                                    gene_name = character(),
                                    rep = character(),
                                    cond = character())
for(i in unique(countMatrixPlot$rep)) {
  for(j in unique(countMatrixPlot$gene_name)) {
    countMatrixPlotTemp <- countMatrixPlot %>% dplyr::filter(rep == i, gene_name == j)
    countMatrixPlotTemp <- countMatrixPlotTemp %>% rowwise() %>% mutate(countNorm = (tpm + 0.1)/(countMatrixPlotTemp %>% filter(cond == "noDrug") %>% .$tpm %>% mean() + 0.1)) %>% ungroup()
    countMatrixPlotFormat <- bind_rows(countMatrixPlotFormat, countMatrixPlotTemp)
  }
}

countMatrixPlotFormat$cond <- factor(countMatrixPlotFormat$cond, levels = c("noDrug", "LSD1i", "DOT1Li"))
countMatrixPlotFormat$gene_name <- factor(countMatrixPlotFormat$gene_name, levels = genesList)

ggplot(countMatrixPlotFormat, aes(x = cond, y = countNorm)) +
  stat_summary(fun = mean, geom = "bar") +
  facet_wrap(~ gene_name, scales = "fixed") +
  geom_hline(yintercept = 1, linetype = "dashed")

ggplot(countMatrixPlotFormat, aes(x = cond, y = countNorm)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2) +
  geom_signif(comparisons = list(c("noDrug", "LSD1i"), c("LSD1i", "DOT1Li"))) +
  geom_signif(comparisons = list(c("noDrug", "DOT1Li")), position = position_nudge(y = 1)) +
  geom_hline(yintercept = 1, linetype = "dashed")
