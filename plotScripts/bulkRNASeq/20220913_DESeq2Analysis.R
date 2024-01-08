library(tidyverse)
library(DESeq2)
library(reshape2)
library(biomaRt)
library(ggrepel)
library(stringr)
library(RColorBrewer)
library(pheatmap)

theme_set(theme_classic())

dataDirectory <- '/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Original Manuscript/extractedData/bulkRNASeq/'
plotDirectory <- '/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Original Manuscript/plots/bulkRNASeq/'

countMatrix <- read.table(file = paste0(dataDirectory, '20220914_RNASeq_ProlifSort_Boosters_finalData.tsv'), sep = '\t', header = TRUE)
countMatrixFormat <- countMatrix %>% dplyr::select(sampleID, counts, gene_id) %>%
  dplyr::filter(!(sampleID %in% c("1_DMSO", "1_LSD1i", "1_DOT1Li", "2_noDrug", "2_LSD1i"))) %>% 
  dcast(., gene_id ~ sampleID, value.var = "counts")
rownames(countMatrixFormat) <- countMatrixFormat$gene_id
countMatrixFormat <- countMatrixFormat %>% dplyr::select(-gene_id)

geneDic <- countMatrix %>% dplyr::select(gene_id, gene_name) %>% unique()
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
proteinCodingGenes <- getBM(attributes = c("ensembl_gene_id"), filters = 'biotype', values = c('protein_coding'), mart = ensembl)
colnames(proteinCodingGenes)[1] <- "gene_id"

condition <- c("control", "high", "low", "midhigh", "midlow",
               "control", "high", "low", "midhigh", "midlow",
               "control", "high", "low", "midhigh", "midlow",
               "control", "high", "low", "midhigh", "midlow")
replicate <- c(rep("1", 5), rep("2", 5), rep("3", 5), rep("4", 5))
colData <- data_frame(cond = condition,
                      rep = replicate)
rownames(colData) <- colnames(countMatrixFormat)

dds <- DESeqDataSetFromMatrix(countData = countMatrixFormat,
                              colData = colData,
                              design = ~ cond + rep)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds <- DESeq(dds)

res <- results(dds, contrast = c("cond", "low", "high")) %>% as.data.frame() %>%
  replace(is.na(.), 0) %>%
  mutate(gene_id = rownames(.)) %>%
  filter(gene_id %in% proteinCodingGenes$gene_id) %>%
  filter(baseMean > 10) %>%
  inner_join(., geneDic, by = "gene_id")

ggplot(res, aes(x = log2FoldChange, y = -log10(padj), label = gene_name)) +
  geom_point() +
  geom_text_repel(data = res %>% filter(abs(log2FoldChange) > 1),
                  aes(x = log2FoldChange, y = -log10(padj), label = gene_name), size = 2)
saveRDS(file = paste0(homeDirectory, "DEProlifSpeed.rds"), object = res)

#### check variance stabilizing distribution ####
vsd <- vst(dds, blind = FALSE)

plotPCA(vsd, intgroup = c("cond"))
ggsave(filename = paste0(plotDirectory, "prolifPCAwoBatchCorrection.pdf"), units = "in", height = 2, width = 5, useDingbats = FALSE)
assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$rep, design = model.matrix(~ vsd$cond))
plotPCA(vsd, intgroup = c("cond"))
ggsave(filename = paste0(plotDirectory, "prolifPCAwBatchCorrection.pdf"), units = "in", height = 2, width = 5, useDingbats = FALSE)

rv <- rowVars(assay(vsd))
select <- order(rv, decreasing=TRUE)[seq_len(min(1000, length(rv)))]
pca <- prcomp(t(assay(vsd)[select,]))
loadings <- as.data.frame(pca$rotation)
loadings <- loadings %>% mutate(gene = rownames(.))
loadings <- left_join(loadings, geneDic, by = c("gene" = "gene_id"))

primedMarkers <- readRDS(file = "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/extractedData/rewind10X/R1/primedMarkersAll.rds")

library(ggrastr)
loadingsMarkers <- full_join(loadings, primedMarkers, by = c("gene_name" = "gene"))
loadingsMarkers$avg_log2FC[is.na(loadingsMarkers$avg_log2FC)] <- 0
loadingsMarkers$PC1[is.na(loadingsMarkers$PC1)] <- 0
ggplot(loadingsMarkers, aes(x = avg_log2FC, y = PC1)) +
  rasterize(geom_point(), dpi = 300) +
  geom_point(data = loadingsMarkers %>% dplyr::filter(gene_name %in% c("SPP1", "GDF15", "CDKN1A", "FTH1",
                                                                       "TOP2A", "MKI67", "CENPF", "SOX21")),
             aes(x = avg_log2FC, y = PC1), color = "blue") +
  geom_text_repel(data = loadingsMarkers %>% dplyr::filter(gene_name %in% c("SPP1", "GDF15", "CDKN1A", "FTH1",
                                                                     "TOP2A", "MKI67", "CENPF", "SOX21")),
                  aes(x = avg_log2FC, y = PC1, label = gene_name), color = "blue") +
  geom_vline(xintercept = 0, linetype = "dashed") + geom_hline(yintercept = 0, linetype = "dashed") +
  xlab("log2(foldchange) in\nprimed cells / nonprimed cells") + ylab("PCA loadings per gene") +
  coord_flip()
ggsave(filename = paste0(plotDirectory, "prolifPCALoadings.pdf"), units = "in", height = 3, width = 5, useDingbats = FALSE)

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
countMatrixPlot <- countMatrix %>% dplyr::select(sampleID, tpm, gene_name) %>%
  dplyr::filter(!(sampleID %in% c("1_DMSO", "1_LSD1i", "1_DOT1Li", "2_noDrug", "2_LSD1i"))) %>%
  dplyr::filter(gene_name %in% c("SPP1", "GDF15", "CDKN1A", "FTH1",
                                 "TOP2A", "MKI67", "CENPF", "SOX21")) %>%
  mutate(rep = word(.$sampleID, 1, sep = "_")) %>%
  mutate(cond = word(.$sampleID, 2, sep = "_"))

countMatrixPlot$cond <- factor(countMatrixPlot$cond, levels = c("high", "midhigh", "control", "midlow", "low"),
                               labels = c("slow", "mid slow", "ungated", "mid fast", "fast"))
countMatrixPlot$gene_name <- factor(countMatrixPlot$gene_name, levels = c("TOP2A", "MKI67", "CENPF", "SOX21",
                                                                          "SPP1", "GDF15", "CDKN1A", "FTH1"
                                                                          ))

ggplot(countMatrixPlot, aes(x = cond, y = tpm)) +
  stat_summary(fun = mean, geom = "bar") +
  geom_jitter(height = 0, width = 0.25, size = 0.25) +
  #stat_summary(fun = mean, geom = "point") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.25) +
  facet_wrap(~ gene_name, ncol = 4) + ylab("TPM") +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave(filename = paste0(plotDirectory, "prolifTPMforMarkerGenes.pdf"), units = "in", height = 3, width = 5, useDingbats = FALSE)

countMatrixPlotFormat <- data.frame(sampleID = character(),
                                    tpm = integer(),
                                    gene_name = character(),
                                    rep = character(),
                                    cond = character())
for(i in unique(countMatrixPlot$rep)) {
  for(j in unique(countMatrixPlot$gene_name)) {
    countMatrixPlotTemp <- countMatrixPlot %>% dplyr::filter(rep == i, gene_name == j)
    countMatrixPlotTemp <- countMatrixPlotTemp %>% rowwise() %>% mutate(countNorm = tpm/(countMatrixPlotTemp %>% filter(cond == "control"))$tpm) %>% ungroup()
    countMatrixPlotFormat <- bind_rows(countMatrixPlotFormat, countMatrixPlotTemp)
  }
}

countMatrixPlotFormat$cond <- factor(countMatrixPlotFormat$cond, levels = c("fast", "midfast", "control", "midslow", "slow"))
countMatrixPlotFormat$gene_name <- factor(countMatrixPlotFormat$gene_name, levels = c("SPP1", "GDF15", "CDKN1A", "TOP2A", "MKI67", "CENPF"))

ggplot(countMatrixPlotFormat, aes(x = cond, y = countNorm)) +
  geom_point() + geom_hline(yintercept = 1, linetype = "dashed") +
  stat_summary(fun.data = mean_se, geom = "errorbar") +
  facet_wrap(~ gene_name) + ylab("TPM normalized to control") +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = paste0(plotDirectory, "prolifTPMforMarkerGenes.pdf"), units = "in", height = 3, width = 6, useDingbats = FALSE)
