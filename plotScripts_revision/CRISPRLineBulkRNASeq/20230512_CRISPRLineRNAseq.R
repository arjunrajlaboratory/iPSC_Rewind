library(tidyverse)
library(tximport) 
library(ensembldb) 
library(EnsDb.Hsapiens.v86)
library(rhdf5)
library(reshape2)
library(DESeq2)
library(biomaRt)
library(gprofiler2)
library(Seurat)
library(ggrepel)

theme_set(theme_classic())

# read in sample info
# at home
targets <- read_tsv("/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Revised Manuscript/extractedData/CRISPRLineBulkRNASeq/samples.tsv") # read in your study design
path <- file.path(paste0("/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Revised Manuscript/extractedData/CRISPRLineBulkRNASeq/", targets$sample, "-1/", "abundance.h5")) # set file paths to your mapped data
all(file.exists(path))

plotDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Revised Manuscript/plots/CRISPRLineBulkRNASeq/"

# at lab
targets <- read_tsv("/Users/naveenj/Dropbox (RajLab)/Shared_Naveen/Revisions/extractedData/CRISPRLineBulkRNASeq/samples.tsv") # read in your study design
path <- file.path(paste0("/Users/naveenj/Dropbox (RajLab)/Shared_Naveen/Revisions/extractedData/CRISPRLineBulkRNASeq/", targets$sample, "-1/", "abundance.h5")) # set file paths to your mapped data
all(file.exists(path))

plotDirectory <- "/Users/naveenj/Dropbox (RajLab)/Shared_Naveen/Revisions/plots/CRISPRLineBulkRNASeq/"

# import abundance 
Tx <- transcripts(EnsDb.Hsapiens.v86, columns=c("tx_id", "gene_name"))
Tx <- as_tibble(Tx)
Tx <- dplyr::rename(Tx, target_id = tx_id)
Tx <- dplyr::select(Tx, "target_id", "gene_name")
Txi_gene <- tximport(path, 
                     type = "kallisto", 
                     tx2gene = Tx, 
                     txOut = FALSE, # determines whether your data represented at transcript or gene level
                     countsFromAbundance = "no",
                     ignoreTxVersion = TRUE)

dataTable <- as_tibble(Txi_gene$counts, rownames = NA)
colnames(dataTable) <- targets$sample
dataTable$gene <- rownames(dataTable)

dataTableMelt <- melt(dataTable, id.vars = c("gene"))
names <- dataTableMelt$variable %>% unique() %>% as.character()
guides <- sapply(names, function(x) {unlist(strsplit(x, "_")[[1]][1:2]) %>% paste(collapse = "_")})
replicates <- sapply(names, function(x) {unlist(strsplit(x, "_")[[1]][3])})

nameTable <- data.frame(names = names, guides = guides, replicates = replicates)
dataTableMeltNames <- inner_join(dataTableMelt, nameTable, c("variable" = "names"))
                                                   
geneList <- c("SPP1", "KDM1A", "UBC", "CDKN1A")
dataSelect <- dataTableMeltNames %>% dplyr::filter(gene %in% geneList)

geneListAverage <- c()
for(i in 1:length(geneList)) {
  dataTemp <- dataSelect %>% dplyr::filter(gene == geneList[i]) %>% dplyr::filter(guides %in% c("backbone_1", "backbone_2"))
  geneListAverage[i] <- mean(dataTemp$value)
}

geneListAverageTable <- data.frame(gene = geneList, average = geneListAverage)
dataSelect <- inner_join(dataSelect, geneListAverageTable, by = "gene")
                                                 
ggplot(dataSelect, aes(x = factor(guides, levels = sort(dataTableMeltNames$guides %>% unique() %>% as.character())), y = value)) +
  stat_summary(fun = mean, geom = "col") +
  #stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  geom_jitter(height = 0, width = 0.25, size = 1) +
  theme_classic() +
  geom_hline(data = dataSelect, aes(yintercept = average)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), axis.title.x = element_blank()) +
  facet_wrap(~gene, scales = "free", ncol = 4)
ggsave(filename = paste0(plotDirectory, "CRISPR_RNASeqKDLevels.pdf"), height = 2, width = 8)

geneList <- c("SMAD2", "SMAD3", "SMAD4", "TGFB1", "TGFB2", "TGFB3", "TGFBR1", "TGFBR2")
dataSelect <- dataTableMeltNames %>% dplyr::filter(gene %in% geneList)

geneListAverage <- c()
for(i in 1:length(geneList)) {
  dataTemp <- dataSelect %>% dplyr::filter(gene == geneList[i]) %>% dplyr::filter(guides %in% c("backbone_1", "backbone_2"))
  geneListAverage[i] <- mean(dataTemp$value)
}

geneListAverageTable <- data.frame(gene = geneList, average = geneListAverage)
dataSelect <- inner_join(dataSelect, geneListAverageTable, by = "gene")
                                                 
ggplot(dataSelect, aes(x = factor(guides, levels = sort(dataTableMeltNames$guides %>% unique() %>% as.character())), y = value)) +
  stat_summary(fun = mean, geom = "col") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  theme_classic() +
  geom_hline(data = dataSelect, aes(yintercept = average)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  facet_wrap(~gene, scales = "free")

nameTable <- data.frame(names = names, guides = guides, replicates = replicates)
dataTableMeltNames <- inner_join(dataTableMelt, nameTable, c("variable" = "names"))
dataTableMeltNames <- dataTableMeltNames %>% dplyr::filter(guides %in% c("backbone_1", "backbone_2", "SPP1_1", "SPP1_2", "SPP1_3"))

geneList <- c("SPP1")
dataSelect <- dataTableMeltNames %>% dplyr::filter(gene %in% geneList)

geneListAverage <- c()
for(i in 1:length(geneList)) {
  dataTemp <- dataSelect %>% dplyr::filter(gene == geneList[i]) %>% dplyr::filter(guides %in% c("backbone_1", "backbone_2"))
  geneListAverage[i] <- mean(dataTemp$value)
}

geneListAverageTable <- data.frame(gene = geneList, average = geneListAverage)
dataSelect <- inner_join(dataSelect, geneListAverageTable, by = "gene")

ggplot(dataSelect, aes(x = factor(variable, levels = sort(dataTableMeltNames$variable %>% unique() %>% as.character())), y = value)) +
  stat_summary(fun = mean, geom = "col") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  theme_classic() +
  geom_hline(data = dataSelect, aes(yintercept = average)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), axis.title.x = element_blank()) +
  facet_wrap(~gene, scales = "free", ncol = 3)

# logFCTable <- dataTableMeltNames %>% dplyr::filter(guides %in% c("backbone_1", "backbone_2", "SPP1_1", "SPP1_2", "SPP1_3"))
# logFCTableCast <- logFCTable %>% dcast(gene ~ guides, value.var = "value", fun.aggregate = mean)
# 
# logFCTableCast <- logFCTableCast %>% rowwise() %>% mutate(fcSPP1_1 = log2((SPP1_1 + 1) / (mean(c(backbone_1, backbone_2)) + 1)))
# logFCTableCast <- logFCTableCast %>% rowwise() %>% mutate(fcSPP1_2 = log2((SPP1_2 + 1) / (mean(c(backbone_1, backbone_2)) + 1)))
# logFCTableCast <- logFCTableCast %>% rowwise() %>% mutate(fcSPP1_3 = log2((SPP1_3 + 1) / (mean(c(backbone_1, backbone_2)) + 1)))
# 
# corrTable <- logFCTableCast %>% dplyr::select(-gene)
# res <- cor(corrTable)
# round(res, 3)
# 
# logFCTableCast <- logFCTableCast %>% rowwise() %>% mutate(fcMean = mean(c(fcSPP1_1, fcSPP1_2, fcSPP1_3)))

samples <- targets$sample
targets <- sapply(samples, function(x) {unlist(strsplit(x, "_")[[1]][1])})
guides <- sapply(samples, function(x) {unlist(strsplit(x, "_")[[1]][1:2]) %>% paste(collapse = "_")})
replicates <- sapply(samples, function(x) {unlist(strsplit(x, "_")[[1]][3])})

sampleTable <- data.frame(samples = samples, targets = targets, guides = guides, replicates = replicates)
rownames(sampleTable) <- colnames(Txi_gene)

dds <- DESeqDataSetFromTximport(Txi_gene, sampleTable, design = ~targets)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

vsd <- vst(dds)
head(assay(vsd), 3)

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vsd$samples
colnames(sampleDistMatrix) <- vsd$targets
library(gplots)
library(RColorBrewer)
library(ggrepel)
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.2( sampleDistMatrix, trace="none", col=colours)
plotPCA(vsd, intgroup = c("targets")) + geom_text_repel(aes(label = samples, color = targets)) + theme_classic() + NoLegend() + xlim(-5, 10) + ylim(-7.5, 7.5)
ggsave(filename = paste0(plotDirectory, "CRISPR_RNASeqPCA.pdf"), height = 4, width = 4)

rv <- rowVars(assay(object))
pca <- prcomp(t(assay(vsd)))
loadings <- as.data.frame(pca$rotation)

pc1_loadings <- loadings %>% dplyr::select(PC1) %>% mutate(gene = rownames(.))
ggplot() +
  geom_point(data = pc1_loadings %>% dplyr::filter(abs(PC1) <= 0.035), aes(x = 0, y = PC1), position = position_jitter(seed = 1234), color = "grey75") +
  geom_point(data = pc1_loadings %>% dplyr::filter(abs(PC1) > 0.035), aes(x = 0, y = PC1), position = position_jitter(seed = 1234), color = "red") +
  geom_text_repel(data = pc1_loadings %>% dplyr::filter(abs(PC1) > 0.035), aes(x = 0, y = PC1, label = gene), position = position_jitter(seed = 1234), color = "red", max.overlaps = 20) + coord_flip()
ggsave(filename = paste0(plotDirectory, "CRISPR_RNASeqPC1.pdf"), height = 5, width = 5)
pc2_loadings <- loadings %>% dplyr::select(PC2) %>% mutate(gene = rownames(.))
ggplot() +
  geom_point(data = pc2_loadings %>% dplyr::filter(abs(PC2) <= 0.035), aes(x = 0, y = PC2), position = position_jitter(seed = 1234), color = "grey75") +
  geom_point(data = pc2_loadings %>% dplyr::filter(abs(PC2) > 0.035), aes(x = 0, y = PC2), position = position_jitter(seed = 1234), color = "red") +
  geom_text_repel(data = pc2_loadings %>% dplyr::filter(abs(PC2) > 0.035), aes(x = 0, y = PC2, label = gene), position = position_jitter(seed = 1234), color = "red")
ggsave(filename = paste0(plotDirectory, "CRISPR_RNASeqPC2.pdf"), height = 5, width = 5)

dds <- DESeq(dds)

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
proteinCodingGenes <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = 'biotype', values = c('protein_coding'), mart = ensembl)
colnames(proteinCodingGenes)[1] <- "gene_id"
colnames(proteinCodingGenes)[2] <- "gene"

res_CDKN1A <- results(dds, contrast = c("targets", "CDKN1A", "backbone")) %>% as.data.frame() %>%
  replace(is.na(.), 0) %>%
  mutate(gene_id = rownames(.)) %>%
  dplyr::filter(baseMean > 20, padj < 0.15, abs(log2FoldChange) > 0.25) %>%
  dplyr::filter(gene_id %in% proteinCodingGenes$gene)

gost(query = res_CDKN1A %>% dplyr::filter(log2FoldChange < -0.25) %>% arrange(log2FoldChange) %>% .$gene,
     organism = "hsapiens", ordered_query = TRUE, as_short_link = TRUE, significant = TRUE,
     sources = c("REAC"))

res_KDM1A <- results(dds, contrast = c("targets", "KDM1A", "backbone")) %>% as.data.frame() %>%
  replace(is.na(.), 0) %>%
  mutate(gene_id = rownames(.)) %>%
  dplyr::filter(baseMean > 20, padj < 0.15, abs(log2FoldChange) > 0.25) %>%
  dplyr::filter(gene_id %in% proteinCodingGenes$gene)

gost(query = res_KDM1A %>% dplyr::filter(log2FoldChange < -0.25) %>% arrange(log2FoldChange) %>% .$gene,
     organism = "hsapiens", ordered_query = TRUE, as_short_link = TRUE, significant = TRUE,
     sources = c("REAC"))

res_SPP1 <- results(dds, contrast = c("targets", "SPP1", "backbone")) %>% as.data.frame() %>%
  replace(is.na(.), 0) %>%
  mutate(gene_id = rownames(.)) %>%
  dplyr::filter(baseMean > 20, padj < 0.15, abs(log2FoldChange) > 0.5) %>%
  dplyr::filter(gene_id %in% proteinCodingGenes$gene)

gost(query = res_SPP1 %>% dplyr::filter(log2FoldChange < -0.25) %>% arrange(log2FoldChange) %>% .$gene,
     organism = "hsapiens", ordered_query = TRUE, as_short_link = TRUE, significant = TRUE,
     sources = c("REAC"))
