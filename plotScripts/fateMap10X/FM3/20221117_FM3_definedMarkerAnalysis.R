rm(list=ls())
gc()

library(Seurat)
library(tidyverse)
library(RANN)
library(tripack)
library(reshape2)
library(svglite)
library(biomaRt)
library(ggpubr)
library(ggrastr)
library(ggsignif)
library(ggrepel)

homeDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/extractedData/fateMap10X/FM3/"
plotDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/plots/fateMap10X/FM3/"

theme_set(theme_classic())

linCountToOverlapsList <- as_tibble(read.table(file = paste0(homeDirectory, "filtered10XCells.txt"), header = TRUE, stringsAsFactors = F, sep = "\t"))

#### UMAP plots for selected barcodes ####
umapCoordinates = as_tibble(read.table(file = paste0(homeDirectory, "umapCoordinates_Scanorama_50pcs_filter.tsv"), header = TRUE, stringsAsFactors = F, sep = "\t"))
jointPCA = inner_join(linCountToOverlapsList %>% dplyr::filter(SampleNum %in% c("S1", "S2", "S5", "S6")), umapCoordinates, by = c("cellID")) %>% dplyr::select(-nLineages)
jointPCA1Barcodes = jointPCA %>% filter(SampleNum=="S1" | SampleNum=="S2") %>% dplyr::select(barcode) %>% unique()
jointPCA2Barcodes = jointPCA %>% filter(SampleNum=="S5" | SampleNum=="S6") %>% dplyr::select(barcode) %>% unique()
  
jointBarcodesOnlyBoth = inner_join(jointPCA2Barcodes, jointPCA1Barcodes, by = "barcode")
jointBarcodesOnlyBoth = jointBarcodesOnlyBoth %>% mutate(barcodeName = c(1:nrow(jointBarcodesOnlyBoth)))
jointBarcodesOnlyBoth$barcodeName = sub("^", "B", jointBarcodesOnlyBoth$barcodeName)

jointBarcodesOnlyBothList = jointBarcodesOnlyBoth$barcode
jointPCAOnlyBoth = jointPCA %>% filter(barcode %in% jointBarcodesOnlyBothList)
finalPCAJoint = inner_join(jointPCAOnlyBoth, jointBarcodesOnlyBoth, by = "barcode") %>% dplyr::select(-barcode,-nUMI)

finalPCAS1 = finalPCAJoint %>% filter(SampleNum == "S1" |
                                        SampleNum == "S2")
finalPCAS2 = finalPCAJoint %>% filter(SampleNum == "S5" |
                                        SampleNum == "S6")

finalPCAS1Big = finalPCAS1 %>% group_by(barcodeName) %>% summarise(nColony = length(barcodeName)) %>% filter(nColony > 0) %>% dplyr::select(-nColony)
finalPCAS2Big = finalPCAS2 %>% group_by(barcodeName) %>% summarise(nColony = length(barcodeName)) %>% filter(nColony > 0) %>% dplyr::select(-nColony)

commonS1S2Big = inner_join(finalPCAS1Big, finalPCAS2Big)
commonS1S2Big = commonS1S2Big$barcodeName
finalPCAJointBig = finalPCAJoint %>% filter(barcodeName %in% commonS1S2Big)
finalPCAJointBig = finalPCAJoint %>% mutate(label = "DMSO")
finalPCAJointBig$label = ifelse(finalPCAJointBig$SampleNum == "S5" | finalPCAJointBig$SampleNum == "S6", "LSD1i", finalPCAJointBig$label)
  
ggplot() +
      rasterise(geom_point(data = umapCoordinates, aes(x = UMAP_1, y = UMAP_2), color = "gray93"), dpi = 100) +
      rasterise(geom_point(data = finalPCAJointBig, aes(x = UMAP_1, y = UMAP_2, color = label), size = 2.5, shape = 16, alpha = 1), dpi = 100) +
      facet_wrap(~label) +
      scale_color_manual(values=c("hotpink3", "turquoise3")) + theme_classic() +
      theme(axis.text = element_blank(),
            axis.ticks = element_blank()) + NoLegend() + NoAxes() 
ggsave(filename = paste0(plotDirectory, "FM3_twinsUMAPDMSOvsLSD1i.pdf"), units = "in", width = 3.75, height = 3.75)

ggplot() +
  rasterise(geom_point(data = umapCoordinates, aes(x = UMAP_1, y = UMAP_2), color = "gray93"), dpi = 100) +
  rasterise(geom_point(data = finalPCAJointBig, aes(x = UMAP_1, y = UMAP_2, color = label), size = 2.5, shape = 16, alpha = 1), dpi = 100) +
  scale_color_manual(values=c("#f39cc2", "#00add0")) + theme_classic() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) + NoLegend() + NoAxes()
ggsave(filename = paste0(plotDirectory, "FM3_twinsUMAPDMSOvsLSD1i.pdf"), units = "in", width = 5, height = 5)

#### adding scoring modules for different markers ####
scanorama_filter <- readRDS(paste0(homeDirectory, "scanorama_filter_subset.rds"))

genesPluri <- list(c("NANOG", "ALPL", "POU5F1", "SOX2", "PODXL", "LIN28A", "UTF1", "TERT", "ZFP42", "DNMT3B"))
# genesPluriFull <- list(c('AASS', 'ACAA2', 'ALDH3B1', 'AMT', 'AP1M2', 'ARFGEF1', 'ASML3B', 'ATP6V0E', 'ATP6V1D', 'BMPR1A', 'CCNB1IP1', 'CDO1', 'CLN5', 'CNTNAP2', 'COBL', 'CRABP1', 'CRIM1', 'CRKL', 'CTSL2', 'CUGBP1', 'CXADR', 'CYP26A1', 'DAPK1', 'DCP2', 'DNMT3B', 'FGFR2', 'FKBP1B', 'FLJ20171', 'FRAT2', 'GCH1', 'GLDC', 'GPM6B', 'H3F3A', 'HEY2', 'HFE', 'HSD17B4', 'INDO', 'JARID2', 'KIF5C', 'LCK', 'LEFTY2', 'LIN28', 'LPP', 'MAP7', 'MBD2', 'MGC8685', 'MRS2L', 'MSH6', 'MYCN', 'NANOG', 'NUP210', 'OVOL2', 'PELI1', 'PHC1', 'PIPOX', 'PODXL', 'POU5F1', 'PPM1E', 'PPP1R16B', 'PRNP', 'PSIP1', 'RAB25', 'RABGAP1L', 'RND2', 'RNF44', 'RNMT', 'SALL2', 'SCNN1A', 'SEMG1', 'SEPHS1', 'SOX2', 'TACSTD1', 'TDGF1', 'TERF1', 'TMEM5', 'TNNI3', 'TNRC9', 'TRIM24', 'TXNRD1', 'UGT8'))
genesFibro <- list(c("LUM", "S100A4", "THY1", "PDGFRA", "COL1A1", "COL5A1", "LOXL1", "FBLN1", "FBLN2", "VTN"))
genesMyo <- list(c("ACTA2", "TPX2", "TAGLN", "MYL9", "CDKN1A", "CDKN2A", "CCN2", "GLI2"))
genesEpi <- list(c("CDH1", "CLDN3", "KRT3", "OCLN", "EPCAM", "ANPEP", "MUC1", "CD24"))
genesMes <- list(c("CDH2", "VIM", "FN1", "ZEB1", "SNAI2", "TWIST1", "TWIST2", "TGFB1"))

DefaultAssay(scanorama_filter) <- "RNA"
scanorama_filter <- NormalizeData(scanorama_filter)

scanorama_filter <- AddModuleScore(scanorama_filter, features = genesFibro, name = 'fibroMarkers')
scanorama_filter <- AddModuleScore(scanorama_filter, features = genesPluri, name = 'pluriMarkers')
scanorama_filter <- AddModuleScore(scanorama_filter, features = genesMyo, name = 'myoMarkers')
scanorama_filter <- AddModuleScore(scanorama_filter, features = genesEpi, name = 'epiMarkers')
scanorama_filter <- AddModuleScore(scanorama_filter, features = genesMes, name = 'mesMarkers')

mid <- c(mean(scanorama_filter$fibroMarkers1),
         median(scanorama_filter$pluriMarkers1),
         median(scanorama_filter$myoMarkers1),
         median(scanorama_filter$mesMarkers1),
         median(scanorama_filter$epiMarkers1))

# scanorama_filter$fibroMarkers1 <- scale(scanorama_filter$fibroMarkers1)
# scanorama_filter$pluriMarkers1 <- scale(scanorama_filter$pluriMarkers1)
# scanorama_filter$myoMarkers1 <- scale(scanorama_filter$myoMarkers1)
# scanorama_filter$mesMarkers1 <- scale(scanorama_filter$mesMarkers1)
# scanorama_filter$epiMarkers1 <- scale(scanorama_filter$epiMarkers1)

scanorama_filter <- FindClusters(scanorama_filter, resolution = 0.3, graph.name = "scanorama_snn_res.")

DimPlot(scanorama_filter)
FeaturePlot(scanorama_filter, features = c("fibroMarkers1", "pluriMarkers1", "myoMarkers1", "mesMarkers1", "epiMarkers1"), repel = TRUE, max.cutoff = "q80", min.cutoff = "q10", raster = TRUE, pt.size = 2, raster.dpi = c(300, 300)) &
  NoAxes() & scale_color_gradient2(low = "#0039a6", mid = "lightgray", high = "red", midpoint = 0.2)
# FeaturePlot(scanorama_filter, features = c("mesMarkers1", "epiMarkers1"), blend = TRUE, blend.threshold = 0.1, max.cutoff = "q80", min.cutoff = "q20") & NoLegend() & NoAxes() + DarkTheme()
# FeaturePlot(scanorama_filter, features = c("pluriMarkers1", "epiMarkers1"), blend = TRUE, blend.threshold = 0.1, max.cutoff = "q80", min.cutoff = "q20") & NoLegend() & NoAxes() + DarkTheme()
# FeaturePlot(scanorama_filter, features = c("mesMarkers1", "fibroMarkers1"), blend = TRUE, blend.threshold = 0.1, max.cutoff = "q80", min.cutoff = "q20") & NoLegend() & NoAxes() + DarkTheme()

FeaturePlot(scanorama_filter, features = "fibroMarkers1", max.cutoff = "q90", min.cutoff = "q2", repel = TRUE, raster = TRUE, pt.size = 2, raster.dpi = c(300, 300)) &
  NoAxes() & scale_color_gradient2(low = "#0039a6", mid = "lightgray", high = "red", midpoint = 0.2) & NoLegend() & labs(title = element_blank())
ggsave(filename = paste0(plotDirectory, "FM3_umap_fibroMarkers.pdf"), units = "in", width = 4, height = 4)
FeaturePlot(scanorama_filter, features = "pluriMarkers1", max.cutoff = "q90", min.cutoff = "q2", repel = TRUE, raster = TRUE, pt.size = 2, raster.dpi = c(300, 300)) &
  NoAxes() & scale_color_gradient2(low = "#0039a6", mid = "lightgray", high = "red", midpoint = 0.15) & NoLegend() & labs(title = element_blank())
ggsave(filename = paste0(plotDirectory, "FM3_umap_pluriMarkers.pdf"), units = "in", width = 4, height = 4)
FeaturePlot(scanorama_filter, features = "epiMarkers1", max.cutoff = "q90", min.cutoff = "q2", repel = TRUE, raster = TRUE, pt.size = 2, raster.dpi = c(300, 300)) &
  NoAxes() & scale_color_gradient2(low = "#0039a6", mid = "lightgray", high = "red", midpoint = 0.1) & NoLegend() & labs(title = element_blank())
ggsave(filename = paste0(plotDirectory, "FM3_umap_epiMarkers.pdf"), units = "in", width = 4, height = 4)
FeaturePlot(scanorama_filter, features = "mesMarkers1", max.cutoff = "q90", min.cutoff = "q2", repel = TRUE, raster = TRUE, pt.size = 2, raster.dpi = c(300, 300)) &
  NoAxes() & scale_color_gradient2(low = "#0039a6", mid = "lightgray", high = "red", midpoint = 0.2) & NoLegend() & labs(title = element_blank())
ggsave(filename = paste0(plotDirectory, "FM3_umap_mesMarkers.pdf"), units = "in", width = 4, height = 4)

FeaturePlot(scanorama_filter, features = c("NANOG", "ALPL", "PODXL", "POU5F1", "SOX2",
                                           "LUM", "TWIST2", "SNAI2", "FOSL1", "COL1A1"), slot = "data", ncol = 5)

cellIDsToJoin <- finalPCAJointBig %>% dplyr::select(cellID, label)
cellIDJoinTable <- left_join(umapCoordinates, finalPCAJointBig, by = "cellID")
cellIDJoinTable$label <- ifelse(is.na(cellIDJoinTable$label), "none", cellIDJoinTable$label)
labelsToAdd <- cellIDJoinTable$label
names(labelsToAdd) <- colnames(scanorama_filter)
scanorama_new <- AddMetaData(scanorama_filter, labelsToAdd, col.name = "condition")

Idents(scanorama_new) <- scanorama_new$condition
VlnPlot(scanorama_new, features = c("fibroMarkers1", "pluriMarkers1"))
scanorama_new_subset <- subset(scanorama_new, subset = condition %in% c("DMSO", "LSD1i"))
Idents(scanorama_new_subset) <- scanorama_new_subset$condition
VlnPlot(scanorama_new_subset, features = c("NANOG")) +
  stat_summary(fun = mean, geom = 'crossbar', size = 0.5, colour = "red") +
  stat_compare_means(comparisons = list(c("DMSO", "LSD1i")), label = "p.signif") + NoLegend()

pluriModuleScore = tibble(cellID = names(scanorama_new$pluriMarkers1), score = scanorama_new$pluriMarkers1)
# pluriModuleScore = tibble(cellID = names(scanorama_new$fibroMarkers1), score = scanorama_new$fibroMarkers1)
# pluriModuleScore = tibble(cellID = names(scanorama_new$epiMarkers1), score = scanorama_new$epiMarkers1)
# pluriModuleScore = tibble(cellID = names(scanorama_new$mesMarkers1), score = scanorama_new$mesMarkers1)
# pluriModuleScore$cellID = sub("-1", "", pluriModuleScore$cellID)
pluriModuleScore$cellID <- sapply(strsplit(sub("S\\d_", "", pluriModuleScore$cellID), "-"), "[", 1)

jointPluriModScore = inner_join(pluriModuleScore, finalPCAJointBig, by = "cellID")

barcodeNameList <- jointPluriModScore$barcodeName %>% unique()
meansDMSOList <- c()
meansLSD1List <- c()
meansDiffList <- c()
for(i in barcodeNameList){
  barcodeTable <- jointPluriModScore %>% filter(barcodeName == i)
  meansDMSO <- barcodeTable %>% filter(label == "DMSO") %>% dplyr::select(score) %>% colMeans()
  meansLSD1 <- barcodeTable %>% filter(label == "LSD1i") %>% dplyr::select(score) %>% colMeans()
  meansDMSOList <- c(meansDMSOList, meansDMSO)
  meansLSD1List <- c(meansLSD1List, meansLSD1)
  meansDiffList <- c(meansDiffList, meansLSD1 - meansDMSO)
}

meansDMSOList <- tibble(meanScores = meansDMSOList) %>% mutate(label = "DMSO")
meansLSD1List <- tibble(meanScores = meansLSD1List) %>% mutate(label = "LSD1i")
meansJointList <- bind_rows(meansDMSOList, meansLSD1List)
meansDiffList <- tibble(meanScoreDiff = meansDiffList)

ggplot(meansJointList, aes(x = label, y = meanScores)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2) +
  geom_signif(comparisons = list(c("DMSO", "LSD1i")))

ggplot(meansDiffList, aes(x = "test", y = meanScoreDiff)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2) + theme_classic()

barcodeNameList <- jointPluriModScore$barcodeName %>% unique()
for(i in 1:length(barcodeNameList)) {
  barcodeTable <- jointPluriModScore %>% filter(barcodeName == barcodeNameList[i])
  meansDMSO <- barcodeTable %>% filter(label == "DMSO") %>% dplyr::select(score) %>% colMeans()
  meansLSD1 <- barcodeTable %>% filter(label == "LSD1i") %>% dplyr::select(score) %>% colMeans()
  barcodeTableTemp <- tibble(barcode = barcodeNameList[i], DMSO = meansDMSO, LSD1 = meansLSD1)
  if(i == 1) {
    barcodeTablePaired <- barcodeTableTemp
  } else {
    barcodeTablePaired <- bind_rows(barcodeTablePaired, barcodeTableTemp)
  }
}

barcodeTablePairedMelt <- barcodeTablePaired %>% melt(., id.vars = c("barcode"))

ggplot(barcodeTablePairedMelt, aes(x = variable, y = value)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_line(aes(group = barcode), color = "gray") +
  geom_point() +
  stat_compare_means(comparisons = list(c("DMSO", "LSD1")), paired = TRUE, method = "wilcox.test")
ggsave(filename = paste0(plotDirectory, "FM3_pluriMarkerScoreBoxplot.pdf"), units = "in", width = 4, height = 4)

ggplot(barcodeTablePaired, aes(x = DMSO, y = LSD1)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  xlim(-1, 1) + ylim(-1, 1)
ggsave(filename = paste0(plotDirectory, "FM3_pluriMarkerScoreScatter.pdf"), units = "in", width = 4, height = 4)

t.test(barcodeTablePaired$DMSO, barcodeTablePaired$LSD1, paired = TRUE, alternative = "two.sided")

#### determining DE between twins in DMSO and LSD1i ####
scanorama_filter <- readRDS(paste0(homeDirectory, "scanorama_filter.rds"))
scanorama_filter <- FindClusters(object = scanorama_filter, resolution = 0.3)
DimPlot(scanorama_filter)

cellIDsToJoin <- finalPCAJointBig %>% dplyr::select(cellID, label)
cellIDJoinTable <- left_join(umapCoordinates, finalPCAJointBig, by = "cellID")
cellIDJoinTable$label <- ifelse(is.na(cellIDJoinTable$label), "none", cellIDJoinTable$label)

umapClusters = as_tibble(read.table(file = paste0(homeDirectory, "umapClusters_Scanorama_50pcs_filter.tsv"), header = TRUE, stringsAsFactors = F, sep = "\t")) %>%
  rename(clusters = 1)

cellIDJoinTable <- full_join(cellIDJoinTable, umapClusters, by = "cellID")

labelsToAdd <- cellIDJoinTable$label
names(labelsToAdd) <- colnames(scanorama_filter)
scanorama_new <- AddMetaData(scanorama_filter, labelsToAdd, col.name = "condition")
labelsToAdd <- cellIDJoinTable$barcodeName
names(labelsToAdd) <- colnames(scanorama_new)
scanorama_new <- AddMetaData(scanorama_new, labelsToAdd, col.name = "barcode")
scanorama_new <- subset(scanorama_new, subset = seurat_clusters %in% c(0, 1, 2, 3))

PrctCellExpringGene <- function(object, genes, group.by = "all"){
  if(group.by == "all"){
    prct = unlist(lapply(genes,calc_helper, object=object))
    result = data.frame(Markers = genes, Cell_proportion = prct)
    return(result)
  }
  
  else{        
    list = SplitObject(object, group.by)
    factors = names(list)
    
    results = lapply(list, PrctCellExpringGene, genes=genes)
    for(i in 1:length(factors)){
      results[[i]]$Feature = factors[i]
    }
    combined = do.call("rbind", results)
    return(combined)
  }
}
calc_helper <- function(object,genes){
  counts = object[['RNA']]@counts
  ncells = ncol(counts)
  if(genes %in% row.names(counts)){
    sum(counts[genes,]>0)/ncells
  }else{return(NA)}
}
percentExpressionPluri <- PrctCellExpringGene(object = scanorama_new, genes = unlist(genesPluri), group.by = "condition")
ggplot(percentExpressionPluri, aes(x = Markers, y = Cell_proportion, fill = Feature)) +
  geom_col(position = position_dodge())

DefaultAssay(scanorama_new) <- "RNA"
scanorama_new <- NormalizeData(scanorama_new)

barcodeNameList <- finalPCAJointBig$barcodeName %>% unique()

#### if using only low mixing coefficient lineages ####
# cellIDJoinTableWBarcode <- left_join(cellIDJoinTable, linCountToOverlapsList, by = "cellID")
# exportTable <- readRDS(file = paste0(homeDirectory, "twinBarcodesDMSOvsLSD1wMixing.rds"))
# lowMixBarcodes <- cellIDJoinTableWBarcode %>% dplyr::filter(barcode %in% (dplyr::filter(exportTable, mean < 0.4) %>% .$barcode)) %>% dplyr::select(label, barcode, barcodeName)
# barcodeNameList <- lowMixBarcodes$barcodeName %>% unique() %>% .[!(is.na(.))]
#######################################################

# barcodeNameList <- barcodeNameList[-22]
# 
# for(i in 1:length(barcodeNameList)) {
#   scanorama_subset_DE <- subset(scanorama_new, subset = barcode == barcodeNameList[i])
#   resultsTableTemp <- FindMarkers(object = scanorama_subset_DE, group.by = "condition", ident.1 = "DMSO", ident.2 = "LSD1i", min.cells.group = 0, logfc.threshold = 0)
#   resultsTableTemp <- resultsTableTemp %>% dplyr::select(avg_log2FC, p_val_adj) %>% mutate(gene = rownames(.), barcode = barcodeNameList[i])
#   if(i == 1) {
#     resultsTableDE <- resultsTableTemp
#   } else {
#     resultsTableDE <- bind_rows(resultsTableDE, resultsTableTemp)
#   }
# }
# 
# saveRDS(object = resultsTableDE, file = paste0(homeDirectory, "resultsDE_DMSOvsLSD1i.rds"))
# saveRDS(object = resultsTableDE, file = paste0(homeDirectory, "resultsDE_DMSOvsLSD1i_lowMix.rds"))

resultsTableDE <- readRDS(paste0(homeDirectory, "resultsDE_DMSOvsLSD1i.rds"))
# resultsTableDE <- readRDS(paste0(homeDirectory, "resultsDE_DMSOvsLSD1i_lowMix.rds"))

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
proteinCodingGenes <- getBM(attributes = c("hgnc_symbol"), filters = 'biotype', values = c('protein_coding'), mart = ensembl)

resultsTableDEFilter <- resultsTableDE %>% dplyr::filter(gene %in% proteinCodingGenes$hgnc_symbol)
resultsTableMean <- resultsTableDEFilter %>% group_by(gene) %>% summarize(meanLog2FC = mean(avg_log2FC), sdev = sd(avg_log2FC, na.rm=TRUE), num = n(), comb.p = -2*sum(log(p_val_adj)))
resultsTableMean <- resultsTableMean %>% rowwise() %>% mutate(chisq = pchisq(comb.p, df = 2*num, lower.tail = FALSE))
resultsTableMeanFilter <- resultsTableMean %>% dplyr::filter(num > 5)
resultsTableMeanPlotL <- resultsTableMeanFilter %>% as.data.frame() %>% dplyr::slice_min(meanLog2FC, n = 25, with_ties = FALSE)
resultsTableMeanPlotH <- resultsTableMeanFilter %>% as.data.frame() %>%dplyr::slice_max(meanLog2FC, n = 25, with_ties = FALSE)
resultsTableMeanPlot <- bind_rows(resultsTableMeanPlotL, resultsTableMeanPlotH)

resultsTableDEFilter$gene <- with(resultsTableDEFilter, reorder(gene, avg_log2FC, mean))
ggplot(resultsTableDEFilter %>% filter(gene %in% resultsTableMeanPlot$gene), aes(x = gene, y = avg_log2FC)) +
  # geom_boxplot() +
  # geom_point() +
  stat_summary(fun = mean, geom = "point", fill = "black", size = 2.5, shape = 16) +
  # stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.25) +
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.25) +
  stat_summary(fun = mean, geom = "text", aes(label = gene), angle = -45, hjust = 0, vjust = 0, size = 3, position = position_nudge(x = 0.1, y = -0.05)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())
ggsave(filename = paste0(plotDirectory, "markerDE_DMSOvsLSD1iTwins.pdf"), height = 2, width = 10, useDingbats = FALSE)

ggplot(resultsTableMeanFilter, aes(x = sdev, y = meanLog2FC)) +
  geom_point() +
  geom_text_repel(data = resultsTableMeanFilter %>% dplyr::filter(gene %in% unlist(genesPluri)), aes(label = gene), color = "red") +
  geom_text_repel(data = resultsTableMeanFilter %>% dplyr::filter(gene %in% unlist(genesFibro)), aes(label = gene), color = "blue") +
  geom_text_repel(data = resultsTableMeanFilter %>% dplyr::filter(gene %in% unlist(genesMes)), aes(label = gene), color = "green") +
  geom_text_repel(data = resultsTableMeanFilter %>% dplyr::filter(gene %in% unlist(genesEpi)), aes(label = gene), color = "orange") +
  geom_text_repel(data = resultsTableMeanFilter %>% dplyr::filter(gene %in% unlist(genesMyo)), aes(label = gene), color = "yellow") +
  geom_hline(yintercept = c(-0.5, 0.5), linetype = "dashed")

ggplot(resultsTableMeanFilter, aes(x = meanLog2FC, y = -log(chisq))) +
  geom_point() +
  geom_hline(yintercept = -log(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  geom_text_repel(data = resultsTableMeanFilter %>% dplyr::filter(abs(meanLog2FC) > 0.5), aes(x = meanLog2FC, y = -log(chisq), label = gene))

#### doing DE just in bulk instead of pairwise for more power ####
# scanorama_subset_DE <- subset(scanorama_new, subset = barcode %in% barcodeNameList)
# resultsTable <- FindMarkers(object = scanorama_subset_DE, group.by = "condition", ident.1 = "DMSO", ident.2 = "LSD1i", logfc.threshold = 0)
# resultsTable <- resultsTable %>% mutate(gene = rownames(.))
# saveRDS(object = resultsTable, file = paste0(homeDirectory, "resultsDE_DMSOvsLSD1i_bulk.rds"))

resultsTableBulk <- readRDS(file = paste0(homeDirectory, "resultsDE_DMSOvsLSD1i_bulk.rds"))

ggplot(resultsTableBulk, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = c(-0.5, 0, 0.5))

geneList <- c("NANOG", "ALPL", "POU5F1", "SOX2", "PODXL", "LIN28A", "UTF1", "TERT", "ZFP42", "DNMT3B", "GAS5", "TERF1", "TRIM28", "ZNF483", "NLRP7", "FGF2", "DPPA5", "EPCAM", "CDH1", "NODAL", "LEFTY1")
ggplot(resultsTableBulk %>% dplyr::filter(gene %in% geneList), aes(x = gene, y = avg_log2FC)) +
  geom_col()

resultsTableFilterBulk <- resultsTableBulk %>% dplyr::filter(p_val_adj <= 0.05)

ggplot(resultsTableFilterBulk, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point() + geom_text_repel(aes(label = gene), max.overlaps = 20)

ggplot(resultsTableBulk %>% dplyr::filter(gene %in% cluster0Markers), aes(x = factor(gene, levels = cluster0Markers), y = avg_log2FC)) +
  geom_col() +
  geom_text(aes(label = signif(p_val_adj, 3)), vjust = 0) + ylim(0, 2)
ggsave(filename = paste0(plotDirectory, "FM3_cluster0MarkerBarGraphs.pdf"), units = "in", height = 1, width = 4)

ggplot(resultsTableBulk %>% dplyr::filter(gene %in% cluster1Markers), aes(x = factor(gene, levels = cluster1Markers), y = -avg_log2FC)) +
  geom_col() +
  geom_text(aes(label = signif(p_val_adj, 3)), vjust = 0) + ylim(0, 2)
ggsave(filename = paste0(plotDirectory, "FM3_cluster1MarkerBarGraphs.pdf"), units = "in", height = 1, width = 4)

ggplot(resultsTableBulk %>% dplyr::filter(gene %in% c(cluster3Markers, cluster2Markers)), aes(x = factor(gene, levels = c(cluster3Markers, cluster2Markers)), y = -avg_log2FC)) +
  geom_col() +
  geom_text(aes(label = signif(p_val_adj, 3)), vjust = 0) + ylim(0, 2)
ggsave(filename = paste0(plotDirectory, "FM3_cluster3and2MarkerBarGraphs.pdf"), units = "in", height = 1, width = 4)

#### analyzing gene expression via clusters ####
scanorama_filter <- readRDS(file = paste0(homeDirectory, "scanorama_filter_iPSC.rds"))
DimPlot(object = scanorama_filter)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
scanorama_filter <- CellCycleScoring(object = scanorama_filter, s.features = s.genes, g2m.features = g2m.genes)
DimPlot(object = scanorama_filter, group.by = "Phase")

# DefaultAssay(scanorama_filter_subset) <- "RNA"
# scanorama_filter_subset <- NormalizeData(scanorama_filter_subset)
# markers <- FindAllMarkers(scanorama_filter_subset)
# markers <- markers %>% mutate(gene = rownames(markers))
# markers$gene <- gsub("\\..*", "", markers$gene)
# saveRDS(markers, file = paste0(homeDirectory, "clusterFindAllMarkers.RDS"))

markersClusters <- readRDS(paste0(homeDirectory, "clusterFindAllMarkers.RDS"))

markersClust0 <- markersClusters %>% filter(cluster == "0")
markersClust1 <- markersClusters %>% filter(cluster == "1")
markersClust2 <- markersClusters %>% filter(cluster == "2")
markersClust3 <- markersClusters %>% filter(cluster == "3")
markersClust4 <- markersClusters %>% filter(cluster == "4")
markersClust5 <- markersClusters %>% filter(cluster == "5")

DefaultAssay(scanorama_filter) <- "RNA"
scanorama_filter <- NormalizeData(scanorama_filter)
scanorama_filter <- ScaleData(scanorama_filter)

# FeaturePlot(scanorama_filter, slot = "data", features = c("OTX2", "SOX1", #ectoderm
#                                                                        "TBXT", "HAND1", #mesoderm
#                                                                        "SOX17", "GATA4"), #endoderm
#             max.cutoff = "q95") & NoAxes() & scale_color_gradient2(low = "blue", mid = "lightgray", high = "red", midpoint = 0.5)

FeaturePlot(scanorama_filter, slot = "data", features = c("GAPDH", "UBC"), max.cutoff = "q98") &
  NoAxes() & scale_color_gradient2(low = "blue", mid = "lightgray", high = "red", midpoint = 0.5)

FeaturePlot(scanorama_filter, slot = "data", features = c("POU5F1", "SOX2", "NANOG", "PODXL", "DNMT3B", "LIN28A", "GAPDH", "UBC"), max.cutoff = "q98", ncol = 2, raster = TRUE, raster.dpi = c(300, 300)) &
  NoAxes() & NoLegend() & scale_color_gradient2(low = "blue", mid = "lightgray", high = "red", midpoint = 0.5)
ggsave(filename = paste0(plotDirectory, "FM3_pluriMarkersUMAPS.pdf"), units = "in", height = 8, width = 4)

# DefaultAssay(scanorama_filter) <- "scanorama"
# FeaturePlot(scanorama_filter, slot = "scale.data", features = c("POU5F1", "SOX2", "NANOG", "PODXL", "DNMT3B", "LIN28A"), max.cutoff = "q98", ncol = 3, raster = TRUE, raster.dpi = c(300, 300)) &
#   NoAxes() & scale_color_gradient2(low = "blue", mid = "lightgray", high = "red", midpoint = -1)

cluster0Markers <- c("TPI1", "LDHA", #glycolysis markers
                     "GAS5", "TERF1", #pluripotency markers
                     "MT1G", "MT1X" #metallothionein
)
FeaturePlot(scanorama_filter, slot = "data", features = c(cluster0Markers), max.cutoff = "q98", ncol = 2, raster = TRUE, raster.dpi = c(300, 300)) &
  NoAxes() & NoLegend() & scale_color_gradient2(low = "blue", mid = "lightgray", high = "red", midpoint = 1)
ggsave(filename = paste0(plotDirectory, "FM3_cluster0MarkerUMAPS.pdf"), units = "in", height = 6, width = 4)

cluster1Markers <- c("TRIM28", "ZNF483",
                     "NLRP7", "FGF2",
                     "LIN28A", "DPPA5" #pluripotency markers
)
FeaturePlot(scanorama_filter, slot = "data", features = c(cluster1Markers), max.cutoff = "q98", ncol = 2, raster = TRUE, raster.dpi = c(300, 300)) &
  NoAxes() & NoLegend() & scale_color_gradient2(low = "blue", mid = "lightgray", high = "red", midpoint = 1)
ggsave(filename = paste0(plotDirectory, "FM3_cluster1MarkerUMAPS.pdf"), units = "in", height = 6, width = 4)

cluster3Markers <- c("EPCAM", "CDH1", #epithelial markers
                     "PODXL", "ALPL" #pluripotency markers
)
FeaturePlot(scanorama_filter, slot = "data", features = c(cluster3Markers), max.cutoff = "q98", ncol = 2, raster = TRUE, raster.dpi = c(300, 300)) &
  NoAxes() & NoLegend() & scale_color_gradient2(low = "blue", mid = "lightgray", high = "red", midpoint = 1)
ggsave(filename = paste0(plotDirectory, "FM3_cluster3MarkerUMAPS.pdf"), units = "in", height = 4, width = 4)

cluster2Markers <- c("NODAL", "LEFTY1" #hematopoietic predisposition
)
FeaturePlot(scanorama_filter, slot = "data", features = c(cluster2Markers), max.cutoff = "q98", ncol = 2, raster = TRUE, raster.dpi = c(300, 300)) &
  NoAxes() & NoLegend() & scale_color_gradient2(low = "blue", mid = "lightgray", high = "red", midpoint = 1)
ggsave(filename = paste0(plotDirectory, "FM3_cluster2MarkerUMAPS.pdf"), units = "in", height = 2, width = 4)

#### check pathways manually using GO ####
pathway_TGFB <- getBM(attributes=c('hgnc_symbol'), filters = 'go', values = 'GO:0007179', mart = ensembl)
pathway_WNT <- getBM(attributes=c('hgnc_symbol'), filters = 'go', values = 'GO:0016055', mart = ensembl)
pathway_NOTCH <- getBM(attributes=c('hgnc_symbol'), filters = 'go', values = 'GO:0007219', mart = ensembl)

DefaultAssay(scanorama_filter) <- "RNA"
scanorama_filter <- NormalizeData(scanorama_filter)

scanorama_filter <- AddModuleScore(scanorama_filter, features = pathway_TGFB, name = 'TGFB_path')
scanorama_filter <- AddModuleScore(scanorama_filter, features = pathway_WNT, name = 'WNT_path')
FeaturePlot(scanorama_filter, features = c("TGFB_path1", "WNT_path1")) &
  NoAxes() & scale_color_gradient2(low = "#0039a6", mid = "lightgray", high = "red", midpoint = 0.2)

#### check log normalized counts for different pluripotency markers ####
logNormCounts = as_tibble(read.table(file = paste0(homeDirectory, "logNormalizedCounts_Scanorama_50pcs_filterRound.tsv"), header = TRUE, stringsAsFactors = F, sep = "\t"))
logNormCountsSubset <- logNormCounts %>% dplyr::select(unlist(genesPluri), cellID, sampleNum)
saveRDS(object = logNormCountsSubset, file = paste0(homeDirectory, "logNormSubsetPluri.rds"))
rm(logNormCounts)

logNormPluriTwins <- inner_join(finalPCAJointBig, logNormCountsSubset, by = "cellID") %>% dplyr::select(-cellID, -SampleNum, -fracUMI, -sampleNum.x, -sampleNum.y)
logNormPluriTwinsMelt <- melt(logNormPluriTwins, id.vars = c("UMAP_1", "UMAP_2", "barcodeName", "label"))
ggplot(logNormPluriTwinsMelt, aes(x = label, y = value)) +
  geom_boxplot() +
  facet_wrap(~variable)

library(ggridges)
ggplot(logNormPluriTwins, aes(y = label, x = numberExp)) +
  geom_density_ridges(quantile_lines = TRUE)

logNormPluriTwins <- logNormPluriTwins %>% rowwise() %>% mutate(genesExp = which(NANOG > 0, ALPL > 0, POU5F1 > 0, SOX2 > 0, PODXL > 0, LIN28A > 0, UTF1 > 0, TERT > 0, ZFP42 > 0, DNMT3B > 0))
pluriGenesCoExp <- logNormPluriTwinsMelt %>% group_by(variable, label) %>% summarise(numberExp = sum(value > 0), n = n()) %>% ungroup()
pluriGenesCoExp <- pluriGenesCoExp %>% group_by(label) %>% mutate(fracExp = numberExp/n)
ggplot(pluriGenesCoExp, aes(x = variable, y = fracExp, fill = label, group = label)) +
  geom_col(position = position_dodge())
