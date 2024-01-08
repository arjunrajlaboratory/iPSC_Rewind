rm(list=ls())
gc()

library(Seurat)
library(tidyverse)
library(RANN)
library(tripack)
library(reshape2)
library(svglite)
library(ggrastr)

theme_set(theme_classic())

homeDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/extractedData/fateMap10X/FM3/"
plotDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/plots/fateMap10X/FM3/"
pcaCoordinates = as_tibble(read.table(file = paste0(homeDirectory, "pcaCoordinates_Scanorama_50pcs_filter.tsv"), header = TRUE, stringsAsFactors = F, sep = "\t"))

#### functions
get_lower_tri<-function(neightborMatrix2){
  neightborMatrix2[upper.tri(neightborMatrix2)] <- NA
  return(neightborMatrix2)
}
get_upper_tri <- function(neightborMatrix){
  neightborMatrix[lower.tri(neightborMatrix)]<- NA
  return(neightborMatrix)
}

linCountToOverlapsList <- as_tibble(read.table(file = paste0(homeDirectory, "filtered10XCells.txt"), header = TRUE, stringsAsFactors = F, sep = "\t"))

################################################################################################################################################################################################################################################
#### calculate pairwise mixing coefficients ####
################################################################################################################################################################################################################################################
plotList1 <- list()
matrixList <- list()
fullMatrixList <- list()
barcodeList <- list()
# twinAList <- c(1, 5, 3, 1, 1, 2, 2, 1, 1, 2, 2)
# twinBList <- c(2, 6, 4, 5, 6, 5, 6, 3, 4, 3, 4)
# comparisonList <- c("DMSO self", "LSD1i self", "DOT1Li self", "DMSO vs. LSD1i", "DMSO vs. LSD1i", "DMSO vs. LSD1i", "DMSO vs. LSD1i", "DMSO vs. DOT1Li", "DMSO vs. DOT1Li", "DMSO vs. DOT1Li", "DMSO vs. DOT1Li")

twinAList <- list(1, 5, c(1, 2))
twinBList <- list(2, 6, c(5, 6))
comparisonList <- c("DMSO self", "LSD1i self", "DMSO vs. LSD1i")

for(n in 1:length(twinAList)){
  if(length(twinAList[[n]] == 1)) {
    jointPCA = inner_join(linCountToOverlapsList %>% dplyr::filter(SampleNum %in% c(paste0("S", twinAList[[n]]), paste0("S", twinBList[[n]]))), pcaCoordinates, by = c("cellID")) %>% dplyr::select(-nLineages)
    jointPCA1Barcodes = jointPCA %>% filter(SampleNum==paste0("S", twinAList[[n]])) %>% dplyr::select(barcode) %>% unique()
    jointPCA2Barcodes = jointPCA %>% filter(SampleNum==paste0("S", twinBList[[n]]))  %>% dplyr::select(barcode) %>% unique()
  } else {
    jointPCA = inner_join(linCountToOverlapsList %>% dplyr::filter(SampleNum %in% c(paste0("S", twinAList[[n]][1]), paste0("S", twinAList[[n]][2]), paste0("S", twinBList[[n]][1]), paste0("S", twinBList[[n]][2]))), pcaCoordinates, by = c("cellID")) %>% dplyr::select(-nLineages)
    jointPCA1Barcodes = jointPCA %>% filter(SampleNum %in% c(paste0("S", twinAList[[n]][1]), paste0("S", twinAList[[n]][2]))) %>% dplyr::select(barcode) %>% unique()
    jointPCA2Barcodes = jointPCA %>% filter(SampleNum %in% c(paste0("S", twinBList[[n]][1]), paste0("S", twinBList[[n]][2]))) %>% dplyr::select(barcode) %>% unique()
  }

  jointBarcodesOnlyBoth = inner_join(jointPCA2Barcodes,jointPCA1Barcodes, by = "barcode") 
  jointBarcodesOnlyBoth = jointBarcodesOnlyBoth %>% mutate(barcodeName = c(1:nrow(jointBarcodesOnlyBoth)))
  jointBarcodesOnlyBoth$barcodeName = sub("^", "B", jointBarcodesOnlyBoth$barcodeName)
  
  jointBarcodesOnlyBothList = jointBarcodesOnlyBoth$barcode
  jointPCAOnlyBoth = jointPCA %>% filter(barcode %in% jointBarcodesOnlyBothList)
  finalPCAJoint = inner_join(jointPCAOnlyBoth,jointBarcodesOnlyBoth, by = "barcode") %>% dplyr::select(-cellID, -barcode, -nUMI)
  
  if(length(twinAList[[n]] == 1)) {
    finalPCAS1 = finalPCAJoint %>% filter(SampleNum==paste0("S", twinAList[[n]]))
    finalPCAS2 = finalPCAJoint %>% filter(SampleNum==paste0("S", twinBList[[n]]))
  } else{
    finalPCAS1 = finalPCAJoint %>% filter(SampleNum %in% c(paste0("S", twinAList[[n]][1]), paste0("S", twinAList[[n]][2])))
    finalPCAS2 = finalPCAJoint %>% filter(SampleNum %in% c(paste0("S", twinBList[[n]][1]), paste0("S", twinBList[[n]][2])))
  }
  
  finalPCAS1Big = finalPCAS1 %>% group_by(barcodeName) %>% summarise(nColony = length(barcodeName)) %>% filter(nColony > 3) %>% dplyr::select(-nColony)
  finalPCAS2Big = finalPCAS2 %>% group_by(barcodeName) %>% summarise(nColony = length(barcodeName)) %>% filter(nColony > 3) %>% dplyr::select(-nColony)
  
  commonS1S2Big = inner_join(finalPCAS1Big, finalPCAS2Big)
  commonS1S2Big = commonS1S2Big$barcodeName
  finalPCAJointBig = finalPCAJoint %>% filter(barcodeName %in% commonS1S2Big)
  barcodeList[[n]] <- unique(finalPCAJointBig$barcodeName)
  
  neightborMatrix = matrix(, nrow = length(unique(finalPCAJointBig$barcodeName)), ncol = length(unique(finalPCAJointBig$barcodeName)))
  for (i in c(1:length(unique(finalPCAJointBig$barcodeName)))) {
    for (j in c(1:length(unique(finalPCAJointBig$barcodeName)))) {
      barcode = unique(finalPCAJointBig$barcodeName)
      iBarcodeBx = finalPCAJointBig %>% filter(SampleNum == paste0("S", twinAList[[n]]), barcodeName == barcode[i]) %>% dplyr::select(-barcodeName)
      jBarcodeBx = finalPCAJointBig %>% filter(SampleNum == paste0("S", twinBList[[n]]), barcodeName == barcode[j]) %>% dplyr::select(-barcodeName)
      ijBarcodeBx = bind_rows(iBarcodeBx,jBarcodeBx)
      BarcodeRef = ijBarcodeBx %>% mutate(num = c(1:nrow(ijBarcodeBx)))
      S1index = BarcodeRef %>% filter(SampleNum == paste0("S", twinAList[[n]])) %>% dplyr::select(num)
      S2index = BarcodeRef %>% filter(SampleNum == paste0("S", twinBList[[n]])) %>% dplyr::select(num)
      
      knnPCA = nn2(ijBarcodeBx[, 3:52], ijBarcodeBx[,3:52])
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
      neightborMatrix[i,j]= mean(withOther$fraction)/mean(withSelf$fraction)
    }
  }
  
  fullMatrixList[[n]] <- neightborMatrix
  
  diagonal = diag(neightborMatrix)
  diagonal = diagonal[1:length(unique(finalPCAJointBig$barcodeName))]
  diagonalMatrix = diag(diagonal)
  diagonalMatrix[col(diagonalMatrix)!=row(diagonalMatrix)] = NA
  
  neightborMatrixRevised = neightborMatrix[1:length(unique(finalPCAJointBig$barcodeName)), 1:length(unique(finalPCAJointBig$barcodeName))]
  
  # rowNames = c(1:length(unique(finalPCAJointBig$barcodeName)))
  # rowNames = sub("^", "Twin ", rowNames)
  rowNames = barcode
  colNames = rowNames
  rownames(diagonalMatrix) = rowNames
  colnames(diagonalMatrix) = colNames
  melted_Matrix <- melt(diagonalMatrix, na.rm = TRUE)
  melted_Matrix$value <- round(melted_Matrix$value, 2)
  melted_Matrix <- melted_Matrix %>% mutate(ident = rep(comparisonList[[n]], nrow(.)))
  matrixList[[n]] <- melted_Matrix
  
  plotList1[[n]] <- ggplot(data = melted_Matrix, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "white", high = "black", 
                         limit = c(0,max(c(1, melted_Matrix$value))), 
                         breaks = c(0,0.5,1),
                         labels = c(0,0.5,1),
                         space = "Lab", 
                         name="Mixing\nCoefficient") +
    theme_minimal() + 
    coord_fixed() +
    geom_text(aes(Var2, Var1, label = value, color = value), size = 5.5) +
    scale_color_gradient(low = "black", high = "black", limits = c(0, 0.7), na.value = "white") +
    theme_classic((base_size = 26)) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.line = element_line(colour = 'black', size = 1.5),
      axis.text.x = element_text(angle = 90),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_line(colour = "black", size = 1.5),
      legend.position = "none",
      text=element_text(family="Helvetica")) +
    guides(fill = guide_colorbar(barwidth = 7, barheight = 1, title.position = "top", title.hjust = 0.5)) +
    scale_y_discrete(position = "right")
  ggsave(plotList1[[n]], file = paste0(plotDirectory, "neighborMatrices/FM3_mixCoeffCorrPlot_", n, ".pdf"), width = 6, height = 6)
}
ggpubr::ggarrange(plotlist = plotList1, ncol = 4)

combinedMixCoeffBoxplot <- bind_rows(matrixList[[1]], matrixList[[2]], matrixList[[3]])
ggplot(combinedMixCoeffBoxplot, aes(x = ident, y = value)) +
  geom_boxplot() +
  geom_jitter(width = 0.25) +
  geom_signif(comparisons = list(c("LSD1i self", "DMSO vs. LSD1i"), c("DMSO self", "DMSO vs. LSD1i")))
ggsave(filename = paste0(plotDirectory, "FM3_mixingCoeffBoxplotsAcrossConditions.pdf"), height = 2, width = 5)

################################################################################################################################################################################################################################################
#### previous analysis (DO NOT RUN) ####
################################################################################################################################################################################################################################################
# matrixDMSOvsLSD1i <- bind_rows(matrixList[[4]], matrixList[[5]], matrixList[[6]], matrixList[[7]])
# diagonalMatrix = diag(matrixDMSOvsLSD1i$value)
# diagonalMatrix[col(diagonalMatrix)!=row(diagonalMatrix)] = NA
# rowNames = rownames(matrixDMSOvsLSD1i)
# rowNames = sub("^", "Twin ", rowNames)
# colNames = rowNames
# rownames(diagonalMatrix) = rowNames
# colnames(diagonalMatrix) = colNames
# matrixDMSOvsLSD1i <- melt(diagonalMatrix, na.rm = TRUE)
# plot <- ggplot(data = matrixDMSOvsLSD1i, aes(Var2, Var1, fill = value)) +
#   geom_tile(color = "white")+
#   scale_fill_gradient2(low = "white", high = "black", 
#                        limit = c(0,max(c(1, melted_Matrix$value))), 
#                        breaks = c(0,0.5,1),
#                        labels = c(0,0.5,1),
#                        space = "Lab", 
#                        name="Mixing\nCoefficient") +
#   theme_minimal()+ 
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) +
#   coord_fixed() +
#   geom_text(aes(Var2, Var1, label = value), color = "black", size = 5.5) +
#   theme_classic((base_size = 26)) +
#   theme(
#     axis.title.x = element_blank(),
#     axis.title.y = element_blank(),
#     axis.line = element_line(colour = 'black', size = 1.5),
#     axis.text.x = element_text(angle = 90),
#     panel.grid.major = element_blank(),
#     panel.border = element_blank(),
#     panel.background = element_blank(),
#     axis.ticks = element_line(colour = "black", size = 1.5),
#     legend.justification = c(1, 0),
#     legend.position = c(0.6, 0.8),
#     legend.direction = "horizontal",
#     text=element_text(family="Helvetica")) +
#   guides(fill = guide_colorbar(barwidth = 7, barheight = 1, title.position = "top", title.hjust = 0.5)) +
#   scale_y_discrete(position = "right") +
#   theme(legend.position = "none")
# ggsave(plot, file = paste0(plotDirectory, "neighborMatrices/FM3_mixCoeffCorrPlot_DMSOvsLSD1i.pdf"), width = 12, height = 12)
# 
# combinedBarcodes <- c(barcodeList[[4]], barcodeList[[5]], barcodeList[[6]], barcodeList[[7]])
# combinedBarcodesMixing <- tibble(barcode = combinedBarcodes, mixing = matrixDMSOvsLSD1i$value) %>% group_by(barcode) %>% summarise(mean = mean(mixing))
# 
# for(n in 4:7) {
#   jointPCA = inner_join(linCountToOverlapsList %>% dplyr::filter(SampleNum %in% c(paste0("S", twinAList[[n]]), paste0("S", twinBList[[n]]))), pcaCoordinates, by = c("cellID")) %>% dplyr::select(-nLineages)
#   jointPCA1Barcodes = jointPCA %>% filter(SampleNum==paste0("S", twinAList[[n]])) %>% dplyr::select(barcode) %>% unique()
#   jointPCA2Barcodes = jointPCA %>% filter(SampleNum==paste0("S", twinBList[[n]]))  %>% dplyr::select(barcode) %>% unique()
#   
#   jointBarcodesOnlyBoth = inner_join(jointPCA2Barcodes,jointPCA1Barcodes, by = "barcode") 
#   jointBarcodesOnlyBoth = jointBarcodesOnlyBoth %>% mutate(barcodeName = c(1:nrow(jointBarcodesOnlyBoth)))
#   jointBarcodesOnlyBoth$barcodeName = sub("^", "B", jointBarcodesOnlyBoth$barcodeName)
#   jointBarcodesOnlyBothFilter = jointBarcodesOnlyBoth %>% dplyr::filter(barcodeName %in% barcodeList[[n]])
#   if(n == 4) {
#     exportTable <- jointBarcodesOnlyBothFilter
#   } else{
#     exportTable <- bind_rows(exportTable, jointBarcodesOnlyBothFilter)
#   }
# }
# 
# exportTable <- exportTable %>% mutate(mixing = matrixDMSOvsLSD1i$value) %>% group_by(barcode) %>% mutate(mean = mean(mixing))
# saveRDS(object = exportTable, file = paste0(homeDirectory, "twinBarcodesDMSOvsLSD1wMixing.rds"))
# 
# matrixCombined <- Reduce(bind_rows, matrixList)
# matrixCombined$ident <- factor(matrixCombined$ident, levels = c("DMSO self", "DMSO vs. LSD1i", "LSD1i self", "DMSO vs. DOT1Li", "DOT1Li self"))
# ggplot(matrixCombined %>% dplyr::filter(ident %in% c("DMSO self", "DMSO vs. LSD1i", "LSD1i self")), aes(x = ident, y = value, fill = ident), show.legend = FALSE) +
#   geom_boxplot(outlier.shape = NA, witdth = 0.5) +
#   geom_jitter(width = 0.25) +
#   stat_summary(fun = mean) +
#   stat_summary(fun.y = mean, geom = "text", aes(label = signif(..y.., 2)), vjust = 2) +
#   ylab("mixing coefficients per twin comparison") + theme(legend.position = "none", axis.title.x = element_blank())
# ggsave(file = paste0(plotDirectory, "FM3_mixCoeffBoxplots_DMSOvsLSD1i.pdf"), width = 3, height = 3)


################################################################################################################################################################################################################################################
#### UMAP plots for selected barcodes ####
################################################################################################################################################################################################################################################
umapCoordinates = as_tibble(read.table(file = paste0(homeDirectory, "umapCoordinates_Scanorama_50pcs_filter.tsv"), header = TRUE, stringsAsFactors = F, sep = "\t"))
umapClusters= as_tibble(read.table(file = paste0(homeDirectory, "umapClusters_Scanorama_50pcs_filter.tsv"), header = TRUE, stringsAsFactors = F, sep = "\t"))
umapCoordinates <- inner_join(umapCoordinates, umapClusters, by = c("cellID", "sampleNum")) %>% dplyr::rename(clusters = 5)

plotList2 <- list()
for(n in 1:length(twinAList)){
  if(length(twinAList[[n]] == 1)) {
    jointPCA = inner_join(linCountToOverlapsList %>% dplyr::filter(SampleNum %in% c(paste0("S", twinAList[[n]]), paste0("S", twinBList[[n]]))), umapCoordinates, by = c("cellID")) %>% dplyr::select(-nLineages)
    jointPCA1Barcodes = jointPCA %>% filter(SampleNum==paste0("S", twinAList[[n]])) %>% dplyr::select(barcode) %>% unique()
    jointPCA2Barcodes = jointPCA %>% filter(SampleNum==paste0("S", twinBList[[n]]))  %>% dplyr::select(barcode) %>% unique()
  } else {
    jointPCA = inner_join(linCountToOverlapsList %>% dplyr::filter(SampleNum %in% c(paste0("S", twinAList[[n]][1]), paste0("S", twinAList[[n]][2]), paste0("S", twinBList[[n]][1]), paste0("S", twinBList[[n]][2]))), umapCoordinates, by = c("cellID")) %>% dplyr::select(-nLineages)
    jointPCA1Barcodes = jointPCA %>% filter(SampleNum %in% c(paste0("S", twinAList[[n]][1]), paste0("S", twinAList[[n]][2]))) %>% dplyr::select(barcode) %>% unique()
    jointPCA2Barcodes = jointPCA %>% filter(SampleNum %in% c(paste0("S", twinBList[[n]][1]), paste0("S", twinBList[[n]][2]))) %>% dplyr::select(barcode) %>% unique()
  }
  
  jointBarcodesOnlyBoth = inner_join(jointPCA2Barcodes,jointPCA1Barcodes, by = "barcode") 
  jointBarcodesOnlyBoth = jointBarcodesOnlyBoth %>% mutate(barcodeName = c(1:nrow(jointBarcodesOnlyBoth)))
  jointBarcodesOnlyBoth$barcodeName = sub("^", "B", jointBarcodesOnlyBoth$barcodeName)
  
  jointBarcodesOnlyBothList = jointBarcodesOnlyBoth$barcode
  jointPCAOnlyBoth = jointPCA %>% filter(barcode %in% jointBarcodesOnlyBothList)
  finalPCAJoint = inner_join(jointPCAOnlyBoth,jointBarcodesOnlyBoth, by = "barcode") %>% dplyr::select(-cellID, -barcode, -nUMI)
  
  if(length(twinAList[[n]] == 1)) {
    finalPCAS1 = finalPCAJoint %>% filter(SampleNum==paste0("S", twinAList[[n]]))
    finalPCAS2 = finalPCAJoint %>% filter(SampleNum==paste0("S", twinBList[[n]]))
  } else{
    finalPCAS1 = finalPCAJoint %>% filter(SampleNum %in% c(paste0("S", twinAList[[n]][1]), paste0("S", twinAList[[n]][2])))
    finalPCAS2 = finalPCAJoint %>% filter(SampleNum %in% c(paste0("S", twinBList[[n]][1]), paste0("S", twinBList[[n]][2])))
  }
  
  finalPCAS1Big = finalPCAS1 %>% group_by(barcodeName) %>% summarise(nColony = length(barcodeName)) %>% filter(nColony > 3) %>% dplyr::select(-nColony)
  finalPCAS2Big = finalPCAS2 %>% group_by(barcodeName) %>% summarise(nColony = length(barcodeName)) %>% filter(nColony > 3) %>% dplyr::select(-nColony)
  
  commonS1S2Big = inner_join(finalPCAS1Big, finalPCAS2Big)
  commonS1S2Big = commonS1S2Big$barcodeName
  
  if(length(twinAList[[n]] == 1)) {
    finalPCAJointBig = finalPCAJoint %>% filter(barcodeName %in% commonS1S2Big) %>% dplyr::filter(SampleNum %in% c(paste0("S", twinAList[[n]]), paste0("S", twinBList[[n]])))
  } else{
    finalPCAJointBig = finalPCAJoint %>% filter(barcodeName %in% commonS1S2Big) %>% dplyr::filter(SampleNum %in% c(paste0("S", twinAList[[n]][1]), paste0("S", twinAList[[n]][2]), paste0("S", twinBList[[n]][1]), paste0("S", twinBList[[n]][2])))
  }
  
  for(i in 1:length(barcodeList[[n]])){
    finalUMAPJointPlot = finalPCAJointBig %>% filter(barcodeName == barcodeList[[n]][i])
    if(n == 3) {
      finalUMAPJointPlot$SampleNum <- ifelse(finalUMAPJointPlot$SampleNum %in% c("S1", "S2"), "DMSO", finalUMAPJointPlot$SampleNum)
      finalUMAPJointPlot$SampleNum <- ifelse(finalUMAPJointPlot$SampleNum %in% c("S5", "S6"), "LSD1", finalUMAPJointPlot$SampleNum)
    }
    plotList2[[length(plotList2) + 1]] <- ggplot() +
      rasterise(geom_point(data = umapCoordinates %>% dplyr::filter(UMAP_2 < 10, clusters %in% c("0", "1", "2", "3", "4", "5")), aes(x = UMAP_1, y = UMAP_2), color = "gray93"), dpi = 100)  +
      geom_point(data = finalUMAPJointPlot %>% dplyr::filter(UMAP_2 < 10, clusters %in% c("0", "1", "2", "3", "4", "5")), aes(x = UMAP_1, y = UMAP_2, color = SampleNum), size = 2.5, shape = 16) +
      scale_color_manual(values=c("hotpink3", "turquoise3")) +
      theme(legend.position = "none") + geom_text(aes(label = finalUMAPJointPlot$barcodeName %>% unique(), x = -Inf, y = Inf), hjust = 0, vjust = 1, size = 5) + NoAxes()
    ggsave(plotList2[[length(plotList2)]], file = paste0(plotDirectory, "neighborUMAPs/FM3_repUMAP_", n, "_", i, ".pdf"), width = 4, height = 4)
    }
}

################################################################################################################################################################################################################################################
#### do analysis for all pairwise twin comparisons #####
################################################################################################################################################################################################################################################
for(n in 1:length(fullMatrixList)) {
  matrixTemp <- fullMatrixList[[n]]
  rownames(matrixTemp) <- barcodeList[[n]]
  colnames(matrixTemp) <- barcodeList[[n]]
  matrixTemp[lower.tri(matrixTemp, diag = FALSE)] <- NA
  
  matrixTempMelt <- melt(matrixTemp)
  matrixTempMelt$value <- round(matrixTempMelt$value, 2)
  
  ggplot(data = matrixTempMelt, aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "white", high = "black", limit = c(0, max(c(1, matrixTempMelt$value))), na.value = "white", breaks = c(0, 0.5, 1)) +
    theme_minimal() + 
    coord_fixed() +
    geom_text(aes(Var2, Var1, label = value, color = value), size = 5.5) +
    scale_color_gradient(low = "black", high = "black", limits = c(0, 0.7), na.value = "white") +
    theme_classic((base_size = 26)) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.line = element_line(colour = 'black', size = 1.5),
      axis.text.x = element_text(angle = 90),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_line(colour = "black", size = 1.5),
      legend.position = "none",
      text=element_text(family="Helvetica")) +
    guides(fill = guide_colorbar(barwidth = 7, barheight = 1, title.position = "top", title.hjust = 0.5)) +
    scale_y_discrete(position = "right")
  ggsave(file = paste0(plotDirectory, "neighborMatrices/FM3_mixCoeffCorrPlotAll_", n, ".pdf"), width = 6, height = 6)
  
  mixingCoeffTwin <- tibble(label = "twin", value = diag(matrixTemp))
  mixingCoeffNonTwin <- tibble(label = "non-twin", value = matrixTemp[upper.tri(matrixTemp, diag = FALSE)])
  mixingAll <- bind_rows(mixingCoeffTwin, mixingCoeffNonTwin)
  
  ggplot(mixingAll, aes(label, y = value)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.25) + ylim(0, 1.5) +
    geom_signif(comparisons = list(c("non-twin", "twin"))) +
    ylab("mixing coefficient") +
    theme(axis.title.x = element_blank())
  ggsave(file = paste0(plotDirectory, "neighborMatrices/FM3_twinVersusNonTwinBoxplot_", n, ".pdf"), width = 6, height = 6)
}

################################################################################################################################################################################################################################################
#### determine examples of moving from clusters ####
################################################################################################################################################################################################################################################
scanorama_filter_subset <- readRDS(file = paste0(homeDirectory, "scanorama_filter_subset.rds"))
DimPlot(object = scanorama_filter_subset)

linCountToOverlapsList <- as_tibble(read.table(file = paste0(homeDirectory, "filtered10XCells.txt"), header = TRUE, stringsAsFactors = F, sep = "\t"))

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
barcodeList <- unique(finalPCAJointBig$barcodeName)

for(i in 1:length(barcodeList[[n]])){
  finalUMAPJointPlot = finalPCAJointBig %>% filter(barcodeName == barcodeList[[n]][i])
  ggplot() +
    rasterise(geom_point(data = umapCoordinates, aes(x = UMAP_1, y = UMAP_2), color = "gray93"), dpi = 100)  +
    geom_point(data = finalUMAPJointPlot, aes(x = UMAP_1, y = UMAP_2, color = SampleNum), size = 2.5, shape = 16) +
    scale_color_manual(values=c("hotpink3", "hotpink3", "turquoise3", "turquoise3")) +
    theme(legend.position = "none") + NoAxes() + geom_text(aes(label = finalUMAPJointPlot$barcodeName %>% unique(), x = -Inf, y = Inf), hjust = 0, vjust = 1, size = 5)
  ggsave(file = paste0(plotDirectory, "neighborUMAPsTransition/FM3_repUMAP_", barcodeList[i], ".pdf"), width = 6, height = 6)
}

barcodePlotList <- c("B1", "B3", "B4", "B5", "B6", "B7", "B9", "B10", "B13", "B14", "B20", "B23")
barcodePlotList <- c("B1", "B6", "B7", "B9", "B10", "B13", "B14", "B20", "B23") #consider removing B3, B4, and B5

switchUMAPJointPlot <- finalPCAJointBig %>% dplyr::filter(barcodeName %in% barcodePlotList) %>% dplyr::filter(UMAP_1 < 6, UMAP_2 < 10)

ggplot() +
  rasterise(geom_point(data = umapCoordinates %>% dplyr::filter(UMAP_1 < 6, UMAP_2 < 10), aes(x = UMAP_1, y = UMAP_2), color = "gray93"), dpi = 100)  +
  geom_point(data = switchUMAPJointPlot, aes(x = UMAP_1, y = UMAP_2, color = label), size = 5, shape = 16) +
  # geom_density_2d(data = switchUMAPJointPlot, aes(x = UMAP_1, y = UMAP_2, color = label), bins = 15) +
  scale_color_manual(values=c("hotpink3", "turquoise3")) +
  theme(legend.position = "none") + NoAxes()
