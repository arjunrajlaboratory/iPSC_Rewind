rm(list=ls())
gc()

library(Seurat)
library(tidyverse)
library(reshape2)
library(ggrepel)
library(ggsignif)
library(egg)
library(pheatmap)
library(monaLisa)
library(JASPAR2020)
library(httr)
library(jsonlite)

options(future.globals.maxSize = 4000 * 1024 ^ 2)
theme_set(theme_classic())

homeDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/extractedData/rewind10X/homerAnalysis/"
plotDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/plots/rewind10X/homerAnalysis/"

markersComb <- readRDS(file = paste0(homeDirectory, "markersComb.rds"))

#### output marker lists for use in ChEA3 for regulators ##############################################################################################
library(mygene)
markers <- markersComb
markersOutput <- markers %>% dplyr::select(gene, avg_log2FC) %>% filter(avg_log2FC < 0)
markersOutput$avg_log2FC <- abs(markersOutput$avg_log2FC)
markersOutput <- markersOutput %>% dplyr::select(gene)
write.table(markersOutput, file = "/Users/naveenjain/Downloads/markers.txt", quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)

test <- queryMany(markersOutput$gene, scopes="symbol", fields="entrezgene", species="human")
list <- test$entrezgene[!(is.na(test$entrezgene))]
write.table(list, file = "/Users/naveenjain/Downloads/HOMERmarkers.txt", quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)

### include R based ChEA3 code here ##################################################################################################################

# for negative markers
genes = markersComb %>% slice_min(., order_by = mean, n = 50) %>% .$gene

url = "https://maayanlab.cloud/chea3/api/enrich/"
encode = "json"
payload = list(query_name = "myQuery", gene_set = genes)

response = POST(url = url, body = payload, encode = encode)
json = content(response, "text")
results = fromJSON(json) %>% .[["Integrated--meanRank"]]

resultsTop <- results %>% dplyr::slice_min(., order_by = as.numeric(Score), n = 12)

ggplot(resultsTop, aes(y = reorder(TF, -as.numeric(Score)), x = as.numeric(Score))) +
  geom_col() + xlab("mean rank") + theme(axis.title.y = element_blank())
ggsave(filename = paste0(plotDirectory, "rankOrderTF_negMarkers.pdf"), height = 4, width = 2)

genes = markersComb %>% slice_min(., order_by = mean, n = 25) %>% .$gene
for(i in 1:nrow(resultsTop)) {
  resultsTopTemp <- resultsTop[i, ]
  overlappingGenes <- resultsTopTemp %>% .$Overlapping_Genes %>% str_split(., pattern = ",") %>% unlist() %>% intersect(., genes)
  overlappingTFTemp <- data.frame(gene = overlappingGenes, TF = resultsTopTemp$TF)
  if(i == 1) {
    overlappingTF <- overlappingTFTemp
  } else {
    overlappingTF <- bind_rows(overlappingTF, overlappingTFTemp)
  }
}
ggplot(overlappingTF, aes(x = TF, y = gene)) +
  geom_tile(aes(fill = TF))

resultsTop <- results %>% dplyr::slice_min(., order_by = as.numeric(Score), n = 200)
genesAssociated <- results %>% .$Overlapping_Genes %>% str_split(., pattern = ",") %>% unlist() %>% table(.) %>% as.data.frame()

genes = markersComb %>% slice_min(., order_by = mean, n = 25) %>% dplyr::filter(gene %in% genesAssociated$.)
geneMatrix <- matrix(nrow = length(genes$gene), ncol = length(genes$gene), data = NA)
for(i in 1:length(genes$gene)) {
  geneA <- genes$gene[i]
  for(j in 1:length(genes$gene)) {
    geneB <- genes$gene[j]
    counter <- 0
    for(n in 1:nrow(resultsTop)) {
      resultsTFTemp <- resultsTop[n, ] %>% .$Overlapping_Genes %>% str_split(., pattern = ",") %>% unlist()
      if(geneA %in% resultsTFTemp) {
        if(geneB %in% resultsTFTemp) {
          counter <- counter + 1
        }
      }
    }
    geneMatrix[i, j] <- counter
  }
}

geneMatrixNorm <- geneMatrix
for(j in 1:length(genes$gene)) {
  for(i in j:length(genes$gene)) {
    geneMatrixNorm[i, j] <- geneMatrix[i, j] / geneMatrix[j, j]
  }
}
for(i in 1:length(genes$gene)) {
  for(j in i:length(genes$gene)) {
    geneMatrixNorm[i, j] <- geneMatrix[i, j] / geneMatrix[i, i]
  }
}

rownames(geneMatrixNorm) <- genes$gene
colnames(geneMatrixNorm) <- genes$gene
pheatmap(geneMatrixNorm)
plot <- pheatmap(geneMatrixNorm, treeheight_row = 0, treeheight_col = 0, legend = FALSE)
ggsave(plot, filename = paste0(plotDirectory, "overlapTF_negMarkers.pdf"), height = 5, width = 5)

# for positive markers
genes = markersComb %>% slice_max(., order_by = mean, n = 50) %>% .$gene

url = "https://maayanlab.cloud/chea3/api/enrich/"
encode = "json"
payload = list(query_name = "myQuery", gene_set = genes)

response = POST(url = url, body = payload, encode = encode)
json = content(response, "text")
results = fromJSON(json) %>% .[["Integrated--meanRank"]]

resultsTop <- results %>% dplyr::slice_min(., order_by = as.numeric(Score), n = 12)

ggplot(resultsTop, aes(y = reorder(TF, -as.numeric(Score)), x = as.numeric(Score))) +
  geom_col() + xlab("mean rank") + theme(axis.title.y = element_blank())
ggsave(filename = paste0(plotDirectory, "rankOrderTF_posMarkers.pdf"), height = 4, width = 2)

genes = markersComb %>% slice_max(., order_by = mean, n = 25) %>% .$gene
for(i in 1:nrow(resultsTop)) {
  resultsTopTemp <- resultsTop[i, ]
  overlappingGenes <- resultsTopTemp %>% .$Overlapping_Genes %>% str_split(., pattern = ",") %>% unlist() %>% intersect(., genes)
  overlappingTFTemp <- data.frame(gene = overlappingGenes, TF = resultsTopTemp$TF)
  if(i == 1) {
    overlappingTF <- overlappingTFTemp
  } else {
    overlappingTF <- bind_rows(overlappingTF, overlappingTFTemp)
  }
}
ggplot(overlappingTF, aes(x = TF, y = gene)) +
  geom_tile(aes(fill = TF))

resultsTop <- results %>% dplyr::slice_min(., order_by = as.numeric(Score), n = 200)
genesAssociated <- results %>% .$Overlapping_Genes %>% str_split(., pattern = ",") %>% unlist() %>% table(.) %>% as.data.frame()

genes = markersComb %>% slice_max(., order_by = mean, n = 25) %>% dplyr::filter(gene %in% genesAssociated$.)
geneMatrix <- matrix(nrow = length(genes$gene), ncol = length(genes$gene), data = NA)
for(i in 1:length(genes$gene)) {
  geneA <- genes$gene[i]
  for(j in 1:length(genes$gene)) {
    geneB <- genes$gene[j]
    counter <- 0
    for(n in 1:nrow(resultsTop)) {
      resultsTFTemp <- resultsTop[n, ] %>% .$Overlapping_Genes %>% str_split(., pattern = ",") %>% unlist()
      if(geneA %in% resultsTFTemp) {
        if(geneB %in% resultsTFTemp) {
          counter <- counter + 1
        }
      }
    }
    geneMatrix[i, j] <- counter
  }
}

geneMatrixNorm <- geneMatrix
for(j in 1:length(genes$gene)) {
  for(i in j:length(genes$gene)) {
    geneMatrixNorm[i, j] <- geneMatrix[i, j] / geneMatrix[j, j]
  }
}
for(i in 1:length(genes$gene)) {
  for(j in i:length(genes$gene)) {
    geneMatrixNorm[i, j] <- geneMatrix[i, j] / geneMatrix[i, i]
  }
}

rownames(geneMatrixNorm) <- genes$gene
colnames(geneMatrixNorm) <- genes$gene
pheatmap(geneMatrixNorm)
plot <- pheatmap(geneMatrixNorm, treeheight_row = 0, treeheight_col = 0, legend = FALSE)
ggsave(plot, filename = paste0(plotDirectory, "overlapTF_posMarkers.pdf"), height = 5, width = 5)

### make a correlation matrix for the top positive and negative markers ##############################################################################
markersPos <- markers %>% slice_max(., order_by = avg_log2FC, n = 50) %>% .$gene
markersNeg <- markers %>% slice_min(., order_by = avg_log2FC, n = 50) %>% .$gene
markersCorrTable <- logNormalizedCounts %>% dplyr::select(all_of(markersPos), all_of(markersNeg))
res <- rcorr(as.matrix(markersCorrTable))
corrplot(res$r, method = "color", order = "hclust", tl.col = "black", tl.srt = 45, p.mat = res$p, sig.level = 0.01, addrect = 6, tl.cex = 0.5, cl.cex = 0.5)

markersPos <- markers %>% slice_max(., order_by = avg_log2FC, n = 100) %>% .$gene
markersNeg <- markers %>% slice_min(., order_by = avg_log2FC, n = 100) %>% .$gene

markersCorrTable <- logNormalizedCounts %>% dplyr::select(all_of(markersPos))
res <- rcorr(as.matrix(markersCorrTable))
corrplot(res$r, method = "color", order = "hclust", tl.col = "black", tl.srt = 45, p.mat = res$p, sig.level = 0.01, addrect = 6, tl.cex = 0.5, cl.cex = 0.5)
#fviz_nbclust(markersCorrTable, kmeans, method = "wss", k.max = 20) + theme_minimal()

markersCorrTable <- logNormalizedCounts %>% dplyr::select(all_of(markersNeg))
res <- rcorr(as.matrix(markersCorrTable))
corrplot(res$r, method = "color", order = "hclust", tl.col = "black", tl.srt = 45, p.mat = res$p, sig.level = 0.01, addrect = 6, tl.cex = 0.5, cl.cex = 0.5)
#fviz_nbclust(markersCorrTable, kmeans, method = "wss", k.max = 20) + theme_minimal()

### reorganize HOMER outputs to make .bed files in IGV ###############################################################################################
homerSNAI2 <- read_tsv(file = paste0(homeDirectory, "SNAI2.txt"))
homerTWIST2 <- read_tsv(file = paste0(homeDirectory, "TWIST2.txt"))
homerOSR1 <- read_tsv(file = paste0(homeDirectory, "OSR1.txt"))
homerPRRX2 <- read_tsv(file = paste0(homeDirectory, "PRXX2.txt"))

homerList <- list(homerSNAI2, homerTWIST2, homerOSR1, homerPRRX2)
nameList <- c("SNAI2", "TWIST2", "OSR1", "PRRX2")

homerOutputFinalCombinedList <- list()
for(j in 1:length(homerList)){
  homerOutput <- homerList[[j]]
  homerOutput <- homerOutput %>% dplyr::rename(PeakID = 1, Info = 22)
  homerOutput <- homerOutput %>% mutate(PeakLoc = 0)
  homerOutput <- homerOutput %>% mutate(PeakStrand = Strand)
  
  for(i in 1:nrow(homerOutput)){
    homerOutputTemp <- homerOutput[i, ]
    peakLocList <- gsub("\\(.*?\\)", "", homerOutputTemp$Info) %>% strsplit(., ",")
    peakStrandList <- gsub(".*\\(.*\\,(.*)\\,.*\\).*?", "\\1", homerOutputTemp$Info) %>% strsplit(., "")
    
    if(!is.na(unlist(peakLocList))){
      homerOutputTemp <- homerOutputTemp %>% dplyr::slice(rep(1:n(), each = length(unlist(peakLocList))))
      homerOutputTemp$PeakLoc <- as.integer(unlist(peakLocList))
      homerOutputTemp$PeakStrand <- unlist(peakStrandList)
      
      if(i == 1){
        homerOutputReformat <- homerOutputTemp
      }
      else{
        homerOutputReformat <- bind_rows(homerOutputReformat, homerOutputTemp)
      }
    }
  }
  
  homerOutputReformat <- homerOutputReformat %>% mutate(PeakCenter = Start + (End - Start)/2 + PeakLoc) %>%
    mutate(PeakStart = PeakCenter - 6) %>%
    mutate(PeakEnd = PeakCenter + 6) %>%
    mutate(Score = 1000) %>%
    mutate(Name = nameList[j])
  
  homerOutputFinal <- homerOutputReformat %>% dplyr::select(Chr, PeakStart, PeakEnd, Name, Score, PeakStrand)
  write_tsv(homerOutputFinal, file = paste0("/Users/naveenjain/Downloads/", nameList[j], ".bed"))
  homerOutputFinalCombinedList[[j]] <- homerOutputFinal
}

homerOutputFinalCombined <- homerOutputFinalCombinedList[[1]]
for(k in 2:length(homerOutputFinalCombinedList)){
  homerOutputFinalCombined <- bind_rows(homerOutputFinalCombined, homerOutputFinalCombinedList[[k]])
}
write_tsv(homerOutputFinalCombined, file = "/Users/naveenjain/Downloads/negativeMarkersTFs.bed")

BiocManager::install("monaLisa")
dumpJaspar(filename = paste0("/Users/naveenjain/Downloads/", "PRRX2.motif"), pkg = "JASPAR2020", opts = list(ID = c("MA0075.3")))
