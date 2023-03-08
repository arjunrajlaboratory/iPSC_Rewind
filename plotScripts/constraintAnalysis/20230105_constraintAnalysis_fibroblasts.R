rm(list=ls())
gc()

library(tidyverse)
library(Seurat)
library(ggsignif)
library(DescTools)
library(RANN)

theme_set(theme_classic())

plotDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/plots/constraintAnalysis/"
homeDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/extractedData/rewind10X/R1/"

#### MIXING CALCULATION FOR ALL CLONES ##############################################################################################################################################
pcaTable <- as_tibble(read.table(file = paste0(homeDirectory, "pcaCoordinates_Scanorama_50pcs_filter.tsv"), header = TRUE, stringsAsFactors = F, sep = "\t"))
filteredBarcodeTable <- as_tibble(read.table(file = paste0(homeDirectory, "filtered10XCells.txt"), header = TRUE, stringsAsFactors = F, sep = "\t"))
jointPCA = inner_join(filteredBarcodeTable, pcaTable, by = c("cellID")) %>% dplyr::select(-SampleNum, -nUMI, -fracUMI, -nLineages)

nCellsPerBarcode <- jointPCA %>% group_by(BC50StarcodeD8) %>% summarise(nCells = n())
nCellsPerBarcodeFilter <- nCellsPerBarcode %>% dplyr::filter(nCells >= 4) %>% .$BC50StarcodeD8

barcodeList <- nCellsPerBarcodeFilter
nReplicates <- 10

for(i in c(1:length(barcodeList))) {
  neighborList = c()
  iBarcodeBx <- jointPCA %>% dplyr::filter(BC50StarcodeD8 == barcodeList[i]) %>% mutate(label = "query")
  for(n in 1:nReplicates) {
    jBarcodeBx <- jointPCA %>% dplyr::filter(BC50StarcodeD8 %in% barcodeList) %>% sample_n(size = nrow(iBarcodeBx)) %>% mutate(label = "rando")
    ijBarcodeBx <- bind_rows(iBarcodeBx, jBarcodeBx)
    BarcodeRef <- ijBarcodeBx %>% mutate(num = 1:nrow(ijBarcodeBx))
    S1index = BarcodeRef %>% dplyr::filter(label == "query") %>% dplyr::select(num)
    S2index = BarcodeRef %>% dplyr::filter(label == "rando") %>% dplyr::select(num)
    
    knnPCA = nn2(ijBarcodeBx[, 3:52], ijBarcodeBx[, 3:52])
    neighborKNN = as_tibble(knnPCA[[1]])
    nS1 = c();
    nS2 = c();
    fraction2 = (nrow(S2index)/nrow(S1index))
    for (k in c(1:nrow(ijBarcodeBx))) {
      nS1[k] =  (sum(neighborKNN[k, 2:ncol(neighborKNN)] %in% S1index$num))*(nrow(S2index)/nrow(S1index))
      nS2[k] =  sum(neighborKNN[k, 2:ncol(neighborKNN)] %in% S2index$num)
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
    neighborList[n] = mean(withOther$fraction)/mean(withSelf$fraction)
  }
  neighborListAllTemp <- tibble(barcode = barcodeList[i], mean = mean(neighborList), sd = sd(neighborList), num = nrow(iBarcodeBx))
  if(i == 1) {
    neighborListAll <-neighborListAllTemp
  } else{
    neighborListAll <- bind_rows(neighborListAll, neighborListAllTemp)
  }
}

ggplot(neighborListAll, aes(x = "", y = mean)) +
  geom_jitter() + ylim(0, 1.5)

sequencedBarcodeTable <- as_tibble(read.table(paste0(homeDirectory, "stepThreeStarcodeShavedReads_BC_10XAndGDNA.txt"), header = TRUE, stringsAsFactors = F, sep = '\t'))
probedBarcodeTable <- sequencedBarcodeTable %>% dplyr::filter(cellID == "dummy") %>% dplyr::select(UMI, BC50StarcodeD8, SampleNum)
probedBarcodeTable$UMI <- as.numeric(probedBarcodeTable$UMI)
probedBarcodeTable <- probedBarcodeTable %>% group_by(BC50StarcodeD8, SampleNum) %>% summarise(nUMI = sum(UMI))
primedBarcodeTable <- probedBarcodeTable %>% ungroup() %>% slice_max(., order_by = nUMI, n = 100)
primedCellTable <- inner_join(primedBarcodeTable, filteredBarcodeTable, by = "BC50StarcodeD8") %>% mutate(label = "primed")

nReplicates <- 100
neighborList = c()
iBarcodeBx <- jointPCA %>% dplyr::filter(BC50StarcodeD8 %in% primedCellTable$BC50StarcodeD8) %>% mutate(label = "query")
for(n in 1:nReplicates) {
  jBarcodeBx <- jointPCA %>% sample_n(size = nrow(iBarcodeBx), seed = n) %>% mutate(label = "rando")
  ijBarcodeBx <- bind_rows(iBarcodeBx, jBarcodeBx)
  BarcodeRef <- ijBarcodeBx %>% mutate(num = 1:nrow(ijBarcodeBx))
  S1index = BarcodeRef %>% dplyr::filter(label == "query") %>% dplyr::select(num)
  S2index = BarcodeRef %>% dplyr::filter(label == "rando") %>% dplyr::select(num)
  
  knnPCA = nn2(ijBarcodeBx[, 3:52], ijBarcodeBx[, 3:52])
  neighborKNN = as_tibble(knnPCA[[1]])
  nS1 = c();
  nS2 = c();
  fraction2 = (nrow(S2index)/nrow(S1index))
  for (k in c(1:nrow(ijBarcodeBx))) {
    nS1[k] =  (sum(neighborKNN[k, 2:ncol(neighborKNN)] %in% S1index$num))*(nrow(S2index)/nrow(S1index))
    nS2[k] =  sum(neighborKNN[k, 2:ncol(neighborKNN)] %in% S2index$num)
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
  neighborList[n] = mean(withOther$fraction)/mean(withSelf$fraction)
}
neighborListQuery <- as_tibble(data.frame(mixingCoeff = neighborList, label = "query"))

ggplot(neighborListQuery, aes(x = "", y = mixingCoeff)) +
  geom_jitter() + ylim(0, 1.5)

nReplicates <- 100
neighborList = c()
for(n in 1:nReplicates) {
  iBarcodeBx <- jointPCA %>% sample_n(size = nrow(jointPCA %>% dplyr::filter(BC50StarcodeD8 %in% primedCellTable$BC50StarcodeD8)), seed = n) %>% mutate(label = "query")
  jBarcodeBx <- jointPCA %>% sample_n(size = nrow(jointPCA %>% dplyr::filter(BC50StarcodeD8 %in% primedCellTable$BC50StarcodeD8)), seed = n+1) %>% mutate(label = "rando")
  ijBarcodeBx <- bind_rows(iBarcodeBx, jBarcodeBx)
  BarcodeRef <- ijBarcodeBx %>% mutate(num = 1:nrow(ijBarcodeBx))
  S1index = BarcodeRef %>% dplyr::filter(label == "query") %>% dplyr::select(num)
  S2index = BarcodeRef %>% dplyr::filter(label == "rando") %>% dplyr::select(num)
  
  knnPCA = nn2(ijBarcodeBx[, 3:52], ijBarcodeBx[, 3:52])
  neighborKNN = as_tibble(knnPCA[[1]])
  nS1 = c();
  nS2 = c();
  fraction2 = (nrow(S2index)/nrow(S1index))
  for (k in c(1:nrow(ijBarcodeBx))) {
    nS1[k] =  (sum(neighborKNN[k, 2:ncol(neighborKNN)] %in% S1index$num))*(nrow(S2index)/nrow(S1index))
    nS2[k] =  sum(neighborKNN[k, 2:ncol(neighborKNN)] %in% S2index$num)
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
  neighborList[n] = mean(withOther$fraction)/mean(withSelf$fraction)
}
neighborListRando <- as_tibble(data.frame(mixingCoeff = neighborList, label = "rando"))

ggplot(neighborListRando, aes(x = "", y = mixingCoeff)) +
  geom_jitter() + ylim(0, 1.5)

neighborListAll <- bind_rows(neighborListQuery, neighborListRando)
ggplot(neighborListAll, aes(x = label, y = mixingCoeff)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.25) + ylim(0, 1.5) +
  geom_signif(comparisons = list(c("query", "rando")))
ggsave(filename = paste0(plotDirectory, "mixingCoeffPrimedVsRandomFibroblasts.pdf"), height = 3, width = 3, useDingbats = FALSE)

#### JSD CALCULATION FOR ALL CLONES #################################################################################################################################################
umapTable <- as_tibble(read.table(file = paste0(homeDirectory, "umapCoordinates_Scanorama_50pcs_filter.tsv"), header = TRUE, stringsAsFactors = F, sep = "\t"))
clusterTable <- as_tibble(read.table(file = paste0(homeDirectory, "umapClusters_Scanorama_50pcs_filter.tsv"), header = TRUE, stringsAsFactors = F, sep = "\t"))
filteredBarcodeTable <- as_tibble(read.table(file = paste0(homeDirectory, "filtered10XCells.txt"), header = TRUE, stringsAsFactors = F, sep = "\t"))

jointUMAP = inner_join(umapTable, clusterTable, by = c("cellID", "sampleNum")) %>% dplyr::rename(cluster = integrated_snn_res.0.45) %>% 
  inner_join(filteredBarcodeTable, ., by = c("cellID")) %>% dplyr::select(-SampleNum, -nUMI, -fracUMI, -nLineages)
jointUMAPFilter <- jointUMAP %>% dplyr::filter(!(cluster %in% c("4", "9", "10"))) %>% dplyr::filter(UMAP_1 > -5)

nCellsPerBarcode <- jointUMAP %>% group_by(BC50StarcodeD8) %>% summarise(nCells = n())
nCellsPerBarcodeFilter <- nCellsPerBarcode %>% dplyr::filter(nCells >= 4) %>% .$BC50StarcodeD8

KLD <- function(x,y) sum(x * log2(x/y))
JSD <- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
mixedrank = function(x) order(gtools::mixedorder(x))

sibACells <- jointUMAP %>% dplyr::filter(!(cluster %in% c("4", "8", "9", "10"))) %>% left_join(., nCellsPerBarcode, by = "BC50StarcodeD8") %>% dplyr::filter(BC50StarcodeD8 %in% nCellsPerBarcodeFilter)
sibACellsDistTable <- sibACells %>% group_by(BC50StarcodeD8, cluster, nCells) %>% summarise(cellsPerCluster = n()) %>% ungroup() %>%
  complete(nesting(BC50StarcodeD8, nCells), cluster) %>% replace(., is.na(.), 0) %>%
  group_by(cluster) %>% mutate(allBarPerClust = sum(cellsPerCluster)) %>% mutate(propBarPerClust = cellsPerCluster / allBarPerClust) %>% ungroup() %>%
  group_by(BC50StarcodeD8) %>% mutate(prop = propBarPerClust / sum(propBarPerClust)) %>% arrange(mixedrank(BC50StarcodeD8))

barPerClustA <- sibACellsDistTable %>% ungroup() %>% dplyr::select(cluster, allBarPerClust) %>% unique()

n = 1000
seeds = 1:n
JSTest = rep(0, n)
JS_A = rep(0, length(nCellsPerBarcodeFilter))
meanJS_A = rep(0, length(nCellsPerBarcodeFilter))
sdJS_A = rep(0, length(nCellsPerBarcodeFilter))
setwd(paste0(plotDirectory, "R1_cloneBarcodes/"))

for (i in 1:length(nCellsPerBarcodeFilter)) {
  barcode = nCellsPerBarcodeFilter[i]
  bartitle = paste0("B", i)
  barcodeclusters = sibACellsDistTable %>%
    filter(BC50StarcodeD8 == barcode)
  allrandomprop = barPerClustA
  #generate random samples
  for (j in seeds) {
    set.seed(j)
    randomsibArows = as.integer(sample(rownames(sibACells), barcodeclusters$nCells[1]))
    randomsibA = suppressMessages(
      sibACells[randomsibArows, ] %>%
        dplyr::select(cluster) %>%
        group_by(cluster) %>%
        summarise(randomcellspercluster = n()) %>%
        full_join(barPerClustA, by = c('cluster')) %>%
        mutate_if(is.numeric, ~ replace(., is.na(.), 0)) %>%
        mutate(randomprop = randomcellspercluster /
                 allBarPerClust) %>%
        mutate(randomprop = randomprop / sum(randomprop)) %>%
        arrange(cluster)
    )
    allrandomprop[, paste0("randomprop", j)] = randomsibA$randomprop
  }
  #for entropy analysis
  allrandomprop = allrandomprop %>% rowwise() %>%
    mutate(avgrandomprop = mean(c_across(starts_with("randomprop"))))
  summaryrandomprop = allrandomprop %>%
    inner_join(barcodeclusters, by = 'cluster') %>%
    dplyr::select(cluster, prop, avgrandomprop)
  #calculate entropy values
  p_observed = summaryrandomprop %>% ungroup() %>%
    arrange(cluster) %>%
    summarise(prob_observed = prop + 0.000001) #add pseudocount to avoid zero in the numerator and/or denominator of KL
  p_random = summaryrandomprop %>% ungroup() %>%
    arrange(cluster) %>%
    summarise(prob_random = avgrandomprop + 0.000001)
  
  observedentropy = Entropy(p_observed, base = 2)
  randomentropy = Entropy(p_random, base = 2)
  JS_A[i] <- JSD(p_random, p_observed)
  
  plotsummaryrandomprop = summaryrandomprop %>% transmute(cluster, observed_prop = prop, randomavg_prop = avgrandomprop) %>%
    pivot_longer(cols = -cluster, names_to = c("sample", ".value"), names_pattern = "(.+)_(.+)")
  sample.labs <- c(paste0("observed ", bartitle, ": entropy = ", round(observedentropy, 3)),
      paste0("random sample: entropy = ", round(randomentropy, 3)))
  names(sample.labs) <- c("observed", "randomavg")
  entropyplot_prop = ggplot(plotsummaryrandomprop, aes(x = as.factor(cluster), y = prop, fill = sample)) +
    geom_bar(stat = 'identity') +
    facet_wrap( ~ sample, ncol = 1, labeller = labeller(sample = sample.labs)) +
    scale_fill_manual(breaks = c("observed", "randomavg"), values = c("#009292", "gray80")) +
    ggtitle(paste(bartitle, ': ', barcodeclusters$nCells[1], ' cells\n', 'Jensen-Shannon distance = ', round(JS_A[i], 3), sep = "")) +
    xlab("Seurat clusters, res=0.5") +
    ylab(paste0("normalized cell proportion")) +
    theme_classic() +
    theme(legend.position = "none")
  pdf(paste(bartitle, "entropyPlot.pdf", sep="_"), width = 4, height = 4)
  print(entropyplot_prop)
  dev.off()
  
  #evaluate JSD significance
  for (k in seeds + n) {
    set.seed(k)
    JSDsibArows = as.integer(sample(rownames(sibACells), barcodeclusters$nCells[1]))
    JSDsibA = sibACells[JSDsibArows, ] %>%
      dplyr::select(cluster) %>%
      full_join(barPerClustA, by = 'cluster') %>%
      group_by(cluster, allBarPerClust) %>%
      summarise(cellspercluster = n(), .groups = 'drop_last') %>%
      mutate(propperbarclust = cellspercluster / allBarPerClust) %>%
      ungroup() %>%
      mutate(prop = propperbarclust / sum(propperbarclust))
    p_test = JSDsibA %>% ungroup() %>%
      arrange(cluster) %>%
      summarise(prob_test = prop + 0.000001)
    JSTest[k - n] <- JSD(p_random, p_test)
  }
  meanJS_A[i] = mean(JSTest)
  sdJS_A[i] = sd(JSTest)
  normp = pnorm(JS_A[i], mean = meanJS_A[i], sd = sdJS_A[i], lower.tail = FALSE)
  if (round(normp, 3) == 0) {
    normp3 <- "< 0.001"
  } else{
    normp3 <- paste("=", round(normp, 3))
  }
  JSplot = ggplot() +
    geom_histogram(aes(x = JSTest), binwidth = 0.01, fill = 'gray80') +
    geom_vline(xintercept = JS_A[i], color = '#009292') +
    ggtitle(paste0(bartitle, ': ', barcodeclusters$nCells[1], ' cells', '\np ', normp3)) +
    xlab("Jensen-Shannon distance") +
    coord_cartesian(xlim = c(0, 1)) + #without coord_cartesian throws an error about rows being removed
    theme_classic()
  pdf(paste(bartitle, "JSDPlot.pdf", sep="_"), width = 4, height = 4)
  print(JSplot)
  dev.off()
  
  #walk through different barcode in split A only
  umap = ggplot() +
    geom_point(data = jointUMAPFilter, aes(x = UMAP_1, y = UMAP_2, color = "other")) +
    geom_point(data = sibACells[which((sibACells$BC50StarcodeD8 == barcode)), ], aes(x = UMAP_1, y = UMAP_2, color = "barcode"), size = 3) +
    scale_color_manual(breaks = c("barcode", "other"), values = c("#009292", "gray80"), labels = c(paste0(barcode), "other")) +
    theme(legend.position = "none", axis.text = element_blank(), axis.ticks = element_blank()) +
    labs(color =  'Split A barcodes\nwith at least 20 cells') +
    ggtitle(paste0(bartitle))
  pdf(paste(bartitle, "umapAlone.pdf", sep="_"), width = 4, height = 4)
  print(umap)
  dev.off()
}
#save JSD values to create table
nCellsPerBarcodeFilterTable <- nCellsPerBarcode %>% dplyr::filter(nCells >= 4)

JSTable_A = tibble(JS_A, meanJS_A, sdJS_A, barNum = paste0("B", rownames(nCellsPerBarcodeFilterTable)), barcode = nCellsPerBarcodeFilterTable$BC50StarcodeD8, nCells = nCellsPerBarcodeFilterTable$nCells)
write_csv(JSTable_A, file = paste0("JSTableInfo.csv"))

ggplot(JSTable_A, aes(x = JS_A, y = meanJS_A)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  xlim(0, 1) + ylim(0, 1)

#### JSD CALCULATION FOR ALL PRIMED CELLS ###########################################################################################################################################
nCellsPerBarcodeFilter <- primedCellTable %>% .$BC50StarcodeD8

clusterList <- as_tibble(data.frame(cluster = c("0", "1", "2", "3", "4", "5", "6", "7")))

sibACells <- jointUMAP %>% dplyr::filter(!(cluster %in% c("4", "9", "10"))) %>% dplyr::filter(BC50StarcodeD8 %in% nCellsPerBarcodeFilter)
barPerClustA <- jointUMAP %>% dplyr::filter(!(cluster %in% c("4", "9", "10"))) %>% dplyr::select(cluster) %>% group_by(cluster) %>% summarise(allBarPerClust = n())
sibACellsDistTable <- sibACells %>% group_by(cluster) %>% summarise(cellsPerCluster = n()) %>% ungroup() %>% right_join(., barPerClustA, by = "cluster") %>%
  replace(., is.na(.), 0) %>%
  mutate(propBarPerClust = cellsPerCluster / allBarPerClust) %>% ungroup() %>%
  mutate(prop = propBarPerClust / sum(propBarPerClust))

n = 1000
seeds = 1:n
JSTest = rep(0, n)
JS_A = rep(0, length(nCellsPerBarcodeFilter))
meanJS_A = rep(0, length(nCellsPerBarcodeFilter))
sdJS_A = rep(0, length(nCellsPerBarcodeFilter))
setwd(paste0(plotDirectory, "R1_cloneBarcodes/"))

bartitle = "primed"
barcodeclusters = sibACellsDistTable
allrandomprop = barPerClustA
#generate random samples
for (j in seeds) {
  set.seed(j)
  randomsibArows = as.integer(sample(rownames(jointUMAPFilter), length(nCellsPerBarcodeFilter)))
  randomsibA = suppressMessages(
    jointUMAPFilter[randomsibArows, ] %>%
      dplyr::select(cluster) %>%
      group_by(cluster) %>%
      summarise(randomcellspercluster = n()) %>%
      full_join(barPerClustA, by = c('cluster')) %>%
      mutate_if(is.numeric, ~ replace(., is.na(.), 0)) %>%
      mutate(randomprop = randomcellspercluster /
               allBarPerClust) %>%
      mutate(randomprop = randomprop / sum(randomprop)) %>%
      arrange(cluster)
  )
  allrandomprop[, paste0("randomprop", j)] = randomsibA$randomprop
}
#for entropy analysis
allrandomprop = allrandomprop %>% rowwise() %>%
  mutate(avgrandomprop = mean(c_across(starts_with("randomprop"))))
summaryrandomprop = allrandomprop %>%
  inner_join(barcodeclusters, by = 'cluster') %>%
  dplyr::select(cluster, prop, avgrandomprop)
#calculate entropy values
p_observed = summaryrandomprop %>% ungroup() %>%
  arrange(cluster) %>%
  summarise(prob_observed = prop + 0.000001) #add pseudocount to avoid zero in the numerator and/or denominator of KL
p_random = summaryrandomprop %>% ungroup() %>%
  arrange(cluster) %>%
  summarise(prob_random = avgrandomprop + 0.000001)

observedentropy = Entropy(p_observed, base = 2)
randomentropy = Entropy(p_random, base = 2)
JS_A[i] <- JSD(p_random, p_observed)

plotsummaryrandomprop = summaryrandomprop %>% transmute(cluster, observed_prop = prop, randomavg_prop = avgrandomprop) %>%
  pivot_longer(cols = -cluster, names_to = c("sample", ".value"), names_pattern = "(.+)_(.+)")
sample.labs <- c(paste0("observed ", bartitle, ": entropy = ", round(observedentropy, 3)),
                 paste0("random sample: entropy = ", round(randomentropy, 3)))
names(sample.labs) <- c("observed", "randomavg")
entropyplot_prop = ggplot(plotsummaryrandomprop, aes(x = as.factor(cluster), y = prop, fill = sample)) +
  geom_bar(stat = 'identity') +
  facet_wrap( ~ sample, ncol = 1, labeller = labeller(sample = sample.labs)) +
  scale_fill_manual(breaks = c("observed", "randomavg"), values = c("#009292", "gray80")) +
  ggtitle(paste(bartitle, ': ', length(nCellsPerBarcodeFilter), ' cells\n', 'Jensen-Shannon distance = ', round(JS_A[i], 3), sep = "")) +
  xlab("Seurat clusters, res=0.5") +
  ylab(paste0("normalized cell proportion")) +
  theme_classic() +
  theme(legend.position = "none")
pdf(paste(bartitle, "entropyPlot.pdf", sep="_"), width = 4, height = 4)
print(entropyplot_prop)
dev.off()

#evaluate JSD significance
for (k in seeds + n) {
  set.seed(k)
  JSDsibArows = as.integer(sample(rownames(jointUMAPFilter), length(nCellsPerBarcodeFilter)))
  JSDsibA = jointUMAPFilter[JSDsibArows, ] %>%
    dplyr::select(cluster) %>%
    full_join(barPerClustA, by = 'cluster') %>%
    group_by(cluster, allBarPerClust) %>%
    summarise(cellspercluster = n(), .groups = 'drop_last') %>%
    mutate(propperbarclust = cellspercluster / allBarPerClust) %>%
    ungroup() %>%
    mutate(prop = propperbarclust / sum(propperbarclust))
  p_test = JSDsibA %>% ungroup() %>%
    arrange(cluster) %>%
    summarise(prob_test = prop + 0.000001)
  JSTest[k - n] <- JSD(p_random, p_test)
}
meanJS_A[i] = mean(JSTest)
sdJS_A[i] = sd(JSTest)
normp = pnorm(JS_A[i], mean = meanJS_A[i], sd = sdJS_A[i], lower.tail = FALSE)
if (round(normp, 3) == 0) {
  normp3 <- "< 0.001"
} else{
  normp3 <- paste("=", round(normp, 3))
}
JSplot = ggplot() +
  geom_histogram(aes(x = JSTest), binwidth = 0.01, fill = 'gray80') +
  geom_vline(xintercept = JS_A[i], color = '#009292') +
  ggtitle(paste0(bartitle, ': ', length(nCellsPerBarcodeFilter), ' cells', '\np ', normp3)) +
  xlab("Jensen-Shannon distance") +
  coord_cartesian(xlim = c(0, 1)) + #without coord_cartesian throws an error about rows being removed
  theme_classic()
pdf(paste(bartitle, "JSDPlot.pdf", sep="_"), width = 4, height = 4)
print(JSplot)
dev.off()

#walk through different barcode in split A only
umap = ggplot() +
  geom_point(data = jointUMAPFilter, aes(x = UMAP_1, y = UMAP_2, color = "other")) +
  geom_point(data = sibACells, aes(x = UMAP_1, y = UMAP_2, color = "barcode"), size = 3) +
  scale_color_manual(breaks = c("barcode", "other"), values = c("#009292", "gray80")) +
  theme(legend.position = "none", axis.text = element_blank(), axis.ticks = element_blank()) +
  labs(color =  'Split A barcodes\nwith at least 20 cells') +
  ggtitle(paste0(bartitle))
pdf(paste(bartitle, "umapAlone.pdf", sep="_"), width = 4, height = 4)
print(umap)
dev.off()

#### PAIRWISE JSD CALCULATION #######################################################################################################################################################


###############pull all barcodes with ≥5 cells per split
plusnegbarcodebreakdown_noNA = plusnegbarcodebreakdown %>% na.omit() #81 of 1022 (2/0.3)
plusnegbarcodebreakdown_noNA_5ormore = plusnegbarcodebreakdown_noNA %>% 
  filter(sibA>4) %>% 
  filter(sibB>4) #18 barcodes

cellsperclustperSampleNum <- plusneglogNormalizedCountsSubsetWBarcodes %>% 
  group_by(BC30StarcodeD6, SampleNum, plusnegcells_clusters) %>% #group by res=0.5 clusters
  summarize(barSamclustCount=n()) %>% #gives number of cells per barcode per sibling per cluster
  inner_join(plusnegbarcodebreakdown_noNA_5ormore, by = "BC30StarcodeD6") %>%
  select(-sibA,-sibB) %>% 
  complete(barNum,plusnegcells_clusters,SampleNum) %>%
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>%
  filter(plusnegcells_clusters!="14") %>% #0 cells in cluster 14 - will throw an error if not removed
  group_by(SampleNum,plusnegcells_clusters) %>%
  mutate(allbarperclust = sum(barSamclustCount)) %>%
  mutate(propperbarclust = barSamclustCount/allbarperclust) %>%
  ungroup() %>%
  group_by(barNum,BC30StarcodeD6,SampleNum) %>%
  mutate(prop = propperbarclust/sum(propperbarclust)) %>%
  arrange(mixedrank(barNum)) #504 rows = 18 barcodes * 14 clusters * 2 siblings

#heatmap of pairwise JSD
propsibA = cellsperclustperSampleNum %>%
  filter(SampleNum == 'sibA') %>%
  ungroup() %>%
  mutate(prob = prop+0.000001) %>% #add to avoid NaNs due to zeroes
  select(barNum,plusnegcells_clusters,prob) %>%
  pivot_wider(id_cols = plusnegcells_clusters,names_prefix = "probA",names_from = barNum,values_from = prob)

propsibB = cellsperclustperSampleNum %>%
  filter(SampleNum == 'sibB') %>%
  ungroup() %>%
  mutate(prob = prop+0.000001) %>% #add to avoid NaNs due to zeroes
  select(barNum,plusnegcells_clusters,prob) %>%
  pivot_wider(id_cols = plusnegcells_clusters,names_prefix = "probB",names_from = barNum,values_from = prob)

JSDMatrix = matrix(, nrow = nrow(plusnegbarcodebreakdown_noNA_5ormore), ncol = nrow(plusnegbarcodebreakdown_noNA_5ormore))
#generates empty matrix that is 18x18 (for the 18 barcodes that are associated with at least 5 cells in both arms)
for (i in c(2:length(propsibA))) {
  for (j in c(2:length(propsibB))) {
    JSDMatrix[i-1,j-1] = JSD(propsibA[,i],propsibB[,j]) #start at 2 so have to index at i-1, j-1
  }
}

usenames = plusnegbarcodebreakdown_noNA_5ormore$barNum
JSDMatrixplot = JSDMatrix
rownames(JSDMatrixplot) = sub("$", "A", usenames)#colnames(propsibA)[2:19]
colnames(JSDMatrixplot) = sub("$", "B", usenames)#colnames(propsibB)[2:19]
melted_JSDMatrix <- melt(JSDMatrixplot, na.rm = TRUE) #list of JSD values
melted_JSDMatrix$value <- round(melted_JSDMatrix$value,1) #round to 1 decimal point
#reverse order per https://stackoverflow.com/questions/40725146/how-do-i-reverse-the-y-axis-in-ggplots-geom-tile

#pull for outlined rectangles
diag_JSDMatrix <- melted_JSDMatrix[which(melted_JSDMatrix$Var2 == sub("A", "B", melted_JSDMatrix$Var1)),]
plots_JSDMatrix <- na.omit(diag_JSDMatrix[which(melted_JSDMatrix$Var1 %in% c("B1A","B2A","B14A","B23A")),])

plot_JSD = ggplot(data = melted_JSDMatrix, aes(x=Var2, y=forcats::fct_rev(Var1), fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient(low = "darkblue", high = "white", 
                      limit = c(0,1), 
                      breaks = seq(0,1,0.2),
                      labels = seq(0,1,0.2),
                      space = "Lab", 
                      name="Jensen-Shannon distance") +
  theme_minimal()+ 
  #https://stackoverflow.com/questions/13258454/marking-specific-tiles-in-geom-tile-geom-raster
  geom_rect(data=diag_JSDMatrix, size=1, fill=NA, colour="black", #used to use #FFFF6D",
            aes(xmin=as.integer(Var2) - 0.5, xmax=as.integer(Var2) + 0.5, ymin=as.integer(forcats::fct_rev(Var1)) - 0.5, ymax=as.integer(forcats::fct_rev(Var1)) + 0.5)) +
  geom_rect(data=plots_JSDMatrix, size=1, fill=NA, colour="#ff068f",
            aes(xmin=as.integer(Var2) - 0.5, xmax=as.integer(Var2) + 0.5, ymin=as.integer(forcats::fct_rev(Var1)) - 0.5, ymax=as.integer(forcats::fct_rev(Var1)) + 0.5)) +
  coord_fixed() +
  geom_text(aes(Var2, forcats::fct_rev(Var1), label = value), color = "black", size = 3) +
  theme_classic((base_size = 10)) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.line = element_line(colour = 'black', size = 1.5),
    axis.text.x = element_text(angle = 90),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_line(colour = "black", size = 1.5),
    legend.position = "bottom",
    legend.direction = "horizontal",
    text=element_text(family="Helvetica"))+
  guides(fill = guide_colorbar(barwidth = 9, barheight = 1,
                               title.position = "top", title.hjust = 0.5)) +
  scale_y_discrete(position = "right")

#save JSD pairwise comparison heatmap
pdf(paste(plotDirectory,"JSDheatmap_1dec.pdf", sep=""), width=7, height=7)
print(plot_JSD)
dev.off()


#barplots and umap
setwd(paste0(plotDirectory))
for (i in c(1:nrow(plusnegbarcodebreakdown_noNA_5ormore))){
  barcode = plusnegbarcodebreakdown_noNA_5ormore$barNum[i]; print(barcode)
  barclusters = cellsperclustperSampleNum %>%
    filter(barNum == barcode) #%>% #should be 14 rows regardless
  bartitle = paste("barcode", barcode)
  
  #create probability distributions
  p_sibA = barclusters %>%
    filter(SampleNum == 'sibA') %>%
    ungroup() %>%
    arrange(plusnegcells_clusters)%>%
    summarise(prob_A = prop+0.000001)
  p_sibB = barclusters %>%
    filter(SampleNum == 'sibB') %>%
    ungroup() %>%
    arrange(plusnegcells_clusters)%>%
    summarise(prob_B = prop+0.000001)
  
  JS <- JSD(p_sibB,p_sibA)
  
  samplenum.labs <- c(paste0("barcode ",barcode,"A",": ",plusnegbarcodebreakdown_noNA_5ormore[i,]$sibA," cells"),
                      paste0("barcode ",barcode,"B",": ",plusnegbarcodebreakdown_noNA_5ormore[i,]$sibB," cells"))
  names(samplenum.labs) <- c("sibA", "sibB")
  plotnormcellsclustsample <- ggplot(barclusters,aes(x=plusnegcells_clusters,y=prop,group=SampleNum,fill=SampleNum)) +
    geom_bar(stat='identity')+
    facet_wrap(~SampleNum,ncol=1,labeller = labeller(SampleNum = samplenum.labs))+
    ylab(paste0("normalized proportion of ",barcode," cells"))+xlab("Seurat clusters, res=0.5")+
    scale_fill_manual(breaks=c("sibA","sibB"),values = c("#FFB6DB","#ff068f")) +
    ggtitle(paste(barcode,"A vs ", barcode,"B\nJensen-Shannon distance = ",round(JS,3),sep=""))+
    theme_classic() +
    theme(legend.position = "none")
  
  #save "probability distribution" bargraphs for each barcode
  pdf(paste(barcode,"normcellsperclusterpersample",".pdf", sep=""), width=4, height=4)
  print(plotnormcellsclustsample)
  dev.off()
  
  umap = ggplot() +
    geom_point(data = plusnegappendumapCoordinates[which((plusnegappendumapCoordinates$barNum %in% plusnegbarcodebreakdown_noNA_5ormore$barNum)),], aes(x = UMAP_1, y = UMAP_2, color = "other"), size = 1)+
    geom_point(data = plusnegappendumapCoordinates[which((plusnegappendumapCoordinates$barNum==barcode) & (plusnegappendumapCoordinates$SampleNum=="sibA")),], aes(x = UMAP_1, y = UMAP_2, color = "sibA"), size = 3) +
    geom_point(data = plusnegappendumapCoordinates[which((plusnegappendumapCoordinates$barNum==barcode) & (plusnegappendumapCoordinates$SampleNum=="sibB")),], aes(x = UMAP_1, y = UMAP_2, color = "sibB"), size = 3) +
    scale_color_manual(breaks = c("sibA","sibB","other"), values = c("sibA"="#FFB6DB","sibB"="#ff068f","other"="gray80"), labels = c(paste("barcode ",barcode,"A",sep=""),paste("barcode ",barcode,"B",sep=""),"other cells")) +
    create_lpr_theme() +
    theme(legend.position = "bottom",
          legend.title = element_text(size = rel(0.3)),
          legend.text = element_text(size = rel(0.3), angle = 20),
          axis.text = element_blank(),
          axis.ticks = element_blank()) +
    labs(color =  "Sample origin of barcoded cells") +
    ggtitle(paste0(bartitle))
  
  #save UMAP marking barcoded cells by which split they were found in
  pdf(paste(barcode,"umap.pdf", sep="_"), width=6, height=6)
  print(umap)
  dev.off()