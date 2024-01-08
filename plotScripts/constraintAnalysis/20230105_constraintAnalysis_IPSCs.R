rm(list=ls())
gc()

library(tidyverse)
library(Seurat)
library(ggsignif)
library(DescTools)
library(reshape2)

plotDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/plots/constraintAnalysis/"

homeDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/extractedData/fateMap10X/FM3/"

#### JSD CALCULATION FOR ALL CLONES #################################################################################################################################################
umapTable <- as_tibble(read.table(file = paste0(homeDirectory, "umapCoordinates_Scanorama_50pcs_filter.tsv"), header = TRUE, stringsAsFactors = F, sep = "\t"))
clusterTable <- as_tibble(read.table(file = paste0(homeDirectory, "umapClusters_Scanorama_50pcs_filter.tsv"), header = TRUE, stringsAsFactors = F, sep = "\t"))
filteredBarcodeTable <- as_tibble(read.table(file = paste0(homeDirectory, "filtered10XCells.txt"), header = TRUE, stringsAsFactors = F, sep = "\t"))

jointUMAP = inner_join(umapTable, clusterTable, by = c("cellID", "sampleNum")) %>% dplyr::rename(cluster = scanorama_snn_res.0.3) %>% 
  inner_join(filteredBarcodeTable, ., by = c("cellID")) %>% dplyr::select(-SampleNum, -nUMI, -fracUMI, -nLineages)
ggplot(jointUMAP, aes(x = UMAP_1, y = UMAP_2, color = factor(cluster))) +
  geom_point() + NoLegend()
jointUMAPFilter <- jointUMAP %>% dplyr::filter(!(cluster %in% c("6", "7", "8"))) %>% dplyr::filter(UMAP_2 < 10) %>% dplyr::filter(sampleNum %in% c("S1", "S2"))
ggplot(jointUMAPFilter, aes(x = UMAP_1, y = UMAP_2, color = factor(cluster))) +
  geom_point() + NoLegend()

nCellsPerBarcode <- jointUMAPFilter  %>% group_by(barcode) %>% summarise(nCells = n())
nCellsPerBarcodeFilter <- nCellsPerBarcode %>% dplyr::filter(nCells >= 4) %>% .$barcode

KLD <- function(x,y) sum(x * log2(x/y))
JSD <- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
mixedrank = function(x) order(gtools::mixedorder(x))

sibACells <- jointUMAPFilter %>% left_join(., nCellsPerBarcode, by = "barcode") %>% dplyr::filter(barcode %in% nCellsPerBarcodeFilter)
sibACellsDistTable <- sibACells %>% group_by(barcode, cluster, nCells) %>% summarise(cellsPerCluster = n()) %>% ungroup() %>%
  complete(nesting(barcode, nCells), cluster) %>% replace(., is.na(.), 0) %>%
  group_by(cluster) %>% mutate(allBarPerClust = sum(cellsPerCluster)) %>% mutate(propBarPerClust = cellsPerCluster / allBarPerClust) %>% ungroup() %>%
  group_by(barcode) %>% mutate(prop = propBarPerClust / sum(propBarPerClust)) %>% arrange(mixedrank(barcode))

barPerClustA <- sibACellsDistTable %>% ungroup() %>% dplyr::select(cluster, allBarPerClust) %>% unique()

n = 1000
seeds = 1:n
JSTest = rep(0, n)
JS_A = rep(0, length(nCellsPerBarcodeFilter))
meanJS_A = rep(0, length(nCellsPerBarcodeFilter))
sdJS_A = rep(0, length(nCellsPerBarcodeFilter))
pvalJS_A = rep(0, length(nCellsPerBarcodeFilter))
setwd(paste0(plotDirectory, "F3_cloneBarcodes/"))

for (i in 1:length(nCellsPerBarcodeFilter)) {
  barcodeTest = nCellsPerBarcodeFilter[i]
  bartitle = paste0("B", i)
  barcodeclusters = sibACellsDistTable %>%
    filter(barcode == barcodeTest)
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
  pvalJS_A[i] = normp
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
    geom_point(data = sibACells[which((sibACells$barcode == barcodeTest)), ], aes(x = UMAP_1, y = UMAP_2, color = "barcode"), size = 3) +
    scale_color_manual(breaks = c("barcode", "other"), values = c("#009292", "gray80"), "other") +
    theme(legend.position = "none", axis.text = element_blank(), axis.ticks = element_blank()) +
    labs(color =  'Split A barcodes\nwith at least 20 cells') +
    ggtitle(paste0(bartitle))
  pdf(paste(bartitle, "umapAlone.pdf", sep="_"), width = 4, height = 4)
  print(umap)
  dev.off()
}
#save JSD values to create table
nCellsPerBarcodeFilterTable <- nCellsPerBarcode %>% dplyr::filter(nCells >= 4)
JSTable_A = tibble(JS_A, meanJS_A, sdJS_A, barNum = paste0("B", rownames(nCellsPerBarcodeFilterTable)), barcode = nCellsPerBarcodeFilterTable$barcode, nCells = nCellsPerBarcodeFilterTable$nCells)
write_csv(JSTable_A, file = paste0("JSTableInfo.csv"))

JSTable_iPSC <- as_tibble(read.table(file = "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/plots/constraintAnalysis/F3_cloneBarcodes/JSTableInfo.csv",
                                  header = TRUE, stringsAsFactors = F, sep = ","))
JSTable_fibroblast <- as_tibble(read.table(file = "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/plots/constraintAnalysis/R1_cloneBarcodes/JSTableInfo.csv",
                                  header = TRUE, stringsAsFactors = F, sep = ","))
JSTable <- bind_rows(JSTable_iPSC %>% mutate(label = "iPSC"), JSTable_fibroblast %>% mutate(label = "fibroblast"))

theme_set(theme_classic())

ggplot(JSTable, aes(x = JS_A, y = meanJS_A)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  xlim(0, 1) + ylim(0, 1) +
  facet_wrap(~label)

ggplot(JSTable %>% bind_rows(., data.frame(JS_A = JSTable$meanJS_A, label = "random")), aes(x = label, y = JS_A)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.25) + ylab("average JSD value compared to random per clone") + theme(axis.title.x = element_blank())

ggplot(JSTable, aes(x = nCells, y = JS_A)) +
  geom_point() +
  geom_point(aes(y = meanJS_A), color = "red") +
  geom_errorbar(aes(ymin = meanJS_A - 1.96*sdJS_A, ymax = meanJS_A + 1.96*sdJS_A), color = "red") +
  facet_wrap(~label) +
  xlim(0, 50) + ylim(0, 1) +
  ylab("average JSD value compared to random per clone")

ggplot(JSTable, aes(x = nCells, y = JS_A)) +
  geom_point(pch = 16) +
  geom_point(aes(y = meanJS_A), color = "red", pch = 16) +
  geom_errorbar(aes(ymin = meanJS_A - 1.96*sdJS_A, ymax = meanJS_A + 1.96*sdJS_A), color = "red") +
  facet_wrap(~label, nrow = 2) +
  ylim(0, 1) + scale_x_continuous(trans = "log", breaks = c(0, 1, 10, 100, 1000)) +
  ylab("average JSD value compared to random per clone")
ggsave(filename = paste0(plotDirectory, "fibroblastVSiPSCJSDValuesCWRandom.pdf"), height = 5.5, width = 4.5, useDingbats = FALSE)

ggplot(JSTable %>% dplyr::filter(label == "iPSC"), aes(x = nCells, y = JS_A)) +
  geom_point(pch = 16) +
  geom_point(aes(y = meanJS_A), color = "red", pch = 16) +
  geom_errorbar(aes(ymin = meanJS_A - 1.96*sdJS_A, ymax = meanJS_A + 1.96*sdJS_A), color = "red") +
  ylim(0, 1) + scale_x_continuous(trans = "log", breaks = c(0, 1, 10, 100, 1000)) +
  ylab("average JSD value compared to random per clone")
ggsave(filename = paste0(plotDirectory, "iPSCJSDValuesCWRandom.pdf"), height = 4, width = 4, useDingbats = FALSE)


ggplot(JSTable, aes(x = nCells, y = meanJS_A)) +
  geom_point() +
  facet_wrap(~label) +
  xlim(0, 50) + ylim(0, 1) +
  ylab("average JSD value compared to random per clone")

#### PAIRWISE JSD CALCULATION #######################################################################################################################################################
umapTable <- as_tibble(read.table(file = paste0(homeDirectory, "umapCoordinates_Scanorama_50pcs_filter.tsv"), header = TRUE, stringsAsFactors = F, sep = "\t"))
clusterTable <- as_tibble(read.table(file = paste0(homeDirectory, "umapClusters_Scanorama_50pcs_filter.tsv"), header = TRUE, stringsAsFactors = F, sep = "\t"))
filteredBarcodeTable <- as_tibble(read.table(file = paste0(homeDirectory, "filtered10XCells.txt"), header = TRUE, stringsAsFactors = F, sep = "\t"))

jointUMAP = inner_join(umapTable, clusterTable, by = c("cellID", "sampleNum")) %>% dplyr::rename(cluster = scanorama_snn_res.0.3) %>% 
  inner_join(filteredBarcodeTable, ., by = c("cellID")) %>% dplyr::select(-SampleNum, -nUMI, -fracUMI, -nLineages)
ggplot(jointUMAP, aes(x = UMAP_1, y = UMAP_2, color = factor(cluster))) +
  geom_point() + NoLegend()
jointUMAPFilter <- jointUMAP %>% dplyr::filter(!(cluster %in% c("6", "7", "8"))) %>% dplyr::filter(UMAP_2 < 10)
ggplot(jointUMAPFilter, aes(x = UMAP_1, y = UMAP_2, color = factor(cluster))) +
  geom_point() + NoLegend()

nCellsPerBarcode <- jointUMAPFilter %>% group_by(barcode, sampleNum) %>% summarise(nCells = n())
nCellsPerBarcodeCast <- nCellsPerBarcode %>% dcast(., barcode ~ sampleNum)
# nCellsPerBarcodeFilter <- nCellsPerBarcodeCast %>% dplyr::filter(S1 > 3 & S2 > 3) %>% .$barcode
# nCellsPerBarcodeFilter <- nCellsPerBarcodeCast %>% dplyr::filter(S5 > 3 & S6 > 3) %>% .$barcode

jointUMAPComb <- jointUMAPFilter
jointUMAPComb$sampleNum <- ifelse(jointUMAPComb$sampleNum %in% c("S1", "S2"), "DMSO", jointUMAPComb$sampleNum)
jointUMAPComb$sampleNum <- ifelse(jointUMAPComb$sampleNum %in% c("S5", "S6"), "LSD1", jointUMAPComb$sampleNum)
nCellsPerBarcode <- jointUMAPComb %>% group_by(barcode, sampleNum) %>% summarise(nCells = n())
nCellsPerBarcodeCast <- nCellsPerBarcode %>% dcast(., barcode ~ sampleNum)
nCellsPerBarcodeFilter <- nCellsPerBarcodeCast %>% dplyr::filter(DMSO > 3 & LSD1 > 3) %>% .$barcode

KLD <- function(x,y) sum(x * log2(x/y))
JSD <- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
mixedrank = function(x) order(gtools::mixedorder(x))

# cellsForAnalysis <- jointUMAPFilter %>% left_join(., nCellsPerBarcode, by = c("barcode", "sampleNum")) %>% dplyr::filter(barcode %in% nCellsPerBarcodeFilter) %>% dplyr::filter(sampleNum %in% c("S1", "S2"))
# cellsForAnalysis <- jointUMAPFilter %>% left_join(., nCellsPerBarcode, by = c("barcode", "sampleNum")) %>% dplyr::filter(barcode %in% nCellsPerBarcodeFilter) %>% dplyr::filter(sampleNum %in% c("S5", "S6"))
cellsForAnalysis <- jointUMAPComb %>% left_join(., nCellsPerBarcode, by = c("barcode", "sampleNum")) %>% dplyr::filter(barcode %in% nCellsPerBarcodeFilter) %>% dplyr::filter(sampleNum %in% c("DMSO", "LSD1"))
cellsperclustperSampleNum <- cellsForAnalysis %>% group_by(barcode, sampleNum, cluster) %>% dplyr::summarise(cellsPerCluster = n()) %>% ungroup() %>%
  complete(nesting(barcode, sampleNum), cluster) %>% replace(., is.na(.), 0) %>%
  group_by(cluster) %>% mutate(allBarPerClust = sum(cellsPerCluster)) %>% mutate(propBarPerClust = cellsPerCluster / allBarPerClust) %>% ungroup() %>%
  group_by(barcode, sampleNum) %>% mutate(prop = propBarPerClust / sum(propBarPerClust)) %>% arrange(mixedrank(barcode))

barcodeNames <- data.frame(barcode = cellsperclustperSampleNum$barcode %>% unique()) %>% mutate(barNum = paste0("B", row.names(.)))
cellsperclustperSampleNum <- inner_join(cellsperclustperSampleNum, barcodeNames, by = "barcode")

#heatmap of pairwise JSD
propsibA = cellsperclustperSampleNum %>%
  # filter(sampleNum == 'S1') %>%
  # filter(sampleNum == 'S5') %>%
  filter(sampleNum == 'DMSO') %>%
  ungroup() %>%
  mutate(prob = prop+0.000001) %>% #add to avoid NaNs due to zeroes
  select(barNum,cluster,prob) %>%
  pivot_wider(id_cols = cluster,names_prefix = "probA",names_from = barNum,values_from = prob)

propsibB = cellsperclustperSampleNum %>%
  # filter(sampleNum == 'S2') %>%
  # filter(sampleNum == 'S6') %>%
  filter(sampleNum == 'LSD1') %>%
  ungroup() %>%
  mutate(prob = prop+0.000001) %>% #add to avoid NaNs due to zeroes
  select(barNum,cluster,prob) %>%
  pivot_wider(id_cols = cluster,names_prefix = "probB",names_from = barNum,values_from = prob)

JSDMatrix = matrix(, nrow = nrow(barcodeNames), ncol = nrow(barcodeNames))
#generates empty matrix that is 18x18 (for the 18 barcodes that are associated with at least 5 cells in both arms)
for (i in c(2:length(propsibA))) {
  for (j in c(2:length(propsibB))) {
    JSDMatrix[i-1,j-1] = JSD(propsibA[,i],propsibB[,j]) #start at 2 so have to index at i-1, j-1
  }
}

usenames = barcodeNames$barNum
JSDMatrix[upper.tri(JSDMatrix, diag = FALSE)] <- NA
JSDMatrixplot = JSDMatrix
rownames(JSDMatrixplot) = sub("$", "A", usenames)#colnames(propsibA)[2:19]
colnames(JSDMatrixplot) = sub("$", "B", usenames)#colnames(propsibB)[2:19]
melted_JSDMatrix <- melt(JSDMatrixplot, na.rm = TRUE) #list of JSD values
melted_JSDMatrix$value <- round(melted_JSDMatrix$value,2) #round to 1 decimal point
#reverse order per https://stackoverflow.com/questions/40725146/how-do-i-reverse-the-y-axis-in-ggplots-geom-tile

#pull for outlined rectangles
diag_JSDMatrix <- melted_JSDMatrix[which(melted_JSDMatrix$Var2 == sub("A", "B", melted_JSDMatrix$Var1)),]
plots_JSDMatrix <- na.omit(diag_JSDMatrix[,])

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

matrixValues <- bind_rows(data.frame(value = diag_JSDMatrix$value, label = "same clone"),
                         data.frame(value = melted_JSDMatrix$value, label = "diff clone"))
ggplot(matrixValues, aes(x = label, y = value)) +
  geom_boxplot() +
  geom_jitter(width = 0.25) +
  geom_signif(comparisons = list(c("diff clone", "same clone"))) +
  ylim(0, 1)

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