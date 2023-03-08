rm(list=ls())
gc()

library(tidyverse)
library(spgs)

theme_set(theme_classic())

homeDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/extractedData/fateMap10X/FM3/"
plotDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/plots/fateMap10X/FM3/"

#### connect 10X lineage barcodes with individual cellIDs ####
#######################################################################################################################################################
barcodes = as_tibble(read.table(paste0(homeDirectory, 'stepThreeStarcodeShavedReads_BC_10XAndGDNA.txt'), stringsAsFactors = F, header = T))
barcodes10X <- barcodes %>% filter(cellID != "dummy") %>% dplyr::select(cellID, UMI, BC50StarcodeD8, SampleNum) %>% rename(barcode = BC50StarcodeD8) %>% unique()
barcodes10X %>% .$cellID %>% unique() %>% length()
barcodes10XFilter <- barcodes10X[grepl("([AT][CG][ATCG]){8,}", barcodes10X$barcode), ]
barcodes10XFilter %>% .$cellID %>% unique() %>% length()

umiCut = 0
fracUMICut = 0.4

upToLineageCounting <- barcodes10XFilter %>% group_by(cellID, barcode, SampleNum) %>% summarise(nUMI = length(cellID)) %>% filter(nUMI >= umiCut) %>%
  group_by(cellID, SampleNum) %>% mutate(fracUMI = nUMI/sum(nUMI)) %>% filter(fracUMI >= fracUMICut) %>%
  group_by(cellID) %>% mutate(nLineages = length(cellID))
upToLineageCounting$SampleNum = sub("^", "S", upToLineageCounting$SampleNum)
upToLineageCounting %>% .$cellID %>% unique() %>% length()
cellsWithManyLineagesSummary = upToLineageCounting %>% group_by(nLineages) %>% summarise(nCells = length(unique(cellID)))
ggplot(cellsWithManyLineagesSummary, aes(x=nLineages, y=nCells/sum(nCells))) +
  geom_bar(stat="identity") + xlab("number of lineages per cell") + ylab("fraction of cells") + ylim(0, 1)
nCellsPerBarcode <- upToLineageCounting %>% group_by(barcode) %>% summarise(nCells = length(cellID))
ggplot(nCellsPerBarcode,aes(x=nCells)) + geom_histogram() + xlab("number of cells per barcode") + ylab("count")

linCut = 1
linCountToOverlaps <- upToLineageCounting %>% ungroup() %>% filter(nLineages <= linCut) %>% unique()
# linCountToOverlaps <- linCountToOverlaps[!(linCountToOverlaps$barcode %in% c("AGCATGCCCGCCAGCGCTACTCCGCACCTGCCACCCTCAGCTCTCGGGCC",
#                                                                              "ATTAGAAGGTCTAGGTCCAGTTCCACTACGACAAGTGAACTACGTGCTCT")), ]
linCountToOverlaps %>% .$cellID %>% length()

nCellsPerBarcodeAll <- linCountToOverlaps %>% group_by(barcode) %>% summarise(nCells = length(cellID))
ggplot(nCellsPerBarcodeAll,aes(x=nCells)) + geom_histogram() + xlab("number of cells per barcode") + ylab("count")
ggplot(nCellsPerBarcodeAll %>% dplyr::filter(nCells > 5),aes(x=nCells)) + geom_histogram(bins = 50) + xlab("number of cells per barcode") + ylab("count") + xlim(0, 100)

write.table(x = linCountToOverlaps, file = paste0(homeDirectory, "filtered10XCells.txt"), sep = "\t", col.names = TRUE)