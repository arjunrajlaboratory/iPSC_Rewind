####LOADING LIBRARIES####
library(tidyverse)
library(msm)
library(reshape2)

####INITIALIZING PARAMETERS####
cellNumber = c(300000, 400000, 300000)
barcodeNumber = 20000000
reprogEffic = 0.01
MOI = 1
divisionNumber = c(3, 6, 3)
splitNumber = c(2, 4, 2)
materialLoss = 0.05
libPrepLoss = 0.05
repetitions = 1000
heritability = 0

####INITIALIZING RESULTS####
results <- data.frame(matrix(, nrow = repetitions, ncol = length(cellNumber)))

####BARCODE SIMULATION####
for(n in 1:length(cellNumber)){
  for (i in 1:repetitions) {
    n = 2
    i = 1
    #SEEDING
    initialCells <- data.frame(
      "Cell ID" = 1:cellNumber[n],
      "Barcode" = rep(0, cellNumber[n]),
      "Daughters" = rep(0, cellNumber[n]),
      "Split ID" = rep(0, cellNumber[n]),
      "Loss1" = rep(FALSE, cellNumber[n]),
      "Loss2" = rep(FALSE, cellNumber[n]),
      "Reprogram" = rep(FALSE, cellNumber[n]))
    
    #BARCODING
    barcodes <- c(1:barcodeNumber)
    sampledBarcodes <- sample(barcodes, cellNumber[n], replace = TRUE)
    
    poissonCutoff <- 1 - dpois(0, MOI, log = FALSE)
    barcodeIndex <- runif(cellNumber[n], 0, 1) < poissonCutoff
    
    sampledBarcodes[barcodeIndex == FALSE] <- 0
    initialCells$Barcode <- sampledBarcodes
    
    #VARIABLE DIVISIONS
    iterativeCells <- initialCells
    
    sampledDaughters <- round(rtnorm(nrow(iterativeCells), 2, 0.1, 0)^divisionNumber[n])
    iterativeCells$Daughters <- sampledDaughters
    
    iterativeCells <- as.data.frame(lapply(iterativeCells, rep, iterativeCells$Daughters))
    
    #SPLITTING
    splits <- c(1:splitNumber[n])
    sampledSplits <- sample(splits, nrow(iterativeCells), replace = TRUE)
    
    iterativeCells$Split.ID <- sampledSplits
    
    #MATERIAL LOSS
    lossIndex1 <- runif(nrow(iterativeCells), 0, 1) < materialLoss
    iterativeCells$Loss1 <- lossIndex1
    
    finalCells <- filter(iterativeCells, Loss1 == FALSE)
    
    #REPROGRAMMING
    plate1 <- filter(finalCells, Split.ID == 1)
    plate2 <- filter(finalCells, Split.ID == 2)
    
    plate1$Reprogram <- runif(nrow(plate1), 0, 1) < reprogEffic
    plate2$Reprogram <- runif(nrow(plate2), 0, 1) < reprogEffic
    
    plate1ReprogSampled <- filter(plate1, Reprogram == TRUE)
    plate2ReprogSampled <- filter(plate2, Reprogram == TRUE)
    
    #LIBPREP LOSS
    lossIndex21 <- runif(nrow(plate1ReprogSampled), 0, 1) < libPrepLoss
    plate1ReprogSampled$Loss2 <- lossIndex21
    
    lossIndex22 <- runif(nrow(plate2ReprogSampled), 0, 1) < libPrepLoss
    plate2ReprogSampled$Loss2 <- lossIndex22
    
    plate1Final <- filter(plate1ReprogSampled, Loss2 == FALSE)
    plate2Final <- filter(plate2ReprogSampled, Loss2 == FALSE)
    
    #OVERLAP ANALYSIS
    plate1ReprogBarcodes <- unique(select(filter(plate1Final, Barcode != 0), Barcode))
    plate2ReprogBarcodes <- unique(select(filter(plate2Final, Barcode != 0), Barcode))
    
    if (nrow(union(plate1ReprogBarcodes, plate2ReprogBarcodes)) == 0) {
      overlap <- 0
    } else{
      overlap <- nrow(intersect(plate1ReprogBarcodes, plate2ReprogBarcodes)) / nrow(union(plate1ReprogBarcodes, plate2ReprogBarcodes))
    }
    
    results[i, n] <- overlap
  }
}

####EXPORTING RESULTS####
saveRDS(results, file = '/project/arjunrajlab/hiFTTM/repo/localScriptsForR/barcodeSimulations/simulationResults1in100.RDS')
