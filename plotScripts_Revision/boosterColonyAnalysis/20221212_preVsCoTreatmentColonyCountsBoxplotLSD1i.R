rm(list=ls())
gc()

library(tidyverse)
library(ggsignif)

homeDirectory = "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Original Manuscript/extractedData/boosterColonyCounts/preVsCoTreatmentLSD1i/"
plotDirectory = "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Original Manuscript/plots/boosterColonyAnalysis/"

sampleTablesList <- list(rep('NA', length(sampleDirList)))
sampleFoldersList <- list(rep('NA', length(sampleDirList)))

sampleFiles <- list.files(path = homeDirectory, pattern = 'results.csv', recursive = TRUE)
for(i in 1:length(sampleFiles)) {
  sampleTableTemp <- read.table(paste0(homeDirectory, sampleFiles[i]), header = TRUE, sep = ",")
  if(i == 1) {
    sampleTableTemp <- sampleTableTemp %>%
      mutate(preCond = gsub("(.*)\\_(.*)\\_(.*).tif", "\\1", .$Dataset)) %>%
      mutate(postCond = gsub("(.*)\\_(.*)\\_(.*).tif", "\\2", .$Dataset)) %>%
      mutate(rep = gsub("(.*)\\_(.*)\\_(.*).tif", "\\3", .$Dataset)) %>%
      mutate(exp = paste0(i)) %>%
      dplyr::filter(preCond %in% c("preControl", "preLSD1i"), postCond %in% c("DMSO", "LSD1i"))
    sampleTableTemp$preCond <- ifelse(sampleTableTemp$preCond == "preControl", "preDMSO", sampleTableTemp$preCond)
    controlMeanTemp <- sampleTableTemp %>% dplyr::filter(preCond == "preDMSO" & postCond == "DMSO") %>% .$NumBlobs %>% mean()
    sampleTableTemp <- sampleTableTemp %>% rowwise() %>% mutate(numBlobNorm = NumBlobs/controlMeanTemp)
    
    sampleTable <- sampleTableTemp
  } else {
    sampleTableTemp <- sampleTableTemp %>%
      mutate(preCond = gsub("(.*)\\_(.*)\\_(.*)\\_(.*).tif", "\\2", .$Dataset)) %>%
      mutate(postCond = gsub("(.*)\\_(.*)\\_(.*)\\_(.*).tif", "\\3", .$Dataset)) %>%
      mutate(rep = gsub("(.*)\\_(.*)\\_(.*)\\_(.*).tif", "\\4", .$Dataset)) %>%
      mutate(exp = paste0(i))
    controlMeanTemp <- sampleTableTemp %>% dplyr::filter(preCond == "preDMSO" & postCond == "DMSO") %>% .$NumBlobs %>% mean()
    sampleTableTemp <- sampleTableTemp %>% rowwise() %>% mutate(numBlobNorm = NumBlobs/controlMeanTemp)
    
    sampleTable <- bind_rows(sampleTable, sampleTableTemp)
  }
}

sampleTableMean <- sampleTable %>% group_by(preCond, postCond, exp) %>% summarize(meanNumBlobNorm = mean(numBlobNorm)) %>% mutate(cond = paste0(preCond, "_", postCond))
sampleTableMean$cond <- factor(sampleTableMean$cond, levels = c("preDMSO_DMSO", "preLSD1i_DMSO", "preDMSO_LSD1i", "preLSD1i_LSD1i"), labels = c("only OKSM", "pre-treatment", "co-treatment", "pre- and\nco-treatment"))

ggplot(sampleTableMean, aes(x = cond, y = meanNumBlobNorm)) +
  stat_summary(fun = mean, geom = "col") +
  geom_jitter(height = 0, width = 0.25, size = 2.5) +
  #stat_summary(fun = mean, geom = "point", size = 5) +
  #stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.25) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_signif(comparisons = list(c("only OKSM", "pre-treatment"), c("co-treatment", "pre- and\nco-treatment")), test = "t.test")

ggsave(file = paste0(plotDirectory, "colonyCountBoxplot_LSD1iOnly.pdf"), height = 4, width = 4)
