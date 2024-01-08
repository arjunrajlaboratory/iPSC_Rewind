library(tidyverse)
library(R.matlab)
library(reshape2)
library(Seurat)
library(ggpubr)

theme_set(theme_classic())

homeDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Revisions/extractedData/prolifAnalysis/"
dataDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/rawData/proliferationAnalysis/CRISPR/"
plotDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Revisions/plots/prolifAnalysis/"

folders1 <- list.dirs(path = paste0(dataDirectory, "R1/"), full.names = FALSE)[-1]
subfolders1 <- folders1[lengths(strsplit(folders1, "/")) == 2]
exp1 <- strsplit(subfolders1, "/") %>% sapply("[", 1) %>% strsplit(., "_") %>% sapply("[", 2)
days1 <- strsplit(subfolders1, "/") %>% sapply("[", 1) %>% strsplit(., "_") %>% sapply("[", 3)
genes1 <- strsplit(subfolders1, "/") %>% sapply("[", 2) %>% strsplit("_") %>% sapply("[", 2)
rep1 <- strsplit(subfolders1, "/") %>% sapply("[", 2) %>% strsplit("_") %>% sapply("[", 3)

folders2 <- list.dirs(path = paste0(dataDirectory, "R2/"), full.names = FALSE)[-1]
subfolders2 <- folders2[lengths(strsplit(folders2, "/")) == 3]
exp2 <- strsplit(subfolders2, "/") %>% sapply("[", 1) %>% strsplit(., "_") %>% sapply("[", 4)
days2 <- strsplit(subfolders2, "/") %>% sapply("[", 1) %>% strsplit(., "_") %>% sapply("[", 2)
genes2 <- strsplit(subfolders2, "/") %>% sapply("[", 2) %>% strsplit("_") %>% sapply("[", 2)
rep2 <- strsplit(subfolders2, "/") %>% sapply("[", 2) %>% strsplit("_") %>% sapply("[", 3)

subfolders <- c(subfolders1, subfolders2)
exp <- c(exp1, exp2)
days <- c(days1, days2)
genes <- c(genes1, genes2)
rep <- c(rep1, rep2)

for(i in 1:length(subfolders)){
  if(i < 29) {
    matTemp <- readMat(paste0(dataDirectory, "R1/", subfolders[i], "/Cell_Info_Scan001.mat")) 
  } else {
    matTemp <- readMat(paste0(dataDirectory, "R2/", subfolders[i], "/Cell_Info_Scan001.mat")) 
  }
  if(i == 1) {
    spotTable <- data.frame(exp = exp[i], day = days[i], gene = genes[i], rep = rep[i], count = as.integer(matTemp$cells[[1]][[1]] %>% nrow(.)))
  } else {
  spotTableTemp <- data.frame(exp = exp[i], day = days[i], gene = genes[i], rep = rep[i], count = as.integer(matTemp$cells[[1]][[1]] %>% nrow(.)))
  spotTable <- bind_rows(spotTable, spotTableTemp)
  }
}

spotTableCast <- dcast(spotTable, exp + genes + rep ~ day, value.var = "count")
spotTableCast <- spotTableCast %>% mutate(prolif = spotTableCast$"day3"/spotTableCast$"day1") %>%
  rowwise() %>% mutate(prolifPerDay = sqrt(prolif)) %>% ungroup()
spotTableCast <- spotTableCast %>% mutate(guide = paste0(genes, "_", rep))
spotTableCast <- spotTableCast %>% dplyr::filter(!(genes %in% c("C2", "C4")))
spotTableCast$guide <- ifelse(spotTableCast$genes %in% c("B1", "B2", "control"), "control", spotTableCast$guide)

prolifPerDayNormList <- c()
for(i in c("r1", "R2")) {
  spotTableCastTemp <- spotTableCast %>% dplyr::filter(exp == i)
  controlMean <- spotTableCastTemp %>% dplyr::filter(guide == "control") %>% .$prolifPerDay %>% mean(.)
  spotTableCastTemp <- spotTableCastTemp %>% rowwise() %>% mutate(prolifPerDayNorm = prolifPerDay/controlMean) %>% ungroup()
  prolifPerDayNormList <- c(prolifPerDayNormList, spotTableCastTemp$prolifPerDayNorm)
}

spotTableCast$prolifPerDayNorm <- prolifPerDayNormList

ggplot(spotTableCast, aes(x = guide, y = prolifPerDayNorm)) +
  stat_summary(fun = "mean", geom = "col") +
  stat_summary(fun = "mean", geom = "point", size = 2.5, color = "black") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.25, color = "black") +
  geom_hline(yintercept = mean(spotTableCast %>% dplyr::filter(guide == "control") %>% .$prolifPerDay), linetype = "dashed") +
  NoLegend() + theme(axis.title.x = element_blank()) + ylab("proliferation rate per day") + ylim(0, 3)

write.csv(spotTableCast, file = paste0(homeDirectory, "prolifAnalysis_CRISPR_noCDKN1A.csv"))
spotTableCast <- read.csv(file = paste0(homeDirectory, "prolifAnalysis_CRISPR_CDKN1A.csv"))

spotTableCast$genesCorr <- factor(spotTableCast$genesCorr, levels = c("control", "MDM2", "CDKN1A", "KDM1A", "SPP1"))

ggplot(spotTableCast, aes(x = genesCorr, y = prolifPerDayNorm)) +
  stat_summary(fun = "mean", geom = "col") +
  stat_summary(fun = "mean", geom = "point", size = 2.5, color = "black") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.25, color = "black") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  stat_compare_means(comparisons = list(c("control", "MDM2"), c("control", "CDKN1A"), c("control", "KDM1A"), c("control", "SPP1"))) +
  NoLegend() + theme(axis.title.x = element_blank()) + ylab("normalized proliferation rate")
ggsave(filename = paste0(plotDirectory, "prolifRatePerCRISPRGuide.pdf"), units = "in", height = 3, width = 4, useDingbats = FALSE)
