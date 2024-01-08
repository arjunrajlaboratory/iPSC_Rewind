library(tidyverse)
library(ggpubr)
library(ggridges)
library(ggsignif)
library(reshape2)

theme_set(theme_classic())

homeDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Original Manuscript/extractedData/inSituCarbonCopies/SQSTM1, SPP1/"
plotDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Original Manuscript/plots/inSituCarbonCopies/"

#### extract data from .csv files ####
#####################################################################################################################
path1 <- '/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Original Manuscript/extractedData/inSituCarbonCopies/SQSTM1, SPP1/TM12_R1_well2_counts.csv'
path2 <- '/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Original Manuscript/extractedData/inSituCarbonCopies/SQSTM1, SPP1/TM12_R2_well2_counts.csv'
path3 <- '/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Original Manuscript/extractedData/inSituCarbonCopies/SQSTM1, SPP1/TM9_well1_counts.csv'
path4 <- '/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Original Manuscript/extractedData/inSituCarbonCopies/SQSTM1, SPP1/TM9_well2_counts.csv'
pathList <- c(path1, path2)

bcTable <- read.csv(path1, header = TRUE) %>% dplyr::filter(isGood != FALSE) %>% mutate(exp = "TM12_1")
bcTableTemp <- read.csv(path2, header = TRUE) %>% dplyr::filter(isGood != FALSE) %>% mutate(exp = "TM12_2")
bcTable <- bind_rows(bcTable, bcTableTemp)
bcTableTemp <- read.csv(path3, header = TRUE) %>% dplyr::filter(isGood != FALSE) %>% mutate(exp = "TM9_1")
bcTable <- bind_rows(bcTable, bcTableTemp)
bcTableTemp <- read.csv(path4, header = TRUE) %>% dplyr::filter(isGood != FALSE) %>% mutate(exp = "TM9_2")
bcTable <- bind_rows(bcTable, bcTableTemp)
bcTable$barcode <- ifelse(bcTable$barcode == "primed", bcTable$barcode, "nonprimed")

give.n <- function(x){
  return(c(y = median(x), label = median(x))) 
}

ggplot(bcTable, aes(x = barcode, y = SPP1)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.25) +
  #stat_summary(fun = mean, geom = "point", size = 5, pch = 23, fill = "white") +
  stat_summary(fun = median, geom = "text", aes(label = round(..y.., 0)), position = position_nudge(x = -0.4, y = 0)) +
  theme(axis.title.x = element_blank()) +
  geom_signif(comparisons = list(c("nonprimed", "primed"))) #+ facet_wrap(~exp)
ggsave(filename = paste0(plotDirectory, "SPP1.pdf"), units = "in", height = 2.5, width = 2.5)

ggplot(bcTable, aes(x = barcode, y = SQSTM1)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.25) +
  #stat_summary(fun = mean, geom = "point", size = 5, pch = 23, fill = "white") +
  stat_summary(fun = median, geom = "text", aes(label = round(..y.., 0)), position = position_nudge(x = -0.4, y = 0)) +
  theme(axis.title.x = element_blank()) +
  geom_signif(comparisons = list(c("nonprimed", "primed"))) #+ facet_wrap(~exp)
ggsave(filename = paste0(plotDirectory, "SQSTM1.pdf"), units = "in", height = 2.5, width = 2.5)
