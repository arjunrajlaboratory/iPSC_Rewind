library(tidyverse)
library(ggpubr)
library(ggridges)
library(ggsignif)
library(reshape2)

theme_set(theme_classic())

homeDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/extractedData/inSituCarbonCopies/R2_TM9/"
plotDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/plots/inSituCarbonCopies/"

#### extract data from .csv files ####
#####################################################################################################################
path1 <- paste0(homeDirectory, "well2_1-25_barcodes,PRMT1,AXL.csv")
path2 <- paste0(homeDirectory, "well2_1-25_SOX21,CENPF,HES1.csv")
path3 <- paste0(homeDirectory, "well2_1-25_TOP2A,ASPM,MKI67.csv")

bcTable1 = read.csv(path1, header = TRUE) %>% dplyr::rename(barcode = cy.RNACounts) %>% dplyr::select(-isGood, -alexa.RNACounts, -tmr.RNACounts)
bcTableTemp = read.csv(path2, header = TRUE) %>% dplyr::rename(SOX21 = tmr.RNACounts, CENPF = alexa.RNACounts, HES1 = cy.RNACounts) %>% dplyr::select(-isGood)
bcTable1 <- inner_join(bcTable1, bcTableTemp, by = c("objArrayNum", "objNum")) %>% mutate(primed = "nonprimed")
bcTableTemp = read.csv(path3, header = TRUE) %>% dplyr::rename(MKI67 = tmr.RNACounts, ASPM = alexa.RNACounts, TOP2A = cy.RNACounts)
bcTable1 <- inner_join(bcTable1, bcTableTemp, by = c("objArrayNum", "objNum")) %>% mutate(primed = "nonprimed")
bcTable1$primed <- ifelse(bcTable1$barcode > 5, "primed", bcTable1$primed)
bcTable1 <- bcTable1 %>% dplyr::filter(isGood == "true")

give.n <- function(x){
  return(c(y = median(x), label = median(x))) 
}

ggplot(bcTable1, aes(x = primed, y = TOP2A)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.25) +
  #stat_summary(fun = mean, geom = "point", size = 5, pch = 23, fill = "white") +
  stat_summary(fun = median, geom = "text", aes(label = round(..y.., 0)), position = position_nudge(x = -0.4, y = 0)) +
  theme(axis.title.x = element_blank()) +
  geom_signif(comparisons = list(c("nonprimed", "primed"))) #+ facet_wrap(~exp)
ggsave(filename = paste0(plotDirectory, "TOP2A.pdf"), units = "in", height = 2.5, width = 2.5)

ggplot(bcTable1, aes(x = primed, y = CENPF)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.25) +
  #stat_summary(fun = mean, geom = "point", size = 5, pch = 23, fill = "white") +
  stat_summary(fun = median, geom = "text", aes(label = round(..y.., 0)), position = position_nudge(x = -0.4, y = 0)) +
  theme(axis.title.x = element_blank()) +
  geom_signif(comparisons = list(c("nonprimed", "primed"))) #+ facet_wrap(~exp)
ggsave(filename = paste0(plotDirectory, "CENPF.pdf"), units = "in", height = 2.5, width = 2.5)
