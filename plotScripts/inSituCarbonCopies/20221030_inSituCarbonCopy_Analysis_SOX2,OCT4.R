library(tidyverse)
library(ggpubr)
library(ggridges)
library(ggsignif)
library(reshape2)

theme_set(theme_classic())

homeDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Original Manuscript/extractedData/inSituCarbonCopies/R1_TM8/"
plotDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Original Manuscript/plots/inSituCarbonCopies/"

#### extract data from .csv files ####
#####################################################################################################################
path1 <- paste0(homeDirectory, "well1_barcodes.csv")
path2 <- paste0(homeDirectory, "well1_SOX2,OCT4,HES1.csv")

pathList <- c(path1, path2)

bcTable = read.csv(pathList[1], header = TRUE) %>% dplyr::rename(barcode = cy.RNACounts) %>% dplyr::select(-isGood)
bcTableTemp = read.csv(pathList[2], header = TRUE) %>% dplyr::rename(SOX2 = tmr.RNACounts, OCT4 = alexa.RNACounts, HES1 = cy.RNACounts)
bcTable <- inner_join(bcTable, bcTableTemp, by = c("objArrayNum", "objNum")) %>% mutate(primed = "nonprimed")
bcTable$primed <- ifelse(bcTable$barcode > 5, "primed", bcTable$primed)
bcTable <- bcTable %>% dplyr::filter(isGood == "true")

#### graph boxplots  ####
#####################################################################################################################
give.n <- function(x){
  return(c(y = median(x), label = median(x))) 
}

ggplot(bcTable, aes(x = OCT4, y = primed, height = stat(density))) + 
  geom_density_ridges2(stat = "binline", bins = 40, scale = 0.95) +
  xlim(-1, 40) +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank())
ggsave(filename = paste0(plotDirectory, "nonprimedVSprimed_OCT4.pdf"), units = "in", height = 2, width = 2.5)

ggplot(bcTable, aes(x = SOX2, y = primed, height = stat(density))) + 
  geom_density_ridges2(stat = "binline", bins = 40, scale = 0.95) +
  xlim(-1, 40) +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank())
ggsave(filename = paste0(plotDirectory, "nonprimedVSprimed_SOX2.pdf"), units = "in", height = 2, width = 2.5)

ggplot(bcTable,aes(x=primed,y=OCT4)) +
  geom_boxplot(aes(fill=primed),outlier.shape=NA) +
  geom_jitter(position = position_jitter(height = 0.25, width = 0.25)) +
  geom_signif(comparisons = list(c("nonprimed", "primed")))

ggplot(bcTable,aes(x=primed,y=SOX2)) +
  geom_boxplot(aes(fill=primed),outlier.shape=NA) +
  geom_jitter(position = position_jitter(height = 0.25, width = 0.25)) +
  geom_signif(comparisons = list(c("nonprimed", "primed")))


### check OKSM expression upon DOX induction ###
homeDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Original Manuscript/extractedData/inSituCarbonCopies/post-DOX_OKSM/"

path1 <- paste0(homeDirectory, "results_1.csv")
path2 <- paste0(homeDirectory, "results_2.csv")
path3 <- paste0(homeDirectory, "results_3.csv")

pathList <- c(path1, path2, path3)

bcTablePostTemp1 = read.csv(pathList[1], header = TRUE) %>% dplyr::rename(KLF4 = alexa.RNACounts) %>% dplyr::filter(isGood == "true")
bcTablePostTemp2 = read.csv(pathList[2], header = TRUE) %>% dplyr::rename(SOX2 = alexa.RNACounts, KLF4 = cy.RNACounts) %>% dplyr::filter(isGood == "true")
bcTablePostTemp3 = read.csv(pathList[3], header = TRUE) %>% dplyr::rename(SOX2 = alexa.RNACounts, POU5F1 = cy.RNACounts) %>% dplyr::filter(isGood == "TRUE")

bcTablePost <- bind_rows(data.frame(group = "KLF4", value = bcTablePostTemp1$KLF4),
                         data.frame(group = "SOX2", value = bcTablePostTemp2$SOX2),
                         data.frame(group = "KLF4", value = bcTablePostTemp2$KLF4),
                         data.frame(group = "SOX2", value = bcTablePostTemp3$SOX2))

nOff1 <- bcTablePost %>% dplyr::filter(group == "KLF4", value < 50) %>% nrow()
nTotal1 <- bcTablePost %>% dplyr::filter(group == "KLF4") %>% nrow()
onRate1 <- (nTotal1 - nOff1) / (nTotal1)

nOff2 <- bcTablePost %>% dplyr::filter(group == "SOX2", value < 50) %>% nrow()
nTotal2 <- bcTablePost %>% dplyr::filter(group == "SOX2") %>% nrow()
onRate2 <- (nTotal2 - nOff2) / (nTotal2)

ggplot(bcTablePost, aes(x = group, y = value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.25) +
  geom_hline(yintercept = 50, linetype = "dashed") +
  stat_summary(fun.data = mean_sdl, geom = "errorbar")

ggplot(bcTablePost, aes(y = group, x = value, height = stat(density))) +
  geom_density_ridges2(stat = "binline", bins = 50, scale = 0.95) +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  geom_vline(xintercept = 50, linetype = "dashed") + xlim(-10, 800)
ggsave(filename = paste0(plotDirectory, "SOX2andKLF4LevelsInduced_histogram.pdf"), height = 3, width = 3, useDingbats = FALSE)

ggplot(bcTablePostTemp2, aes(x = SOX2, y = KLF4)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_smooth(method = "lm", formula = y~x-1, se = FALSE) +
  stat_cor(aes(label = paste0(..rr.label.., "~", ..p.label..)), geom = "text")
ggsave(filename = paste0(plotDirectory, "SOX2andKLF4LevelsInduced_scatterplot.pdf"), height = 3, width = 3, useDingbats = FALSE)

bcTableComb <- bind_rows(bcTable, bcTablePostTemp2, bcTablePostTemp3 %>% dplyr::select(-isGood))
bcTableComb$label <- "-DOX"
bcTableComb$label <- ifelse(is.na(bcTableComb$primed), "+DOX", bcTableComb$label)

ggplot(bcTableComb, aes(x = label, y = SOX2)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.25) +
  geom_hline(yintercept = 50, linetype = "dashed") +
  geom_signif(comparisons = list(c("-DOX", "+DOX"))) +
  theme(axis.title.x = element_blank())
ggsave(filename = paste0(plotDirectory, "SOX2LevelsUninducedvsInduced_boxplot.pdf"), height = 3, width = 3, useDingbats = FALSE)

ggplot(bcTableComb, aes(y = label, x = SOX2, height = stat(density))) +
  geom_density_ridges2(stat = "binline", bins = 50, scale = 0.95) +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  geom_vline(xintercept = 50, linetype = "dashed") + xlim(-10, 800)

ggplot(bcTableComb, aes(y = label, x = SOX2)) +
  geom_density_ridges(bins = 50, scale = 0.95, quantile_lines = TRUE) +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank()) +
  geom_vline(xintercept = 50, linetype = "dashed")
