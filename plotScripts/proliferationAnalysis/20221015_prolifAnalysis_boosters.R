library(tidyverse)
library(R.matlab)
library(reshape2)
library(Seurat)
library(ggsignif)

theme_set(theme_classic())

dataDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/rawData/proliferationAnalysis/boosters/"
plotDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/plots/proliferationAnalysis/"

folders <- list.dirs(path = dataDirectory, full.names = FALSE)[-1]
subfolders <- folders[lengths(strsplit(folders, "/")) == 3] %>% .[c(9:10, 21:23, 33:56)]
exp <- strsplit(subfolders, "/") %>% sapply("[", 1)
days <- strsplit(subfolders, "/") %>% sapply("[", 2)
names <- strsplit(subfolders, "/") %>% sapply("[", 3) %>% strsplit("_") %>% sapply("[", 2)
rep <- strsplit(subfolders, "/") %>% sapply("[", 3) %>% strsplit("_") %>% sapply("[", 3)
names[1:5] = rep[1:5]
rep[1:5] = "1"

for(i in 1:length(subfolders)){
  if(i %in% 14:30) {
    matTemp <- readMat(paste0(dataDirectory, subfolders[i], "/splitScan/Cell_Info_Scan001.mat"))
  } else {
    matTemp <- readMat(paste0(dataDirectory, subfolders[i], "/Cell_Info_Scan001.mat"))
  }
  if(i == 1) {
    spotTable <- data.frame(exp = exp[i], day = days[i], condition = names[i], count = as.integer(matTemp$cells[[1]][[1]] %>% nrow(.)), rep = rep[i])
  } else {
  spotTableTemp <- data.frame(exp = exp[i], day = days[i], condition = names[i], count = as.integer(matTemp$cells[[1]][[1]] %>% nrow(.)), rep = rep[i])
  spotTable <- bind_rows(spotTable, spotTableTemp)
  }
}

spotTable$exp <- spotTable$exp %>% strsplit("_") %>% sapply("[", 1)
spotTable$day[6:29] <- spotTable$day[6:29] %>% strsplit("_") %>% sapply("[", 2)

spotTableCast <- dcast(spotTable, exp + condition + rep ~ day, value.var = "count")
spotTableCast$'day 3'[2] <- mean(c(spotTableCast$'day 3'[1], spotTableCast$'day 3'[3]))
spotTableCast <- spotTableCast %>% mutate(prolif = NA, prolifPerDay = NA)
spotTableCast$prolif <- ifelse(!(is.na(spotTableCast$'day 3')) & !(is.na(spotTableCast$'day 6')), spotTableCast$'day 6'/spotTableCast$'day 3', spotTableCast$prolif)
spotTableCast$prolif <- ifelse(!(is.na(spotTableCast$day1)) & !(is.na(spotTableCast$day4)), spotTableCast$day4/spotTableCast$day1, spotTableCast$prolif)
spotTableCast$prolif <- ifelse(!(is.na(spotTableCast$day1)) & !(is.na(spotTableCast$day3)), spotTableCast$day3/spotTableCast$day1, spotTableCast$prolif)
spotTableCast$prolifPerDay <- ifelse(!(is.na(spotTableCast$'day 3')) & !(is.na(spotTableCast$'day 6')), (spotTableCast$prolif)^1/3, spotTableCast$prolifPerDay)
spotTableCast$prolifPerDay <- ifelse(!(is.na(spotTableCast$day1)) & !(is.na(spotTableCast$day4)), (spotTableCast$prolif)^1/2.5, spotTableCast$prolifPerDay)
spotTableCast$prolifPerDay <- ifelse(!(is.na(spotTableCast$day1)) & !(is.na(spotTableCast$day3)), (spotTableCast$prolif)^1/2, spotTableCast$prolifPerDay)

spotTableCast$condition <- ifelse(spotTableCast$condition %in% c("control", "DMSO"), "control", spotTableCast$condition)
spotTableCast$condition <- factor(spotTableCast$condition, levels = c("control", "LSD1i", "DOT1Li"))

means <- spotTableCast %>% group_by(condition) %>% summarise(mean = mean(prolifPerDay))

ggplot(spotTableCast, aes(x = condition, y = prolifPerDay, fill = condition)) +
  stat_summary(fun = "mean", geom = "col") +
  geom_hline(yintercept = means %>% dplyr::filter(condition == "control") %>% .$mean, linetype = "dashed") +
  stat_summary(fun = "mean", geom = "point", size = 2.5, color = "black") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.25, color = "black") +
  geom_signif(comparisons = list(c("control", "LSD1i"), c("control", "DOT1Li"), c("LSD1i", "DOT1Li")), test = "t.test") +
  NoLegend() + theme(axis.title.x = element_blank()) + ylab("normalized proliferation rate")
ggsave(filename = paste0(plotDirectory, "prolifRatesPerBooster_onlyLSD1i.pdf"), units = "in", height = 3, width = 1.5, useDingbats = FALSE)

ggplot(spotTableCast %>% dplyr::filter(condition != "DOT1Li"), aes(x = condition, y = prolifPerDay, fill = condition)) +
  stat_summary(fun = "mean", geom = "col") +
  geom_hline(yintercept = means %>% dplyr::filter(condition == "control") %>% .$mean, linetype = "dashed") +
  stat_summary(fun = "mean", geom = "point", size = 2.5, color = "black") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.25, color = "black") +
  geom_signif(comparisons = list(c("control", "LSD1i")), test = "t.test") +
  NoLegend() + theme(axis.title.x = element_blank()) + ylab("normalized proliferation rate")
ggsave(filename = paste0(plotDirectory, "prolifRatesPerBooster_LSD1iAndDOT1Li.pdf"), units = "in", height = 3, width = 2, useDingbats = FALSE)

ggplot(spotTableCast %>% dplyr::filter(condition != "control"), aes(x = condition, y = prolifPerDay, fill = condition)) +
  stat_summary(fun = "mean", geom = "col") +
  geom_hline(yintercept = means %>% dplyr::filter(condition == "control") %>% .$mean, linetype = "dashed") +
  stat_summary(fun = "mean", geom = "point", size = 2.5, color = "black") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.25, color = "black") +
  NoLegend() + theme(axis.title.x = element_blank()) + ylab("normalized proliferation rate")
