library(tidyverse)
library(R.matlab)
library(reshape2)
library(Seurat)
library(ggsignif)

theme_set(theme_classic())

dataDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Original Manuscript/rawData/proliferationAnalysis/prolifSort/"
plotDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Original Manuscript/plots/proliferationAnalysis/"

folders <- list.dirs(path = dataDirectory, full.names = FALSE)[-1]
subfolders <- folders[lengths(strsplit(folders, "/")) == 4]
exp <- strsplit(subfolders, "/") %>% sapply("[", 1)
days <- strsplit(subfolders, "/") %>% sapply("[", 2)
names <- strsplit(subfolders, "/") %>% sapply("[", 3) %>% strsplit("_") %>% sapply("[", 2)
rep <- strsplit(subfolders, "/") %>% sapply("[", 3) %>% strsplit("_") %>% sapply("[", 3)

for(i in 1:length(subfolders)){
  matTemp <- readMat(paste0(dataDirectory, subfolders[i], "/Cell_Info_Scan001.mat"))
  if(i == 1) {
    spotTable <- data.frame(exp = exp[i], day = days[i], condition = names[i], count = as.integer(matTemp$cells[[1]][[1]] %>% nrow(.)))
  } else {
  spotTableTemp <- data.frame(exp = exp[i], day = days[i], condition = names[i], count = as.integer(matTemp$cells[[1]][[1]] %>% nrow(.)))
  spotTable <- bind_rows(spotTable, spotTableTemp)
  }
}

spotTableCast <- dcast(spotTable, exp + condition + rep ~ day, value.var = "count")
spotTableCast <- spotTableCast %>% mutate(prolif = NA, prolifPerDay = NA)
spotTableCast$prolif <- ifelse(!(is.na(spotTableCast$day2)) & !(is.na(spotTableCast$day4)), spotTableCast$day4/spotTableCast$day2, spotTableCast$prolif)
spotTableCast$prolif <- ifelse(!(is.na(spotTableCast$day1)) & !(is.na(spotTableCast$day3)), spotTableCast$day3/spotTableCast$day1, spotTableCast$prolif)
spotTableCast$prolif <- ifelse(!(is.na(spotTableCast$day1)) & !(is.na(spotTableCast$day4)), spotTableCast$day4/spotTableCast$day1, spotTableCast$prolif)
spotTableCast$prolifPerDay <- ifelse(!(is.na(spotTableCast$day2)) & !(is.na(spotTableCast$day4)), (spotTableCast$prolif)^(1/2), spotTableCast$prolifPerDay)
spotTableCast$prolifPerDay <- ifelse(!(is.na(spotTableCast$day1)) & !(is.na(spotTableCast$day3)), (spotTableCast$prolif)^(1/2), spotTableCast$prolifPerDay)
spotTableCast$prolifPerDay <- ifelse(!(is.na(spotTableCast$day1)) & !(is.na(spotTableCast$day4)), (spotTableCast$prolif)^(1/3), spotTableCast$prolifPerDay)

spotTableCast$condition <- factor(spotTableCast$condition, levels = c("high", "midhigh", "control", "midlow", "low"))

means <- spotTableCast %>% group_by(condition) %>% summarise(mean = mean(prolifPerDay))

ggplot(spotTableCast, aes(x = condition, y = prolifPerDay, fill = condition)) +
  stat_summary(fun = "mean", geom = "col") +
  geom_hline(yintercept = means %>% dplyr::filter(condition == "control") %>% .$mean, linetype = "dashed") +
  stat_summary(fun = "mean", geom = "point", size = 2.5, color = "black") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.25, color = "black") +
  geom_signif(comparisons = list(c("control", "low"), c("control", "high"), c("high", "low")), test = "t.test")
  NoLegend() + theme(axis.title.x = element_blank()) + ylab("normalized proliferation rate")

spotTableCast$prolifPerDayNorm <- spotTableCast$prolifPerDay/(means %>% dplyr::filter(condition == "control") %>% .$mean)
ggplot(spotTableCast %>% dplyr::filter(condition != "control"), aes(x = condition, y = prolifPerDayNorm, fill = condition)) +
  stat_summary(fun = "mean", geom = "col") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_jitter(height = 0, width = 0.25, size = 2) +
  #stat_summary(fun = "mean", geom = "point", size = 2.5, color = "black") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.25, color = "black") +
  NoLegend() + theme(axis.title.x = element_blank()) + ylab("normalized proliferation rate")
ggsave(filename = paste0(plotDirectory, "prolifRatesPerCondition.pdf"), units = "in", height = 1.5, width = 3, useDingbats = FALSE)

dataDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Original Manuscript/extractedData/proliferationSort/"
plotDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Original Manuscript/plots/proliferationSort/"
data <- read.csv(file = paste0(dataDirectory, "proliferationSortAggregatecColonyCounts.csv"))
dataFilter <- data %>% filter(!is.na(AverageNormalizedCounts)) %>% filter(!(Condition == "control"))
dataFilter$Condition <- factor(dataFilter$Condition, levels = c("high", "midhigh", "midlow", "low"), labels = c("slow", "mid slow", "mid fast", "fast"))

ggplot(dataFilter, aes(x = Condition, y = AverageNormalizedCounts)) +
  stat_summary(fun = "mean", geom = "col") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_jitter(height = 0, width = 0.25, size = 2) +
  #stat_summary(fun = "mean", geom = "point", size = 2.5, color = "black") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.25, color = "black") +
  NoLegend() + theme(axis.title.x = element_blank()) + ylab("normalized reprogramming rate")
ggsave(paste0(plotDirectory, "cyclingReprogrammingRatePlot.pdf"), height = 1.5, width = 3)

prolifRateTable <- spotTableCast %>% dplyr::group_by(condition) %>% summarize(meanProlif = mean(prolifPerDay), sdProlif = sd(prolifPerDay)) %>%
  dplyr::filter(condition %in% c("high", "low"))
prolifRateTable$condition <- factor(prolifRateTable$condition, levels = c("high", "low"), labels = c("slow", "fast"))
reprogRateTable <- dataFilter %>% dplyr::group_by(Condition) %>% summarize(meanReprog = mean(AverageNormalizedCounts), sdReprog = sd(AverageNormalizedCounts)) %>%
  dplyr::filter(Condition %in% c("slow", "fast"))

combRateTable <- inner_join(prolifRateTable, reprogRateTable, by = c("condition" = "Condition"))
combRateTable <- combRateTable %>% rowwise() %>%
  mutate(day0 = meanReprog / meanProlif^0) %>%
  mutate(day2 = meanReprog / meanProlif^2) %>%
  mutate(day4 = meanReprog / meanProlif^4) %>%
  mutate(day6 = meanReprog / meanProlif^6) %>%
  mutate(day8 = meanReprog / meanProlif^8)
combRateTableMelt <- combRateTable %>% melt(id.vars = c("condition", "meanProlif", "sdProlif", "meanReprog", "sdReprog"))

ggplot(combRateTableMelt, aes(x = variable, y = value, fill = condition)) +
  geom_col(position = "dodge") + theme(legend.position = "none")
ggsave(filename = paste0(plotDirectory, "correctedReprogRateByProlifRate.pdf"), units = "in", height = 3, width = 6, useDingbats = FALSE)
