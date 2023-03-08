library(tidyverse)
library(ggsignif)
library(Seurat)

dataDirectory = "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/extractedData/CRISPR/"
plotDirectory = "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/plots/CRISPR/"

#### ROUND 1 ####
data1 <- read.csv(file = paste0(dataDirectory, "results_round1.csv"))

data1 <- data1 %>% mutate(cond = gsub("(.*)\\_(.*)\\_(.*).tif", "\\1", .$Dataset)) %>%
  mutate(guide = gsub("(.*)\\_(.*)\\_(.*).tif", "\\2", .$Dataset)) %>%
  mutate(rep = gsub("(.*)\\_(.*)\\_(.*) - Stitched.tif", "\\3", .$Dataset)) %>%
  mutate(round = "1")

ggplot(data1, aes(x = guide, y = NumBlobs)) +
  geom_jitter(height = 0.1, width = 0.1) + facet_grid(~cond, scales = "free_x", space = "free") +
  geom_hline(yintercept = mean((data1 %>% filter(cond == "backbone"))$NumBlobs), linetype = "dashed") +
  stat_summary(fun = mean, geom = "bar", fun.min = mean, fun.max = mean, alpha = 0.5) +
  stat_summary(fun.data = "mean_se", fun.args = list(mult=1), geom = "errorbar", width = 0.5) +
  theme_classic(base_size = 10) + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("number of colonies") + NoLegend()

controlMean1 <- data1 %>% filter(cond == "backbone") %>% .$NumBlobs %>% mean()

data1Norm <- data1 %>% rowwise() %>% mutate(NumBlobsNorm = NumBlobs/controlMean1) %>% ungroup()
data1Norm $cond <- ifelse(data1Norm $cond == "backbone", "B1", data1Norm$cond)

ggplot(data1Norm, aes(x = guide, y = NumBlobsNorm)) +
  geom_jitter(height = 0.1, width = 0.1) + facet_grid(~cond, scales = "free_x", space = "free") +
  geom_hline(yintercept = mean((data1Norm %>% filter(cond == "B1"))$NumBlobsNorm), linetype = "dashed") +
  stat_summary(fun = mean, geom = "bar", fun.min = mean, fun.max = mean, alpha = 0.5) +
  stat_summary(fun.data = "mean_se", fun.args = list(mult=1), geom = "errorbar", width = 0.5) +
  theme_classic(base_size = 10) + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("number of colonies normalized to backbone control") + NoLegend()

#### ROUND 2 ####
data2.1 <- read.csv(file = paste0(dataDirectory, "results_round2.1.csv"))
data2.1 <- data2.1 %>% mutate(cond = gsub("(.*)\\_(.*)\\_(.*).tif", "\\1", .$Dataset)) %>%
  mutate(guide = gsub("(.*)\\_(.*)\\_(.*).tif", "\\2", .$Dataset)) %>%
  mutate(rep = gsub("(.*)\\_(.*)\\_(.*) - Stitched.tif", "\\3", .$Dataset)) %>%
  mutate(round = "2.1")
controlMean2.1 <- data2.1 %>% filter(cond %in% c("B1", "B2")) %>% .$NumBlobs %>% mean()
dataNorm2.1 <- data2.1 %>% rowwise() %>% mutate(NumBlobsNorm = NumBlobs/controlMean2.1) %>% ungroup()
ggplot(dataNorm2.1, aes(x = guide, y = NumBlobsNorm)) +
  geom_jitter(height = 0.1, width = 0.1) + facet_grid(~cond, scales = "free_x", space = "free") +
  geom_hline(yintercept = mean(dataNorm2.1 %>% filter(cond %in% c("B1", "B2")) %>% .$NumBlobsNorm), linetype = "dashed") +
  stat_summary(fun = mean, geom = "bar", fun.min = mean, fun.max = mean, alpha = 0.5) +
  stat_summary(fun.data = "mean_se", fun.args = list(mult=1), geom = "errorbar", width = 0.5) +
  theme_classic(base_size = 10) + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("number of colonies normalized to backbone control") + NoLegend()

data2.2 <- read.csv(file = paste0(dataDirectory, "results_round2.2.csv"))
data2.2 <- data2.2 %>% mutate(cond = gsub("(.*)\\_(.*)\\_(.*).tif", "\\1", .$Dataset)) %>%
  mutate(guide = gsub("(.*)\\_(.*)\\_(.*).tif", "\\2", .$Dataset)) %>%
  mutate(rep = gsub("(.*)\\_(.*)\\_(.*) - Stitched.tif", "\\3", .$Dataset)) %>%
  mutate(round = "2.2")
controlMean2.2 <- data2.2 %>% filter(cond %in% c("B1", "B2")) %>% .$NumBlobs %>% mean()
dataNorm2.2 <- data2.2 %>% rowwise() %>% mutate(NumBlobsNorm = NumBlobs/controlMean2.2) %>% ungroup()
ggplot(dataNorm2.2, aes(x = guide, y = NumBlobsNorm)) +
  geom_jitter(height = 0.1, width = 0.1) + facet_grid(~cond, scales = "free_x", space = "free") +
  geom_hline(yintercept = mean(dataNorm2.2 %>% filter(cond %in% c("B1", "B2")) %>% .$NumBlobsNorm), linetype = "dashed") +
  stat_summary(fun = mean, geom = "bar", fun.min = mean, fun.max = mean, alpha = 0.5) +
  stat_summary(fun.data = "mean_se", fun.args = list(mult=1), geom = "errorbar", width = 0.5) +
  theme_classic(base_size = 10) + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("number of colonies normalized to backbone control") + NoLegend()
dataNorm2.2$IntensityThresh <- as.integer(dataNorm2.2$IntensityThresh)

data2.3 <- read.csv(file = paste0(dataDirectory, "results_round2.3.csv"))
data2.3 <- data2.3 %>% mutate(cond = gsub("(.*)\\_(.*)\\_(.*).tif", "\\1", .$Dataset)) %>%
  mutate(guide = gsub("(.*)\\_(.*)\\_(.*).tif", "\\2", .$Dataset)) %>%
  mutate(rep = gsub("(.*)\\_(.*)\\_(.*) - Stitched.tif", "\\3", .$Dataset)) %>%
  mutate(round = "2.3")
controlMean2.3 <- data2.3 %>% filter(cond %in% c("B1", "B2")) %>% .$NumBlobs %>% mean()
dataNorm2.3 <- data2.3 %>% rowwise() %>% mutate(NumBlobsNorm = NumBlobs/controlMean2.3) %>% ungroup()
ggplot(dataNorm2.3, aes(x = guide, y = NumBlobsNorm)) +
  geom_jitter(height = 0.1, width = 0.1) + facet_grid(~cond, scales = "free_x", space = "free") +
  geom_hline(yintercept = mean(dataNorm2.3 %>% filter(cond %in% c("B1", "B2")) %>% .$NumBlobsNorm), linetype = "dashed") +
  stat_summary(fun = mean, geom = "bar", fun.min = mean, fun.max = mean, alpha = 0.5) +
  stat_summary(fun.data = "mean_se", fun.args = list(mult=1), geom = "errorbar", width = 0.5) +
  theme_classic(base_size = 10) + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("number of colonies normalized to backbone control") + NoLegend()

data2Norm <- bind_rows(dataNorm2.1, dataNorm2.2, dataNorm2.3)
ggplot(data2Norm, aes(x = guide, y = NumBlobsNorm)) +
  geom_jitter(height = 0.1, width = 0.1) + facet_grid(~cond, scales = "free_x", space = "free") +
  geom_hline(yintercept = mean(data2Norm %>% filter(cond %in% c("B1", "B2")) %>% .$NumBlobsNorm), linetype = "dashed") +
  stat_summary(fun = mean, geom = "bar", fun.min = mean, fun.max = mean, alpha = 0.5) +
  stat_summary(fun.data = "mean_se", fun.args = list(mult=1), geom = "errorbar", width = 0.5) +
  theme_classic(base_size = 10) + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("number of colonies normalized to backbone control") + NoLegend()

ggplot(data2Norm, aes(x = guide, y = NumBlobs)) +
  geom_jitter(height = 0.1, width = 0.1) + facet_grid(~cond, scales = "free_x", space = "free") +
  geom_hline(yintercept = mean(data2Norm %>% filter(cond %in% c("B1", "B2")) %>% .$NumBlobs), linetype = "dashed") +
  stat_summary(fun = mean, geom = "bar", fun.min = mean, fun.max = mean, alpha = 0.5) +
  stat_summary(fun.data = "mean_se", fun.args = list(mult=1), geom = "errorbar", width = 0.5) +
  theme_classic(base_size = 10) + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("number of colonies normalized to backbone control") + NoLegend()

data2Norm$cond <- factor(data2Norm$cond, levels = c("B1", "B2", "C2", "C4", "KDM1A", "DOT1L", "SPP1", "FTH1", "CDKN1A", "MDM2", "MYBL2", "NFE2L2"))
ggplot(data2Norm %>% dplyr::filter(!(cond %in% c("C2", "C4"))), aes(x = guide, y = NumBlobsNorm)) +
  stat_summary(fun = mean, geom = "bar", fun.min = mean, fun.max = mean) +
  stat_summary(fun.data = "mean_se", fun.args = list(mult=1), geom = "errorbar", width = 0.5) +
  geom_jitter(height = 0.1, width = 0.1) + facet_grid(~cond, scales = "free_x", space = "free") +
  geom_hline(yintercept = mean(data2Norm %>% filter(cond %in% c("B1", "B2")) %>% .$NumBlobsNorm), linetype = "dashed") +
  theme_classic(base_size = 10) + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("number of colonies normalized to backbone control") + NoLegend()
ggsave(filename = paste0(plotDirectory, "R2_CRISPR_AllGuides.pdf"), height = 2.5, width = 10, useDingbats = FALSE)

#### combined ####
dataNormComb <- bind_rows(data1Norm %>% group_by(cond, guide) %>% summarize(NumBlobsNormMean = mean(NumBlobsNorm)) %>% mutate(rep = "1"),
                          data2Norm %>% group_by(cond, guide) %>% summarize(NumBlobsNormMean = mean(NumBlobsNorm)) %>% mutate(rep = "2"))
dataNormComb$cond <- ifelse(dataNormComb$cond == "backbone", "B1", dataNormComb$cond)

dataNormComb <- bind_rows(data1Norm %>% mutate(exp = "1"), data2Norm %>% mutate(exp = "2")) %>% filter(cond %in% c("B1", "KDM1A", "MDM2", "SPP1"))
ggplot(dataNormComb %>% filter(exp == "1"), aes(x = guide, y = NumBlobsNorm)) +
  geom_jitter(height = 0.1, width = 0.1) + facet_grid(~cond, scales = "free_x", space = "free") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  stat_summary(fun = mean, geom = "bar", fun.min = mean, fun.max = mean, alpha = 0.5) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.5) + theme_classic() +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 0, hjust = 1)) +
  ylab("foldchange in number of colonies compared to backbone control") + NoLegend() + ylim(0,15)
ggsave(paste0(plotDirectory, "R1_CRISPRKO.pdf"), units = "in", height = 2, width = 6, useDingbats = FALSE)
ggplot(dataNormComb %>% filter(exp == "2"), aes(x = guide, y = NumBlobsNorm)) +
  geom_jitter(height = 0.1, width = 0.1) + facet_grid(~cond, scales = "free_x", space = "free") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  stat_summary(fun = mean, geom = "bar", fun.min = mean, fun.max = mean, alpha = 0.5) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.5) + theme_classic() +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 0, hjust = 1)) +
  ylab("foldchange in number of colonies compared to backbone control") + NoLegend() + ylim(0,15)
ggsave(paste0(plotDirectory, "R2_CRISPRKO.pdf"), units = "in", height = 2, width = 6, useDingbats = FALSE)
