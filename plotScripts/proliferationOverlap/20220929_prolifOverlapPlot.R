library(tidyverse)
library(ggsignif)
library(ggpubr)
library(Seurat)

theme_set(theme_classic())

dataDirectory <- '/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Original Manuscript/extractedData/prolifOverlap/'
plotDirectory <- '/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Original Manuscript/plots/prolifOverlap/'

exp1Table <- readRDS(file = paste0(dataDirectory, "overlapMatrixTable_exp1.rds")) %>% mutate(exp = "1")
exp2Table <- readRDS(file = paste0(dataDirectory, "overlapMatrixTable_exp2.rds")) %>% mutate(exp = "2")
exp3Table <- readRDS(file = paste0(dataDirectory, "overlapMatrixTable_exp3.rds")) %>% mutate(exp = "3")

finalTable <- bind_rows(exp1Table, exp2Table, exp3Table)

finalTable$name <- factor(finalTable$name, levels = c("slow", "control", "fast"), labels = c("slow", "control", "fast"))

ggplot(finalTable %>% dplyr::filter(exp %in% c("1", "2", "3")), aes(x = name, y = value, fill = name)) +
  stat_summary(fun = "mean", geom = "bar") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.25) +
  # geom_point() +
  # stat_compare_means(comparisons = list(c("slow", "control"), c("control", "fast"), c("slow", "fast"))) +
  facet_grid(variable ~ colonyThresh, scales = "fixed") + ylab("normalized ratio of overlap in lineage barcodes\nbetween each sorted population and iPSCs") +
  theme(legend.position = "none", axis.title.x = element_blank())
ggsave(filename = paste0(plotDirectory, "prolifOverlapAggregateAll.pdf"), units = "in", height = 5, width = 5, useDingbats = FALSE)

ggplot(finalTable %>% filter(exp %in% c("1", "2", "3"), !(colonyThresh == 10), !(name == "control")), aes(x = colonyThresh, y = value, color = name, group = name)) +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.25) +
  facet_wrap(~variable, ncol = 4) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  ylab("normalized ratio of overlap in lineage barcodes\nbetween each sorted population and iPSCs") +
  xlab("cutoff for number of barcoded cells to\nto be called an iPSC colony") +
  theme(legend.position = "none", axis.title.x = element_blank(), panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  scale_color_manual(values = c("#ca6702", "#e9d8a6"))
ggsave(filename = paste0(plotDirectory, "prolifOverlapAggregateAllLine.pdf"), units = "in", height = 1.75, width = 6, useDingbats = FALSE)

ggplot(finalTable %>% filter(exp %in% c("1", "2", "3")), aes(x = colonyThresh, y = value, color = name, group = name)) +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.25) +
  facet_wrap(~variable+exp, ncol = 3) +
  ylab("normalized ratio of overlap in lineage barcodes\nbetween each sorted population and iPSCs") +
  xlab("cutoff for number of barcoded cells to\nto be called an iPSC colony") +
  theme(legend.position = "none", axis.title.x = element_blank())

ggplot(finalTable %>% filter(exp %in% c("1", "2", "3")), aes(x = name, y = value, fill = name, group = name)) +
  stat_summary(fun = "mean", geom = "bar") +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.25) +
  facet_wrap(~variable, ncol = 2) + ylab("normalized ratio of overlap in lineage barcodes\nbetween each sorted population and iPSCs") +
  theme(legend.position = "none", axis.title.x = element_blank())

finalTableFilter <- finalTable %>% dplyr::filter(variable == 0.5, colonyThresh == 25)

ggplot(finalTableFilter %>% dplyr::filter(exp %in% c("1", "2", "3"), !(name == "control")), aes(x = name, y = value, fill = name)) +
  stat_summary(fun = "mean", geom = "bar") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.25) +
  geom_jitter(height = 0, width = 0.25, size = 1.5) +
  ylab("normalized ratio of overlap in lineage barcodes\nbetween each sorted population and iPSCs") +
  theme(legend.position = "none", axis.title.x = element_blank()) +
  geom_hline(yintercept = 1, linetype = "dashed")
ggsave(filename = paste0(plotDirectory, "prolifOverlapAggregateMain.pdf"), units = "in", height = 2, width = 2, useDingbats = FALSE)

#### plot overlap as a function of cutoff ####
sampleOverlapPlot1 <- readRDS(file = paste0(dataDirectory, "sampleOverlapPlot_exp1.rds")) %>% mutate(exp = "1")
sampleOverlapPlot2 <- readRDS(file = paste0(dataDirectory, "sampleOverlapPlot_exp2.rds")) %>% mutate(exp = "2")
sampleOverlapPlot3 <- readRDS(file = paste0(dataDirectory, "sampleOverlapPlot_exp3.rds")) %>% mutate(exp = "3")

finalOverlap <- bind_rows(sampleOverlapPlot1, sampleOverlapPlot2, sampleOverlapPlot3)

ggplot(finalOverlap, aes(x = cutoff, y = overlap, group = sample, color = sample)) +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun = mean, geom = "point") +
  # stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.25) +
  ylim(0, 1) + theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  scale_color_manual(values = c("#ca6702", "#a7a9ac", "#e9d8a6")) + NoLegend()
ggsave(filename = paste0(plotDirectory, "prolifSampleOverlapAggregate.pdf"), units = "in", height = 1.75, width = 1.75, useDingbats = FALSE)
