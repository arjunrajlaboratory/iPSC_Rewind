rm(list=ls())
gc()

library(tidyverse)
library(spgs)
library(concaveman)
library(egg)
library(reshape2)
library(ggrepel)
library(ggforce)
library(ggsignif)
library(ggridges)
library(ggpmisc)
library(viridis)
library(VennDiagram)
library(ggalluvial)

theme_set(theme_classic())

homeDirectory <- "/Users/naveenj/Dropbox (RajLab)/Shared_Naveen/Revisions/extractedData/continuousCulture/"
plotDirectory <- "/Users/naveenj/Dropbox (RajLab)/Shared_Naveen/Revisions/plots/continuousCulture/"

overlapData_1 <- readRDS(file = paste0(homeDirectory, "exp1/overlapTableNorm_observed.rds")) %>% mutate(rep = "1")
overlapData_2 <- readRDS(file = paste0(homeDirectory, "exp3/overlapTableNorm_observed.rds")) %>% mutate(rep = "2")
overlapData <- bind_rows(overlapData_1, overlapData_2)

ggplot(overlapData, aes(x = cond, y = overlap, fill = rep)) +
  geom_col(stat = "identity", position = "dodge") +
  theme(axis.title.x = element_blank(), legend.position = "none")
ggsave(filename = paste0(plotDirectory, "normalizedOverlap.pdf"), height = 2, width = 4)

abundanceData_1 <- readRDS(file = paste0(homeDirectory, "exp1/primedAbundanceCorrelation_observed.rds")) %>% mutate(rep = "1")
abundanceData_2 <- readRDS(file = paste0(homeDirectory, "exp3/primedAbundanceCorrelation_observed.rds")) %>% mutate(rep = "2")
abundanceData_2$SampleNum <- factor(abundanceData_2$SampleNum, levels = c("rep1", "rep2", "rep3", "rep4"), labels = c("1", "2", "3", "4"))
abudnanceData <- bind_rows(abundanceData_1, abundanceData_2)

ggplot(abudnanceData, aes(x = SampleNum, y = nUMIFracMean, group = rep, color = rep)) +
  geom_line() + facet_grid(~label) + theme(legend.position = "none")
ggsave(filename = paste0(plotDirectory, "cloneAbundancePlot.pdf"), height = 2, width = 5)

overlapTablePrimed_observed_1 <- readRDS(file = paste0(homeDirectory, "exp1/overlapTablePrimed_observed.rds")) %>% mutate(rep = "1")
overlapTablePrimed_observed_2 <- readRDS(file = paste0(homeDirectory, "exp3/overlapTablePrimed_observed.rds")) %>% mutate(rep = "2")
overlapTablePrimed_observed <- bind_rows(overlapTablePrimed_observed_1, overlapTablePrimed_observed_2)
ggplot(overlapTablePrimed_observed, aes(x = cond, y = overlap, fill = rep)) +
  geom_col(stat = "identity", position = "dodge") +
  theme(axis.title.x = element_blank(), legend.position = "none") + ylim(0, 0.75)
ggsave(filename = paste0(plotDirectory, "overlap_primed_observed.pdf"), height = 2, width = 4)

overlapTablePrimed_same_1 <- readRDS(file = paste0(homeDirectory, "exp1/overlapTablePrimed_same.rds")) %>% mutate(rep = "1")
overlapTablePrimed_same_2 <- readRDS(file = paste0(homeDirectory, "exp3/overlapTablePrimed_same.rds")) %>% mutate(rep = "2")
overlapTablePrimed_same <- bind_rows(overlapTablePrimed_same_1, overlapTablePrimed_same_2)
ggplot(overlapTablePrimed_same, aes(x = cond, y = overlap, fill = rep)) +
  geom_col(stat = "identity", position = "dodge") +
  theme(axis.title.x = element_blank(), legend.position = "none") + ylim(0, 0.75)
ggsave(filename = paste0(plotDirectory, "overlap_primed_same.pdf"), height = 2, width = 4)

overlapTableBase_1 <- readRDS(file = paste0(homeDirectory, "exp1/overlapTableBase_observed.rds")) %>% mutate(rep = "1")
overlapTableBase_2 <- readRDS(file = paste0(homeDirectory, "exp3/overlapTableBase_observed.rds")) %>% mutate(rep = "2")
overlapTableBase <- bind_rows(overlapTableBase_1, overlapTableBase_2)
ggplot(overlapTableBase, aes(x = cond, y = overlap, fill = rep)) +
  geom_col(stat = "identity", position = "dodge") +
  theme(axis.title.x = element_blank(), legend.position = "none") + ylim(0, 0.75)
ggsave(filename = paste0(plotDirectory, "overlap_base.pdf"), height = 2, width = 4)

ggplot(overlapData, aes(x = cond, y = overlap, fill = rep)) +
  geom_col(stat = "identity", position = "dodge") +
  theme(axis.title.x = element_blank(), legend.position = "none") + ylim(0, 0.75)
ggsave(filename = paste0(plotDirectory, "overlap_norm_observed.pdf"), height = 2, width = 4)

overlapTableNorm_same_1 <- readRDS(file = paste0(homeDirectory, "exp1/overlapTableNorm_same.rds")) %>% mutate(rep = "1")
overlapTableNorm_same_2 <- readRDS(file = paste0(homeDirectory, "exp3/overlapTableNorm_same.rds")) %>% mutate(rep = "2")
overlapTableNorm <- bind_rows(overlapTableNorm_same_1, overlapTableNorm_same_2)
ggplot(overlapTableNorm, aes(x = cond, y = overlap, fill = rep)) +
  geom_col(stat = "identity", position = "dodge") +
  theme(axis.title.x = element_blank(), legend.position = "none") + ylim(0, 0.75)
ggsave(filename = paste0(plotDirectory, "overlap_norm_same.pdf"), height = 2, width = 4)
