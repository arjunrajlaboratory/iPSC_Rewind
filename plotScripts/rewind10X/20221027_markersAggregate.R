rm(list=ls())
gc()

library(tidyverse)

theme_set(theme_classic())

homeDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/extractedData/rewind10X/"
plotDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/plots/rewind10X/"

markers_R1 <- readRDS("/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/extractedData/rewind10X/R1/primedMarkersAll.rds") %>% mutate(rep = "1")
markers_R2 <- readRDS("/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/extractedData/rewind10X/R2/primedIndMarkersAll.rds") %>% mutate(rep = "2")
markers_R3 <- readRDS("/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/extractedData/rewind10X/R3/primedMarkersAll.rds") %>% mutate(rep = "3")

# markersComb <- bind_rows(markers_R1, markers_R2, markers_R3) %>% group_by(gene) %>% summarise(mean = mean(avg_log2FC))
overlap1 <- inner_join(markers_R1, markers_R2, by = "gene")$gene
overlap2 <- inner_join(markers_R1, markers_R3, by = "gene")$gene
overlap3 <- inner_join(markers_R2, markers_R3, by = "gene")$gene
overlapAll <- bind_rows(markers_R1, markers_R2, markers_R3) %>% filter(gene %in% overlap1 | gene %in% overlap2 | gene %in% overlap3)
markersComb <- overlapAll %>% group_by(gene) %>% summarise(mean = mean(avg_log2FC))
saveRDS(markersComb, file = paste0(homeDirectory, "markersComb.rds"))

posInclude <- markersComb %>% slice_max(., order_by = mean, n = 25) %>% .$gene
negInclude <- markersComb %>% slice_min(., order_by = mean, n = 25) %>% .$gene
controlInclude <- c("ACTB", "GAPDH", "PGK1", "UBC")

ggplot(overlapAll %>% dplyr::filter(gene %in% posInclude), aes(x = reorder(gene, avg_log2FC, mean), y = avg_log2FC)) +
  # geom_point(aes(color = rep)) +
  geom_text(aes(label = gene, y = 0), angle = 90, size = 2.5) +
  stat_summary(fun = mean, geom = "crossbar", fill = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.25) +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

ggplot(overlapAll %>% dplyr::filter(gene %in% negInclude), aes(x = reorder(gene, avg_log2FC, mean), y = avg_log2FC)) +
  # geom_point(aes(color = rep)) +
  geom_text(aes(label = gene, y = 0), angle = 90, size = 2.5) +
  stat_summary(fun = mean, geom = "crossbar", fill = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.25) +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

ggplot(overlapAll %>% dplyr::filter(gene %in% controlInclude), aes(x = reorder(gene, avg_log2FC, mean), y = avg_log2FC)) +
  # geom_point(aes(color = rep)) +
  geom_text(aes(label = gene, y = 0), angle = 90, size = 2.5) +
  stat_summary(fun = mean, geom = "crossbar", fill = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.25) +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

ggplot(overlapAll %>% dplyr::filter(gene %in% c(negInclude, controlInclude, posInclude)), aes(x = reorder(gene, avg_log2FC, mean), y = avg_log2FC)) +
  # geom_point(aes(color = rep)) +
  # geom_text(aes(label = gene, y = 0), angle = 90, size = 2.5, hjust = 1) +
  stat_summary(fun = mean, geom = "point", fill = "black", size = 2.5) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.25) +
  stat_summary(fun = mean, geom = "text", aes(label = gene), angle = -45, hjust = 0, vjust = 0, size = 3, position = position_nudge(x = 0.1, y = -0.05)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())
ggsave(filename = paste0(plotDirectory, "markersAggregatePlot.pdf"), units = "in", height = 2, width = 10)
