library(tidyverse)
library(ggsignif)
library(ggforce)
library(ggridges)
library(egg)
library(Seurat)
library(PupillometryR)

theme_set(theme_classic())

dataDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/rawData/CRISPR/R2_smFISH/UBC A594, KDM1A CY3/"
plotDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/plots/CRISPR/R2/"

folders <- list.dirs(path = dataDirectory, full.names = FALSE)[-1]

for(i in 1:length(folders)){
  if(i == 1) {
    spotTable <- read_csv(file = paste0(dataDirectory, folders[i], "/spotsSummary.csv")) %>% mutate(condition = rep(folders[i] %>% sub("^[^_]*_", "", .), nrow(.)))
  }
  spotTableTemp <- read_csv(file = paste0(dataDirectory, folders[i], "/spotsSummary.csv")) %>% mutate(condition = rep(folders[i] %>% sub("^[^_]*_", "", .), nrow(.)))
  spotTable <- bind_rows(spotTable, spotTableTemp)
}

spotTable <- spotTable %>% unique()
spotTable %>% group_by(condition) %>% summarize(number = n())
spotTable <- spotTable %>% mutate(label = condition)
spotTable$label <- ifelse(spotTable$condition %in% c("B1_1", "B1_2", "B2_1", "B2_2"), "control", spotTable$label)

spotTablePlot <- spotTable %>% filter(channel == "tmr") %>% filter(!(condition %in% c("B2_1", "B2_2", "C2", "C4")))
means <- aggregate(GroupCount ~ label, spotTablePlot, mean)
means$GroupCount <- round(means$GroupCount, digits = 2)
ggplot(spotTablePlot, aes(x = label, y = GroupCount)) +
  geom_flat_violin(aes(fill = label), position = position_nudge(x = 0.2, y = 0)) +
  geom_jitter(aes(fill = label), width = 0.15, alpha = 0.1, pch = 21, color = "black") +
  geom_text(data = means, aes(x = label, y = GroupCount, label = round(GroupCount, 2)), vjust = -2, hjust = -0.25, size = 3, angle = 90, color = "blue") +
  geom_errorbar(data = means, aes(x = label, ymax = GroupCount, ymin = GroupCount), color = "blue") + NoLegend() +
  geom_signif(comparisons = list(c("control", "KDM1A_1"), c("control", "KDM1A_2"), c("control", "KDM1A_3"), c("control", "KDM1A_4")), step_increase = 0.075, map_signif_level = TRUE) +
  theme(axis.title.x = element_blank()) + ylab("RNA counts per cell") +
  scale_fill_manual(values = c(rep("grey", 1), rep("red", 4)))

plot1 <- ggplot(spotTablePlot, aes(x = label, y = GroupCount, fill = label)) +
  geom_hline(yintercept = means %>% filter(label == "control") %>% .$GroupCount, linetype = "dashed") +
  geom_violin(draw_quantiles = TRUE) + NoLegend() +
  geom_signif(comparisons = list(c("control", "KDM1A_1"), c("control", "KDM1A_2"), c("control", "KDM1A_3"), c("control", "KDM1A_4")), step_increase = 0.075, map_signif_level = TRUE) +
  theme(axis.title.x = element_blank()) + ylab("RNA counts per cell") +
  scale_fill_manual(values = c(rep("#D8D8D8", 1), rep("#ED1C24", 4)))

plot2 <- ggplot(spotTablePlot, aes(x = label, y = GroupCount, fill = label)) +
  stat_summary(fun = "mean", geom = "col") +
  geom_hline(yintercept = means %>% filter(label == "control") %>% .$GroupCount, linetype = "dashed") +
  stat_summary(fun = "mean", geom = "point", size = 2.5, color = "black") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.25, color = "black") +
  # geom_text(data = means, aes(x = label, y = GroupCount, label = round(GroupCount, 2)), vjust = 0.5, hjust = 0, angle = 45, color = "blue", nudge_x = 0.1, nudge_y = 1) +
  NoLegend() + theme(axis.title.x = element_blank()) + ylab("mean RNA counts per cell") +
  scale_fill_manual(values = c(rep("#D8D8D8", 1), rep("#ED1C24", 4)))
ggsave(filename = paste0(plotDirectory, "smFISHCheckGraphs_barOnly_R2_LSD1.pdf"), plot = plot2, units = "in", height = 3, width = 4.5, useDingbats = FALSE)

plot <- egg::ggarrange(plots = list(plot2, plot1), ncol = 2)
ggsave(filename = paste0(plotDirectory, "smFISHCheckGraphs_R2_LSD1.pdf"), plot = plot, units = "in", height = 3, width = 6, useDingbats = FALSE)

spotTablePlot <- spotTable %>% filter(channel == "tmr")
means <- aggregate(GroupCount ~ condition, spotTablePlot, mean)
means$GroupCount <- round(means$GroupCount, digits = 2)
ggplot(spotTablePlot, aes(x = condition, y = GroupCount)) +
  geom_flat_violin(aes(fill = condition), position = position_nudge(x = 0.2, y = 0)) +
  geom_jitter(aes(fill = condition), width = 0.15, alpha = 0.1, pch = 21, color = "black") +
  geom_text(data = means, aes(x = condition, y = GroupCount, label = round(GroupCount, 2)), vjust = -2, hjust = -0.25, size = 3, angle = 90, color = "blue") +
  geom_errorbar(data = means, aes(x = condition, ymax = GroupCount, ymin = GroupCount), color = "blue") + NoLegend() +
  theme(axis.title.x = element_blank()) + ylab("RNA counts per cell") +
  scale_fill_manual(values = c(rep("grey", 6), rep("seagreen3", 4)))

spotTablePlot <- spotTable %>% filter(channel == "alexa")
means <- aggregate(GroupCount ~ condition, spotTablePlot, mean)
means$GroupCount <- round(means$GroupCount, digits = 2)
ggplot(spotTablePlot, aes(x = condition, y = GroupCount)) +
  geom_flat_violin(aes(fill = condition), position = position_nudge(x = 0.2, y = 0)) +
  geom_jitter(aes(fill = condition), width = 0.15, alpha = 0.1, pch = 21, color = "black") +
  geom_text(data = means, aes(x = condition, y = GroupCount, label = round(GroupCount, 2)), vjust = -2, hjust = -0.25, size = 3, angle = 90, color = "blue") +
  geom_errorbar(data = means, aes(x = condition, ymax = GroupCount, ymin = GroupCount), color = "blue") + NoLegend() +
  theme(axis.title.x = element_blank()) + ylab("RNA counts per cell") +
  scale_fill_manual(values = c(rep("grey", 4), rep("seagreen3", 4)))

library(reshape2)
spotTableNorm <- spotTable %>% dcast(value.var = "GroupCount", formula = condition + nucID + x + y ~ channel, fun.aggregate = mean) %>%
  rowwise() %>% mutate(normCount = (cy + 1) / (alexa + 1)) %>% ungroup()
means <- aggregate(normCount ~ condition, spotTableNorm, mean)
means$normCount <- round(means$normCount, digits = 2)
ggplot(spotTableNorm, aes(x = condition, y = normCount)) +
  geom_flat_violin(aes(fill = condition), position = position_nudge(x = 0.2, y = 0)) +
  geom_jitter(aes(fill = condition), width = 0.15, alpha = 0.1, pch = 21, color = "black") +
  geom_text(data = means, aes(x = condition, y = normCount, label = round(normCount, 2)), vjust = -2, hjust = -0.25, size = 3, angle = 90, color = "blue") +
  geom_errorbar(data = means, aes(x = condition, ymax = normCount, ymin = normCount), color = "blue") + NoLegend() +
  theme(axis.title.x = element_blank()) + ylab("RNA counts per cell") +
  scale_fill_manual(values = c(rep("grey", 4), rep("seagreen3", 4)))
