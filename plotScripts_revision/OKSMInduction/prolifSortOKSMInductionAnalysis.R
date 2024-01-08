rm(list=ls())
gc()

library(tidyverse)
library(reshape2)
library(ggrepel)
library(ggsignif)

theme_set(theme_classic())

homeDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Revised Manuscript/rawData/OKSMInduction/smFISH/"
plotDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Revised Manuscript/plots/OKSMInduction/"

filePaths <- list.files(homeDirectory, recursive = TRUE, pattern = "spotsSummary.csv", full.names = TRUE)

for(i in 1:length(filePaths)) {
  spotTableTemp <- read.csv(file = filePaths[i], header = TRUE)
  if (i == 1) {
    spotTable <- spotTableTemp
  } else {
    spotTable <- bind_rows(spotTable, spotTableTemp)
  }
}

spotTableCast <- dcast(data = spotTable, formula = nucID + x + y + Well + Condition + Time + Rep ~ channel, value.var = "GroupCount") %>%
  dplyr::rename(UBC = alexa, SOX2 = tmr, OCT4 = cy) %>% dplyr::filter(Rep != 3)

spotTableCast$Condition <- factor(spotTableCast$Condition, levels = c("slow", "ungated", "fast"))

spotTableCast <- spotTableCast %>% dplyr::filter(!(Time != 0 & (UBC < 10)))

ggplot(spotTableCast, aes(x = Condition, y = OCT4)) +
  geom_violin(draw_quantiles = c(0.50), scale = "width", aes(fill = Condition), trim = FALSE) +
  geom_hline(yintercept = 10, linetype = "dashed") +
  facet_grid(Rep ~ Time) +
  stat_summary(fun.y = mean, geom = "point", size = 2.5) +
  stat_summary(fun.y = median, geom = "text", aes(label = round(..y.., 0)), position = position_nudge(y = Inf), vjust = 1) +
  # geom_signif(comparisons = list(c("slow", "ungated"), c("ungated", "fast"), c("slow", "fast")), step_increase = 0.075, y_position = 500, map_signif_level = TRUE) +
  theme(axis.title.x = element_blank(), legend.position = "none")
ggsave(filename = paste0(plotDirectory, "OKSMInduction_lowThreshold.pdf"), height = 6, width = 8)

spotTableCast <- spotTableCast %>% dplyr::filter(!(Time != 0 & (UBC < 100)))

ggplot(spotTableCast, aes(x = Condition, y = OCT4)) +
  geom_violin(draw_quantiles = c(0.50), scale = "width", aes(fill = Condition), trim = FALSE) +
  geom_hline(yintercept = 10, linetype = "dashed") +
  facet_grid(Rep ~ Time) +
  stat_summary(fun.y = mean, geom = "point", size = 2.5) +
  stat_summary(fun.y = median, geom = "text", aes(label = round(..y.., 0)), position = position_nudge(y = Inf), vjust = 1) +
  # geom_signif(comparisons = list(c("slow", "ungated"), c("ungated", "fast"), c("slow", "fast")), step_increase = 0.075, y_position = 500, map_signif_level = TRUE) +
  theme(axis.title.x = element_blank(), legend.position = "none")
ggsave(filename = paste0(plotDirectory, "OKSMInduction_highThreshold.pdf"), height = 6, width = 8)
