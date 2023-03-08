library(tidyverse)
library(ggsignif)

dataDirectory = "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/extractedData/boosterColonyCounts/"
plotDirectory = "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/plots/boosterColonyAnalysis/"

data <- read.csv(file = paste0(dataDirectory, "manualColonYCounts.csv")) %>%
  dplyr::select(experiment, condition, fold.change.average) %>%
  filter(!is.na(fold.change.average))

data$condition <- factor(data$condition, levels = c("DMSO", "LSD1i", "DOT1Li"))

plot <- ggplot(data %>% filter(condition %in% c("LSD1i")), aes(x = condition, y = fold.change.average)) +
  geom_jitter(height = 0.1, width = 0.1) + geom_hline(yintercept = 1, linetype = "dashed") +
  stat_summary(fun = mean, geom = "point", width = 0.1, color = "blue", size = 5) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.1, color = "blue", size = 1) +
  theme_classic(base_size = 15) + theme(axis.title.x = element_blank()) + ylab("fold change in number of colonies (LSD1i/DMSO)") + ylim(0, 8)

ggsave(plot, file = paste0(plotDirectory, "log2FCColonyCountBoxplot_boosters.pdf"), height = 4, width = 2)