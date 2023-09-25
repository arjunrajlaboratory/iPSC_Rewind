library(tidyverse)
library(reshape2)

theme_set(theme_classic())

plotDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Revisions/plots/prolifReversion/"

prolifPath <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Revisions/extractedData/prolifRevesion/cyclingReversionAggregated_prolifRate.csv"
prolifTable <- read.csv(prolifPath, header = TRUE)
prolifTableMelt <- melt(prolifTable, id.vars = c("exp", "collection"))

ggplot(prolifTableMelt, aes(x = collection, y = value, group = variable, color = variable)) +
  geom_point() +
  geom_line() +
  ylim(0, 2.1) +
  facet_grid(~exp) +
  xlab("days after sorting by cycling speed") + ylab("normalized proliferation rate") +
  theme(legend.position = "none")
ggsave(filename = paste0(plotDirectory, "prolifReversion_prolifRate.pdf"), height = 2, width = 5)

reprogPath <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Revisions/extractedData/prolifRevesion/cyclingReversionAggregated_reprogRate.csv"
reprogTable <- read.csv(reprogPath, header = TRUE)
reprogTableMelt <- melt(reprogTable, id.vars = c("exp", "collection"))

ggplot(reprogTableMelt, aes(x = collection, y = value, group = variable, color = variable)) +
  geom_point() +
  geom_line() +
  ylim(0, 2.1) +
  facet_grid(~exp) +
  xlab("days after sorting by cycling speed") + ylab("normalized reprogramming rate") +
  theme(legend.position = "none")
ggsave(filename = paste0(plotDirectory, "prolifReversion_reprogRate.pdf"), height = 2, width = 5)
