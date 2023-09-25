library(tidyverse)
theme_set(theme_classic())

homeDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Revisions/extractedData/prolifSort/"
plotDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Revisions/plots/prolifSort/"

homeDirectory <- "/Users/naveenj/Dropbox (RajLab)/Shared_Naveen/Revisions/extractedData/prolifSort/"
plotDirectory <- "/Users/naveenj/Dropbox (RajLab)/Shared_Naveen/Revisions/plots/prolifSort/"

prolifTable <- read.csv(file = paste0(homeDirectory, "prolifSort_MEF_prolif_aggregate.csv"))
prolifTable$condition <- factor(prolifTable$condition, levels = c("slow", "ungated", "fast"))

ggplot(prolifTable %>% dplyr::filter(!(condition == "ungated")), aes(x = condition, y = normfoldchange)) +
  stat_summary(fun = mean, geom = "bar") +
  stat_summary(fun = mean, geom = "point", size = 5) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.25) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  ylab("normalized proliferation rate") + theme(axis.title.x = element_blank())
ggsave(filename = paste0(plotDirectory, "prolifSort_MEF_prolifRate.pdf"), height = 3, width = 4.5)

colonyTable <- read.csv(file = paste0(homeDirectory, "prolifSort_MEF_colonies_aggregate.csv"))
colonyTableGroup <- colonyTable %>% group_by(exp, condition) %>% summarise(meanColonies = mean(colonies))
colonyTableGroup$normColonies <- colonyTableGroup$meanColonies/prolifTable$d1

for(i in 1:length(colonyTableGroup$exp %>% unique())) {
  expColonyTable <- colonyTableGroup %>% dplyr::filter(exp == (colonyTableGroup$exp %>% unique())[i])
  controlMean <- dplyr::filter(expColonyTable, condition == "ungated") %>% .$meanColonies %>% mean()
  expColonyTable$relColonies <- expColonyTable$meanColonies / controlMean
  
  if(i == 1) {
    colonyTableRel <- expColonyTable
  } else {
    colonyTableRel <- bind_rows(colonyTableRel, expColonyTable)
  }
}
colonyTableRel$condition <- factor(colonyTableRel$condition, levels = c("slow", "ungated", "fast"))

ggplot(colonyTableRel %>% dplyr::filter(!(condition == "ungated")), aes(x = condition, y = relColonies)) +
  stat_summary(fun = mean, geom = "bar") +
  stat_summary(fun = mean, geom = "point", size = 5) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.25) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  ylab("normalized reprogramming rate") + theme(axis.title.x = element_blank())
ggsave(filename = paste0(plotDirectory, "prolifSort_MEF_reprogRate.pdf"), height = 3, width = 4.5)
