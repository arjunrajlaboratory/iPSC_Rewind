library(tidyverse)
theme_set(theme_classic())

plotDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Revisions/plots/CRISPRLinesColonies/"

data <- read.csv(file = "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Revisions/extractedData/CRISPRLineColonies/CRISPRKDColoniesAggregated.csv")

for(i in 1:length(unique(data$round))){
  data_round <- data %>% dplyr::filter(round == unique(data$round)[i])
  for(j in 1:length(unique(data_round$plate))){
   data_plate <- data_round %>% dplyr::filter(plate == unique(data_round$plate)[j])
   control_value <- data_plate %>% dplyr::filter(gene == "control") %>% .$colonies %>% mean()
   data_plate$colonies<- data_plate$colonies/control_value
   data_plate_comb <- data_plate %>% group_by(gene, guide) %>% summarize(mean_colonies = mean(colonies)) %>% mutate(round = i)
   if(i == 1 & j == 1) {
     data_norm <- data_plate_comb
   } else{
     data_norm <- bind_rows(data_norm, data_plate_comb)
   }
  }
}

data_norm$name <- paste0(data_norm$gene, "_", data_norm$guide)
data_norm <- data_norm %>% group_by(gene, guide, round) %>% summarize(mean_colonies = mean(mean_colonies))
data_norm$gene <- factor(data_norm$gene, levels = c("control", "KDM1A", "DOT1L", "CDKN1A", "MDM2", "FTH1", "GDF15", "SPP1", "SQSTM1", "NFE2L2", "MYBL2"))
data_norm$guide <- factor(data_norm$guide)

ggplot(data_norm, aes(x = guide, y = mean_colonies)) +
  stat_summary(fun = mean, geom = "col", fill = "gray75") +
  #stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.25) +
  geom_jitter(width = 0.1, height = 0) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  facet_grid(~gene, scales = "free_x", space = "free")

ggsave(filename = paste0(plotDirectory, "CRISPR_hiFT_normColonyCounts.pdf"), height = 3, width = 15)
