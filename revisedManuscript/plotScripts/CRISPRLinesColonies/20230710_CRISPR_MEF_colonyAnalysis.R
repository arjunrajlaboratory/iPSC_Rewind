library(tidyverse)
theme_set(theme_classic())

plotDirectory <- "/Users/naveenj/Dropbox (RajLab)/Shared_Naveen/Revisions/plots/CRISPRLinesColonies/"

r1p1_data <- read.csv(file = "/Users/naveenj/Dropbox (RajLab)/Shared_Naveen/Revisions/rawData/CRISPRKnockdownReprog/MEF/20230709_CRISPR_MEF_round1/plate1/CRISPR_MEF_r1_p1_keep.csv")
r1p2_data <- read.csv(file = "/Users/naveenj/Dropbox (RajLab)/Shared_Naveen/Revisions/rawData/CRISPRKnockdownReprog/MEF/20230709_CRISPR_MEF_round1/plate2/CRISPR_MEF_r1_p2_keep.csv")

r1_data <- bind_rows(r1p1_data, r1p2_data)
r2_data <- read.csv(file = "/Users/naveenj/Dropbox (RajLab)/Shared_Naveen/Revisions/rawData/CRISPRKnockdownReprog/MEF/20230730_CRISPR_MEF_round2/results.csv")
r3_data <- read.csv(file = "/Users/naveenj/Dropbox (RajLab)/Shared_Naveen/Revisions/rawData/CRISPRKnockdownReprog/MEF/20230821_CRISPR_MEF_round3/CRISPR_MEF_round3.csv")
r4_data <- read.csv(file = "/Users/naveenj/Dropbox (RajLab)/Shared_Naveen/Revisions/rawData/CRISPRKnockdownReprog/MEF/20230821_CRISPR_MEF_round4/CRISPR_MEF_round4.csv")
r5_data <- read.csv(file = "/Users/naveenj/Dropbox (RajLab)/Shared_Naveen/Revisions/rawData/CRISPRKnockdownReprog/MEF/20230821_CRISPR_MEF_round5/CRISPR_MEF_round5.csv")

data <- bind_rows(r1_data, r2_data, r3_data, r4_data, r5_data)

data_metadata <- data$Dataset %>% strsplit("_")
data_well <- lapply(data_metadata, function(lst) lst[[1]]) %>% unlist()
data_round <- lapply(data_metadata, function(lst) lst[[2]]) %>% unlist()
data_plate <- lapply(data_metadata, function(lst) lst[[3]]) %>% unlist()
data_gene <- lapply(data_metadata, function(lst) lst[[4]]) %>% unlist()
data_guide <- lapply(data_metadata, function(lst) lst[[5]]) %>% unlist()

data <- data %>% mutate(Round = data_round,
                        Well = data_well,
                        Plate = data_plate,
                        Gene = data_gene,
                        Guide = data_guide)

for(i in 1:length(unique(data$Round))){
  data_temp <- data %>% dplyr::filter(Round == unique(data$Round)[i])
  if(length(unique(data_temp$Plate)) > 1) {
    for(j in 1:length(unique(data_temp$Plate))) {
      data_plate_temp <- data_temp %>% dplyr::filter(Plate == unique(data_temp$Plate)[j])
      control_value <- data_plate_temp %>% dplyr::filter(Gene == "control") %>% .$NumBlobs %>% mean()
      data_plate_temp$NumBlobsNorm <- data_plate_temp$NumBlobs/control_value
      
      if(j == 1) {
        data_norm <- data_plate_temp 
      } else{
        data_norm <- bind_rows(data_norm, data_plate_temp)
      } 
    }
    
    data_norm <- data_norm %>% group_by(Round, Gene, Guide) %>% summarize(NumBlobsNorm = mean(NumBlobsNorm))
  }
  
  else{
    control_value <- data_temp %>% dplyr::filter(Gene == "control") %>% .$NumBlobs %>% mean()
    data_temp$NumBlobsNorm <- data_temp$NumBlobs/control_value
    data_temp <- data_temp %>% group_by(Round, Gene, Guide) %>% summarize(NumBlobsNorm = mean(NumBlobsNorm))
    
    if(i == 1) {
      data_norm <- data_temp 
    } else{
      data_norm <- bind_rows(data_norm, data_temp)
    }
  }
}

data_norm$Target <- paste0(data_norm$Gene, "_", data_norm$Guide)

data_norm$Gene <- factor(data_norm$Gene, levels = c("control", "KDM1A", "SPP1", "FTH1", "SQSTM1"))

ggplot(data_norm, aes(x = Target, y = NumBlobsNorm)) +
  stat_summary(fun = mean, geom = "col") +
  # stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.25) +
  geom_jitter(width = 0.1, height = 0) +
  facet_grid(~Gene, scales = "free_x", space = "free") +
  geom_hline(yintercept = 1, linetype = "dashed") + ylab("normalized colony counts") +
  theme(axis.title.x = element_blank())

ggsave(filename = paste0(plotDirectory, "CRISPR_MEF_normColonyCounts.pdf"), height = 2.5, width = 8)
