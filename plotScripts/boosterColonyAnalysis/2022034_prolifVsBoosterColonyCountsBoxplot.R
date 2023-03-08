library(tidyverse)
library(ggsignif)

dataDirectory = "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/extractedData/boosterColonyCounts/"
plotDirectory = "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/plots/boosterColonyAnalysis/"

data_R1 <- read.csv(file = "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/rawData/proliferationSort/colonyCounts/prolifSortwLSD1i/20220329_prolifVsBoosters/results.csv")
data_R2 <- read.csv(file = "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/rawData/proliferationSort/colonyCounts/prolifSortwLSD1i/20220808_prolifVsBoosters/results.csv")

data_R1 <- data_R1 %>% mutate(cond = gsub("(.*)\\_(.*)\\_(.*).tif", "\\1", .$Dataset)) %>%
  mutate(prolif = gsub("(.*)\\_(.*)\\_(.*).tif", "\\2", .$Dataset)) %>%
  mutate(rep = gsub("(.*)\\_(.*)\\_(.*).tif", "\\3", .$Dataset))
baseline_R1 <- mean(data_R1 %>% dplyr::filter(cond == "control", prolif == "C") %>% .$NumBlobs)
data_R1 <- data_R1 %>% rowwise() %>% mutate(normBlob = NumBlobs / baseline_R1) %>% ungroup()

data_R2 <- data_R2 %>% mutate(cond = gsub("(.*)\\_(.*)\\_(.*)\\_(.*).tif", "\\2", .$Dataset)) %>%
  mutate(prolif = gsub("(.*)\\_(.*)\\_(.*)\\_(.*).tif", "\\3", .$Dataset)) %>%
  mutate(rep = gsub("(.*)\\_(.*)\\_(.*)\\_(.*).tif", "\\4", .$Dataset))
baseline_R2 <- mean(data_R2 %>% dplyr::filter(cond == "control", prolif == "C") %>% .$NumBlobs)
data_R2 <- data_R2 %>% rowwise() %>% mutate(normBlob = NumBlobs / baseline_R2) %>% ungroup()

data <- bind_rows(data_R2 %>% dplyr::select(-Dataset) %>% dplyr::mutate(exp = "1"),
                  data_R1 %>% dplyr::select(-Dataset) %>% dplyr::mutate(exp = "2"))
data$cond <- ifelse(data$cond == "LSD1", "LSD1i", data$cond)

data$prolif <- factor(data$prolif, levels = c("H", "MH", "C", "ML", "L"), labels = c("slow", "mid slow", "control", "mid fast", "fast"))
data <- data %>% dplyr::filter(prolif %in% c("fast", "control", "slow"))

ggplot(data, aes(x = cond, y = normBlob)) +
  facet_wrap(~prolif) +
  stat_summary(fun = mean, geom = "col") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.25) +
  geom_signif(comparisons = list(c("control", "LSD1i")), test = "t.test", y_position = 12.5) +
  geom_hline(yintercept = 1, linetype = "dashed")
ggsave(filename = paste0(plotDirectory, "prolifVsBoostersReprogrammingRate.pdf"), unit = "in", height = 3, width = 4)

foldchangeTable <- data %>% dplyr::select(cond, prolif, rep, normBlob, exp) %>%
  dcast(prolif + exp ~ cond, value.var = "normBlob", fun.aggregate = mean) %>%
  rowwise() %>% mutate(foldchange = (LSD1i + 0.5) / (control + 0.5)) %>% ungroup()

ggplot(foldchangeTable, aes(x = prolif, y = foldchange)) +
  stat_summary(fun = mean, geom = "col") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.25) +
  stat_signif(comparisons = list(c("slow", "control"), c("control", "fast"), c("slow", "fast")))
