library(tidyverse)
library(reshape2)
library(VennDiagram)

theme_set(theme_classic())

dataDirectory <- '/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/extractedData/barcodeOverlap/'
plotDirectory <- '/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/plots/barcodeOverlap/'

overlapValListSim1 <- readRDS(file = paste0(dataDirectory, 'simulationResults1in10000.RDS'))
overlapValListSim2 <- readRDS(file = paste0(dataDirectory, 'simulationResults1in1000.RDS'))
overlapValListSim3 <- readRDS(file = paste0(dataDirectory, 'simulationResults1in100.RDS'))

overlapValPlotTable <- bind_rows(overlapValListSim1 %>% mutate(label = "0.01% efficiency"),
                                 overlapValListSim2 %>% mutate(label = "0.1% efficiency"),
                                 overlapValListSim3 %>% mutate(label = "1% efficiency")) %>%
  rename("(1)" = X3, "(2)" = X1, "(3)" = X2) %>% melt(id.vars = c("label"))
overlapValPlotTable$variable <- factor(overlapValPlotTable$variable, levels = c("(1)", "(2)", "(3)"))

ggplot(overlapValPlotTable, aes(x = label, y = value)) +
  geom_boxplot(fill = NA) +
  geom_jitter(shape = 16) +
  facet_wrap(~variable) +
  #stat_summary(fun = 'mean', geom = 'crossbar') +
  ylim(0, 1) + ylab('percent of barcode overlap') + theme(axis.title.x = element_blank()) + theme_classic()

interceptDummy <- data.frame(variable = c("(1)", "(2)", "(3)"), intercept = c(7.56, 26.1, 37.1))

plot <- ggplot(overlapValPlotTable, aes(x = value * 100, fill = label)) +
  geom_histogram(alpha = 0.5, position = "identity", binwidth = 0.2) + facet_wrap(~variable, nrow = 3, scales = "fixed") +
  geom_vline(data = interceptDummy, aes(xintercept = intercept)) +
  xlab("percent of overlapping barcodes between split A and split B") + ylab("frequency")
ggsave(plot, file = paste0(plotDirectory, 'simulationDistributionPlots.pdf'), width = 6, height = 9)

meanOverlap1 <- overlapValPlotTable %>% filter(label == "1% efficiency") %>% summarise(mean = mean(value)) %>% .$mean
grid.newpage()
venndiagram <- draw.pairwise.venn(area1 = 0.5, area2 = 0.5, cross.area = meanOverlap1/(1+meanOverlap1),
                                  euler.d = TRUE, fill = c(alpha("gray5", 0.5), alpha("gray5", 0.5)),
                                  fontfamily = "sans", cex = 2.5, print.mode = "percent", sigdigs = 3)
ggsave(plot = venndiagram, filename = paste0(plotDirectory, 'overlapVennDiagram_1in100.pdf'), width = 5, height = 5, useDingbats = F)

meanOverlap2 <- overlapValPlotTable %>% filter(label == "0.1% efficiency") %>% summarise(mean = mean(value)) %>% .$mean
grid.newpage()
venndiagram <- draw.pairwise.venn(area1 = 0.5, area2 = 0.5, cross.area = meanOverlap2/(1+meanOverlap2),
                                  euler.d = TRUE, fill = c(alpha("gray5", 0.5), alpha("gray5", 0.5)),
                                  fontfamily = "sans", cex = 2.5, print.mode = "percent", sigdigs = 3)
ggsave(plot = venndiagram, filename = paste0(plotDirectory, 'overlapVennDiagram_1in1000.pdf'), width = 5, height = 5, useDingbats = F)

meanOverlap3 <- overlapValPlotTable %>% filter(label == "0.01% efficiency") %>% summarise(mean = mean(value)) %>% .$mean
grid.newpage()
venndiagram <- draw.pairwise.venn(area1 = 0.5, area2 = 0.5, cross.area = meanOverlap3/(1+meanOverlap3),
                                  euler.d = TRUE, fill = c(alpha("gray5", 0.5), alpha("gray5", 0.5)),
                                  fontfamily = "sans", cex = 2.5, print.mode = "percent", sigdigs = 3)
ggsave(plot = venndiagram, filename = paste0(plotDirectory, 'overlapVennDiagram_1in10000.pdf'), width = 5, height = 5, useDingbats = F)