rm(list=ls())
gc()

library(tidyverse)
library(metafor)

theme_set(theme_classic())

homeDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/extractedData/rewind10X/"
plotDirectory <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/plots/rewind10X/"

oddsTable_R1 <- readRDS("/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/extractedData/rewind10X/R1/R1_oddsRatioListFilter.rds")
oddsTable_R2 <- readRDS("/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/extractedData/rewind10X/R2/R2_oddsRatioListFilter.rds")
oddsTable_R3 <- readRDS("/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/extractedData/rewind10X/R3/R3_oddsRatioListFilter.rds")

#### combine samples by using a random effects model for OR ####
#######################################################################################################################################################
geneList <- c("CENPF", "MKI67", "TOP2A", "SOX21", "SPP1", "GDF15", "CDKN1A", "FTH1", "GAPDH", "UBC")
propList <- c(0.01, 0.05, 0.10, 0.15, 0.25)

for(i in 1:length(geneList)) {
  for(j in 1:length(propList)) {
    orTemp1 <- oddsTable_R1 %>% dplyr::filter(gene == geneList[i], prop == propList[j])
    orTemp2 <- oddsTable_R2 %>% dplyr::filter(gene == geneList[i], prop == propList[j])
    orTemp3 <- oddsTable_R3 %>% dplyr::filter(gene == geneList[i], prop == propList[j])
    
    orTempAll <- bind_rows(orTemp1, orTemp2, orTemp3)
    
    combORTemp <- as_tibble(predict(rma(yi = log_estimate, sei = se, data = orTempAll)))
    combORTemp <- combORTemp %>% mutate(gene = geneList[i], prop = propList[j])
    
    if(i == 1 & j == 1) {
      combORList <- combORTemp
    } else{
      combORList <- bind_rows(combORList, combORTemp)
    }
  }
}

combORList$prop <- factor(combORList$prop, levels = propList)
combORList$gene <- factor(combORList$gene, levels = geneList)

ggplot(combORList, aes(x = prop, y = pred)) +
  geom_point(size = 5) +
  geom_errorbar(aes(ymin = pred - se, ymax = pred + se), width = 0.25) +
  # geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub), width = 0.25) +
  facet_wrap(~gene, ncol = 4) +
  geom_hline(yintercept = 0, linetype = "dashed")

combORListFilter <- combORList %>% filter(prop == 0.1) %>%
  rowwise() %>% mutate(pval1 = pnorm(abs(abs(pred) - abs(0.3109228))/sqrt(se^2 + 0.1581213^2),
                                     mean = 0, sd = 1, lower.tail = FALSE)* 2) %>% # comparing with fast
  rowwise() %>% mutate(pval2 = pnorm(abs(abs(pred) - abs(-0.3909865))/sqrt(se^2 + 0.1816956^2),
                                     mean = 0, sd = 1, lower.tail = FALSE)* 2) # comparing with slow

prolifORList <- oddsTable_R3 %>% filter(gene %in% c("S2", "S3"))
prolifORList$gene <- c("fast", "slow")
prolifORList <- prolifORList %>% dplyr::rename(pred = log_estimate, se = se) %>% mutate(ci.lb = NA, ci.ub = NA, pi.lb = NA, pi.ub = NA, pval1 = NA, pval2 = NA)

plotORList <- bind_rows(combORListFilter, prolifORList)
plotGeneList <- c("fast", "CENPF", "MKI67", "TOP2A", "SOX21", "slow", "SPP1", "GDF15", "CDKN1A", "FTH1", "GAPDH", "UBC")
plotORList$gene <- factor(plotORList$gene, levels = plotGeneList)
plotORList$pval <- ifelse(plotORList$gene %in% plotGeneList[1:5], plotORList$pval1, plotORList$pval2)

ggplot(plotORList %>% dplyr::filter(gene %in% plotGeneList[1:5]), aes(x = gene, y = pred)) +
  geom_errorbar(aes(ymin = pred - se, ymax = pred + se), width = 0.25) +
  geom_text(aes(label = round(pval, 2)), hjust = -1) +
  geom_point(size = 5) + ylim(-1.5, 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("log(odds ratio) of being primed with high expression of each gene") + theme(axis.title.x = element_blank())
ggsave(filename = paste0(plotDirectory, "oddsRatioMarkers_positive.pdf"), unit = "in", height = 2, width = 5)

ggplot(plotORList %>% dplyr::filter(gene %in% plotGeneList[6:10]), aes(x = gene, y = pred)) +
  geom_errorbar(aes(ymin = pred - se, ymax = pred + se), width = 0.25) +
  geom_text(aes(label = round(pval, 2)), hjust = -1) +
  geom_point(size = 5) + ylim(-1.5, 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("log(odds ratio) of being primed with high expression of each gene") + theme(axis.title.x = element_blank())
ggsave(filename = paste0(plotDirectory, "oddsRatioMarkers_negative.pdf"), unit = "in", height = 2, width = 5)

ggplot(plotORList %>% dplyr::filter(gene %in% plotGeneList[11:12]), aes(x = gene, y = pred)) +
  geom_errorbar(aes(ymin = pred - se, ymax = pred + se), width = 0.25) +
  geom_text(aes(label = round(pval, 2)), hjust = -1) +
  geom_point(size = 5) + ylim(-1.5, 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("log(odds ratio) of being primed with high expression of each gene") + theme(axis.title.x = element_blank())
ggsave(filename = paste0(plotDirectory, "oddsRatioMarkers_control.pdf"), unit = "in", height = 2, width = 2)

#### combine samples by simply pooling means and standard deviations ####
#######################################################################################################################################################
plotGeneList <- c("S2", "CENPF", "MKI67", "TOP2A", "SOX21", "S3", "SPP1", "GDF15", "CDKN1A", "FTH1", "GAPDH", "UBC")
simpORTable <- bind_rows(oddsTable_R1, oddsTable_R2, oddsTable_R3)
simpORTable <- simpORTable %>% dplyr::filter(prop %in% c(NA, 0.1))
for(i in 1:length(plotGeneList)){
  simpORTableTemp <- simpORTable %>% dplyr::filter(gene == plotGeneList[i])
  simpORTableCombTemp <- tibble(gene = plotGeneList[i], mean = mean(simpORTableTemp$log_estimate), se = sqrt(sum(simpORTableTemp$se^2)/nrow(simpORTableTemp)))
  if(i == 1) {
    simpORTableComb <- simpORTableCombTemp
  } else{
    simpORTableComb <- bind_rows(simpORTableComb, simpORTableCombTemp)
  }
}

simpORTable$gene <- factor(simpORTable$gene, levels = plotGeneList)
ggplot(simpORTable, aes(x = gene, y = log_estimate)) +
  geom_errorbar(aes(ymin = log_estimate - se, ymax = log_estimate + se), width = 0.25) +
  geom_point(size = 5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("log(odds ratio) of being primed with high expression of each gene") + theme(axis.title.x = element_blank())
