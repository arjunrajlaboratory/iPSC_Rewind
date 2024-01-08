# load required libraries and functions
library(tidyverse)
library(reshape2)
library(ggridges)
library(stringdist)
library(stringr)

# load dataset into data frames
sampleDir <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/extractedData/embryoidBodyAssay/analyzed/C4/"
plotDir <- "/Users/naveenjain/Dropbox (RajLab)/Shared_Naveen/Paper/plots/embryoidBodyAssay/"
sampleFolders <- list.files(path=sampleDir,pattern='*d8.txt',recursive=TRUE)

sampleTables <- list(rep("NA",length(sampleFolders)))
for (i in 1:length(sampleFolders)) {
  sampleTable = read.table(paste0(sampleDir, sampleFolders[i]), header = FALSE, sep = "\t")
  sampleTables[[i]] <- sampleTable
}

# check heritability control
filterThreshold <- 10
sampleA <- sampleTables[[1]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2/sum(V2)*(10^6))
sampleB <- sampleTables[[2]] %>% filter(.,V2>filterThreshold) %>% mutate(V2=V2/sum(V2)*(10^6))
sampleEB <- sampleTables[[3]] %>% filter(.,V2>filterThreshold*10) %>% mutate(V2=V2/sum(V2)*(10^6))

heritability <- full_join(sampleA,sampleB, by="V1") %>% replace(is.na(.), 0)
heritability <- heritability %>% mutate(EB = rep("no", nrow(heritability)))
heritability$EB <- ifelse(heritability$V1 %in% sampleEB$V1 ,"yes", heritability$EB)

heritabilityEB <- inner_join(heritability, sampleEB,by="V1")
plot <- ggplot() +
  geom_point(aes(x=heritabilityEB$V2.x,y=heritabilityEB$V2.y,fill=heritabilityEB$V2),pch=21,size=3) +
  scale_fill_gradient(low="white",high="blue", name = "reads per million\nin embryoid body split") +
  geom_abline(slope=1,intercept=0, linetype = "dashed") +
  xlab("barcodes reads per million\nin reprogramming split A") + ylab("barcode reads per million\nin reprogramming split B") + theme_classic() +
  theme(legend.position="bottom") +
  annotate(geom = "text", x = 9500, y = 8000, label = "large groups\n(i.e. iPSC colonies)") +
  annotate(geom = "text", x = 3000, y = 1250, label = "small groups\n(i.e. singletons)")

ggsave(filename = paste0(plotDir, "embryoidBodyAssay.pdf"), plot = plot, height = 5.5, width = 5, units = "in")
