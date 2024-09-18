library(tidyverse)
library(DESeq2)
library(gtools)

## need normalization factors for autosomal genes only
## for dosage compensation analyses
rnaCountsClean <- read_csv("data/rnaCountsClean.csv.gz")
readsMatAll<-rnaCountsClean %>% column_to_rownames("Geneid") %>% select("39dF":"43bM") # convert to matrix of all genes and inds

readsMatAuto<-rnaCountsClean %>% filter(grepl("A",Chr)) %>% # matrix with just autosomes, for generating norm factors
  column_to_rownames("Geneid") %>% select("39dF":"43bM") 

sampleNames<- colnames(readsMatAll)
#
sampleInfo<-data.frame(sampleNames) %>% 
  mutate(.,Sex = if_else(grepl("*M", sampleNames),"M", "F"))

# DEseq all counts
df_all<-DESeqDataSetFromMatrix(readsMatAll,sampleInfo,design = ~ 1)

# DESeq autosomal gene counts only
df_auto<-DESeqDataSetFromMatrix(readsMatAuto,sampleInfo,design = ~ 1)

## run deseq on both
dds_all <- DESeq(df_all)
dds_auto <- DESeq(df_auto)

## get size factors for all genes
allchrFactors<-sizeFactors(dds_all)
## get size factors for autosomal genes only
autoFactors<-sizeFactors(dds_auto)

## compare factors in plot
compareFactors<-data.frame(allchrFactors,autoFactors) %>% rownames_to_column("sample") 
compareFactorsLong <- compareFactors %>% pivot_longer(cols = -sample, names_to = "normType", values_to = "value")
ggplot(compareFactorsLong) + geom_point(aes(x = sample, y = value, color = normType)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#ggsave("compare_norm_factors.pdf", h = 5, w = 18)
## very similar
## use autosome factors anyway *for dc results

## save autosomal factors to file
autoFactorsdf<- as.data.frame(autoFactors) %>% rownames_to_column("Sample")
write_csv(autoFactorsdf, "data/autosomalSizeFactors.csv")
## save deseq output incl auto size factors
dds_all$sizeFactor<-autoFactors
dds_all$sizeFactor


