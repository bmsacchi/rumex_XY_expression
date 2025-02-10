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

autoFactors<- read_csv("data/autosomalSizeFactors.csv") %>% column_to_rownames("Sample")
## save deseq output incl auto size factors
dds_all$sizeFactor<-autoFactors
dds_all$sizeFactor


#normalized_counts<-counts(dds_all, normalized=TRUE)


####
### compare model where male expression is halved for X

### need to extract just X genes
## divede female expression by 2

X_only_all<-rnaCountsClean %>% column_to_rownames("Geneid") %>% filter(Chr =="X") %>% select("39dF":"43bM") 
OtherChr_all <- rnaCountsClean %>% column_to_rownames("Geneid") %>% filter(Chr != "X") %>% select("39dF":"43bM") 

################ filter by hemiz and y/x genes ##########################
## for later
hemiz_blast_confirm<-read_csv("data/tx_yloss_blastconfirmed_Aug2024.csv.gz") %>% 
  select(tx_mat)
reads_pg_pvals<-read_csv("data/rna_ase_results_eqtl_sept12.csv.gz")
x_overexpressed<-filter(reads_pg_pvals,
                        res.allele.padj < 0.1 & 
                          res.allele.log2FoldChange < 1) %>%
  select(id_tx_mat, id_tx_pat)

y_overexpressed<-filter(reads_pg_pvals, res.allele.padj < 0.1 &
                          res.allele.log2FoldChange > 1) %>%
  select(id_tx_mat,id_tx_pat)

all_reads_pg<-reads_pg_pvals %>% select(id_tx_mat, id_tx_pat)


#### divide female X chr in 2 
X_only_F <- X_only_all %>% select(contains("F", ignore.case = "F")) 
X_only_M <- X_only_all %>% select(contains("M", ignore.case = "F"))
X_only_F_half <- round(X_only_F/2)
X_only_M_2x<- X_only_M*2
###
### sample level normalizing? how does it change? keep the same as above?
# auto factors only for all?
X_only_2x<-bind_cols(X_only_F,X_only_M_2x)
X_only_half<-bind_cols(X_only_F_half,X_only_M)
################
all_halfX<- bind_rows(X_only_half, OtherChr_all)
all_2x<- bind_rows(X_only_2x, OtherChr_all)

####### run DESEq on both ######################

df_half<-DESeqDataSetFromMatrix(all_halfX, sampleInfo, design = ~Sex)
df_2x<-DESeqDataSetFromMatrix(all_2x, sampleInfo, design = ~Sex)

### half female X expression, rounded to nearest integer. 0.5 rounded to closest even
dds_half<-DESeq(df_half)
res_half<-results(dds_half)
summary(res_half)
resultsNames(dds_half)
## tidy results, get chr info
genechr<-rnaCountsClean %>% select(Geneid, Chr)
res_half_tidy<-results(dds_half, tidy = TRUE) %>% dplyr::rename(Geneid = row)

res_half_chr<-full_join(genechr,res_half_tidy)
res_half_X<- res_half_chr %>% filter(Chr == "X")  %>%
  mutate(genetype = case_when(Geneid%in%x_overexpressed$id_tx_mat ~ "xOverexp",
                              Geneid%in%y_overexpressed$id_tx_mat ~ "yOverexp",                                             Geneid%in%hemiz_blast_confirm$tx_mat ~ "hemizygous", 
                              TRUE ~ "other")) %>% 
  filter(Geneid%in%all_reads_pg$id_tx_mat | 
           Geneid%in%all_reads_pg$id_tx_pat | 
           Geneid%in%hemiz_blast_confirm$tx_mat) 



res_half_X_summary<- res_half_X %>% group_by(sig = padj < 0.1, LFC_pos = log2FoldChange >1, LFC_neg = log2FoldChange < -1) %>% summarise(n=n()) %>%
  mutate(expType = case_when(sig ==TRUE & LFC_pos == TRUE ~ "maleUpreg", 
                             sig ==TRUE & LFC_neg == TRUE ~ "femaleUpreg", TRUE ~ NA)) %>%
 group_by(expType) %>%
  summarise(n=sum(n))
res_half_X_summary 

## Take two
############### doubling X expression in males to avoid rounding ##########

dds_2x<-DESeq(df_2x)
res_2x<-results(dds_2x)
summary(res_2x)
resultsNames(dds_2x)
## tidy results, get chr info
genechr<-rnaCountsClean %>% select(Geneid, Chr)
res_2x_tidy<-results(dds_2x, tidy = TRUE) %>% dplyr::rename(Geneid = row)

res_2x_chr<-full_join(genechr,res_2x_tidy)
res_2x_X<- res_2x_chr %>% filter(Chr == "X") %>%
  mutate(genetype = case_when(Geneid%in%x_overexpressed$id_tx_mat ~ "xOverexp",
                              Geneid%in%y_overexpressed$id_tx_mat ~ "yOverexp",                                             Geneid%in%hemiz_blast_confirm$tx_mat ~ "hemizygous", 
                              TRUE ~ "other")) %>% 
  filter(Geneid%in%all_reads_pg$id_tx_mat | 
         Geneid%in%all_reads_pg$id_tx_pat | 
         Geneid%in%hemiz_blast_confirm$tx_mat) 

res_2x_X_summary<- res_2x_X %>% 
  group_by(sig = padj < 0.1, LFC_pos = log2FoldChange >1, LFC_neg = log2FoldChange < -1) %>%
  summarise(n=n())  %>%
  mutate(expType = case_when(sig ==TRUE & LFC_pos == TRUE ~ "maleUpreg", 
                             sig ==TRUE & LFC_neg == TRUE ~ "femaleUpreg", TRUE ~ NA)) %>%
  group_by(expType) %>%
 summarise(n=sum(n))
res_2x_X_summary


#### all genes SBGE
#readsMatAll
df_sex<-DESeqDataSetFromMatrix(readsMatAll, sampleInfo, design = ~Sex)
dds_sex<-DESeq(df_sex)
res_sex <-results(dds_sex)
summary_sex <- summary(res_sex)
res_sex_tidy <- results(dds_sex, tidy = TRUE) %>% dplyr::rename(Geneid = row)
res_sex_chr<- full_join(genechr,res_sex_tidy)
res_sex_X<- res_sex_chr %>% filter(Chr == "X")  %>%
  mutate(genetype = case_when(Geneid%in%x_overexpressed$id_tx_mat ~ "xOverexp",
                              Geneid%in%y_overexpressed$id_tx_mat ~ "yOverexp",                                             Geneid%in%hemiz_blast_confirm$tx_mat ~ "hemizygous", TRUE ~ "other")) %>% 
  filter(Geneid%in%all_reads_pg$id_tx_mat | Geneid%in%all_reads_pg$id_tx_pat | Geneid%in%hemiz_blast_confirm$tx_mat)

res_sex_X_summary<-res_sex_X %>% 
  group_by(sig = padj < 0.1, LFC_pos = log2FoldChange >1, LFC_neg = log2FoldChange < -1) %>%     summarise(n=n()) %>% 
  mutate(expType = case_when(sig ==TRUE & LFC_pos == TRUE ~ "maleUpreg",
                             sig ==TRUE & LFC_neg == TRUE ~ "femaleUpreg", TRUE ~ NA)) %>% group_by(expType) %>% summarise(n=sum(n))
res_sex_X_summary
#res_sex_X_summary%>%filter(sig == TRUE)

######## redo-FDR correction ############
## i may be mistaken, probably not adviseable?
p.adjust(p, method = p.adjust.methods, n = length(p))

res_sex_X_fdr<-res_sex_X %>% mutate(newP = p.adjust(pvalue, method = "fdr", n = 612))
res_2x_X_fdr<-res_2x_X %>% mutate(newP = p.adjust(pvalue, method = "fdr", n = 612))

res_sex_X_summary_fdr<-res_sex_X_fdr %>% 
  group_by(sig = newP < 0.1, LFC_pos = log2FoldChange >1, LFC_neg = log2FoldChange < -1) %>%     summarise(n=n()) %>% 
  mutate(expType = case_when(sig ==TRUE & LFC_pos == TRUE ~ "maleUpreg",
                             sig ==TRUE & LFC_neg == TRUE ~ "femaleUpreg", TRUE ~ NA)) %>% group_by(expType) %>% summarise(n=sum(n))
res_sex_X_summary_fdr

res_2x_X_summary_fdr<-res_2x_X_fdr %>% 
  group_by(sig = newP < 0.1, LFC_pos = log2FoldChange >1, LFC_neg = log2FoldChange < -1) %>%     summarise(n=n()) %>% 
  mutate(expType = case_when(sig ==TRUE & LFC_pos == TRUE ~ "maleUpreg",
                             sig ==TRUE & LFC_neg == TRUE ~ "femaleUpreg", TRUE ~ NA)) %>% group_by(expType) %>% summarise(n=sum(n))
res_2x_X_summary_fdr
### barely different. not worth the trouble!

