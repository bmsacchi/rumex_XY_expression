library(tidyverse)
library(DESeq2)
library(gtools)

## load in reads_pg
reads_pg<-read_csv("data/reads_pg_Sept2024.csv.gz")


### ASE test
ase_data<-reads_pg %>% 
  select(c(pairID,contains("M", ignore.case = FALSE ))) %>% 
  column_to_rownames(var ="pairID")
## ase matrix
IDs<-colnames(ase_data)
sorted_col_names <- mixedsort(IDs)
ase_data <- ase_data[, sorted_col_names]
#colnames(ase_data)
# yay!!
sample<-str_remove(sorted_col_names,"..at") # tm1 tm1 etc.
allele<-str_extract(sorted_col_names,".at$") # mat pat mat pat

colInfo<-data.frame(sorted_col_names,sample,allele) # wooo

ase_matrix<-as.matrix(ase_data)
design<- ~sample+allele

dds <- DESeqDataSetFromMatrix(ase_data, colInfo, design)
# load size factors

library(dplyr)

autoFactorsdf<-read_csv("data/autosomalSizeFactors.csv") %>% # filter by samples containing M, case sensitive
  filter(str_detect(Sample, "M")) 


df.mat<-autoFactorsdf %>% mutate(Sample = paste0(Sample, ".mat"))  # append .mat to each sample in new df using dplyr
df.pat<-autoFactorsdf %>% mutate(Sample = paste0(Sample, ".pat"))  # append .pat to each sample in new df using dplyr

autoFactorsdf <- rbind(df.mat, df.pat) # combine the two dfs
## sort based on order in sorted_col_names
autoFactorsdf<-autoFactorsdf[match(sorted_col_names, autoFactorsdf$Sample),] # match the order of the samples in the deseq object

sizeFactors(dds) <- autoFactorsdf$autoFactors # set size factors to unity
dds <- DESeq(dds, fitType="local") # yay!
#sizeFactors(dds)
resultsNames(dds)

res.allele <- results(dds, name="allele_pat_vs_mat")
head(res.allele$log2FoldChange)

adj_pvals<-res.allele$padj


pvals_genes_norm<-data.frame(res.allele@rownames,res.allele$padj,
                        res.allele$log2FoldChange,
                        res.allele$lfcSE) 
write_csv(pvals_genes_norm, "data/norm_rna_ase_resultsOnlySept2024.csv.gz")


##### NEED TO PLOT NORMALIZED COUNTS

## normalized counts

normCounts<-as.data.frame(counts(dds, normalized=T)) %>%
  rownames_to_column(var = "pairID") %>%
  select(!contains("F", ignore.case = FALSE )) %>%
  rowwise() %>%
  mutate(
    matCountsTotal = sum(c_across(contains("mat", ignore.case = FALSE))), 
    matMean = mean(c_across(contains("mat", ignore.case = FALSE))),
    matSE = plotrix::std.error(c_across(contains("mat", ignore.case = FALSE))),
    patCountsTotal = sum(c_across(contains("pat", ignore.case = FALSE))), 
    patMean = mean(c_across(contains("pat", ignore.case = FALSE))),
    patSE = plotrix::std.error(c_across(contains("pat", ignore.case = FALSE)))) %>%
  ungroup()

#write_csv(normCounts, "data/normalizedCounts_xyExp.csv")

## unite norm counts with pvals
reads_pg_pvals_norm<-inner_join(normCounts, pvals_genes_norm, 
                                by = c("pairID" ="res.allele.rownames"), keep=TRUE) %>%
  mutate(
    sigTest = # just do significance, leave the log2foldchane to be filtered later
      case_when(res.allele.padj < 0.1 ~ "sig", 
                TRUE ~"nonsig" ))
### summarize
get_obs_table(reads_pg_pvals_norm, 0.5)
fold_change_cutoffs <- c(0.5, 1, 1.5, 2)
obs_tables <- lapply(fold_change_cutoffs, function(cutoff) {
  get_obs_table((filter(reads_pg_pvals_norm, res.allele.padj < 0.1)), cutoff)
})

combined_table <- do.call(rbind, obs_tables)
print(combined_table) ## same results, not surprisingly
## test occurs at the sample level, normalized counts won't change that

reads_pg_pvals2_norm <- reads_pg_pvals_norm %>% mutate(sigLFC = case_when(
  res.allele.padj < 0.1 &
    abs(res.allele.log2FoldChange) > 1 ~ "True",
  TRUE ~"False" )) #%>%
  mutate(log2male_mean.mat = log2(male_mean.mat)) %>%
  mutate(log2male_mean.pat = log2(male_mean.pat))

reads_pg_pvals2$sigLFC<-factor(reads_pg_pvals2_norm$sigLFC, levels = c("True","False"))
#mycolors<-RColorBrewer::brewer.pal(2,"Set1")[2:3]

figure1a_norm <- 
  ggplot(reads_pg_pvals2_norm,aes(x=matMean,y=patMean, color = sigLFC)) + 
  #geom_point(aes(size = res.allele.lfcSE),alpha = 0.60) +
  geom_point(alpha = 0.5, size =1) +
  geom_errorbar(aes(ymin = patMean-patSE, ymax = patMean+patSE)) +
  geom_errorbarh(aes(xmin = matMean-matSE, xmax =matMean+matSE)) +
  geom_abline(slope = 1, intercept = 0, linetype=1) +
  theme_bw()+
  geom_smooth(method ="lm") +
  xlab("mean normalized expression of X\n gametolog in males") +
  ylab("mean normalized expression of Y\n gametolog in males") +
  scale_x_continuous(transform = "log2",limits = c(NA,4000),breaks = c(0.0625,2,64,2048)) +
  scale_y_continuous(transform = "log2",limits = c(NA,4000),breaks = c(0.0625,2,64,2048)) +
  labs(color= "p-val < 0.1 & |log2FC| >1 ", size = "log2FC std. error") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.text = element_text(size =12),
        legend.title = element_text(size =12)) 

ggsave("figures/Figure1_gametologExpressionNormalized.png", h =6, w =10)

figure1a_norm
figure1a
