# if running in isolation
library(tidyverse)
library(DESeq2)
library(gtools)

`%nin%` <- Negate(`%in%`)
## source prev script
#rnaCountsClean<-read_csv("data/rnaCountsClean.csv.gz")
## contains pollen 
raw_counts <- read_delim("../rhastYchromosome/data/deseq2/rhastTXeQTLrna_pollenleaf.txt")
my_cols <- colnames(raw_counts)#my_cols_multi <- colnames(raw_counts_multi)
my_cols_repaired <- gsub(
    ".*NEBNext_dual_i._.+\\.([0-9]+[a-z]+[FM][PL])RAligned.+",
    "\\1",
    my_cols
) # extract sample names
colnames(raw_counts) <- my_cols_repaired
cols_to_remove <- c(
  "24fMP","35aMP","7bMP","27eMP","53bMP","5aMP")
## pollen RNA counts


rnaCountsClean <- raw_counts %>%
  tidyr::separate(Chr, c("chr"), sep = ";", extra = "drop") %>%
  tidyr::separate(Start, c("start"), sep = ";", extra = "drop") %>%
  tidyr::separate(End, c("end"), sep = ";", extra = "drop") %>%
  mutate(
    start = as.integer(start),
    end = as.integer(end)
  ) %>%
  select(-all_of(cols_to_remove)) %>%
  select(Geneid, chr, start, end, contains("P", ignore.case = "FALSE")) %>%
  mutate(
      # Vectorized row-wise sums
      pollen_total = base::rowSums(
          select(., "78eMP":"70bMP"),
          na.rm = TRUE
      ),
      pollen_mean = base::rowMeans(
        select(.,  "78eMP":"70bMP"),
        na.rm = TRUE
      ))

#pollen_cols <- grep("P$", names(data_cleaned), value = TRUE)


pgMatOrths<-read_csv("data/pgMatOrths.csv.gz")
source("scripts/gametologExpressionDNA.R")
#### RNA counts

## formatting data for deseq ASE approach
eqtl_x<- filter(rnaCountsClean, chr =="X") %>% 
  select(c(Geneid,pollen_total,pollen_mean,"78eMP":"70bMP")) ## extract X gene counts
eqtl_y<- filter(rnaCountsClean, chr =="Y") %>%
  select(c(Geneid,pollen_total,pollen_mean,"78eMP":"70bMP"))# extract Y gene counts

# no filtering
reads_pg_nofilter <- left_join(
  pgMatOrths,
  eqtl_x,
  by = c("id_tx_mat" = "Geneid")
) %>%
  left_join(
    .,
    eqtl_y,
    by = c("id_tx_pat" = "Geneid"),
    suffix = c(".mat", ".pat")
  ) %>%
  drop_na() %>%
  summarise(n = n())

# filtering out pseudoautosomal region
reads_pg_PARfilter<-left_join(pgMatOrths,eqtl_x, by = c("id_tx_mat"="Geneid")) %>%
  left_join(.,eqtl_y, by = c("id_tx_pat" = "Geneid"),suffix = c(".mat",".pat")) %>%
  filter(start_tx_pat > 45000000) %>% summarise(n = n())
#filtering out genes with female Y mismapping in the genomic DNA samples
reads_pg_badDNAgenesFilter<-left_join(pgMatOrths,eqtl_x, by = c("id_tx_mat"="Geneid")) %>%
  left_join(.,eqtl_y, by = c("id_tx_pat" = "Geneid"),suffix = c(".mat",".pat")) %>%
  #filter(id_tx_mat%nin%badgenes_x$id_tx_mat) %>% # weird expression filter
  #filter(id_tx_pat%nin%badgenes$Geneid) %>%
  filter(start_tx_pat > 45000000) %>%
  filter(pairID%nin%biasedgenes_dna$pairID) %>% # biased gene filter
  filter(pairID%nin%badgenes_dna$pairID) %>%
  summarise(n = n())

# filtering out genes with with evidence of mapping bias in DNA
# significant LFC >1 or < -1
reads_pg_biasedDNAgenesFilter<-left_join(pgMatOrths,eqtl_x, by = c("id_tx_mat"="Geneid")) %>%
  left_join(.,eqtl_y, by = c("id_tx_pat" = "Geneid"),suffix = c(".mat",".pat")) %>%
  #filter(id_tx_mat%nin%badgenes_x$id_tx_mat) %>% # weird expression filter
  #filter(id_tx_pat%nin%badgenes$Geneid) %>%
  filter(start_tx_pat > 45000000) %>%
  filter(pairID%nin%biasedgenes_dna$pairID) %>% # biased gene filter
  filter(pairID%nin%badgenes_dna$pairID) %>%
  summarise(n = n())

# last filter: all previous, and total read counts across all samples must be greater than 20 
reads_pg<-left_join(pgMatOrths,eqtl_x, by = c("id_tx_mat"="Geneid")) %>%
  left_join(.,eqtl_y, by = c("id_tx_pat" = "Geneid"),suffix = c(".mat",".pat")) %>% 
  filter(start_tx_pat > 45000000) %>% # par filter
  filter(pairID%nin%biasedgenes_dna$pairID) %>% # biased gene filter
  filter(pairID%nin%badgenes_dna$pairID) %>% # filter weird Y mapping
  drop_na() %>% select(-flag_tx_pat,-flag_tx_mat) %>% 
  unique() %>%
  rowwise() %>% 
  mutate(totalCountAll=(pollen_total.mat+pollen_total.pat)) %>%
  ungroup() %>% filter(totalCountAll >20)

write_csv(reads_pg, "data/reads_pg_pollen_April_2025.csv.gz")

#reads_pg_allfilters<-reads_pg %>% summarise(n=n())

# # table of read numbers
# category<-c("reads_pg_nofilter","reads_pg_PARfilter","reads_pg_badRNAgenesFilter","reads_pg_badDNAgenesFilter","reads_pg_biasedDNAgenesFilter","reads_pg_allfilters")
# counts<-c(reads_pg_nofilter$n,reads_pg_PARfilter$n,reads_pg_badRNAgenesFilter$n,reads_pg_badDNAgenesFilter$n,reads_pg_biasedDNAgenesFilter$n,reads_pg_allfilters$n)

# filtering_summary<-data.frame(category,counts)


#### deseq ase test

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

m <- 69 
#length(IDs)

dds <- DESeqDataSetFromMatrix(ase_data, colInfo, design)
sizeFactors(dds) <- rep(1, 2*m)
dds <- DESeq(dds, fitType="local")

resultsNames(dds)

res.allele <- results(dds, name="allele_pat_vs_mat")
head(res.allele$log2FoldChange)
adj_pvals<-res.allele$padj


pvals_genes<-data.frame(res.allele@rownames,res.allele$padj,
                         res.allele$log2FoldChange,
                         res.allele$lfcSE) 
write_csv(pvals_genes, "data/rna_ase_results_pollen_April2025.csv")
# } else {
#   print(paste("Skipping command because", file_name, "exists."))
# # write pvals to file
#   pvals_genes<-read_csv("data/rna_AsE_resultsOnlySept2024.csv")
#   }# read in data

reads_pg_pvals <-
  inner_join(
    reads_pg,
    pvals_genes,
    by = c("pairID" = "res.allele.rownames")
  ) %>%
  mutate(
    # just do significance, leave the log2foldchane to be filtered later
    sigTest = case_when(res.allele.padj < 0.1 ~ "sig", TRUE ~ "nonsig")
  ) %>%
  select(
    pairID:pollen_total.mat,
    pollen_total.pat,
    pollen_mean.mat,
    pollen_mean.pat,
    totalCountAll,
    res.allele.padj,
    res.allele.log2FoldChange,
    res.allele.lfcSE,
    sigTest
  )

write_csv(reads_pg_pvals,"data/rna_ase_results_summary_pollen_April2025.csv")

# } else { 
  
#   reads_pg_pvals<-read_csv("data/rna_ase_results_eqtl_sept12.csv.gz")# read in data
#   }# read in data



# fold_change_cutoffs <- c(0.5, 1, 1.5, 2)
# obs_tables <- lapply(fold_change_cutoffs, function(cutoff) {
#   get_obs_table((filter(reads_pg_pvals, res.allele.padj < 0.1)), cutoff)
# })

# combined_table <- do.call(rbind, obs_tables)
# print(combined_table)
# #write_csv(combined_table, "data/n_xy_oe_RNA_eQTL_lfcRanges.csv")
# #267 genes total


# Call:
# lm(formula = log2male_mean.pat ~ log2male_mean.mat, data = reads_pg_pvals2)

# Residuals:
#     Min      1Q  Median      3Q     Max 
# -9.9247 -0.7468  0.5183  1.3206  8.5569 

# Coefficients:
#                   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)       -0.73865    0.22803  -3.239  0.00131 ** 
# log2male_mean.mat  0.97108    0.04266  22.764  < 2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 2.437 on 349 degrees of freedom
# Multiple R-squared:  0.5976,	Adjusted R-squared:  0.5964 
# F-statistic: 518.2 on 1 and 349 DF,  p-value: < 2.2e-16
# > summary(lm_sig)

# Call:
# lm(formula = log2male_mean.pat ~ log2male_mean.mat, data = reads_pg_pvals2 %>% 
#     filter(sigLFC == "True"))

# Residuals:
#     Min      1Q  Median      3Q     Max 
# -8.5810 -1.8735  0.2087  2.6839  7.0879 

# Coefficients:
#                   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        -0.1766     0.5632  -0.314    0.754    
# log2male_mean.mat   0.6963     0.1097   6.346 2.78e-09 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 3.23 on 142 degrees of freedom
# Multiple R-squared:  0.221,	Adjusted R-squared:  0.2155 
# F-statistic: 40.28 on 1 and 142 DF,  p-value: 2.776e-09