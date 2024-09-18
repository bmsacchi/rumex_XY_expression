# if running in isolation
library(tidyverse)
library(DESeq2)
library(gtools)
## source prev script
rnaCountsClean<-read_csv("data/rnaCountsClean.csv.gz")

source("scripts/gametologExpressionDNA.R")
pgMatOrths<-read_csv("data/pgMatOrths.csv")
#### RNA counts
### most up-to date results
## main pipeline won't run without it
file_name<-"data/rna_ase_resultsOnlySept2024.csv" ## ase stat resutls
file_name2 = "data/rna_ase_results_eqtl_sept12.csv.gz" ## results merged with pangene info

# filtering out problematic genes from RNA
badgenes<-rnaCountsClean %>% 
  filter(female_mean>= 20 & Chr =="Y") %>% select(Geneid) # remove genes with improper female Y mapping
# some may be on Y PAR - which is fromm 0-45MB. Disregard and remove all for now!
badgenes_x <- pgMatOrths %>% select(id_tx_mat,id_tx_pat) %>% filter(id_tx_pat%in%badgenes$Geneid) %>% select(id_tx_mat)
## find corresponding X gene for every "bad" Y gene

# how does this compare to DNA
shared_bad_biased<-intersect(badgenes$Geneid,biasedgenes_dna$id_tx_pat) # 150 genes
# a fair number of genes with mismapping in the RNA are biased in the DNA
shared_bad<-intersect(badgenes$Geneid,badgenes_dna$id_tx_pat) # 394 genes
# a high number of genes with mismapping in the RNA are also mismapping in the DNA

## formatting data for deseq ASE approach
eqtl_x<- filter(rnaCountsClean, Chr =="X") %>% select(c(Geneid,totalCountsGene:"43bM")) ## extract X gene counts
eqtl_y<- filter(rnaCountsClean, Chr =="Y") %>% select(c(Geneid,totalCountsGene:"43bM")) ## extract Y gene counts

# no filtering
reads_pg_nofilter<-left_join(pgMatOrths,eqtl_x, by = c("id_tx_mat"="Geneid")) %>%
  left_join(.,eqtl_y, by = c("id_tx_pat" = "Geneid"),suffix = c(".mat",".pat")) %>%
  drop_na()%>%  summarise(n = n())

# filtering out pseudoautosomal region
reads_pg_PARfilter<-left_join(pgMatOrths,eqtl_x, by = c("id_tx_mat"="Geneid")) %>%
  left_join(.,eqtl_y, by = c("id_tx_pat" = "Geneid"),suffix = c(".mat",".pat")) %>%
  filter(start_tx_pat > 45000000) %>% summarise(n = n())

# filtering out genes with female Y mismapping in the rnaseq samples
reads_pg_badRNAgenesFilter<-left_join(pgMatOrths,eqtl_x, by = c("id_tx_mat"="Geneid")) %>%
  left_join(.,eqtl_y, by = c("id_tx_pat" = "Geneid"),suffix = c(".mat",".pat")) %>%
  filter(id_tx_mat%nin%badgenes_x$id_tx_mat) %>% # weird expression filter
  filter(id_tx_pat%nin%badgenes$Geneid) %>%
  filter(start_tx_pat > 45000000) %>%
  summarise(n = n())

# filtering out genes with female Y mismapping in the genomic DNA samples
reads_pg_badDNAgenesFilter<-left_join(pgMatOrths,eqtl_x, by = c("id_tx_mat"="Geneid")) %>%
  left_join(.,eqtl_y, by = c("id_tx_pat" = "Geneid"),suffix = c(".mat",".pat")) %>%
  filter(id_tx_mat%nin%badgenes_x$id_tx_mat) %>% # weird expression filter
  filter(id_tx_pat%nin%badgenes$Geneid) %>%
  filter(start_tx_pat > 45000000) %>%
  #filter(pairID%nin%biasedgenes_dna$pairID) %>% # biased gene filter
  filter(pairID%nin%badgenes_dna$pairID) %>%
  summarise(n = n())

# filtering out genes with with evidence of mapping bias in DNA
# significant LFC >1 or < -1
reads_pg_biasedDNAgenesFilter<-left_join(pgMatOrths,eqtl_x, by = c("id_tx_mat"="Geneid")) %>%
  left_join(.,eqtl_y, by = c("id_tx_pat" = "Geneid"),suffix = c(".mat",".pat")) %>%
  filter(id_tx_mat%nin%badgenes_x$id_tx_mat) %>% # weird expression filter
  filter(id_tx_pat%nin%badgenes$Geneid) %>%
  filter(start_tx_pat > 45000000) %>%
  filter(pairID%nin%biasedgenes_dna$pairID) %>% # biased gene filter
  filter(pairID%nin%badgenes_dna$pairID) %>%
  summarise(n = n())

# last filter: all previous, and total read counts across all samples must be greater than 20 
reads_pg<-left_join(pgMatOrths,eqtl_x, by = c("id_tx_mat"="Geneid")) %>%
  left_join(.,eqtl_y, by = c("id_tx_pat" = "Geneid"),suffix = c(".mat",".pat")) %>% 
  filter(id_tx_mat%nin%badgenes_x$id_tx_mat) %>% # weird expression filter
  filter(id_tx_pat%nin%badgenes$Geneid) %>%
  filter(start_tx_pat > 45000000) %>% # par filter
  filter(pairID%nin%biasedgenes_dna$pairID) %>% # biased gene filter
  filter(pairID%nin%badgenes_dna$pairID) %>% # filter weird Y mapping
  drop_na() %>% select(-flag_tx_pat,-flag_tx_mat) %>% 
  unique() %>%
  rowwise() %>% 
  mutate(totalCountAll=(totalCountsGene.mat+totalCountsGene.pat)) %>%
  ungroup() %>% filter(totalCountAll >20)

#write_csv(reads_pg, "data/reads_pg_Sept2024.csv.gz")

reads_pg_allfilters<-reads_pg %>% summarise(n=n())

# table of read numbers
category<-c("reads_pg_nofilter","reads_pg_PARfilter","reads_pg_badRNAgenesFilter","reads_pg_badDNAgenesFilter","reads_pg_biasedDNAgenesFilter","reads_pg_allfilters")
counts<-c(reads_pg_nofilter$n,reads_pg_PARfilter$n,reads_pg_badRNAgenesFilter$n,reads_pg_badDNAgenesFilter$n,reads_pg_biasedDNAgenesFilter$n,reads_pg_allfilters$n)

filtering_summary<-data.frame(category,counts)


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

m <- 75



if (!file.exists(file_name)) {
  # Your R command goes here
  print("Running the command because the file does not exist.")

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
  #write_csv(pvals_genes, "data/rna_ase_resultsOnlySept2024.csv")
} else { 
  print(paste("Skipping command because", file_name, "exists."))
# write pvals to file
  pvals_genes<-read_csv("data/rna_AsE_resultsOnlySept2024.csv") 
  }# read in data


if (!file.exists(file_name2)) {
reads_pg_pvals<-inner_join(reads_pg, pvals_genes, 
                           by = c("pairID" ="res.allele.rownames")) %>%
  mutate(
    sigTest = # just do significance, leave the log2foldchane to be filtered later
    case_when(res.allele.padj < 0.1 ~ "sig", 
              TRUE ~"nonsig" )) %>% 
  relocate(c(pairID:male_mean.mat,
             meanCounts.pat:male_mean.pat,
             res.allele.padj:sigTest))
   write_csv(reads_pg_pvals,"data/rna_ase_results_eqtl_sept12.csv.gz")
} else { 
  
  reads_pg_pvals<-read_csv("data/rna_ase_results_eqtl_sept12.csv.gz")# read in data
  }# read in data



fold_change_cutoffs <- c(0.5, 1, 1.5, 2)
obs_tables <- lapply(fold_change_cutoffs, function(cutoff) {
  get_obs_table((filter(reads_pg_pvals, res.allele.padj < 0.1)), cutoff)
})

combined_table <- do.call(rbind, obs_tables)
print(combined_table)
#write_csv(combined_table, "data/n_xy_oe_RNA_eQTL_lfcRanges.csv")
#267 genes total






