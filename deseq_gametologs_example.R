### DESEQ gametolog ASE pipeline
#
# files needed
# - file containing orthologs shared between X and Y
# - raw read counts per sample (e.g. featurecounts or htseq output)
library(tidyverse)
library(DESeq2)
# load file containing gametologs
pgMatOrths<-read_csv("data/pgMatOrths.csv")
# load readcounts table for all samples
# or load a separate file for each sample and merge together
rnaCounts<-read_table("readCountsRNAeQTL/rhastTXeQTLrna.txt")
# fix the long weird sample names
rnaColnames<-colnames(rnaCounts)
# i used chatgpt to get the correct regular expression for my samples
# it is very useful for this
rnaColnamesRepaired <- gsub(".*NEBNext_dual_i._.+\\.([0-9]+[a-z]+[FM])LRAligned.+", "\\1", rnaColnames)
colnames(rnaCounts)<-rnaColnamesRepaired 

# make sure the readcounts file is organized. 
# you won't need to do all of these steps 
# the male_mean column will be useful for plotting
rnaCountsClean <- rnaCounts %>% 
  relocate(1:6,ends_with("F"), ends_with("M")) %>% # sort columns for ease
  mutate(Chr = ifelse(grepl("X",Chr), "X", # relabel Chrs for convenience
                      ifelse(grepl("Y",Chr), "Y",
                             ifelse(grepl("A1",Chr), "A1",
                                    ifelse(grepl("A2",Chr), "A2",
                                           ifelse(grepl("A3",Chr), "A3",
                                                  ifelse(grepl("A4",Chr), "A4", NA))))))) %>%
  rowwise() %>% 
  mutate(totalCountsGene = sum(c_across(7:155))) %>% # summarising across all samples (columns 7 to 155)
  mutate(meanCounts = totalCountsGene/149) %>%
  mutate(femaleCountsTotal = sum(c_across(contains("F", ignore.case = FALSE)))) %>%
  mutate(female_mean = femaleCountsTotal/74) %>%
  mutate(maleCountsTotal = sum(c_across(contains("M", ignore.case = FALSE)))) %>%
  mutate(male_mean = maleCountsTotal/75) 


## extract X gene counts
eqtl_x<- filter(rnaCountsClean, Chr =="X") %>% select(c(Geneid,totalCountsGene:"43bM")) 
## extract Y gene counts
eqtl_y<- filter(rnaCountsClean, Chr =="Y") %>% select(c(Geneid,totalCountsGene:"43bM")) 


# filter out weird genes

# get genes where females have too much Y coverage
# some may be on Y PAR - which is fromm 0-45MB. Disregard and remove all for now!
badgenes<-rnaCountsClean %>% 
  filter(female_mean>= 20 & Chr =="Y") %>% select(Geneid) 
# there may be corresponding X genes to the "bad genes" on Y
badgenes_x <- pgMatOrths %>% select(id_tx_mat,id_tx_pat) %>% filter(id_tx_pat%in%badgenes$Geneid) %>% select(id_tx_mat)

reads_pg<-left_join(pgMatOrths,eqtl_x, by = c("id_tx_mat"="Geneid")) %>%
  left_join(.,eqtl_y, by = c("id_tx_pat" = "Geneid"),suffix = c(".mat",".pat")) %>% 
  filter(id_tx_mat%nin%badgenes_x$id_tx_mat) %>% # weird expression filter
  filter(id_tx_pat%nin%badgenes$Geneid) %>%
  filter(start_tx_pat > 45000000) %>% #  filter out pseudoautosomal region
  drop_na() %>% select(-flag_tx_pat,-flag_tx_mat) %>% 
  unique() %>%
  rowwise() %>% 
  mutate(totalCountAll=(totalCountsGene.mat+totalCountsGene.pat)) %>%
  ungroup() %>% filter(totalCountAll >20) # filter out gametologs where totalcount <20 reads


### ASE test

# select males only!
ase_data<-reads_pg %>% select(c(pairID,contains("M", ignore.case = FALSE ))) %>%
  column_to_rownames(var ="pairID") # turn pairID column to rowname
## ase matrix
IDs<-colnames(ase_data) # column names are sample IDs
library(gtools) # optional
sorted_col_names <- mixedsort(IDs) #helpful for sorting. 
# you'll want the columns to be like sample1.pat sample1.mat sample2.pat sample2.mat etc...
ase_data <- ase_data[, sorted_col_names] # create dataframe

sample<-str_remove(sorted_col_names,"..at") # get sample name without .pat or .mat suffix
allele<-str_extract(sorted_col_names,".at$") # extract and keep only .mat or .pat suffix

colInfo<-data.frame(sorted_col_names,sample,allele) # create column info dataframe

# save ase data as a matrix
ase_matrix<-as.matrix(ase_data)
# study design
design<- ~sample+allele
# 

m <- ####number of males

dds <- DESeqDataSetFromMatrix(ase_data, colInfo, design)
sizeFactors(dds) <- rep(1, 2*m) # size factors for all samples should be 1
dds <- DESeq(dds, fitType="local") # run deseq

resultsNames(dds) # get the results name. choose the last one.

# get results table for the pat vs mat allele comparison
res.allele <- results(dds, name="allele_pat_vs_mat")

# create table with rownames, adjusted p values, log2FC and log2FC Standard error
pvals_genes<-data.frame(res.allele@rownames,res.allele$padj,res.allele$log2FoldChange,res.allele$lfcSE) 

# merge with reads_pg dataframe
# choose significance filters.
reads_pg_pvals<-inner_join(reads_pg, pvals_genes, by =c("pairID" ="res.allele.rownames"))%>%
  mutate(sigTest = 
           case_when(res.allele.padj < 0.1 & abs(res.allele.log2FoldChange) > 1 ~ "sig",
                     TRUE ~"nonsig" )) %>% 
  relocate(c(pairID:male_mean.mat, meanCounts.pat:male_mean.pat, res.allele.padj:sigTest))
# save your results
# write_csv(reads_pg_pvals, "rna_ase_results_eqtl_may17.csv")


# plot male x mean vs male y expression mean

ggplot(reads_pg_pvals,aes(x=male_mean.mat,y=male_mean.pat, color = sigTest)) + geom_point(aes(size = res.allele.lfcSE),alpha =0.5) +
  geom_abline(slope = 1, intercept = 0, linetype=1) +
  theme_bw()+
  geom_smooth(method ="lm") +
  xlab("mean exp of X gametolog") +
  ylab("mean exp of Y gametolog") +
  scale_x_continuous(transform = "log2") +
  scale_y_continuous(transform = "log2")


ggsave("filename.pdf",h =5, w=7)



