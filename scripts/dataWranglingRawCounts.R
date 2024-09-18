library(tidyverse)
library(DESeq2)
library(gtools)
source("scripts/generalFunctions.R")


#### Raw read count data
dnaCounts<-read_table("data/readCountsDNAeQTL/rhastTXeQTLdna.txt")
rnaCounts<-read_table("data/readCountsRNAeQTL/rhastTXeQTLrna.txt")


# rna data
rnaColnames<-colnames(rnaCounts) # save colnames for later
rnaColnamesRepaired <- gsub(".*NEBNext_dual_i._.+\\.([0-9]+[a-z]+[FM])LRAligned.+", "\\1", rnaColnames) # extract sample names
colnames(rnaCounts)<-rnaColnamesRepaired # replace colnames

# dna data
dnaColnames<-colnames(dnaCounts) # save colnames for later
dnaColnamesRepaired<- gsub("NS\\.[0-9]+\\.[0-9]+.IDT_i._\\d+---IDT_i._\\d+\\.(\\d+[a-z]+[FM])LD_sortMarkDups.*", "\\1", dnaColnames) # extract sample names
colnames(dnaCounts)<-dnaColnamesRepaired # replace colnames
setdiff(rnaColnamesRepaired,dnaColnamesRepaired) # check if all samples are present in both datasets
# [1] "67eF" "14bF" "45aF" are not shared. Need to be removed from rna data
# also need to remove 40aF and 24fM from both rna and dna
# 40aF and 35aM have aberrant coverage patterns

## clean data, remove those samples
rnacolsremove<-c("67eF","14bF","45aF","40aF","24fM","35aM")
dnacolsremove<-c("40aF","24fM","35aM")
rnaCountsClean <- clean_counts(rnaCounts,rnacolsremove)
write_csv(rnaCountsClean, "data/rnaCountsClean.csv.gz")# get total gene-level counts for filtering
#rnaCountsClean<-read_csv("data/rnaCountsClean.csv.gz")

dnaCountsClean <- clean_counts(dnaCounts,dnacolsremove)
write_csv(dnaCountsClean, "data/dnaCountsClean.csv.gz")# get total gene-level counts for filtering
#dnaCountsClean<-read_csv("data/dnaCountsClean.csv.gz")
