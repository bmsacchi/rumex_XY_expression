library(DESeq2)
library(tidyverse)

# calc log-transformed RPKM

eQTL_counts<-read_delim("data/readCountsRNAeQTL/rhastTXeQTLrna.txt")

sag_counts<-read_delim("data/readCountsSagittatus/sagittatusRnaCounts.txt") 
