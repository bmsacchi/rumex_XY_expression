library(tidyverse)
#dnaSamples<-read_table("DNA.samples_full.txt", col_names = "sampleName")

head(dnaSample)
# format:

# Family ID
# Individual ID
# Paternal ID
# Maternal ID
# Sex (1=male; 2=female; other=unknown)
# Phenotype


#remove 24fM
# low coverage
samples_short<-dnaColnamesRepaired[7:158]

samplesPED<-as.data.frame(samples_short) %>% filter(samples_short != "24fM") %>%
  dplyr::rename(indivID = samples_short) %>%
  mutate(.,Sex = if_else(grepl("*M", indivID),1, 2)) %>% # 
  mutate(familyID = as.numeric(gsub("[^0-9]", "", indivID))) %>% 
  mutate(patID = 0) %>%
  mutate(matID = 0) %>%
  mutate(pheno = 0) %>%
  relocate(familyID, indivID,patID,matID,Sex,pheno)
 
write_csv(samplesPED, "PEDfile_DNA_eqtl.csv")
