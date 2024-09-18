#source("scripts/dataWranglingRawCounts.R")
#source("scripts/gametologExpressionDNA.R")
#source("scripts/gametologExpressionRNA.R")
#source("scripts/autosomalNormFactors.R")
#source("scripts/gametologNormExpressionRNA.R")
source("scripts/xyExpressionPlots.R")
source("scripts/geneloss.R")

library(kableExtra)

#### filtering steps #### 
filtering_summary %>% kable() %>% kable_styling(full_width = F) %>% 
  save_kable("tables/filtering_summary.html")

#### dna bias ####
combined_table_dna %>% kable() %>% kable_styling(full_width = F) %>% 
  save_kable("tables/dna_foldchange.html")
#### rna - ase ####
combined_table %>% kable() %>% kable_styling(full_width = F) %>% 
  save_kable("tables/rna_foldchange.html")
## figures

figure1a
#ggsave("figures/Figure1_gametologExpression.png", h =6, w =10)

figure1a_norm
#ggsave("figures/Figure1_gametologExpressionNormalized.png", h =6, w =10)
figure1b
## we have plots with normalized reads
## same overall pattern, but the error bars are somewhat different


##### windows #####


#### gene loss ####
## chisquare test for different gene categories
resTable
testTable
test1
## summary table
geneloss_summary %>% kable() %>% kable_styling(full_width = F) %>% 
  save_kable("tables/geneloss_summary.html")

## gene loss and expression



## dosage compensation
