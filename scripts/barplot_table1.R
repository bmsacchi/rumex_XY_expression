library(ggplot2)
library(tidyverse)
library(khroma)


reads_pg_pvals<-read_csv("data/rna_ase_results_eqtl_sept12.csv.gz")# read in data

fold_change_cutoffs <- c(0.5, 1, 1.5, 2)
obs_tables <- lapply(fold_change_cutoffs, function(cutoff) {
  get_obs_table((filter(reads_pg_pvals, res.allele.padj < 0.1)), cutoff)
})

combined_table <- do.call(rbind, obs_tables)
print(combined_table)

## make combined_table long
df<-pivot_longer(combined_table, 
             cols = c(n_Y_oe, n_X_oe), 
             names_to = "Expression_Type", 
             values_to = "Count") %>%
  mutate(foldChangeCutoff = as.character(foldChangeCutoff)) #%>%
  #mutate(Expression_Type = fct_relevel(Expression_Type,"n_X_oe","n_Y_oe"))

ggplot(df, aes(x = foldChangeCutoff, 
                           y = Count, 
                           fill = Expression_Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Log2 Fold Change Cutoff", 
       y = "Gene Count") +
  #scale_fill_discrete(breaks = "n_Y_oe","n_X_oe") +
 scale_fill_manual(name = "", 
                  values=alpha(c("darkblue","darkorange"),
                  c(0.75, 0.75)), 
                  na.value = alpha("grey",0.3), 
                  labels = c("Y/X < 1", "Y/X > 1")) +
  #scale_fill_discrete(name = "Dose", 
                      #labels = c("Y/X < 1", "Y/X > 1") +
  theme_classic()
ggsave("figures/barplot_table1.png", width = 6, height = 4, dpi = 300)

