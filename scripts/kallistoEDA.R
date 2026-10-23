library(tidyverse)

## read in kallisto output
kallisto_df <- read_delim("expression_03_2026/kallisto_out/merged_est_counts.tsv") %>%
  mutate(target_id = str_replace(target_id, "-RA$", ""))


### clean up colnames
kallisto_colnames <- colnames(kallisto_df)

## regex to fix colnames
kallisto_colnames_repaired <- str_extract(kallisto_colnames, "target_id|\\d+[a-z][FM](?=LR)")

### change colnames
colnames(kallisto_df) <- kallisto_colnames_repaired
## prev samples that did not pass QC
rnacolsremove<-c("67eF","14bF","45aF","40aF","24fM","35aM")

# remove col if name matches rnacolsremove
kallisto_df_clean <- kallisto_df %>% select(-(all_of(rnacolsremove))) 



## number of sample columns per sex (target_id is excluded by the case-sensitive match)
sample_cols <- setdiff(colnames(kallisto_df_clean), "target_id")
num_female_cols <- sum(str_detect(sample_cols, "F"))
num_male_cols <- sum(str_detect(sample_cols, "M"))

# sort sample_cols alphanumerically
sample_cols <- gtools::mixedsort(sample_cols)

kallisto_df_sum<-kallisto_df_clean %>%rowwise() %>% mutate( # Calculate mean counts
  femaleCountsTotal = sum(c_across(contains("F", ignore.case = FALSE))),  # Calculate total female counts
  female_mean = (femaleCountsTotal / num_female_cols), # Calculate mean female counts dynamically
  female_SD = stats::sd(c_across(contains("F", ignore.case = FALSE))),  # Calculate standard deviation
  female_SE = plotrix::std.error(c_across(contains("F", ignore.case = FALSE))),  # Calculate standard error
  maleCountsTotal = sum(c_across(contains("M", ignore.case = FALSE))),  # Calculate total male counts
  male_mean = (maleCountsTotal / num_male_cols),  # Calculate mean male counts dynamically
  male_SE = plotrix::std.error(c_across(contains("M", ignore.case = FALSE))),
  totalCountsGene = (femaleCountsTotal+maleCountsTotal)) %>%
  relocate(totalCountsGene, 
           femaleCountsTotal, female_mean, female_SE,
           maleCountsTotal, male_mean, male_SE)

## test: filter just to see summary columns and the indiv column containing "3dM"
kallisto_df_sum %>% select(totalCountsGene, femaleCountsTotal, female_mean, female_SE, maleCountsTotal, male_mean, male_SE, contains("3dM"))


#dplyr remove ### read in orthology information
pgMatOrthsRaw<- read_csv("data/pgMatOrths.csv.gz")  

pgMatOrths<-pgMatOrthsRaw %>% 
  select(pairID:chr_tx_pat) %>%
  unique()
  # filter(flag_tx_mat == "PASS" & flag_tx_pat == "PASS") 
### read in annotation information

## merge pgmatorths and kallisto

#Y genes

#ncol(x_exp)


y_exp <- kallisto_df_sum %>% filter(target_id %in% pgMatOrths$id_tx_pat)
x_exp <- kallisto_df_sum %>% filter(target_id %in% pgMatOrths$id_tx_mat)
df1<-left_join(pgMatOrths,x_exp, 
                             by = c("id_tx_mat"="target_id")) %>%
  left_join(.,y_exp, by = c("id_tx_pat" = "target_id"),
          suffix = c(".mat",".pat")) 

ase_data<-df1 %>% 
  select(c(pairID,contains("M", ignore.case = FALSE ))) %>% 
  column_to_rownames(var ="pairID")

# filter to find where df1 has NA for 3dM.pat
#df1_na_3dM_mat <- df1 %>% filter(is.na(`3dM.mat`))


IDs<-colnames(ase_data)
sorted_col_names <- gtools::mixedsort(IDs)
ase_data <- ase_data[, sorted_col_names]
#colnames(ase_data)
# yay!!
sample<-str_remove(sorted_col_names,"..at") # tm1 tm1 etc.
allele<-str_extract(sorted_col_names,".at$") # mat pat mat pat

colInfo<-data.frame(sorted_col_names,sample,allele) # wooo

ase_matrix<-as.matrix(ase_data)
# drop pairs missing a mat or pat gene in kallisto (NAs introduced by left_join)


# select 3dM in matrix 
ase_matrix <- ase_matrix[, grepl("3dM", colnames(ase_matrix))]

n_before <- nrow(ase_data)
ase_data <- ase_data[complete.cases(ase_data), ]
message("Dropped ", n_before - nrow(ase_data), " pairs with NA counts; ", nrow(ase_data), " remain")
ase_data <- round(ase_data)
colInfo$sample <- factor(colInfo$sample)
colInfo$allele <- factor(colInfo$allele, levels = c("mat", "pat"))
rownames(colInfo) <- colInfo$sorted_col_names
design<- ~sample+allele

library(DESeq2)

dds <- DESeqDataSetFromMatrix(ase_data, colInfo, design)
sizeFactors(dds) <- rep(1, ncol(ase_data))
dds <- DESeq(dds, fitType="local")

resultsNames(dds)

res.allele <- results(dds, name="allele_pat_vs_mat")
head(res.allele$log2FoldChange)
adj_pvals<-res.allele$padj


pvals_genes<-data.frame(res.allele@rownames,res.allele$padj,
                        res.allele$log2FoldChange,
                        res.allele$lfcSE) 


#### compare crossmapping with kallisto vs. previous tools

## how many female reads map to Y?
badgenes<-df1 %>% 
  filter(female_mean.pat > 20) #%>% select(ta) # remove genes with improper female Y mapping
# some may be on Y PAR - which is fromm 0-45MB. Disregard and remove all for now!
badgenes_x <- pgMatOrths %>% select(id_tx_mat,id_tx_pat) %>% 
  filter(id_tx_pat%in%badgenes$id_tx_pat) %>% select(id_tx_mat)

pargenes <- pgMatOrthsRaw %>% filter(start_tx_mat < 450000)
#1500... thats a lot oh no  

`%nin%` <- purrr::negate(`%in%`)

df1_nomismap<-df1 %>% filter(id_tx_mat%nin%badgenes_x$id_tx_mat)  %>% 
  filter(id_tx_pat%nin%badgenes) %>%
  filter(id_tx_mat%nin%pargenes$id_tx_mat)
## 1267 from 2778 genesbadgenes_x## 1267 from 2778 genes
reads_pg_pvals<-inner_join(df1_nomismap, pvals_genes, 
                           by = c("pairID" ="res.allele.rownames")) %>%
  mutate(
    sigTest = # just do significance, leave the log2foldchane to be filtered later
      case_when(res.allele.padj < 0.05 ~ "sig", 
                TRUE ~"nonsig" )) 

  
write_csv(reads_pg_pvals,"data/rna_ase_results_eqtl_kallisto.csv")

library(ggplot2)
library(tidyverse)
library(cowplot)
#library(RColorBrewer)
library(khroma)
df2 <- reads_pg_pvals %>% mutate(sigLFC = case_when(
  res.allele.padj < 0.05 &
    abs(res.allele.log2FoldChange) > 1 ~ "True",
  TRUE ~ "False")) %>%
  mutate(log2male_mean.mat = log2(male_mean.mat)) %>%
  mutate(log2male_mean.pat = log2(male_mean.pat)) %>%
  #filter(male_mean.mat > 10 & male_mean.pat >10) #%>% 
  filter(log2male_mean.pat != "-Inf" & log2male_mean.mat != "-Inf") 


df2$sigLFC<-factor(df2$sigLFC, levels = c("True","False"))

colnames(df2)


figure1a <- ggplot(df2,aes(x=male_mean.mat,y=male_mean.pat)) + 
  geom_point(alpha = 0.5, size =1,aes(color = sigLFC)) +
  geom_errorbar(aes(color = sigLFC,
                    ymin = ((male_mean.pat)-(male_SE.pat)), 
                    ymax = ((male_mean.pat)+(male_SE.pat)))) +
  geom_errorbarh(aes(color = sigLFC,
                     xmin = male_mean.mat-male_SE.mat,
                     xmax = male_mean.mat+male_SE.mat)) +
  geom_abline(slope = 1, intercept = 0, linetype=1) +
  theme_classic() +
  geom_smooth(method ="lm") + 
    xlab("mean expression of X\n gametolog in males") +
  ylab("mean expression of Y\n gametolog in males")# all points regression line
  #geom_smooth(method = "lm", 
#              data = function(x) { filter(x, sigLFC == "True") },
 #             se = T, fullrange = TRUE) +
 

figure1a
## we've just removed X genes that mismap to Y only, so the other situation is not covered...
## keeping the cutoff of 20 read count or greater

  #scale_x_continuous(transform = "log2",limits = c(NA,4000),breaks = c(0.0625,2,64,2048)) +
  #scale_y_continuous(transform = "log2",limits = c(NA,4000),breaks = c(0.0625,2,64,2048)) +
  #labs(color= "p-val < 0.1 & |log2FC| >1 ") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.position = "none") +
  scale_color_manual(values= colors) # from the vibran

ggsave("figures/figure1a_kallisto.png") 
#### compare biases?
reads_pg_pvals_old<-read_csv("data/rna_ase_results_eqtl_sept12.csv.gz")

comparison <- full_join(
  reads_pg_pvals %>%
    transmute(pairID,
              new_lfc = res.allele.log2FoldChange,
              new_padj = res.allele.padj),
  reads_pg_pvals_old %>%
    transmute(pairID,
              old_lfc = res.allele.log2FoldChange,
              old_padj = res.allele.padj),
  by = "pairID"
) %>%
  mutate(
    new_lfc_available = !is.na(new_lfc),
    old_lfc_available = !is.na(old_lfc),
    new_positive = coalesce(
      !is.na(new_padj) & new_padj < 0.05 & !is.na(new_lfc) & new_lfc > 1,
      FALSE
    ),
    old_positive = coalesce(
      !is.na(old_padj) & old_padj < 0.05 & !is.na(old_lfc) & old_lfc > 1,
      FALSE
    ),
    new_negative = coalesce(
      !is.na(new_padj) & new_padj < 0.05 & !is.na(new_lfc) & new_lfc < -1,
      FALSE
    ),
    old_negative = coalesce(
      !is.na(old_padj) & old_padj < 0.05 & !is.na(old_lfc) & old_lfc < -1,
      FALSE
    )
  )

common_results <- comparison %>%
  filter(new_lfc_available, old_lfc_available)

summarise_direction <- function(direction) {
  new_call <- paste0("new_", direction)
  old_call <- paste0("old_", direction)
  new_all <- sum(comparison[[new_call]])
  old_all <- sum(comparison[[old_call]])
  new_shared <- sum(common_results[[new_call]])
  old_shared <- sum(common_results[[old_call]])
  both <- sum(common_results[[new_call]] & common_results[[old_call]])

  tibble(
    direction = direction,
    new_pairs_in_results = n_distinct(reads_pg_pvals$pairID),
    old_pairs_in_results = n_distinct(reads_pg_pvals_old$pairID),
    new_pairs_with_lfc = sum(comparison$new_lfc_available),
    old_pairs_with_lfc = sum(comparison$old_lfc_available),
    shared_pairs_with_lfc = nrow(common_results),
    new_significant_all = new_all,
    old_significant_all = old_all,
    new_significant_shared = new_shared,
    old_significant_shared = old_shared,
    significant_in_both = both,
    new_only_significant_shared = new_shared - both,
    old_only_significant_shared = old_shared - both,
    fraction_new_also_old = if (new_shared > 0) both / new_shared else NA_real_,
    fraction_old_also_new = if (old_shared > 0) both / old_shared else NA_real_,
    call_jaccard = if (new_shared + old_shared - both > 0) {
      both / (new_shared + old_shared - both)
    } else {
      NA_real_
    },
    lfc_spearman_shared = if (nrow(common_results) > 1) {
      cor(common_results$new_lfc, common_results$old_lfc,
          method = "spearman", use = "complete.obs")
    } else {
      NA_real_
    },
    lfc_sign_agreement_shared = if (nrow(common_results) > 0) {
      mean(sign(common_results$new_lfc) == sign(common_results$old_lfc),
           na.rm = TRUE)
    } else {
      NA_real_
    }
  )
}

directional_overlap_summary <- bind_rows(
  summarise_direction("positive"),
  summarise_direction("negative")
)

print(directional_overlap_summary)
write_csv(directional_overlap_summary, "data/kallisto_vs_old_method_summary.csv")

common_results %>%
  mutate(positive_call = case_when(
    new_positive & old_positive ~ "positive_in_both",
    new_positive ~ "positive_new_only",
    old_positive ~ "positive_old_only",
    TRUE ~ "not_positive_in_either"
  ),
  negative_call = case_when(
    new_negative & old_negative ~ "negative_in_both",
    new_negative ~ "negative_new_only",
    old_negative ~ "negative_old_only",
    TRUE ~ "not_negative_in_either"
  )) %>%
  select(pairID, new_lfc, new_padj, old_lfc, old_padj,
         positive_call, negative_call) %>%
  write_csv("data/kallisto_vs_old_method_shared_pairs.csv")



#### new plotts
comparison_plot_data <- common_results %>%
  mutate(
    new_call = case_when(
      new_negative ~ "X higher",
      new_positive ~ "Y higher",
      TRUE ~ "Not called"
    ),
    old_call = case_when(
      old_negative ~ "X higher",
      old_positive ~ "Y higher",
      TRUE ~ "Not called"
    ),
    agreement = case_when(
      new_call == "X higher" & old_call == "X higher" ~ "Both call X higher",
      new_call == "Y higher" & old_call == "Y higher" ~ "Both call Y higher",
      new_call == "Not called" & old_call == "Not called" ~ "Neither calls significant",
      new_call == "Not called" | old_call == "Not called" ~ "Only one method calls significant",
      TRUE ~ "Opposite significant calls"
    )
  )

lfc_comparison_plot <- ggplot(
  comparison_plot_data,
  aes(x = old_lfc, y = new_lfc, color = agreement)
) +
  geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "grey50") +
  geom_hline(yintercept = c(-1, 1), linetype = "dotted", color = "grey50") +
  geom_abline(slope = 1, intercept = 0, color = "black") +
  geom_point(alpha = 0.75, size = 1.5) +
  coord_equal() +
  theme_classic() +
  labs(
    x = "featureCounts log2 fold change (pat vs mat)",
    y = "kallisto log2 fold change (pat vs mat)",
    color = "Significant call",
    title = "Gametolog expression: kallisto vs featureCounts"
  )

lfc_comparison_plot
ggsave("figures/kallisto_vs_featureCounts_lfc.png",
       lfc_comparison_plot, width = 7, height = 6)


df_ratio <- df2 %>%
  mutate(
    log2_yx_ratio = log2male_mean.pat - log2male_mean.mat,
    ratio_call = case_when(
      sigLFC == "True" & res.allele.log2FoldChange < -1 ~
        "Significant: X higher",
      sigLFC == "True" & res.allele.log2FoldChange > 1 ~
        "Significant: Y higher",
      TRUE ~ "Not significant / below LFC cutoff"
    )
  )

ratio_plot <- ggplot(df_ratio, aes(x = log2_yx_ratio, fill = ratio_call)) +
  geom_histogram(bins = 40, color = "white") +
  geom_vline(xintercept = 0, linewidth = 0.6) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dotted",
    color = "grey35"
  ) +
  scale_fill_manual(values = c(
    "Significant: X higher" = "#4C78A8",
    "Significant: Y higher" = "#E07B39",
    "Not significant / below LFC cutoff" = "grey75"
  )) +
  theme_classic() +
  labs(
    x = "log2(mean Y expression / mean X expression)",
    y = "Number of gametolog pairs",
    fill = "DESeq2 result",
    title = "Distribution of Y:X expression ratios",
    subtitle = "Solid line: equal expression; dotted lines: 2-fold expression ratio"
  )

ratio_plot
ggsave("figures/yx_expression_ratio_kallisto.png",
       ratio_plot, width = 7, height = 6)

lfc_plot_data <- reads_pg_pvals %>%
  mutate(
    lfc_call = case_when(
      !is.na(res.allele.padj) & res.allele.padj < 0.05 &
        !is.na(res.allele.log2FoldChange) &
        res.allele.log2FoldChange < -1 ~ "Significant: X higher",
      !is.na(res.allele.padj) & res.allele.padj < 0.05 &
        !is.na(res.allele.log2FoldChange) &
        res.allele.log2FoldChange > 1 ~ "Significant: Y higher",
      TRUE ~ "Not significant / below LFC cutoff"
    )
  ) %>%
  filter(!is.na(res.allele.log2FoldChange))

lfc_distribution_plot <- ggplot(
  lfc_plot_data,
  aes(x = res.allele.log2FoldChange, fill = lfc_call)
) +
  geom_histogram(bins = 40, color = "white") +
  geom_vline(xintercept = 0, linewidth = 0.6) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dotted",
    color = "grey35"
  ) +
  scale_fill_manual(values = c(
    "Significant: X higher" = "#4C78A8",
    "Significant: Y higher" = "#E07B39",
    "Not significant / below LFC cutoff" = "grey75"
  )) +
  theme_classic() +
  labs(
    x = "DESeq2 log2 fold change (pat vs mat)",
    y = "Number of gametolog pairs",
    fill = "DESeq2 result",
    title = "Distribution of kallisto gametolog log2 fold changes",
    subtitle = "Solid line: no allelic difference; dotted lines: LFC cutoffs at -1 and 1"
  )

lfc_distribution_plot
ggsave("figures/kallisto_lfc_distribution.png",
       lfc_distribution_plot, width = 7, height = 6)


## to do next time: explore simulated reads as a way to establish better cutoffs and compare kallisto and featurecounts results
