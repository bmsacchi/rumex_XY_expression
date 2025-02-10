## functions relevant for rumex xy expression project

# invers of %in%
`%nin%` <- purrr:negate(`%in%`)

# create table of log fold change cutoffs
get_obs_table <- 
  function(data, cutoff) {
    data %>%
      filter(res.allele.padj < 0.1, abs(res.allele.log2FoldChange) > cutoff) %>%
      summarise(n_Y_oe = sum(res.allele.log2FoldChange > 0),
                n_X_oe = sum(res.allele.log2FoldChange <= 0)) %>%
      mutate(foldChangeCutoff = cutoff)
  }

# clean raw read counts file
clean_counts <- function(data, columns_to_remove) {
  # Calculate number of columns that contain "F" and "M"
  num_female_cols <- ncol(data %>% select(-all_of(columns_to_remove)) %>%
                            select(contains("F", ignore.case = FALSE)))
  num_male_cols <- ncol(data %>% select(-all_of(columns_to_remove)) %>%
                          select(contains("M", ignore.case = FALSE)))
  num_cols<-num_female_cols+num_male_cols
  
  data %>%
    select(-all_of(columns_to_remove)) %>%  # Remove specified columns
    relocate(1:6, ends_with("F"), ends_with("M")) %>%  # Sort columns
    mutate(Chr = case_when(  # Relabel chromosomes
      grepl("X", Chr) ~ "X",
      grepl("Y", Chr) ~ "Y",
      grepl("A1", Chr) ~ "A1",
      grepl("A2", Chr) ~ "A2",
      grepl("A3", Chr) ~ "A3",
      grepl("A4", Chr) ~ "A4",
      TRUE ~ NA_character_
    )) %>%
    rowwise() %>%
    mutate( # Calculate mean counts
      femaleCountsTotal = sum(c_across(contains("F", ignore.case = FALSE))),  # Calculate total female counts
      female_mean = (femaleCountsTotal / num_female_cols), # Calculate mean female counts dynamically
      female_SD = stats::sd(c_across(contains("F", ignore.case = FALSE))),  # Calculate standard deviation
      female_SE = plotrix::std.error(c_across(contains("F", ignore.case = FALSE))),  # Calculate standard error
      maleCountsTotal = sum(c_across(contains("M", ignore.case = FALSE))),  # Calculate total male counts
      male_mean = (maleCountsTotal / num_male_cols),  # Calculate mean male counts dynamically
      male_SE = plotrix::std.error(c_across(contains("M", ignore.case = FALSE))),
      totalCountsGene = (femaleCountsTotal+maleCountsTotal),
      meanCounts = totalCountsGene / num_cols) %>%
    relocate(1:6, totalCountsGene, meanCounts, 
             femaleCountsTotal, female_mean, female_SE,
             maleCountsTotal, male_mean, male_SE)  # Reorder columns
}

