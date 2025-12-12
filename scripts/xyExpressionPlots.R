library(ggplot2)
library(tidyverse)
library(cowplot)
#library(RColorBrewer)
library(khroma)

## use palette
vibrant<- color("vibrant")
#vibrant
#vibrant(2)
vibrant()
reads_pg_pvals<-read_csv("data/rna_ase_results_eqtl_sept12.csv.gz")

reads_pg_pvals2 <- reads_pg_pvals %>% mutate(sigLFC = case_when(
    res.allele.padj < 0.1 &
      abs(res.allele.log2FoldChange) > 1 ~ "True",
    TRUE ~ "False")) %>%
  mutate(log2male_mean.mat = log2(male_mean.mat)) %>%
  mutate(log2male_mean.pat = log2(male_mean.pat)) %>%
  #filter(male_mean.mat > 10 & male_mean.pat >10) #%>% 
  filter(log2male_mean.pat != "-Inf" & log2male_mean.mat != "-Inf") 

  
  
reads_pg_pvals2$sigLFC<-factor(reads_pg_pvals2$sigLFC, levels = c("True","False"))
### retransforming the data bc I hate myself!
colnames(reads_pg_pvals)

df<-reads_pg_pvals %>% select(c(1:9,23:26,28:176,184:332))  %>% mutate(sigLFC = case_when(
  res.allele.padj < 0.1 &
    abs(res.allele.log2FoldChange) > 1 ~ "True",
  TRUE ~ "False"))

readcounts_long<- df %>% pivot_longer(cols = 14:311, names_to = c("sample",".value"),
                    names_sep = "\\.") %>% 
  mutate(sex = case_when(grepl("M",sample) ~ "Male", 
                       grepl("F", sample) ~ "Female", TRUE ~ NA_character_)) 


sum_filt_readcounts_long <- readcounts_long %>% 
  filter(sex == "Male") %>%
  #filter(pat >=0 & mat >=0) %>%
  group_by(pairID,sigTest,sigLFC) %>% 
  summarise(n = n(), totalReadCount = sum(pat+mat), male_mean.pat = mean(pat),male_mean.mat = mean(mat), 
            male_SE.pat = plotrix::std.error(pat),male_SE.mat = plotrix::std.error(mat))# %>%
  filter(male_mean.pat !="-Inf")
readcounts_long$sigLFC<-factor(readcounts_long$sigLFC, levels = c("True","False"))
  

colors <- vibrant(7)[c(4,7)]

figure1a <- ggplot(reads_pg_pvals2,aes(x=male_mean.mat,y=male_mean.pat)) + 
  #geom_point(aes(size = res.allele.lfcSE),alpha = 0.60) +
  geom_point(alpha = 0.5, size =1,aes(color = sigLFC)) +
  geom_errorbar(aes(color = sigLFC,ymin = ((male_mean.pat)-(male_SE.pat)), ymax = ((male_mean.pat)+(male_SE.pat)))) +
  geom_errorbarh(aes(color = sigLFC,xmin = male_mean.mat-male_SE.mat, xmax = male_mean.mat+male_SE.mat)) +
  geom_abline(slope = 1, intercept = 0, linetype=1) +
  theme_classic()+
  geom_smooth(method ="lm", color =  vibrant(7)[7]) + # all points regression line
  geom_smooth(method = "lm", data = function(x) { filter(x, sigLFC == "True") }, se = T, fullrange = TRUE, color =  vibrant(7)[4]) +
  # regression sig only
  xlab("mean expression of X\n gametolog in males") +
  ylab("mean expression of Y\n gametolog in males") +
  scale_x_continuous(transform = "log2",limits = c(NA,4000),breaks = c(0.0625,2,64,2048)) +
  scale_y_continuous(transform = "log2",limits = c(NA,4000),breaks = c(0.0625,2,64,2048)) +
  #labs(color= "p-val < 0.1 & |log2FC| >1 ") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.position = "none") +
  scale_color_manual(values= colors) # from the vibrant palette, see vignette("tol")
## not the long version
## no filter
figure1a
ggsave("figures/Figure1a_2024Jan.png", h =6, w =10)
## this is the final figure

#######
# an old attempt to remake figure with long format data for some reason
#######
figure1a_long <- ggplot(sum_filt_readcounts_long,aes(x=male_mean.mat,y=male_mean.pat)) +
  geom_point(alpha = 0.5, size =1,aes(color = sigLFC)) +
  geom_errorbar(aes(color = sigLFC,ymin = ((male_mean.pat)-(male_SE.pat)), ymax = ((male_mean.pat)+(male_SE.pat)))) +
  geom_errorbarh(aes(color = sigLFC,xmin = male_mean.mat-male_SE.mat, xmax = male_mean.mat+male_SE.mat)) +
  geom_abline(slope = 1, intercept = 0, linetype=1) +
  theme_bw()+
  geom_smooth(method ="lm", color =  vibrant(7)[7]) + # all points regression line
  geom_smooth(method = "lm", data = function(x) { filter(x, sigLFC == "True") }, se = T, fullrange = TRUE, color =  vibrant(7)[3]) +
  xlab("mean expression of X\n gametolog in males") +
  ylab("mean expression of Y\n gametolog in males") +
  scale_x_continuous(transform = "log2",limits = c(0.001,4000),breaks = c(0.0625,2,64,2048)) +
  scale_y_continuous(transform = "log2",limits = c(0.001,4000),breaks = c(0.0625,2,64,2048)) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.position = "none") +
  scale_color_manual(values= colors) # from the vibrant palette, see vignette("tol")
  
  
figure1a_long
cowplot::plot_grid(figure1a,figure1a_long)

#### plot with filtering samples

sum_filt_readcounts_long_filt <- readcounts_long %>% 
  filter(sex == "Male") %>%
  filter(pat >=10 & mat >=10) %>%
  group_by(pairID,sigTest,sigLFC) %>% 
  summarise(n = n(), totalReadCount = sum(pat+mat), male_mean.pat = mean(pat),male_mean.mat = mean(mat), 
            male_SE.pat = plotrix::std.error(pat),male_SE.mat = plotrix::std.error(mat)) %>%
  filter(n >=10)

figure1a_long_filt <- ggplot(sum_filt_readcounts_long_filt,aes(x=male_mean.mat,y=male_mean.pat)) +
  geom_point(alpha = 0.5, size =1,aes(color = sigLFC)) +
  geom_errorbar(aes(color = sigLFC,ymin = ((male_mean.pat)-(male_SE.pat)), ymax = ((male_mean.pat)+(male_SE.pat)))) +
  geom_errorbarh(aes(color = sigLFC,xmin = male_mean.mat-male_SE.mat, xmax = male_mean.mat+male_SE.mat)) +
  geom_abline(slope = 1, intercept = 0, linetype=1) +
  theme_bw()+
  geom_smooth(method ="lm", color =  vibrant(7)[7]) + # all points regression line
  geom_smooth(method = "lm", data = function(x) { filter(x, sigLFC == "True") }, se = T, fullrange = TRUE, color =  vibrant(7)[3]) +
  xlab("mean expression of X\n gametolog in males") +
  ylab("mean expression of Y\n gametolog in males") +
  scale_x_continuous(transform = "log2",limits = c(9,4000),breaks = c(0.0625,2,64,2048)) +
  scale_y_continuous(transform = "log2",limits = c(9,4000),breaks = c(0.0625,2,64,2048)) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.position = "none") +
  scale_color_manual(values= colors) 


cowplot::plot_grid(figure1a_long,figure1a_long_filt)
figure1a_long_filt

# what are the R^2 etc
### ok just keep what you're originally doing, and describe properly in the results.
## save copies for supmat


##### regression of all points, and sig points #####


colnames(reads_pg_pvals2) # rm NaN

lm_all<-lm(formula = log2male_mean.pat ~ log2male_mean.mat, data = reads_pg_pvals2) 
#plot(ks_lm)
lm_sig<-lm(formula = log2male_mean.pat ~ log2male_mean.mat, data = reads_pg_pvals2 %>% filter(sigLFC == "True"))
summary(lm_all)
summary(lm_sig)


########## sanity check - normalized counts #############
##### using the norm factors from the dosage analyses (each sample has a norm factor)
# doesn't change anything but the standard error bars
# love et al. indicate norm not needed for this method

#### norm counts using size factors from dosage analysis
pvals_genes_norm <- read_csv("data/norm_rna_ase_resultsOnlySept2024.csv.gz")

## normalized counts

normCounts<-read_csv("data/normalizedCounts_xyExp.csv.gz")

## unite norm counts with pvals
reads_pg_pvals_norm<-inner_join(normCounts, pvals_genes_norm, 
                                by = c("pairID" ="res.allele.rownames"), keep=TRUE) %>%
  mutate(
    sigTest = # just do significance, leave the log2foldchane to be filtered later
      case_when(res.allele.padj < 0.1 ~ "sig", 
                TRUE ~"nonsig" ))
### summarize
get_obs_table(reads_pg_pvals_norm, 0.5)
fold_change_cutoffs <- c(0.5, 1, 1.5, 2)
obs_tables <- lapply(fold_change_cutoffs, function(cutoff) {
  get_obs_table((filter(reads_pg_pvals_norm, res.allele.padj < 0.1)), cutoff)
})

combined_table <- do.call(rbind, obs_tables)
print(combined_table) ## same results, not surprisingly
## test occurs at the sample level, normalized counts won't change that

reads_pg_pvals2_norm <- reads_pg_pvals_norm %>% mutate(sigLFC = case_when(
  res.allele.padj < 0.1 &
    abs(res.allele.log2FoldChange) > 1 ~ "True",
  TRUE ~"False" )) #%>%
mutate(log2male_mean.mat = log2(male_mean.mat)) %>%
  mutate(log2male_mean.pat = log2(male_mean.pat))

reads_pg_pvals2_norm$sigLFC<-factor(reads_pg_pvals2_norm$sigLFC, levels = c("True","False"))

figure1a_norm <- 
  ggplot(reads_pg_pvals2_norm,aes(x=matMean,y=patMean, color = sigLFC)) + 
  #geom_point(aes(size = res.allele.lfcSE),alpha = 0.60) +
  geom_point(alpha = 0.5, size =1) +
  geom_errorbar(aes(ymin = patMean-patSE, ymax = patMean+patSE)) +
  geom_errorbarh(aes(xmin = matMean-matSE, xmax =matMean+matSE)) +
  geom_abline(slope = 1, intercept = 0, linetype=1) +
  theme_bw()+
  geom_smooth(method ="lm") +
  xlab("mean normalized expression of X\n gametolog in males") +
  ylab("mean normalized expression of Y\n gametolog in males") +
  scale_x_continuous(transform = "log2",limits = c(NA,4000),breaks = c(0.0625,2,64,2048)) +
  scale_y_continuous(transform = "log2",limits = c(NA,4000),breaks = c(0.0625,2,64,2048)) +
  labs(color= "p-val < 0.1 & |log2FC| >1 ", size = "log2FC std. error") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.text = element_text(size =12),
        legend.title = element_text(size =12)) 
figure1a_norm
#ggsave("figures/Figure1_gametologExpressionNormalized.png", h =6, w =10)

figure1a_norm
figure1a
######

