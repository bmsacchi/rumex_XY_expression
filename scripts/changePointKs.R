library(tidyverse)
library(mcp)
########################### Ks step change ##################################
ks_mb<-read_csv("data/ks_mb.csv.gz")

ks_mb %>% ggplot() +geom_point(aes(x = start_tx_mat_mb, y = XYdS)) +
#  geom_step(aes(x = start_tx_mat_mb, y = XYdS)) +
  xlab("X-bearing Haplotype\nposition(Mb)") +
  ylab("Ks") +
  theme_classic()
ks_mb %>% ggplot() +geom_point(aes(x = start_tx_pat_mb, y = XYdS)) +
  #  geom_step(aes(x = start_tx_mat_mb, y = XYdS)) +
  xlab("X-bearing Haplotype\nposition(Mb)") +
  ylab("Ks") +
  theme_classic()
