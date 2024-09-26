library(tidyverse)
library(mcp)
library(rjags)
########################### Ks step change ##################################
hyphy_anno_reads <- read_csv("data/hyphy_anno_reads.csv.gz")

hyphy_anno_reads %>% filter(seqname.pat =="Y" & XYdS <0.5) %>% 
  ggplot() + geom_point(aes(x = start.pat.mb, y = XYdS)) +
#  geom_step(aes(x = start_tx_mat_mb, y = XYdS)) +
  xlab("X-bearing Haplotype\nposition(Mb)") +
  ylab("Ks") +
  theme_classic()

## define the model
model = list(
  response ~ 1,  # plateau (int_1)
  ~ 0 + time,    # joined slope (time_2) at cp_1
  ~ 1 + time     # disjoined slope (int_3, time_3) at cp_2
)
#get data and fit it
ex = mcp_example("demo")
fit = mcp(model, ex$data)
## run this on a server my god!
plot(fit)

