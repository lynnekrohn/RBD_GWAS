# GTEx_Comparisons
Comparing the functional associations of RBD and PD GWAS loci in GTEx v8.  
  
All expression data was extracted from the Genotype-Tissue Expression (GTEx) project (version 8) portal: (https://www.gtexportal.org/home/).  
The top SCARB2 and SNCA/SNCA-AS1 GWAS risk variants for RBD (current project) and PD (Nalls et. al. 2019) were examined. These loci were chosen because out of the 6 GWAS-nominated risk variants for RBD, the SCARB2 and SNCA/SNCA-AS1 were the only ones independent of the top risk loci for PD. 

## Set up:
```R
require(tidyr)
require(dplyr)
require(ggplot2)

setwd("~/Desktop/Projects/GWAS/MR/QTL/")
data = read.csv("GTEx_browser_results.csv", header = T)
data = subset(data, Tissue == "Brain")
attach(data)
```

## Heatmap:
```R
hm <- ggplot(data = data, mapping = aes(x = Concat,
                                                y = Subtissue,
                                                fill = NES)) + geom_tile(color = "black") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_gradient2(low="red", mid="white", high="blue", midpoint=0) +
  geom_text(aes(label = GTEx_Flair))
  
hm = hm + xlab("") + ylab("Brain Tissue") + ggtitle("GTEx v8 Expression") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  labs(fill = "eQTL effect size") 

hm

hm + facet_wrap(~Gene, scales = "free_x")

ggsave("gtex_v8_PDvRBD.png", dpi=300)
```
![GTEx_Heatmap](gtex_v8_PDvRBD.png)
