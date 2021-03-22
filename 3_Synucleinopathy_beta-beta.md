# Synucleinopathy_comparisons

PD is the most genetically understood of the synucleinopathies (PD, DLB, MSA, RBD) because of the large sample sizes available for GWAS. For this reason, I've compared the results for the PD top GWAS hits (Nalls et. al. 2019) in both RBD (current) and LBD (Chia et. al. 2021) GWAS.  

## Beta-Beta Plot

**Prepare:**
```R
setwd("~/Desktop/Projects/GWAS/beta-beta/")
data <- read.csv("RBD_GWAS_meta-5-comparisons.csv")
names(data)
attach(data)

require("ggplot2")
require("reshape2")
require("corrplot")
require("ggrepel")
```

### Create plots
**Idiopathic RBD:**
```R
png('iRBD_GWAS_beta-beta_tricolor.png', units="in", width=6, height=5, res=300, compression = 'lzw')

irbd <- ggplot(data, aes(x=Beta_Meta5, y=meta5_beta_adjusted_iRBD)) + geom_vline(xintercept = 0, colour = "grey") + 
  geom_hline(yintercept = 0, colour = "grey") +
 geom_point(shape = 16, size = 4, aes(color = Direction_compare_iRBD), alpha=0.7) + 
  xlim(-0.8, 0.8) 

irbd2 <- irbd + scale_color_manual(values = c("red2", "blue", "grey65"))

irbd3 <- irbd2 + geom_vline(xintercept = 0, colour = "grey") + geom_hline(yintercept = 0, colour = "grey") + theme_light()

irbd4 <- irbd3 +
  geom_label_repel(aes(label = NAME),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') 

irbd5 <- irbd4 + ggtitle("Beta-Beta Plot: PD vs iRBD") +
  xlab("PD Betas (Nalls et. al. 2019)") + ylab("RBD Betas (iRBD GWAS)") + theme(plot.title = element_text(face="bold"))

irbd5 + theme(legend.position = "none") + theme(plot.title = element_text(hjust = 0.5)) 

dev.off()
```
![iRBD beta-beta](iRBD_GWAS_beta-beta_tricolor.png)

**23andMe cohort: PD+RBD**
```R
png('23andMe_GWAS_beta-beta_noLRRK2_tricolor.png', units="in", width=6, height=5, res=300)

pdwrbd <- ggplot(data, aes(x=Beta_Meta5, y=meta5_beta_adjusted_23andMe)) + geom_vline(xintercept = 0, colour = "grey") + 
  geom_hline(yintercept = 0, colour = "grey") +
  geom_point(shape = 16, size = 4, aes(color = direction_23andMe), alpha=0.7) + xlim(-0.8, 0.8) 


pdwrbd2 <- pdwrbd + scale_color_manual(values = c("blue", "grey65")) + theme_light()

pdwrbd3 <- pdwrbd2 + theme(legend.position = "none") + geom_vline(xintercept = 0, colour = "grey") + geom_hline(yintercept = 0, colour = "grey")

pdwrbd4 <- pdwrbd3 +
  geom_label_repel(aes(label = NAME),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') 

pdwrbd5 <- pdwrbd4 + ggtitle("Beta-Beta Plot: PD vs PD+RBD") +
  xlab("PD Betas (Nalls et. al. 2019)") + ylab("PD+RBD Betas (23andMe)") + theme(plot.title = element_text(face="bold"))

pdwrbd5 + theme(plot.title = element_text(hjust = 0.5)) 


dev.off()
````
![PDwRBD beta-beta](23andMe_PDwRBD_GWAS_beta-beta_tricolor.png)

**LBD:**
```R
png('DLB_GWAS_beta-beta.png', units="in", width=6, height=5, res=300)

dlb <- ggplot(data, aes(x=Beta_Meta5, y=Beta_meta5_adjusted_DLB)) + geom_vline(xintercept = 0, colour = "grey") + 
  geom_hline(yintercept = 0, colour = "grey") +
  geom_point(shape = 16, size = 4, aes(color = Direction_compare_DLB), alpha=0.7) + 
  xlim(-0.8, 0.8)

dlb2 <- dlb + scale_color_manual(values = c("red2", "blue", "grey65")) + theme_light()

dlb3 <- dlb2 + theme(legend.position = "none") + geom_vline(xintercept = 0, colour = "grey") + geom_hline(yintercept = 0, colour = "grey") 

dlb4 <- dlb3 +
  geom_label_repel(aes(label = NAME),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') 

dlb5 <- dlb4 + ggtitle("Beta-Beta Plot: PD vs DLB") +
  xlab("PD Betas (Nalls et. al. 2019)") + ylab("DLB Betas (FILL IN et. al. 2020") + theme(plot.title = element_text(face="bold"))

dlb5 + theme(plot.title = element_text(hjust = 0.5)) 

dev.off()
```
![DLB beta-beta](LBD_beta-beta.jpg)


