############################

## Retinol - SNP h2: LDSR vs LDAK plots

## William Reay (2022)

###########################

library(ggplot2)
library(ggpubr)
library(readxl)
library(RColorBrewer)

DF <- read_excel("~/Desktop/Retinol_GWAS/h2_est/h2_plot_LDSR_LDAK.xlsx")

h2 <- ggplot(DF, aes(x=GWAS, y=h2, fill=LD)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=h2-SE, ymax=h2+SE), width=.2,
                position=position_dodge(.9)) +
  facet_wrap(~Model,strip.position="top",ncol=3,scales = "free_y") +
  theme_bw() +
  labs(y=expression(paste(h^2, " (+/- SE)"))) +
  theme( axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  ylim(0,0.2) +
  ggtitle("Circulating all-trans retinol")

h2 + scale_fill_manual(values=c("#999999", "#E69F00"))

                        