##################################################

## Testing associations between PGS  with supra and infra normal individuals 

## Visit 1,2,3

## RapidoPGS

## No BMI in normative model (GAMLSS - BCT family)

## Dr William Reay (2023)

##################################################

library(data.table)
library(dplyr)
library(ggplot2)

setwd("~/Desktop/Retinol_GWAS/TwinsUK_app/E1205_20.01.2023/")

## Inferred SD

PGS <- fread("RapidoPGS_approach/Retinol_inferred_SD_prior.all_score")
PGS$publicid <- PGS$IID


## PCs

PCs <- fread("PCA_calculation/PCA_Subset_1_all_chr_PCA_input.eigenvec")

Val_subset <- merge(PGS, PCs, by = "IID")



Long_ret <- fread("Subset_1_females_half_of_twin_pairs_TwinsUK.txt")

Long_ret$IID <- Long_ret$PublicID

Long_ret <- Long_ret %>% filter(!is.na(Retinol_visit1)) %>% filter(!is.na(Retinol_visit2)) %>% filter(!is.na(Retinol_visit3))  

Val_subset <- merge(Val_subset, Long_ret, by = "IID")

Val_subset$Scaled_PGS <- as.numeric(scale(Val_subset$Pt_1))

Val_subset$Scaled_measured_retinol_1 <- as.numeric(scale(Val_subset$Retinol_visit1))
Val_subset$Scaled_measured_retinol_2 <- as.numeric(scale(Val_subset$Retinol_visit2))
Val_subset$Scaled_measured_retinol_3 <- as.numeric(scale(Val_subset$Retinol_visit3))
## Read in deviations

Dev_1 <- fread("Maria_normative_modelling/Visit_1_deviations_BCT.txt")
Dev_1$IID <- Dev_1$publicid

Merged_final_1 <- merge(Val_subset, Dev_1, by = "IID")

## Infra-normal and supra-normal PGS

PGS_infra <- glm(band_infra ~
                   PC1 + PC2 + PC3 + PC4 + PC5 + Scaled_PGS,
                 data = Merged_final_1, family = "binomial")

PGS_supra <- glm(band_supra ~
                   PC1 + PC2 + PC3 + PC4 + PC5 + Scaled_PGS,
                 data = Merged_final_1, family = "binomial")

## Not great performance for visit 1

## Visit 2

Dev_2 <- fread("Maria_normative_modelling/Visit_2_deviations_BCT.txt")
Dev_2$IID <- Dev_2$publicid

Merged_final_2 <- merge(Val_subset, Dev_2, by = "IID")

## Infra-normal and supra-normal PGS

PGS_infra_2 <- glm(band_infra ~
                     PC1 + PC2 + PC3 + PC4 + PC5 + Scaled_PGS,
                   data = Merged_final_2, family = "binomial")

PGS_supra_2 <- glm(band_supra ~
                     PC1 + PC2 + PC3 + PC4 + PC5 + Scaled_PGS,
                   data = Merged_final_2, family = "binomial")

## Visit 3

Dev_3 <- fread("Maria_normative_modelling/Visit_3_deviations_BCT.txt")
Dev_3$IID <- Dev_3$publicid

Merged_final_3 <- merge(Val_subset, Dev_3, by = "IID")

## Infra-normal and supra-normal PGS

PGS_infra_3 <- glm(band_infra ~
                     PC1 + PC2 + PC3 + PC4 + PC5 + Scaled_PGS,
                   data = Merged_final_3, family = "binomial")

PGS_supra_3 <- glm(band_supra ~
                     PC1 + PC2 + PC3 + PC4 + PC5 + Scaled_PGS,
                   data = Merged_final_3, family = "binomial")

## Make plot of deviations OR in validation subset 

Plot_df <- as.data.frame(matrix(nrow = 6, ncol = 0))
Plot_df$Name <- c("Supra-normal", "Infra-normal","Supra-normal", "Infra-normal",
                  "Supra-normal", "Infra-normal")
Plot_df$Visit <- c("First visit", "First visit", "Second visit", "Second visit",
                   "Third visit", "Third visit")
Plot_df$Beta <- c(0.4932, -0.1558,  0.5853,-0.3680,
                  0.4322, -0.08099)
Plot_df$SE <- c(0.1797, 0.1561, 0.1763, 0.1533,
                0.2048, 0.14764)

Plot_df$OR <- exp(Plot_df$Beta)

Plot_df$LOR <- exp(Plot_df$Beta - (1.96*Plot_df$SE))

Plot_df$UOR <- exp(Plot_df$Beta + (1.96*Plot_df$SE))

write.table(Plot_df, file="Normative_modelling_PGS_integration/Half_subset_RapidoPGS_validation_prediction_of_infra_and_supra_normal_for_age.txt",
            sep = "\t", row.names = F, quote = F)
ggplot(Plot_df,aes(x=Name, y=OR, ymin=LOR, ymax=UOR, colour = Name)) +
  geom_pointrange() +
  coord_flip() +
  facet_wrap(~Visit,strip.position="top",nrow=3,scales = "free_y") +
  theme_bw() +
  theme(legend.position = "none") +
  geom_hline(yintercept = 1, lty=2) +
  ylab("OR [95% CI] of falling outside normative retinol for age per SD in PGS") +
  xlab(" ") +
  scale_colour_brewer(palette = "Dark2")

## Retinol + beta-carotene

FFQ <- fread("Processed_FFQ_data_within_a_year/Females_only_within_a_year_visit.txt")

Merged_FFQ_1 <- merge(FFQ, Dev_1, by = "publicid")
Merged_FFQ_1$Retinol_carotene <- Merged_FFQ_1$retinol_ug+Merged_FFQ_1$carotene_ug

## Infra-normal and supra-normal PGS

FFQ_infra <- glm(band_infra ~
                   scale(Retinol_carotene),
                 data = Merged_FFQ_1, family = "binomial")

FFQ_supra <- glm(band_supra ~
                   scale(Retinol_carotene),
                 data = Merged_FFQ_1, family = "binomial")

## Not great performance for visit 1


## Vist 2

Merged_FFQ_2 <- merge(FFQ, Dev_2, by = "publicid")

## Infra-normal and supra-normal PGS

FFQ_infra_2 <- glm(band_infra ~
                     scale(Retinol_carotene),
                   data = Merged_FFQ_2, family = "binomial")

FFQ_supra_2 <- glm(band_supra ~
                     scale(retinol_ug),
                   data = Merged_FFQ_2, family = "binomial")

## Vist 3

Merged_FFQ_3 <- merge(FFQ, Dev_3, by = "publicid")

## Infra-normal and supra-normal PGS

FFQ_infra_3 <- glm(band_infra ~
                     scale(retinol_ug),
                   data = Merged_FFQ_3, family = "binomial")

FFQ_supra_3 <- glm(band_supra ~
                     scale(retinol_ug),
                   data = Merged_FFQ_3, family = "binomial")