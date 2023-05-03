###########################

## Retinol PGS tuning and association

## William Reay (2023)

###########################

library(dplyr)
library(data.table)
library(ggplot2)
library(lmvar)
library(ComplexHeatmap)
library(circlize)
library(ggpubr)
library(caret)

setwd("~/Desktop/Retinol_GWAS/TwinsUK_app/E1205_20.01.2023/PGS_calc/")

## Read in subset 2 and PCs - merge: larger tuning set

Subset_2 <- fread("../Subset_2_females_half_of_twin_pairs_TwinsUK.txt")
PCs_2 <- fread("../PCA_calculation/PCA_Subset_2_all_chr_PCA_input.eigenvec")
PCs_2$PublicID <- PCs_2$IID

Merged_ret_PCs_Subset_2 <- merge(Subset_2, PCs_2, by = "PublicID")

## Remove the four individuals with missing retinol at visit 1

Merged_ret_PCs_Subset_2 <- Merged_ret_PCs_Subset_2 %>% filter(!is.na(Retinol_visit1))

Merged_ret_PCs_Subset_2_output_1 <- Merged_ret_PCs_Subset_2 %>% select(PublicID,
                                                                     PC1, PC2,
                                                                     PC3, PC4,
                                                                     PC5, Age_visit1,
                                                                     Batch_visit1,
                                                                     Age_visit2, Batch_visit2,
                                                                     Age_visit3, Batch_visit3)

write.table(Merged_ret_PCs_Subset_2_output_1, file="../Covariates_TwinsUK_subset_2.txt",
            sep = "\t", row.names = F, quote = F)

## Read in PGS info before merging

PGS_ret <- fread("Retinol_IVW_MET_INT_PGS_TwinsUK.all_score")
PGS_ret <- rename(PGS_ret, "Pt_5e_08"="Pt_5e-08",
                  "Pt_1e_05"="Pt_1e-05", "Pt_0_0001"="Pt_0.0001",
                  "Pt_0_01"="Pt_0.01", "Pt_0_05"="Pt_0.05",
                  "Pt_0_1"="Pt_0.1", "PublicID"="IID")

## Merge

PGS_merge_subset_2 <- merge(PGS_ret, Merged_ret_PCs_Subset_2, by = "PublicID")

## Function to define rsq of PGS

PGS_tuning_r2_func <- function(Outcome, Score, Age, Batch, df) {
  Input_df <- df
  Input_df$scaled_score <- as.numeric(scale(Input_df[[Score]]))
  Input_df$scaled_outcome <- as.numeric(scale(Input_df[[Outcome]]))
  baseline_fmla <- as.formula(paste("scaled_outcome ~ + ",  Age, " + ", Batch, " + PC1 + PC2 + PC3 + PC4 + PC5"))
  PGS_fmla <- as.formula(paste("scaled_outcome ~ + ",  Age, " + ", Batch, " + PC1 + PC2 + PC3 + PC4 + PC5 + scaled_score"))
  mod1 <- lm(baseline_fmla, data = Input_df)
  mod2 <- lm(PGS_fmla, data=Input_df)
  mod1_r2 <- summary(mod1)$adj.r.squared
  mod2_r2 <- summary(mod2)$adj.r.squared
  Delta_r2 <- mod2_r2 - mod1_r2
  return(Delta_r2)
}

## Define list of scores to test

Scores_to_test <- as.list(colnames(PGS_merge_subset_2 %>% select(starts_with("Pt"))))

Visit_1 <- sapply(Scores_to_test, PGS_tuning_r2_func, Outcome = "Retinol_visit1", Age = "Age_visit1", Batch = "Batch_visit1", df = PGS_merge_subset_2)
Visit_1 <- as.data.frame(Visit_1)
Visit_2 <- sapply(Scores_to_test, PGS_tuning_r2_func, Outcome = "Retinol_visit2", Age = "Age_visit2", Batch = "Batch_visit2", df = PGS_merge_subset_2)
Visit_2 <- as.data.frame(Visit_2)
Visit_3 <- sapply(Scores_to_test, PGS_tuning_r2_func, Outcome = "Retinol_visit3", Age = "Age_visit3", Batch = "Batch_visit3", df = PGS_merge_subset_2)
Visit_3 <- as.data.frame(Visit_3)

Three_visit_concat <- bind_cols(Visit_1, Visit_2, Visit_3)
Three_visit_concat$Score <- unlist(Scores_to_test)
Three_visit_concat$Mean_r2 <- ((Three_visit_concat$Visit_1 + Three_visit_concat$Visit_2 +
                                 Three_visit_concat$Visit_3)/3)

write.table(Three_visit_concat, file="Tuning_subset_r2_per_visit.txt",
            sep = "\t", row.names = F, quote = F)


# Summarize the results

## Tuning plots

Plot_melt <- melt(Three_visit_concat, id.vars=c("Score"), 
                      measure.vars = c("Visit_1", "Visit_2", "Visit_3", "Mean_r2"))

Plot_melt$P <- c("5e-08", "1e-05", "1e-03", "0.01", "0.05", "0.1", "0.5", "1",
                          "5e-08", "1e-05", "1e-03", "0.01", "0.05", "0.1", "0.5", "1",
                          "5e-08", "1e-05", "1e-03", "0.01", "0.05", "0.1", "0.5", "1",
                          "5e-08", "1e-05", "1e-03", "0.01", "0.05", "0.1", "0.5", "1")

Names <-c("Visit_1"="First visit", "Visit_2"="Second visit",
          "Visit_3"="Third visit", "Mean_r2"="Mean R2 across visits")

Plot_melt$P <- factor(Plot_melt$P, levels = c("5e-08",
                                                      "1e-05", 
                                                      "1e-03", 
                                                      "0.01", 
                                                      "0.05", 
                                                      "0.1", 
                                                      "0.5", 
                                                      "1"))

ggplot(Plot_melt, aes(x = P, y = value, colour = as.factor(Score), group = 1)) +
  facet_wrap(~variable,strip.position="top",nrow=2, scales = "free_y", labeller = as_labeller(Names)) +
  theme_bw() + geom_point(cex = 3) +
  geom_line(colour = "black", linetype = "dashed") +
  theme(legend.position = "none") +
  ylim(0.00001, 0.04) +
  ylab(expression(Delta~R^{"2"})) +
  xlab("P-value threshold") +
  scale_colour_viridis(discrete = T) +
  ggtitle("Retinol PGS tuning")

## Plot mean age per visit

Age_per_visit <- read_excel("../Metabolon_v4_HLI_RETINOL/RETINOL_metabolon_v4_hli_scaleRunDayMedian.xlsx")
Merge_age <- merge(Age_per_visit, PGS_merge_subset_2, by ="PublicID")

ggplot(Merge_age, aes(x=visit, y=age.x, fill = visit)) +
  geom_violin() +
  geom_boxplot(width=0.1, fill = "white") +
  theme_bw() +
  ylab("Age at visit") +
  xlab("Visit") +
  theme(legend.position = "none")

## Validation - subset 1 ##

## Read in subset 1 and PCs - merge: 

Subset_1 <- fread("../Subset_1_females_half_of_twin_pairs_TwinsUK.txt")
PCs_1 <- fread("../PCA_calculation/PCA_Subset_1_all_chr_PCA_input.eigenvec")
PCs_1$PublicID <- PCs_1$IID

Merged_ret_PCs_Subset_1 <- merge(Subset_1, PCs_1, by = "PublicID")

## Remove the 8 IDs with a missing val

Merged_ret_PCs_Subset_1 <- Merged_ret_PCs_Subset_1 %>% filter(!is.na(Retinol_visit1) & !is.na(Retinol_visit2) &
                                                                !is.na(Retinol_visit3))

Merged_ret_PCs_Subset_1_output_1 <- Merged_ret_PCs_Subset_1 %>% select(PublicID,
                                                                       PC1, PC2,
                                                                       PC3, PC4,
                                                                       PC5, Age_visit1,
                                                                       Batch_visit1,
                                                                       Age_visit2, Batch_visit2,
                                                                       Age_visit3, Batch_visit3)

write.table(Merged_ret_PCs_Subset_1_output_1, file="../Covariates_TwinsUK_subset_1.txt",
            sep = "\t", row.names = F, quote = F)

## Merge

PGS_merge_Subset_1 <- merge(PGS_ret, Merged_ret_PCs_Subset_1, by = "PublicID")



## Define list of scores to test

Val_Visit_1 <- sapply(Scores_to_test, PGS_tuning_r2_func, Outcome = "Retinol_visit1", Age = "Age_visit1", Batch = "Batch_visit1", df = PGS_merge_Subset_1)
Val_Visit_1 <- as.data.frame(Val_Visit_1)
Val_Visit_2 <- sapply(Scores_to_test, PGS_tuning_r2_func, Outcome = "Retinol_visit2", Age = "Age_visit2", Batch = "Batch_visit2", df = PGS_merge_Subset_1)
Val_Visit_2 <- as.data.frame(Val_Visit_2)
Val_Visit_3 <- sapply(Scores_to_test, PGS_tuning_r2_func, Outcome = "Retinol_visit3", Age = "Age_visit3", Batch = "Batch_visit3", df = PGS_merge_Subset_1)
Val_Visit_3 <- as.data.frame(Val_Visit_3)

Val_Three_visit_concat <- bind_cols(Val_Visit_1, Val_Visit_2, Val_Visit_3)
Val_Three_visit_concat$Score <- unlist(Scores_to_test)
Val_Three_visit_concat$Mean_r2 <- ((Val_Three_visit_concat$Val_Visit_1 + Val_Three_visit_concat$Val_Visit_2 +
                                      Val_Three_visit_concat$Val_Visit_3)/3)

write.table(Val_Three_visit_concat, file="Validation_subset_r2_per_visit.txt",
            sep = "\t", row.names = F, quote = F)

## Tuning plots

Val_Plot_melt <- melt(Val_Three_visit_concat, id.vars=c("Score"), 
                  measure.vars = c("Val_Visit_1", "Val_Visit_2", "Val_Visit_3", "Mean_r2"))

Val_Plot_melt$P <- c("5e-08", "1e-05", "1e-03", "0.01", "0.05", "0.1", "0.5", "1",
                 "5e-08", "1e-05", "1e-03", "0.01", "0.05", "0.1", "0.5", "1",
                 "5e-08", "1e-05", "1e-03", "0.01", "0.05", "0.1", "0.5", "1",
                 "5e-08", "1e-05", "1e-03", "0.01", "0.05", "0.1", "0.5", "1")

Names <-c("Val_Visit_1"="First visit", "Val_Visit_2"="Second visit",
          "Val_Visit_3"="Third visit", "Mean_r2"="Mean R2 across visits")

Val_Plot_melt$P <- factor(Val_Plot_melt$P, levels = c("5e-08",
                                              "1e-05", 
                                              "1e-03", 
                                              "0.01", 
                                              "0.05", 
                                              "0.1", 
                                              "0.5", 
                                              "1"))

ggplot(Val_Plot_melt, aes(x = P, y = value, colour = as.factor(Score), group = 1)) +
  facet_wrap(~variable,strip.position="top",nrow=2, scales = "free_y", labeller = as_labeller(Names)) +
  theme_bw() + geom_point(cex = 3) +
  geom_line(colour = "black", linetype = "dashed") +
  theme(legend.position = "none") +
  ylim(0.00001, 0.04) +
  ylab(expression(Delta~R^{"2"})) +
  xlab("P-value threshold") +
  scale_colour_viridis(discrete = T) +
  ggtitle("Retinol PGS validation")

## No interaction
## Show corr of retinol PGS with retinol at each visit in validation

Linear_1 <- ggplot(PGS_merge_Subset_1, aes(x=scale(Pt_0_1), y=scale(Retinol_visit1))) +
  geom_point(cex = 1) +
  stat_smooth(method = "lm", formula = y ~ x, size = 1) +
  theme_bw() +
  xlab("Retinol PGS (SD)") +
  ylab("Retinol (SD)") +
  ggtitle("First visit") +
  ylim(-3.5, 6)

Linear_2 <- ggplot(PGS_merge_Subset_1, aes(x=scale(Pt_0_1), y=scale(Retinol_visit2))) +
  geom_point(cex = 1) +
  stat_smooth(method = "lm", formula = y ~ x, size = 1) +
  theme_bw() +
  xlab("Retinol PGS (SD)") +
  ylab("") +
  ggtitle("Second visit") +
  ylim(-3.5, 6)

Linear_3 <- ggplot(PGS_merge_Subset_1, aes(x=scale(Pt_0_1), y=scale(Retinol_visit3))) +
  geom_point(cex = 1) +
  stat_smooth(method = "lm", formula = y ~ x, size = 1) +
  theme_bw() +
  xlab("Retinol PGS (SD)") +
  ylab("") +
  ggtitle("Third visit") +
  ylim(-3.5, 6)

ggarrange(Linear_1, Linear_2, Linear_3, nrow = 1)

## Split retinol into quintiles

PGS_merge_Subset_1$Ret_quint_1 <- ntile(PGS_merge_Subset_1$Retinol_visit1, 5)
PGS_merge_Subset_1$Ret_quint_2 <- ntile(PGS_merge_Subset_1$Retinol_visit2, 5)
PGS_merge_Subset_1$Ret_quint_3 <- ntile(PGS_merge_Subset_1$Retinol_visit3, 5)

Q1 <- ggplot(data = PGS_merge_Subset_1, aes(x=as.factor(Ret_quint_1), y = scale(Pt_0_1), fill = as.factor(Ret_quint_1))) +
  geom_boxplot(width=0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
  theme_bw() +
  theme(legend.position = "none") +
  scale_fill_viridis(discrete = T, option = "magma") +
  ylim(-3.5, 6) +
  ylab("Retinol (SD)") +
  xlab("Retinol PGS quintiles") +
  ggtitle("First visit")

Q2 <- ggplot(data = PGS_merge_Subset_1, aes(x=as.factor(Ret_quint_2), y = scale(Pt_0_1), fill = as.factor(Ret_quint_2))) +
  geom_boxplot(width=0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
  theme_bw() +
  theme(legend.position = "none") +
  scale_fill_viridis(discrete = T, option = "magma") +
  ylim(-3.5, 6) +
  ylab("") +
  xlab("Retinol PGS quintiles") +
  ggtitle("Second visit")

Q3 <- ggplot(data = PGS_merge_Subset_1, aes(x=as.factor(Ret_quint_3), y = scale(Pt_0_1), fill = as.factor(Ret_quint_3))) +
  geom_boxplot(width=0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
  theme_bw() +
  theme(legend.position = "none") +
  scale_fill_viridis(discrete = T, option = "magma") +
  ylim(-3.5, 6) +
  ylab("") +
  xlab("Retinol PGS quintiles") +
  ggtitle("Third visit")

ggarrange(Q1, Q2, Q3, nrow = 1)

## Test effect of retinol PGS per quintile of age

Assoc_per_quintile <- function(quintile, name, df) {
  new_df <- df[df[[name]] == quintile,]
  mod <- lm(scale(Retinol_visit1) ~ Batch_visit1 + scale(Pt_0_1), data = new_df)
  return(summary(mod))
}

Vec_quintile <- c(1,2,3,4,5)

Age_strat <- sapply(Vec_quintile, Assoc_per_quintile, name = "Ret_quint_1",
                    df = PGS_merge_Subset_1)

Age_extract <- apply(Age_strat, 2, function(x) return(as.data.frame(x$coefficients)[6, 1:4]))

## Underpowered to test this

## Test association with coefficient of variation - no assoc

CoV <- lm(Retinol_coefficient_of_variation ~ scale(Pt_0_05) +
            year.birth, data = PGS_merge_Subset_1)


