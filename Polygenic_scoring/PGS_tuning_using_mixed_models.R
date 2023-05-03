#############################

## PGS tuning comparison

## Mixed-model with family level random intercepts and modelling increased correlation of MZ vs DZ twins

## Dr William Reay (2023)

#############################

library(dplyr)
library(data.table)
library(ggplot2)
library(lmvar)
library(ComplexHeatmap)
library(circlize)
library(ggpubr)
library(caret)
library(lme4)
library(sjstats)

setwd("~/Desktop/Retinol_GWAS/TwinsUK_app/E1205_20.01.2023/PGS_calc/")

set.seed(47373)

## Read in PGS info before merging

PGS_ret <- fread("Retinol_IVW_MET_INT_PGS_TwinsUK.all_score")
PGS_ret <- rename(PGS_ret, "Pt_5e_08"="Pt_5e-08",
                  "Pt_1e_05"="Pt_1e-05", "Pt_0_0001"="Pt_0.0001",
                  "Pt_0_01"="Pt_0.01", "Pt_0_05"="Pt_0.05",
                  "Pt_0_1"="Pt_0.1", "PublicID"="IID")

## Read in Twin data

Twin_df <- fread("../Long_format_full_retinol_TwinsUK.txt")
Twin_df <- Twin_df %>% filter(sex == "F")

PGS_merge_subset <- merge(PGS_ret, Twin_df, by = "PublicID")

## Split into 70/30 training/test split

smp_size <- floor(0.70 * nrow(PGS_merge_subset))

train_ind <- sample(seq_len(nrow(PGS_merge_subset)), size = smp_size)

PGS_train <- PGS_merge_subset[train_ind, ]
PGS_test <- PGS_merge_subset[-train_ind, ]

## Function to define marginal rsq of PGS

PGS_tuning_r2_func_MM <- function(Outcome, Score, Age, Batch, df) {
  Input_df <- df
  Input_df$scaled_score <- as.numeric(scale(Input_df[[Score]]))
  Input_df$scaled_outcome <- as.numeric(scale(Input_df[[Outcome]]))
  baseline_fmla <- as.formula(paste("scaled_outcome ~ + ",  Age, " + ", Batch, " + (1|FamilyID)"))
  PGS_fmla <- as.formula(paste("scaled_outcome ~ + ",  Age, " + ", Batch, " + (1|FamilyID) + scaled_score"))
  mod1 <- lmer(baseline_fmla, data = Input_df)
  mod2 <- lmer(PGS_fmla, data=Input_df)
  mod1_r2 <- r2(mod1)$R2_marginal
  mod2_r2 <- r2(mod2)$R2_marginal
  Delta_r2 <- mod2_r2 - mod1_r2
  return(Delta_r2)
}

## Training 

Scores_to_test <- as.list(colnames(PGS_merge_subset %>% dplyr::select(starts_with("Pt_"))))

Visit_1 <- sapply(Scores_to_test, PGS_tuning_r2_func_MM, Outcome = "Retinol_visit1", Age = "Age_visit1", Batch = "Batch_visit1", df = PGS_train)
Visit_1 <- as.data.frame(Visit_1)
Visit_2 <- sapply(Scores_to_test, PGS_tuning_r2_func_MM, Outcome = "Retinol_visit2", Age = "Age_visit2", Batch = "Batch_visit2", df =  PGS_train)
Visit_2 <- as.data.frame(Visit_2)
Visit_3 <- sapply(Scores_to_test, PGS_tuning_r2_func_MM, Outcome = "Retinol_visit3", Age = "Age_visit3", Batch = "Batch_visit3", df =  PGS_train)
Visit_3 <- as.data.frame(Visit_3)

Three_visit_concat <- bind_cols(Visit_1, Visit_2, Visit_3)
Three_visit_concat$Score <- unlist(Scores_to_test)
Three_visit_concat$Mean_r2 <- ((Three_visit_concat$Visit_1 + Three_visit_concat$Visit_2 +
                                  Three_visit_concat$Visit_3)/3)

write.table(Three_visit_concat, file="MIXED_MODEL_Tuning_subset_r2_per_visit.txt",
            sep = "\t", row.names = F, quote = F)

## Testing/validation

## Define list of scores to test

Val_Visit_1 <- sapply(Scores_to_test, PGS_tuning_r2_func_MM, Outcome = "Retinol_visit1", Age = "Age_visit1", Batch = "Batch_visit1", df = PGS_test)
Val_Visit_1 <- as.data.frame(Val_Visit_1)
Val_Visit_2 <- sapply(Scores_to_test, PGS_tuning_r2_func_MM, Outcome = "Retinol_visit2", Age = "Age_visit2", Batch = "Batch_visit2", df = PGS_test)
Val_Visit_2 <- as.data.frame(Val_Visit_2)
Val_Visit_3 <- sapply(Scores_to_test, PGS_tuning_r2_func_MM, Outcome = "Retinol_visit3", Age = "Age_visit3", Batch = "Batch_visit3", df = PGS_test)
Val_Visit_3 <- as.data.frame(Val_Visit_3)

Val_Three_visit_concat <- bind_cols(Val_Visit_1, Val_Visit_2, Val_Visit_3)
Val_Three_visit_concat$Score <- unlist(Scores_to_test)
Val_Three_visit_concat$Mean_r2 <- ((Val_Three_visit_concat$Val_Visit_1 + Val_Three_visit_concat$Val_Visit_2 +
                                      Val_Three_visit_concat$Val_Visit_3)/3)

write.table(Val_Three_visit_concat, file="MM_Validation_subset_r2_per_visit.txt",
            sep = "\t", row.names = F, quote = F)

## Make mixed model tuning plots

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
  ggtitle("Retinol PGS validation (mixed model - 30% split)")
