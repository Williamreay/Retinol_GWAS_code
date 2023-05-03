#############################

## PGS tuning comparison - RapidoPGS via finemapping

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

RapidoPGS_1 <- fread("../RapidoPGS_approach/Retinol_default_SD_prior.all_score")
RapidoPGS_2 <- fread("../RapidoPGS_approach/Retinol_inferred_SD_prior.all_score")

Merge_Rapido <- merge(RapidoPGS_1, RapidoPGS_2, by = "IID")

Merge_Rapido <- rename(Merge_Rapido, "Default_prior_SD_PGS"="Pt_1.x",
                       "Inferred_prior_SD_PGS"="Pt_1.y", "PublicID"="IID")

## Read in Twin data

Twin_df <- fread("../Long_format_full_retinol_TwinsUK.txt")
Twin_df <- Twin_df %>% filter(sex == "F")

PGS_merge_subset <- merge(Merge_Rapido, Twin_df, by = "PublicID")

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

Scores_to_test <- as.list(colnames(PGS_merge_subset %>% dplyr::select(starts_with("Default") | starts_with("Inferred"))))

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

write.table(Three_visit_concat, file = "../RapidoPGS_approach/Training_Mixed_model_RapidoPGS_results.txt",
            sep = "\t", row.names = T, quote = F)

## Testing/validation

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

write.table(Val_Three_visit_concat, file="../RapidoPGS_approach/MM_Validation_subset_r2_per_visit.txt",
            sep = "\t", row.names = F, quote = F)



TE <- lmer(scale(Retinol_visit1) ~ Age_visit1 + Batch_visit1 + (1|FamilyID) +
             scale(Inferred_prior_SD_PGS), data = PGS_merge_subset)

TE_2 <- lmer(scale(Retinol_visit2) ~ Age_visit2 + Batch_visit2 + (1|FamilyID) +
             scale(Inferred_prior_SD_PGS), data =PGS_merge_subset)

TE_3 <- lmer(scale(Retinol_visit3) ~ Age_visit3 + Batch_visit3 + (1|FamilyID) +
             scale(Inferred_prior_SD_PGS), data =PGS_merge_subset)
