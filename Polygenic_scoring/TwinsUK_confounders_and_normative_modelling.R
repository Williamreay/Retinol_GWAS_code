#############################

## TwinsUK retinol data - effect of confounders and normative modelling

## Dr William Reay (2023)

#########################

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(cowplot)
library(gridExtra)
library(ggpubr)

setwd("~/Desktop/Retinol_GWAS/TwinsUK_app/E1205_20.01.2023/")

Ret_dat <- fread("Long_format_full_retinol_TwinsUK.txt")

## Retain only females

Ret_dat <- Ret_dat %>% filter(sex == "F")

## Split cohort into two subsets based on public ID last digit so each twin pair is separated - use code below to confirm this will work as anticipated

Ret_dat %>% filter(Last_digit == "1") %>% select(FamilyID) %>% group_by_all() %>% summarise(COUNT = n()) %>% arrange(desc(COUNT)) %>% filter(COUNT > 1)
Ret_dat %>% filter(Last_digit == "2") %>% select(FamilyID) %>% group_by_all() %>% summarise(COUNT = n()) %>% arrange(desc(COUNT)) %>% filter(COUNT > 1)

## Split into subsets

Subset_1 <- Ret_dat %>% filter(Last_digit == "1")
Subset_2 <- Ret_dat %>% filter(Last_digit == "2")

write.table(Subset_1, "Subset_1_females_half_of_twin_pairs_TwinsUK.txt",
            sep = "\t", row.names = F, quote = F)

write.table(Subset_2, "Subset_2_females_half_of_twin_pairs_TwinsUK.txt",
            sep = "\t", row.names = F, quote = F)
## Model effects of age - linear, quadratic and cubic - scale retinol outcome to SD units

## Test linear first

Linear_age_effect <- function(Retinol, Age, Batch, df) {
  ## Remove missing variables
  cleaned_df <- df[!is.na(df[[Retinol]]),]
  cleaned_df$scaled_retinol <- as.numeric(scale(cleaned_df[[Retinol]]))
  fmla  <- as.formula(paste("scaled_retinol ~ ", Age, " + ", Batch, sep = ""))
  mod <- lm(fmla, data = cleaned_df)
  Mod_sum <- summary(mod)
  return(Mod_sum)
}

List_age <- as.list(c("Age_visit1", "Age_visit2", "Age_visit3"))
List_retinol <- as.list(c("Retinol_visit1", "Retinol_visit2", "Retinol_visit3"))
List_batch <- as.list(c("Batch_visit1", "Batch_visit2", "Batch_visit3"))

Age_linear_subset1 <- mapply(Linear_age_effect, Retinol = List_retinol, Age = List_age, Batch = List_batch,
                             MoreArgs = list(Subset_1),
                             SIMPLIFY = T)

Age_linear_subset2 <- mapply(Linear_age_effect, Retinol = List_retinol, Age = List_age, Batch = List_batch,
                             MoreArgs = list(Subset_2),
                             SIMPLIFY = T)

Process_age_subset1 <-  apply(Age_linear_subset1, 2, function(x) return(as.data.frame(x$coefficients)[2, 1:4]))
Process_age_subset2 <-  apply(Age_linear_subset2, 2, function(x) return(as.data.frame(x$coefficients)[2, 1:4]))

TE <- rbind(Process_age_subset1, Process_age_subset2)

Output <- data.frame()
for (i in 1:length(TE)) {
  Output <- rbind(Output, TE[[i]])
}

Output$Visit <- c("Visit1", "Visit1", "Visit2", "Visit2", "Visit3", "Visit3")
Output$Subset <- c("Subset1", "Subset2", "Subset1", "Subset2","Subset1", "Subset2")

write.table(Output, file="Effect_of_age_normative_modelling/Age_linear_effects_full_cohort.txt",
            sep="\t", row.names = F, quote = F)

## Look at visit 3 in more detail as age effect seems to weaken here

Linear_1 <- ggplot(Subset_1, aes(x=Age_visit1, y=scale(Retinol_visit1), colour = Batch_visit1)) +
  geom_point(cex = 0.75) +
  stat_smooth(method = "lm", formula = y ~ x, size = 1) +
  theme_bw() +
  xlab("Age visit 1") +
  ylab("Retinol visit 1 (SD)") +
  ggtitle("Linear") +
  scale_colour_brewer(palette =)

Linear_3 <- ggplot(Subset_1, aes(x=Age_visit3, y=scale(Retinol_visit3), colour = Batch_visit3)) +
  geom_point(cex = 0.75) +
  stat_smooth(method = "lm", formula = y ~ x, size = 1) +
  theme_bw() +
  xlab("Age visit 3") +
  ylab("Retinol visit 3 (SD)") +
  ggtitle("Linear")

Quadratic_1 <- ggplot(Subset_1, aes(x=Age_visit1, y=scale(Retinol_visit1), colour = Batch_visit1)) +
  geom_point(cex = 0.75) +
  stat_smooth(method = "lm", formula = y ~ poly(x, 2, raw=TRUE), size = 1, colour = "red") +
  theme_bw() +
  xlab("Age visit 1") +
  ylab(" ") +
  ggtitle("Quadratic")

Quadratic_3 <- ggplot(Subset_1, aes(x=Age_visit3, y=scale(Retinol_visit3), colour = Batch_visit3)) +
  geom_point(cex = 0.75) +
  stat_smooth(method = "lm", formula = y ~ poly(x, 2, raw=TRUE), size = 1, colour = "red") +
  theme_bw() +
  xlab("Age visit 3") +
  ylab(" ") +
  ggtitle("Quadratic")

Cubic_1 <- ggplot(Subset_1, aes(x=Age_visit1, y=scale(Retinol_visit1), colour = Batch_visit1)) +
  geom_point(cex = 0.75) +
  stat_smooth(method = "lm", formula = y ~ poly(x, 3, raw=TRUE), size = 1, colour = "purple") +
  theme_bw() +
  xlab("Age visit 1") +
  ylab(" ") +
  ggtitle("Cubic")

Cubic_3 <- ggplot(Subset_1, aes(x=Age_visit3, y=scale(Retinol_visit3), colour = Batch_visit3)) +
  geom_point(cex = 0.75) +
  stat_smooth(method = "lm", formula = y ~ poly(x, 3, raw=TRUE), size = 1, colour = "purple") +
  theme_bw() +
  xlab("Age visit 3") +
  ylab(" ") +
  ggtitle("Cubic")

V1 <- ggarrange(Linear_1, Quadratic_1, Cubic_1, nrow=1)

V1_out <- annotate_figure(V1, top = text_grob("First visit", size = 16, face = "bold"))


V3 <- ggarrange(Linear_3, Quadratic_3, Cubic_3, nrow=1)

V3_out <- annotate_figure(V3, top = text_grob("Third visit", size = 16, face = "bold"))

ggarrange(V1_out, V3_out, nrow = 2)

## Date of birth vs coefficient of variation

CV_S1 <- ggplot(Subset_1, aes(x=year.birth, y=Retinol_coefficient_of_variation)) +
  geom_point(cex = 0.75) +
  stat_smooth(method = "lm", formula = y ~ x, size = 1, colour = "red") +
  theme_bw() +
  xlab("Year of birth") +
  ylab("Coefficient of variation (retinol)") +
  ggtitle("Subset 1")

CV_S2 <- ggplot(Subset_2, aes(x=year.birth, y=Retinol_coefficient_of_variation)) +
  geom_point(cex = 0.75) +
  stat_smooth(method = "lm", formula = y ~ x, size = 1, colour = "red") +
  theme_bw() +
  xlab("Year of birth") +
  ylab("Coefficient of variation (retinol)") +
  ggtitle("Subset 2")

ggarrange(CV_S1, CV_S2)

## Effect of BMI

List_BMI <- as.list(c("BMI_visit1", "BMI_visit2", "BMI_visit3"))

BMI_linear_subset1 <- mapply(Linear_age_effect, Retinol = List_retinol, Age = List_BMI, Batch = List_batch,
                             MoreArgs = list(Subset_1),
                             SIMPLIFY = T)

BMI_linear_subset2 <- mapply(Linear_age_effect, Retinol = List_retinol, Age = List_BMI, Batch = List_batch,
                             MoreArgs = list(Subset_2),
                             SIMPLIFY = T)

Process_BMI_subset1 <-  apply(BMI_linear_subset1, 2, function(x) return(as.data.frame(x$coefficients)[2, 1:4]))
Process_BMI_subset2 <-  apply(BMI_linear_subset2, 2, function(x) return(as.data.frame(x$coefficients)[2, 1:4]))

## No large effect of BMI that is statistically significant