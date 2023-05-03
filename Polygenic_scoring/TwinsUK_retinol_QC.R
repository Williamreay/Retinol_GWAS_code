#############################

## TwinsUK retinol data QC 

## Dr William Reay (2023)

#########################

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(haven)
library(readxl)
library(table1)
library(reshape2)
library(cowplot)
library(gridExtra)

setwd("~/Desktop/Retinol_GWAS/TwinsUK_app/E1205_20.01.2023/")

## Read in retinol data

Retinol_raw <- read_excel("Metabolon_v4_HLI_RETINOL/RETINOL_metabolon_v4_hli_scaleRunDayMedian.xlsx")

Demo_ID_config <- Retinol_raw %>% separate(PublicID, c("FamilyID", "Last_digit"),
                                        sep = -1)
Demo_ID_config$PublicID <- paste(Demo_ID_config$FamilyID, "", Demo_ID_config$Last_digit, sep = "")

Demo_ID_config <- rename(Demo_ID_config, "Retinol"="M01806")
Demo_ID_config$date <- as.character(Demo_ID_config$date)

write.table(Demo_ID_config, "Retinol_data_text_format.txt",
            sep = "\t", row.names = F, quote = F)

## Reformat data so retinol from the three visits are in separate columns

Wide <- dcast(Demo_ID_config, PublicID  ~ visit, value.var = "Retinol")
Wide <- rename(Wide, "Retinol_visit1"="v1", "Retinol_visit2"="v2", "Retinol_visit3"="v3")

Wide1 <- dcast(Demo_ID_config, PublicID  ~ visit, value.var = "bmi")
Wide1 <- rename(Wide1, "BMI_visit1"="v1", "BMI_visit2"="v2", "BMI_visit3"="v3")

Wide2 <- dcast(Demo_ID_config, PublicID  ~ visit, value.var = "age")
Wide2 <- rename(Wide2, "Age_visit1"="v1", "Age_visit2"="v2", "Age_visit3"="v3")

Wide3 <- dcast(Demo_ID_config, PublicID  ~ visit, value.var = "date")
Wide3 <- rename(Wide3, "Date_visit1"="v1", "Date_visit2"="v2", "Date_visit3"="v3")

Wide4 <- dcast(Demo_ID_config, PublicID  ~ visit, value.var = "batch")
Wide4 <- rename(Wide4, "Batch_visit1"="v1", "Batch_visit2"="v2", "Batch_visit3"="v3")

Wide_merged <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "PublicID"),
                      list(Demo_ID_config, Wide, Wide1, Wide2, Wide3, Wide4))

## Select columns of interest

Long_format_retinol <- Wide_merged %>% filter(!duplicated(PublicID)) %>%
  select(-Retinol, -batch, -date, -year, -visit, -bmi)

## Mean retinol and BMI

Long_format_retinol$Mean_retinol_across_visits <- (Long_format_retinol$Retinol_visit1 + Long_format_retinol$Retinol_visit2 +
                                                     Long_format_retinol$Retinol_visit3)/3

Long_format_retinol$Mean_BMI_across_visits <- (Long_format_retinol$BMI_visit1 + Long_format_retinol$BMI_visit2 +
                                                     Long_format_retinol$BMI_visit3)/3

## Coefficient of variation across three visits

Long_format_retinol$Retinol_coefficient_of_variation <- apply(Long_format_retinol[ ,8:10],1,sd)/apply(Long_format_retinol[ ,8:10],1,mean)
Long_format_retinol$Retinol_coefficient_of_variation <- Long_format_retinol$Retinol_coefficient_of_variation*100

Long_format_retinol$BMI_coefficient_of_variation <- apply(Long_format_retinol[ ,11:13],1,sd)/apply(Long_format_retinol[ ,11:13],1,mean)
Long_format_retinol$BMI_coefficient_of_variation <- Long_format_retinol$BMI_coefficient_of_variation*100

write.table(Long_format_retinol, file="Long_format_full_retinol_TwinsUK.txt",
            sep = "\t", row.names = F, quote = F)

## Basic descriptive stats

table1( ~ Age_visit1 + Age_visit2 + Age_visit3 +
          factor(zygosity) + Retinol_visit1 + Retinol_visit2 +
          Retinol_visit3 + Retinol_coefficient_of_variation + BMI_visit1 +
          BMI_visit2 + BMI_visit3 +
          BMI_coefficient_of_variation | factor(sex), data = Long_format_retinol)

## Descriptive statistics (kernel density estimation plots)

Three_visits <- ggplot(Demo_ID_config, aes(x=Retinol, fill=visit)) +
  geom_density(alpha = .25) +
  theme_bw() +
  xlab("Retinol (Metabolon)") +
  ylab(" ") +
  labs(fill = "Visit") +
  ggtitle("Retinol - per visit")

Mean_retinol <- ggplot(Long_format_retinol, aes(x=Mean_retinol_across_visits)) +
  geom_density(alpha = .25, fill = "orange") +
  theme_bw() +
  ylab(" ") +
  xlab("Mean retinol") +
  ggtitle("Retinol - mean across three visits")


CV <- ggplot(Long_format_retinol, aes(x=Retinol_coefficient_of_variation)) +
  geom_density(alpha = .25, fill = "violet") +
  theme_bw() +
  ylab(" ") +
  xlab("Retinol CoV%") +
  ggtitle("Retinol - coefficient of variation (CoV)")

PLT <- ggarrange(Mean_retinol, CV, nrow = 1)

PLT2 <- ggarrange(Three_visits, PLT, nrow = 2)

## Split twins based on whether last digit is one or two
