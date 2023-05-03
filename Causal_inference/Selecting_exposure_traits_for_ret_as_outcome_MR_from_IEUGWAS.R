###########################

## Selecting exposures to test with retinol as outcome

## Dr William Reay (2023)

###########################

library(dplyr)
library(ieugwasr)
library(data.table)

setwd("~/Desktop/Retinol_GWAS/MR_retinol_as_outcome/")

List_of_traits <- fread("../IEUGWAS_available_outcomes.txt")

## Steps to identify continous traits

## Non-missing unit that is not logOR

Units <- List_of_traits %>% filter(!is.na(unit)) %>% filter(unit != "logOR" | unit != "log odds")

## Remove traits with number of cases and number of controls

Units <- Units %>% filter(is.na(ncase)) %>% filter(!is.na(category == "Disease"))

## Impose an additional filter of continuous (category) to get the Nightingale UKBB metaboltes and brain imaging phenotypes

Continous <- List_of_traits %>% filter(category == "Continuous") %>% filter(is.na(ncase)) %>% filter(!is.na(category == "Disease"))

## Concatenate and remove duplicates

Concat <- rbind(Units, Continous)

De_duplicated <- Concat %>% filter(!duplicated(id))

write.table(De_duplicated, file="IEUGWAS_ids_to_use_as_exposures.txt",
            sep= "\t", quote = F, row.names = F)
