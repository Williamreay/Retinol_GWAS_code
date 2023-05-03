###############################

## Cauchy aggregation at gene-level

## Retinol rare variants - Stouffer's meta-analysis (INTERVAL+METSIM)

## William Reay (2022)

###############################

library(dplyr)
library(data.table)

setwd("~/Desktop/Retinol_GWAS/Cauchy_rare_var_gene/")

## Load ACAT function

source("~/cloudstor/23andMe_pneumonia/2022_FinnGen_r6_meta_analysis/Scripts/ACAT_function.R")

## Load all rare variants and the P value from the meta-analysis

Raw_meta <- fread("../Meta_analyses/Rare_var/RARE_VAR_Stouffers_METSIM_INTERVAL1.tbl.gz")

## Read in annotations for all variants annotated

All_annot_rare <- fread("Annotated_retinol_rare_var_Stouffer.variant_function", header = F)

All_annot_rare <- All_annot_rare %>% filter(V1 != "intergenic")

## Retain gene name and SNP ID

Filt_rare <- All_annot_rare %>% select(V8, V2)

Filt_rare <- rename(Filt_rare, "MarkerName"="V8", "Gene"="V2")

## Merge with P values

Annot_meta <- merge(Raw_meta, Filt_rare, by="MarkerName")

## Load NCBI 38 genes

Gene_anno <- read.table("NCBI38.gene.loc", header = F)

Gene_anno <- rename(Gene_anno, "Gene"="V6")

## Merge with NCBI protein coding genes

Coord_annot_meta <- merge(Annot_meta, Gene_anno, by ="Gene")

## Retain genes with at least 2 variants

Genes <- Coord_annot_meta %>% select(Gene)

Retained_genes <- Genes %>% group_by_all() %>% summarise(COUNT = n()) %>% arrange(desc(COUNT)) %>% filter(COUNT > 1)

Coord_annot_meta <- merge(Coord_annot_meta, Retained_genes, by = "Gene")

## Process data for gene-wise Cauchy aggregation, followed by FWER and FDR calculation

genes <- sort(unique(Coord_annot_meta$Gene))
gene_ps <- list()

for (i in genes) { gene_ps[[i]] <- c(Coord_annot_meta[Coord_annot_meta$Gene==i,][,7])}

genes_acat <- list()
for (i in genes) { genes_acat[[i]] <- ACAT(unlist(gene_ps[[i]], use.names=FALSE))}

af_acat <- data.frame(cbind(as.vector(genes), as.vector(unlist(genes_acat))))

af_acat <- rename(af_acat, c("Gene"="X1", "P"="X2"))

af_acat$FDR <- p.adjust(af_acat$P, method = "fdr")
af_acat$FWER <- p.adjust(af_acat$P, method="bonferroni")

write.table(af_acat, file="ALL_genic_var_Cauchy_2022_Stouffers_retinol.txt",
            sep = "\t", row.names = F, quote = F)

## Repeat the above with just exonic annotated variants

Exonic_annot_rare <- All_annot_rare %>% filter(V1 == "exonic" | V1 == "exonic;splicing")

## Retain gene name and SNP ID

Exonic_Filt_rare <- Exonic_annot_rare %>% select(V8, V2)

Exonic_Filt_rare <- rename(Exonic_Filt_rare, "MarkerName"="V8", "Gene"="V2")

## Merge with P values

Exonic_Annot_meta <- merge(Raw_meta, Exonic_Filt_rare, by="MarkerName")

## Merge with NCBI protein coding genes

Exonic_coord_annot_meta <- merge(Exonic_Annot_meta, Gene_anno, by ="Gene")

## Retain genes with at least 2 variants

Genes_ex <- Exonic_coord_annot_meta %>% select(Gene)

Ex_Retained_genes <- Genes_ex %>% group_by_all() %>% summarise(COUNT = n()) %>% arrange(desc(COUNT)) %>% filter(COUNT > 1)

Exonic_coord_annot_meta <- merge(Exonic_coord_annot_meta, Ex_Retained_genes, by = "Gene")

## Process data for gene-wise Cauchy aggregation, followed by FWER and FDR calculation

ex_genes <- sort(unique(Exonic_coord_annot_meta$Gene))
ex_gene_ps <- list()

for (i in ex_genes) { ex_gene_ps[[i]] <- c(Exonic_coord_annot_meta[Exonic_coord_annot_meta$Gene==i,][,7])}

ex_genes_acat <- list()
for (i in ex_genes) { ex_genes_acat[[i]] <- ACAT(unlist(ex_gene_ps[[i]], use.names=FALSE))}

ex_af_acat <- data.frame(cbind(as.vector(ex_genes), as.vector(unlist(ex_genes_acat))))

ex_af_acat <- rename(ex_af_acat, c("Gene"="X1", "P"="X2"))

ex_af_acat$FDR <- p.adjust(ex_af_acat$P, method = "fdr")
ex_af_acat$FWER <- p.adjust(ex_af_acat$P, method="bonferroni")

write.table(ex_af_acat, file="Exonic_var_Cauchy_2022_Stouffers_retinol.txt",
            sep = "\t", row.names = F, quote = F)

## CHD1 and FREM2 assoc driven by single signal
