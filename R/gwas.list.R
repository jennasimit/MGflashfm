#' List of GWAS results for two traits in two groups
#'
#' Group 1 GWAS results are in the first component, which has a list of two data frames, one for each trait. Each of these data frames has
#' 201 rows; each row is a variant. Likewise for the Group2 structure, but these data frames have 329 rows.
#' These data are for a region in chromosome 19:45386029-45439498 (GRCh37/hg19), containing APOE.
#' We used HapGen2 to generate two groups of individuals based on cohorts from the 1000 Genomes Project:
#' N1=90,000 based on CEU+TSI and N2=10,000 based on YRI+LWK.
#' Within each group, we retained variants with MAF > 0.001, regardless of whether the variant satisfies our MAF threshold in the other group.
#' This resulted in 201 variants in Group 1 and 329 variants in Group 2; 132 variants appear in both groups and there are 398 unique variants present $
#'
"gwas.list"



