#' Covariance matrices for two simulated traits in two groups
#' 
#' Simulated data were generated using the mvrnorm function in R and causal variants are specified in the cvs object. 
#' The two traits have correlation 0.4 and each has two causal variants, of which one, "chr19_45387034_A_C" is shared between traits; 
#' trait 1 has second causal variant "chr19_45425175_T_G" and trait 2 has second causal variant "chr19_45427353_C_G". 
#' The shared variant has MAF < 0.05 in Group 1 and MAF > 0.05 in Group 2, and the second causal variants of each trait have MAF > 0.05 in both groups; 
#' r2 < 0.5 for any pair of causal variants.
#' The Group 1 traits are from a sample of 90,000 and Group 2 traits are from a sample of size 10,000; there are no missing measurements.
#' 
#' The matrix has row and column names coinciding with the trait names.
#'
"covY.list"

