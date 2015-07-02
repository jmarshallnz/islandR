#' Genotype attribution dataset
#'
#' A dataset containing posterior samples of attribution to sources for a list
#' of genotypes.
#'
#' @format A data frame containing attribution information
#' \describe{
#'   \item{ST}{The sequence type identifier}
#'   \item{ASP}{The ASP housekeeping gene}
#'   \item{GLN}{The GLN housekeeping gene}
#'   \item{GLT}{The GLT housekeeping gene}
#'   \item{GLY}{The GLY housekeeping gene}
#'   \item{PGM}{The PGM housekeeping gene}
#'   \item{TKT}{The TKT housekeeping gene}
#'   \item{UNC}{The UNC housekeeping gene}
#'   \item{Poultry}{The probability the sample is from poultry}
#'   \item{Ruminants}{The probability the sample is from ruminants}
#'   \item{Water}{The probability the sample is from environmental water}
#'   \item{Other}{The probability the sample is from other sources}
#'   \item{Iteration}{The posterior sample iteration}
#' }
"genotype_attribution"
