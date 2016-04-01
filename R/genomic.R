#' Get the number of posterior samples available for attribution
get_num_samples <- function() {
  return(max(genotype_attribution$Iteration))
}

#' Get the genotypes available for attribution
get_genotypes <- function() {
  return(unique(genotype_attribution[,1:8])) # TODO: Make indicies a function of data
}

#' Retrieve the i-th posterior sample for attribution for a given genotype
#' @param genotype the genotype to retrieve.
#' @param sample the sample to retrieve. Defaults to NULL, where all samples will be retrieved.
#' @return A data frame of sequence types, their allelic profile, and clonal complex.
#' @seealso pubmlst
get_source_probability_sample <- function(genotype, sample = NULL) {
  wch <- genotype_attribution$ST == genotype;
  if (!is.null(sample))
    wch <- wch & genotype_attribution$Iteration == sample
  return(genotype_attribution[wch, 9:12]) # TODO: Make indicies a function of data
}
