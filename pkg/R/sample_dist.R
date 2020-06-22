#' @importFrom graphics barplot
#' @importFrom stats na.omit
NULL


prep_genotype_dist <- function(phi) {
  # for each row, factor out the largest exponent
  max_exponent <- apply(phi, 1, max)
  remainder <- exp(sweep(phi, 1, max_exponent, '-'))
  list(max_exponent, remainder)
}

#' Fit the sampling distribution of genotypes to sources using different models.
#' @export
#' @param formula A formula of the form Source ~ Genotype
#' @param non_primary one or more sources that should be considered 'output' sources. No genotype distribution is computed for these,
#'        but P(ST | source) will be computed for these STs in addition to those observed on other sources.
#' @param method the method to use to fit the genotype distribution. \code{"island"} and \code{"dirichlet"} are supported currently.
#' @param data optional data frame from which to take variables in \code{formula}.
#' @param ... further paramters to pass to the method-specific fitting algorithms.
#' @return an object derived from class sampling_dist.
#' @seealso \code{\link{st_fit_island}}, \code{\link{print.sample_dist}}, \code{\link{plot.sample_dist}}, \code{\link{summary.sample_dist}}
st_fit <- function(formula, non_primary = "Human", method=c("island", "dirichlet"), data, ...) {
  type <- match.arg(method)
  switch(type,
         island = st_fit_island(formula=formula, non_primary=non_primary, data=data, ...),
         dirichlet = st_fit_dirichlet(formula=formula, non_primary=non_primary, data=data, ...))
}

#' Print a sample_dist object
#' @export
#' @param x sample_dist object to print
#' @param ... further parameters supplied to the print function
print.sample_dist <- function(x, ...) {
  cat("Sampling distribution of genotypes\n")
  cat("----------------------------------\n")
  cat("Model:      ", x$model, "\n")
  cat("Genotypes:  ", length(x$types), "\n")
  cat("Sources:    ", length(x$sources), "\n")
  cat("Iterations: ", length(x$sampling_distribution), "\n")
}

#' Generic for retrieving the number of iterations from an object
#' @param x an object
#' @return the number of iterations
#' @export
iterations <- function(x) {
  UseMethod("iterations", x)
}

#' Retrieve the number of iterations from a sample_dist object
#' @export
#' @param x sample_dist object
#' @return the number of iterations
iterations.sample_dist <- function(x) {
  length(x$sampling_distribution)
}

#' Print a summary for a sample_dist object
#' @export
#' @param x sample_dist object to summarise
#' @param ... further parameters supplied to summary
#' @return posterior means of the ST distribution on each source
summary.sample_dist <- function(x, ...) {
  # compute posterior mean, assuming apriori all sources equally likely
  post_mean <- Reduce('+', x$sampling_distribution)/length(x$sampling_distribution)
  post_mean <- prep_genotype_dist(post_mean)[[2]]
  post_mean <- post_mean / rowSums(post_mean)
  print(post_mean)
  invisible(post_mean)
}

#' Produce a summary plot (stacked barplot) of posterior means
#' @export
#' @param x sample_dist object to plot
#' @param ... further parameters supplied to plot
plot.sample_dist <- function(x, ...) {
  # compute posterior mean and plot that maybe?
  summ <- summary(x)
  barplot(t(summ), beside=FALSE)
}

#' Convert sampling distribution to a data frame
#' @export
#' @param x sample_dist object to convert
#' @param ... further parameters supplied to as.data.frame
#' @param types the sequence types to extract. Defaults to NULL (all types)
#' @return a data.frame with columns type, iteration, source and log_p
as.data.frame.sample_dist <- function(x, ..., types = NULL) {
  which_types <- seq_along(x$types)
  if (!is.null(types)) {
    which_types <- na.omit(match(types, x$types))
  }
  p_types <- unlist(lapply(x$sampling_distribution, function(x) { x[which_types,] }))
  data.frame(type=rep(x$types[which_types], by=iterations(x)*length(x$sources)),
             iteration = rep(seq_len(iterations(x)), each=length(x$sources)*length(which_types)),
             source = rep(rep(x$sources, each=length(which_types)), by=iterations(x)*length(which_types)),
             log_p = p_types,
             stringsAsFactors = FALSE)
}

