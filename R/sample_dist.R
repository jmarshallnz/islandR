#' @importFrom graphics barplot
NULL

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
#' @param object sample_dist object to summarise
#' @param ... further parameters supplied to summary
#' @return posterior means of the ST distribution on each source
summary.sample_dist <- function(object, ...) {
  # compute posterior mean, assuming apriori all sources equally likely
  # TODO FIXME!
  print("Not currently implemented")
}

#' Produce a summary plot (stacked barplot) of posterior means
#' @export
#' @param x sample_dist object to plot
#' @param ... further parameters supplied to plot
plot.sample_dist <- function(x, ...) {
  # compute posterior mean and plot that maybe?
  # TODO: FIXME!
  print("Not currently implemented")
}
