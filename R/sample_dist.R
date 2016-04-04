
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
  cat("Iterations: ", dim(x$sampling_distribution)[3], "\n")
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
  dim(x$sampling_distribution)[3]
}

#' Print a summary for a sample_dist object
#' @export
#' @param object sample_dist object to summarise
#' @param ... further parameters supplied to summary
#' @return posterior means of the ST distribution on each source
summary.sample_dist <- function(object, ...) {
  # compute posterior mean, assuming apriori all sources equally likely
  post_mean = apply(object$sampling_distribution, 1:2, mean)
  post_mean = post_mean / rowSums(post_mean)
  print(round(post_mean,2))
  invisible(post_mean)
}

#' Produce a summary plot (stacked barplot) of posterior means
#' @export
#' @param x sample_dist object to plot
#' @param ... further parameters supplied to plot
plot.sample_dist <- function(x, ...) {
  # compute posterior mean and plot that maybe?
  post_mean = apply(x$sampling_distribution, 1:2, mean)
  post_mean = post_mean / rowSums(post_mean)
  barplot(t(post_mean), beside=FALSE)
}
