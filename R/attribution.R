#' Attribute cases to sources via MCMC
#' @export
#' @param formula A formula of the form `GenoType ~ Covariates` for cases
#' @param sampling_dist A sampling_dist object previously fitted to source genotypes
#' @param data An optional data frame to take the variables in `formula` from
#' @param iterations the number of iterations to use in the MCMC
#' @param burnin the number of iterations to eliminate due to burnin
#' @param thinning how many iterations to perform before taking a sample
#' @return an object of type attribution
#' @seealso print.attribution, summary.attribution, plot.attribution
attribution <- function(formula, sampling_dist, data, iterations=10000, burnin=1000, thinning=100) {

  # check inputs
  if("sampling_dist" %in% class(sampling_dist)) {
    stop("sampling_dist must be of class 'sampling_dist'")
  }

  # pull out the formula terms etc in order to compute the model matrix
  mod.terms = terms(formula, data=data)
  mod.frame = model.frame(formula, data=data)

  if (attr(mod.terms, "response") == 0) {
    stop("formula needs a left hand side")
  }

  mod.vars = attr(mod.terms, "variables")
  response = all.vars(mod.terms[[2]]) # response variable

  mod.matrix = model.matrix(formula, data=data)

  # run through the model matrix and find the unique entries
  reduced.matrix = list()
  reduced.response = list()
  for (i in 1:nrow(mod.matrix)) {
    # check if this is similar to one of the rows we already have
    row = mod.matrix[i,]
    found = FALSE
    for (j in seq_along(reduced.matrix)) {
      if (isTRUE(all.equal(row, reduced.matrix[[j]]))) {
        # yes! Accumulate up the response
        reduced.response[[j]] = c(reduced.response[[j]], mod.frame[i,response])
        found = TRUE;
        break;
      }
    }
    if (!found) {
      # add to the end
      reduced.matrix[[length(reduced.matrix)+1]] = row
      reduced.response[[length(reduced.response)+1]] = mod.frame[i,response]
    }
  }
  # now accumulate up the response variable and matrix
  reduced.matrix = t(as.matrix(simplify2array(reduced.matrix)))
  colnames(reduced.matrix) = colnames(mod.matrix)
  y = lapply(reduced.response, function(x) { as.data.frame(table(Type=x), responseName = "Number") })

  # now map types in y to types from our sampling distribution on the sources
  y = lapply(y, function(x, map) { x$Type = match(x$Type, map); x }, sampling_dist$types)
  # check that there are no mismatches
  if (any(lapply(y, function(x) { sum(is.na(x$Type)) })>0)) {
    stop("Have types that aren't found in the sampling distribution of the sources")
  }
  # convert y to a matrix so it's happy in c++-land
  y = lapply(y, as.matrix)

  # now fit our MCMC
  post = NULL
  ar = c(0,0)
  for (i in seq_len(iterations(sampling_dist))) {
    cat("Performing iteration", i, "of", iters, "\n")
    iter = mcmc_no_ar1(y, reduced.matrix, sampling_dist$sampling_distribution[,,i], iterations, burnin, thinning)
    post = c(post, iter$post)
    ar = ar + iter$ar
  }
  # now assemble the final object
  x = list(posterior = post,
           acceptance_rate = ar[1] / sum(ar),
           formula = formula,
           data = data,
           cases = y,
           model_matrix = reduced.matrix,
           genotype_data = sampling_dist)
  class(x) = "attribution"
  x
}

#' Print an overview of an attribution object
#' @export
#' @param x An object of class `attribution`
print.attribution <- function(x) {
  cat("Attribution of cases\n")
  cat("====================\n")
  cat("Formula:", deparse(x$formula), "\n")
  cat("Cases:", sum(simplify2array(lapply(y, function(x) { sum(x[,2]) }))), "\n")
  cat("Iterations:", length(x$posterior), "\n")
  cat("Sources:", paste(x$genotype_data$sources, collapse=", "), "\n")
}

.attribution.probabilities <- function(x) {

  X = x$model_matrix
  post_theta = get_var(x$posterior, "theta")

  p_pred = list()
  for (j in seq_along(post_theta)) {
    p_pred[[j]] = t(exp(X %*% t(post_theta[[j]])))
    colnames(p_pred[[j]]) = apply(X, 1, function(x, y) { paste(y[x == 1], collapse=":")}, colnames(X))
  }
  p_pred[[length(p_pred)+1]] = matrix(1, nrow(p_pred[[1]]), ncol(p_pred[[1]]))

  p_sum = matrix(0, nrow(p_pred[[1]]), ncol(p_pred[[1]]))
  for (j in seq_along(p_pred)) {
    p_sum = p_sum + p_pred[[j]]
  }
  for (j in seq_along(p_pred)) {
    p_pred[[j]] = p_pred[[j]] / p_sum
  }
  names(p_pred) = x$genotype_data$sources

  p_pred
}

#' Retrieve a summary of an attribution object
#' @export
#' @param x An object of class `attribution`
summary.attribution <- function(x) {

  theta = get_var(x$posterior, "theta")

  j = 1
  theta = simplify2array(lapply(x$posterior, function(x) { x$theta[,j] }))

  med = apply(theta, 1:2, median)
  lci = apply(theta, 1:2, quantile, 0.05)
  uci = apply(theta, 1:2, quantile, 0.95)

  # length of post_theta is the number of sources-1
  med = lapply(theta, function(x) { apply(x, 2, median) })
  lci = lapply(theta, function(x) { apply(x, 2, quantile, 0.05) })
  uci = lapply(theta, function(x) { apply(x, 2, quantile, 0.95) })

  s = sprintf("%02.2f (%02.2f, %02.2f)", med, lci, uci)
  names(s) = names(med)
  sh = matrix(s, ncol=3)
  colnames(sh) = names(med)
  print(s)

}

#' Plot an attribution object
#' @export
#' @param x An object of class `attribution`
plot.attribution <- function(x) {
  # hmm, what is a useful plot?
}

#' Retrieve the number of MCMC samples from an attribution object
#' @export
#' @param x An object of class `attribution`
#' @return The number of iterations (integer)
iterations.attribution <- function(x) {
  length(x$posterior)
}
