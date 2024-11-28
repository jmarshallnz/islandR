#' Attribute cases to sources via an intermediate source with MCMC
#' @export
#' @param formula A formula of the form `GenoType ~ Covariates` for human 'sink' cases
#' @param int_genotypes A vector of intermediate genotypes (e.g. from water or similar)
#' @param sampling_dist A sampling_dist object previously fitted to source genotypes (not intermediate)
#' @param data An optional data frame to take the variables in `formula` from
#' @param iterations the number of iterations to use in the MCMC
#' @param burnin the number of iterations to eliminate due to burnin
#' @param thinning how many iterations to perform before taking a sample
#' @param priors a list of priors to use of the form `priors=list(theta=list(mean=0,prec=0.1))`
#' @return an object of type intermediate_attribution
#' @seealso print.attribution, summary.attribution, plot.attribution
attribution_with_intermediate <- function(formula, int_genotypes, sampling_dist, data, iterations=10000, burnin=1000, thinning=100, priors=NULL) {

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

  # count our intermediate genotypes
  w = as.data.frame(table(Type=int_genotypes), responseName = "Number") |> list()

  # map types in w to types from our sampling distribution on the sources
  w = lapply(w, function(x, map) { x$Type = match(x$Type, map); x }, sampling_dist$types)

  # check that there are no mismatches
  if (any(lapply(w, function(x) { sum(is.na(x$Type)) })>0)) {
    stop("Have types in `int_genotypes` aren't found in the sampling distribution of the sources")
  }
  # convert w to a matrix so it's happy in c++-land
  w = lapply(w, as.matrix)

  # now fit our MCMC
  post = NULL
  ar = c(0,0)
  for (i in seq_len(iterations(sampling_dist))) {
    cat("Performing iteration", i, "of", iterations(sampling_dist), "\n")
    iter = mcmc_intermediate(y, reduced.matrix, w, sampling_dist$sampling_distribution[[i]], iterations, burnin, thinning, priors)
    post = c(post, iter$post)
    ar = ar + iter$ar
  }
  # now assemble the final object
  x = list(posterior = post,
           acceptance_rate = ar[1] / sum(ar),
           formula = formula,
           data = data,
           cases = y,
           int_cases = w,
           model_matrix = reduced.matrix,
           genotype_data = sampling_dist)
  class(x) = c("intermediate_attribution", "attribution")
  # TODO: extend this class to attribution.intermediate or similar so summaries work for int_cases as well
  x
}
