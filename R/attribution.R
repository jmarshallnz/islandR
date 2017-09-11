#' @importFrom reshape2 melt
#' @importFrom stats terms model.frame model.matrix median delete.response quantile setNames
NULL

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
    cat("Performing iteration", i, "of", iterations(sampling_dist), "\n")
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

#' Attribute cases to sources via MCMC with an AR(1) model on time
#' @export
#' @param formula A formula of the form `GenoType ~ Covariates` for cases
#' @param time The name of the time column
#' @param sampling_dist A sampling_dist object previously fitted to source genotypes
#' @param data An optional data frame to take the variables in `formula` from
#' @param iterations the number of iterations to use in the MCMC
#' @param burnin the number of iterations to eliminate due to burnin
#' @param thinning how many iterations to perform before taking a sample
#' @param priors A list specifying the priors to use.
#' @return an object of type attribution
#' @seealso print.attribution, summary.attribution, plot.attribution
attribution_ar1 <- function(formula, time, sampling_dist, data, iterations=10000, burnin=1000, thinning=100, priors=NULL) {

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

  # pull out the times
  time_range = range(data[,time])
  times = 1:time_range[2]

  # hmm, I think we'll need to be lots smarter here, but for now this is a hack
  # you'd basically run expand.grid for any other covariates here...
  # we kind of need times expanded (so no missing ones). Perhaps that can be a requirement?
  # and then need to reduce things in the usual fashion, then sort by time and then
  # remove time from the matrix.
  x0 = expand.grid(Times=times, UR_bool=factor(c("Rural", "Urban")))
  x0$Intervention=factor(ifelse(x0$Times <= 36, "Before", "After"))
  mod.matrix = model.matrix(~Intervention*UR_bool, data=x0)
#  mod.matrix = matrix(1, nrow=length(times), ncol=1)
#  colnames(mod.matrix) = "(Intercept)"

  # run through the model matrix and find the unique entries
#   reduced.matrix = list()
#   reduced.response = list()
#   for (i in 1:nrow(mod.matrix)) {
#     # check if this is similar to one of the rows we already have
#     row = mod.matrix[i,]
#     found = FALSE
#     for (j in seq_along(reduced.matrix)) {
#       if (isTRUE(all.equal(row, reduced.matrix[[j]]))) {
#         # yes! Accumulate up the response
#         reduced.response[[j]] = c(reduced.response[[j]], mod.frame[i,response])
#         found = TRUE;
#         break;
#       }
#     }
#     if (!found) {
#       # add to the end
#       reduced.matrix[[length(reduced.matrix)+1]] = row
#       reduced.response[[length(reduced.response)+1]] = mod.frame[i,response]
#     }
#   }
  reduced.matrix = list()
  reduced.response = list()
  for (i in 1:nrow(mod.matrix)) {
    row = mod.matrix[i,]
    # find all rows in the data that match this time
    wch = data[,time] == x0$Times[i] & data$UR_bool == x0$UR_bool[i] # fuck this is bad...
    reduced.matrix[[length(reduced.matrix)+1]] = row
    reduced.response[[length(reduced.response)+1]] = data[wch,response]
  }
  # now accumulate up the response variable and matrix
  # TODO: this seems to screw up if we don't actually have a matrix??
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
    cat("Performing iteration", i, "of", iterations(sampling_dist), "\n")
    iter = mcmc(y,
                t=times,
                X=reduced.matrix,
                phi=sampling_dist$sampling_distribution[,,i],
                iterations=iterations,
                burnin=burnin,
                thinning=thinning,
                priors=priors)
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
#' @param x an object of class `attribution`
#' @param ... further parameters passed to base print functions
print.attribution <- function(x, ...) {
  cat("Attribution of cases\n")
  cat("====================\n")
  cat("Formula:", deparse(x$formula), "\n")
  cat("Cases:", sum(simplify2array(lapply(x$cases, function(x) { sum(x[,2]) }))), "\n")
  cat("Iterations:", length(x$posterior), "\n")
  cat("Sources:", paste(x$genotype_data$sources, collapse=", "), "\n")
}

.attribution.probabilities <- function(x, model_matrix) {

  X = model_matrix
  post_theta = get_var(x$posterior, "theta")

  p_pred = list()
  for (j in seq_along(post_theta)) {
    p_pred[[j]] = t(exp(X %*% t(post_theta[[j]])))
    # TODO: This only works for factor variables...
    #       What do we do for others?
    colnames(p_pred[[j]]) = apply(X, 1, function(x, y) { paste(paste(y,x,sep='='), collapse=":")}, colnames(X))
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
#' @param object an object of class `attribution`, usually, a result of a call to `attribution`
#' @param ... further parameters passed to summary functions
summary.attribution <- function(object, ...) {

  theta = simplify2array(lapply(object$posterior, function(x) { x$theta }))

  # NOTE: These likely have too large a range to be able to say
  #       with any confidence that we have differences.
  #       This is because they're capturing the variation of the
  #       baseline source. i.e. we can add a constant to each row
  #       of theta, and then add that constant to 0 and we'll get the same
  #       probabilities... So can we maybe try and find an a for each
  #       iteration that makes sum(1 + e^theta) constant?
  # e^a e^0, e^a e^t1, e^a e^t2, e^a e^t3
 # etheta = apply(theta, c(1,3), function(x) { log(sum(exp(x)+1)) })
#  etheta_med = apply(etheta, 1, mean)
#  etheta_scale = sweep(etheta, 1, etheta_med, "-")
  # now get rid of this out of theta
#  thetahat = sweep(theta, c(1,3), etheta_scale, "-")
#  med = apply(thetahat, 1:2, quantile, c(0.50, 0.05, 0.95))
  # doesn't seem to really help much though!

  med = apply(theta, 1:2, quantile, c(0.50, 0.05, 0.95))

  split.along.dim <- function(a, n)
    setNames(lapply(split(a, arrayInd(seq_along(a), dim(a))[, n]),
                    array, dim = dim(a)[-n], dimnames(a)[-n]),
             dimnames(a)[[n]])

  summ = list(summ = lapply(split.along.dim(med, 3), t),
              baseline = object$genotype_data$sources[length(object$genotype_data$sources)],
              n = sum(simplify2array(lapply(object$cases, function(x) { sum(x[,2]) }))))
  class(summ) = "summary.attribution"
  summ
}

#' Print a summary of an attribution object
#' @export
#' @param x an object of class `summary.attribution`, usually, a result of a call to `summary.attribution`
#' @param ... further parameters passed to base print functions
print.summary.attribution = function(x, ...) {
  cat("Attribution model fit: Posterior quantiles\n")
  cat("==========================================\n")
  cat("n =", x$n, "\n\n")
  cat("Baseline:", x$baseline, "\n\n")
  snames = names(x$summ)
  lapply(seq_along(snames), function(i, y) { cat(snames[i], ":\n", sep=""); print(y[[i]]); cat("\n") }, x$summ)
}

#' Plot an attribution object
#' @export
#' @param x An object of class `attribution`
#' @param ... further parameters passed to lower-level plot commands.
plot.attribution <- function(x, ...) {
  # hmm, what is a useful plot? I guess posterior predictions for each
  # combination of covariates?
  print("Not currently implemented")
}

#' Predict attribution on a new data set
#' @export
#' @param object an object of class `attribution`, usually a result of a call to `attribution`.
#' @param newdata a `data.frame` to predict attribution values. Set to NULL to use the reduced model matrix from the fitted model.
#' @param FUN a function to operate on the posterior attribution, defaults to `median`, use `identity` to retrieve all posterior samples.
#' @param ... further parameters to pass to FUN
#' @return a data frame containing the posterior attribution (p), the source (Source), the covariate
#'         combination (X), and if multiple results are returned from FUN, an additional column named
#'         after the function call.
predict.attribution <- function(object, newdata=NULL, FUN=median, ...) {

  # generate a model matrix for the new data
  if (missing(newdata) || is.null(newdata)) {
    model_matrix = object$model_matrix
  } else {
    Terms = delete.response(terms(object$formula))
    model_matrix = model.matrix(Terms, newdata)
  }

  # get prediction
  pred = .attribution.probabilities(object, model_matrix)

  # operate on the prediction(s) using FUN
  to_matrix <- function(x) {
    if (is.matrix(x)) {
      x
    } else {
      t(as.matrix(x))
    }
  }
  m_pred = lapply(pred, function(x) { to_matrix(apply(x, 2, FUN, ...)) })

  func_name = deparse(substitute(FUN))

  m_pred_df = reshape2::melt(m_pred, varnames=c(func_name, "X"), value.name = "p")
  names(m_pred_df)[names(m_pred_df) == "L1"] = "Source"
  # drop func_name column if it is identical
  if (length(unique(m_pred_df[,func_name])) == 1) {
    m_pred_df[,func_name] = NULL
  }

  m_pred_df
}

#' Retrieve the number of MCMC samples from an attribution object
#' @export
#' @param x An object of class `attribution`
#' @return The number of iterations (integer)
iterations.attribution <- function(x) {
  length(x$posterior)
}
