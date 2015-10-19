#' Island model, metropolis-hastings fit on logit scale
#'
#' Eventually this will include AR1() model
library(MASS)

#' Log-likelihood function
log_lik_R = function(humans, phi, p) {
  loglik = 0
  # run through human isolates

  for (h in seq_len(nrow(humans))) {
    # calculate the likelihood for this human isolate
    lik_h = 0;
    for (j in seq_along(p)) {
      lik_h = lik_h + p[j] * phi[humans[h,1], j]
    }
    loglik = loglik + humans[h,2] * log(lik_h)
  }
  return(loglik)
}

# inverse logit function
inverse_logit = function(logit_p) {
  # Basically, P[i] = exp(logit_p[i]) / (1 + sum(exp(logit_p[i]))
  #            P[n] = 1 / (1 + sum(exp(logit_p[i]))
  pn = exp(c(logit_p, 0)) # as logit_p[n] = 0
  pn / sum(pn)
}

update_priors = function(curr) {
  # hyper-priors
  theta_0    = rep(0, length(curr$theta))
  theta_prec = diag(0.1, length(curr$theta))
  tau_shape  = 0.1
  tau_rate   = 0.1

  # 1. Sample theta
  #
  # prior theta ~ Normal(theta_0, theta_prec)
  #
  # and p ~iid Normal(X theta, tau)
  #
  # the posterior is
  #
  # theta ~ Normal(theta_hat, prec_hat)
  #
  # where
  #
  # prec_hat = theta_prec + tau X'X
  # theta_hat = prec_hat^{-1} (theta_prec*theta_0 + tau X'y)

  var_hat   = solve(theta_prec + curr$tau * t(curr$X) %*% curr$X)
  theta_hat = var_hat %*% (theta_prec %*% theta_0 + curr$tau * t(curr$X) %*% curr$p)

  curr$theta = mvrnorm(1, theta_hat, var_hat)

  # 2. Sample tau
  #
  # prior tau ~ Gamma(tau_shape, tau_rate)
  #
  # posterior will be (tau_shape + n/2, tau_rate + (p - X theta)'(p - X theta)/2)
  rss = sum((curr$p - curr$X %*% curr$theta)^2)

  n     = length(curr$p)
  shape = tau_shape + n / 2
  rate  = tau_rate  + rss / 2
  curr$tau   = rgamma(1, shape, rate=rate)

  curr
}

update_p = function(curr, humans, phi) {

  # proposal distribution
  p_proposal_sigma = 1

  # first up, compute the fitted values
  mu = curr$X %*% curr$theta

  # TODO: Could do this at random instead
  # for each covariate patterns
  for (i in 1:length(curr$p)) {

    # we can optimise this quite a bit, as we need only update
    # the log likelihood corresponding to this covariate pattern
    # as every other param is staying the same

    # find the covariate pattern corresponding to this observation
    t = curr$t[i]
    h = humans[[t]]

    # update the corresponding p
    p     = curr$p
    p[i]  = rnorm(1, curr$p[i], p_proposal_sigma);

    # compute prior-hastings ratio
    # Prior-Hastings ratio = Proposal(f,f')/Proposal(f',f) * Prior(f')/Prior(f)
    # Proposal is normal distribution so is symmetric, so this drops down to the prior.

    # Prior in our hierarchical model is determined by mu + tau
    # exp(((p-mu)^2-(p'-mu)^2)/2*tau)
    log_hastings_ratio = ((curr$p[i] - mu[i])^2 - (p[i] - mu[i])^2)*0.5*curr$tau;

    # compute likelihood ratio
    log_likelihood = log_lik(h, phi, p[curr$t == t])
    log_likelihood_ratio = log_likelihood - curr$log_likelihood[t]

    # accept/reject
    log_alpha = log_likelihood_ratio + log_hastings_ratio
    if (log_alpha > 0 || runif(1) < exp(log_alpha)) {
      # accept new values, copy to old
      curr$p[i] = p[i]
      curr$log_likelihood[t] = log_likelihood
      curr$accept = curr$accept + 1
    } else {
      curr$reject = curr$reject + 1
    }
  }
  return(curr)
}

mcmc = function(humans, phi) {

  #' MCMC parameters
  iterations = 10000
  burnin     = 1000
  thinning   = 100

  # priors
  logit_p_sigma = 1

  # accept/reject
  accept_reject = numeric(2)

  n_sources = ncol(phi)
  n_times   = length(humans)

  # posterior
  posterior = list()

  #' sample p from the prior
  p = rnorm(n_times * (n_sources-1), 0, logit_p_sigma)

  #' time column (in case time is important)
  t = rep(1:n_times, 3)

  #' design matrix
  X = matrix(0, n_times * (n_sources-1), n_sources-1)
  X[            1:n_times, 1] = 1
  X[  n_times + 1:n_times, 2] = 1
  X[2*n_times + 1:n_times, 3] = 1

  #' parameter vector
  theta   = rep(0, n_sources - 1)

  #' precision
  tau     = 1

  #' compute log-likelihood for each covariate pattern
  log_likelihood = numeric(n_times)
  for (x in 1:n_times) {
    log_likelihood[x] = log_lik(humans[[x]], phi, p[t == x])
  }

  # storage for the current iteration
  curr = list(p = p, X = X, t = t, theta = theta, tau = tau, log_likelihood = log_likelihood, accept=0, reject=0)

  # main MCMC loop
  post_i = 0;
  for (i in seq_len(iterations+burnin)) {

    # update priors
    curr = update_priors(curr)

    # update the p's
    curr = update_p(curr, humans, phi)

    # sample
    if (i %% 1000 == 0) {
      cat("Up to iteration", i, "of", burnin + iterations, "\n")
    }
    if (i > burnin && i %% thinning == 0) {
      post_i = post_i + 1;
      posterior[[post_i]] = list(p     = curr$p,
                                 theta = curr$theta,
                                 tau   = curr$tau)
    }
  }

  list(post = posterior, ar = curr$accept / (curr$accept + curr$reject))
}
