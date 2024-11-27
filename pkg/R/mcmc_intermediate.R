
# update intermediate logits eta
update_eta = function(curr, humans, intermediate, phi, priors) {

  # proposal distribution
  eta_proposal_sigma = 2

  # prior distribution
  eta_0    = matrix(priors$eta$mean, nrow(curr$eta), ncol(curr$eta))
  eta_prec = diag(priors$eta$prec, nrow(curr$eta)) # same prior across all sources

  # for each source, sampled at random
  i = 1 # just one row
  cols = sample(ncol(curr$eta))
  for (j in cols) {
    # update eta[i,j] with RWMH
    eta = curr$eta
    eta[i,j] = rnorm(1, curr$eta[i,j], eta_proposal_sigma)

    # proposal ratio is symmetric, so need only the prior ratio
    t0_prop = eta[,j]      - eta_0[,j]
    t0_curr = curr$eta[,j] - eta_0[,j]
    # exp(-0.5*(t(t0_prop) %*% theta_prec %*% t0_prop)) / exp(-0.5*(t(t0_curr) %*% theta_prec %*% t0_curr))

    log_hastings_ratio = -0.5*(t(t0_prop) %*% eta_prec %*% t0_prop -
                                 t(t0_curr) %*% eta_prec %*% t0_curr)

    # compute likelihood ratio.

    # need both the human ratio and the intermediate ratio
    log_likelihood_int = log_lik(intermediate[[1]], phi, eta)

    # compute the probability of each ST on the water source
    h_phi <- append_log_sum_exp(phi, eta)

    # compute log-likelihood for each covariate pattern
    log_likelihood = log_lik_full(humans, h_phi, curr$p)

    # compute log-likelihood for intermediate source
    log_likelihood_ratio = sum(log_likelihood - curr$log_likelihood) +
      log_likelihood_int - curr$log_likelihood_int

    # accept/reject
    log_alpha = log_likelihood_ratio + log_hastings_ratio
    if (log_alpha > 0 || runif(1) < exp(log_alpha)) {
      # accept new values, copy to old
      curr$eta[i,j] = eta[i,j]
      curr$log_likelihood = log_likelihood
      curr$log_likelihood_int = log_likelihood_int
      curr$accept = curr$accept + 1
    } else {
      curr$reject = curr$reject + 1
    }
  }
  return(curr)
}

append_log_sum_exp <- function(phi, eta) {
  # convert p to exp(p)
  max_exponent <- phi[[1]]
  remainder <- phi[[2]]

  # convert eta from logit
  p <- c(exp(eta), 1)
  p <- p / sum(p)

  # product up
  remainder <- cbind(remainder, remainder %*% p)
  colnames(remainder)[ncol(remainder)] <- "intermediate"

  list(max_exponent, remainder)
}

mcmc_intermediate = function(humans, X, intermediate, genotype_dist, iterations = 10000, burnin = 1000, thinning = 100, priors = NULL) {

  if (is.null(priors)) {
    priors <- list()
    priors$theta <- list(mean = 0, prec = 0.1)
    priors$eta   <- list(mean = 0, prec = 0.1)
  }

  phi <- prep_genotype_dist(genotype_dist)

  # accept/reject
  accept_reject = numeric(2)

  # posterior
  posterior = list()

  # parameter vector
  n_sources = ncol(genotype_dist)
  theta   = matrix(0, ncol(X), n_sources)
  rownames(theta) <- colnames(X)
  colnames(theta) <- colnames(genotype_dist)

  eta     = matrix(0, 1, n_sources-1)
  rownames(eta) <- "(Intercept)"
  colnames(eta) <- colnames(genotype_dist)[-n_sources]

  # initialise p, the per-covariate pattern logit predictors
  p = X %*% theta

  # compute log-likelihood for intermediate source
  log_likelihood_int = log_lik(intermediate[[1]], phi, eta)

  # compute log-likelihood for each covariate pattern
  h_phi <- append_log_sum_exp(phi, eta)
  log_likelihood = log_lik_full(humans, h_phi, p)

  # storage for the current iteration
  curr = list(p = p, X = X, t = NULL, theta = theta, eta = eta, tau = NULL, rho = NULL, log_likelihood = log_likelihood, log_likelihood_int = log_likelihood_int, accept=0, reject=0)

  # main MCMC loop
  post_i = 0;
  for (i in seq_len(iterations+burnin)) {

    # update the eta (needs full likelihood on intermediate and humans)
    curr = update_eta(curr, humans, intermediate, phi, priors)

    # update theta (needs only partial likelihood)
    # recompute h_phi if eta has changed
    if (any(curr$eta != eta)) {
      eta <- curr$eta
      h_phi <- append_log_sum_exp(phi, curr$eta)
    }
    curr = update_theta(curr, humans, h_phi, priors)

    # sample
    if (i %% 1000 == 0) {
      cat("Up to iteration", i, "of", burnin + iterations, "\n")
    }
    if (i > burnin && i %% thinning == 0) {
      post_i = post_i + 1;
      posterior[[post_i]] = list(p      = curr$p,
                                 eta    = curr$eta,
                                 theta  = curr$theta,
                                 loglik = sum(curr$log_likelihood),
                                 loglik_int = curr$log_likelihood_int)
    }
  }

  list(post = posterior, ar = c(curr$accept, curr$reject))
}
