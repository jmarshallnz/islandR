#' Island model, metropolis-hastings fit on logit scale
#'
#' Eventually this will include AR1() model

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

update_p = function(curr, humans, phi) {

  # proposal distribution
  logit_p_proposal_sigma_2 = 1

  # update each of the logit_p's
  # TODO: Could do this at random instead
  for (id in 1:length(curr$logit_p)) {
    # update one of the logit_p's
    logit_p = curr$logit_p;
    logit_p[id] = rnorm(1, curr$logit_p[id], logit_p_proposal_sigma);

    # inverse logit to new P
    p = inverse_logit(logit_p)

    # compute prior-hastings ratio
    # Prior-Hastings ratio = Proposal(f,f')/Proposal(f',f) * Prior(f')/Prior(f)
    # Proposal is normal distribution so is symmetric, so this drops down to the prior.
    # Prior ratio exp((p^2-p'^2)/(2*sigma^2))
    #
    log_hastings_ratio = (curr$logit_p[id]^2-logit_p[id]^2)/(2*logit_p_proposal_sigma);

    # compute likelihood ratio
    log_likelihood = log_lik(humans, phi, p)
    log_likelihood_ratio = log_likelihood - curr$log_likelihood

    # accept/reject
    log_alpha = log_likelihood_ratio + log_hastings_ratio
    if (log_alpha > 0 || runif(1) < exp(log_alpha)) {
      # accept new values, copy to old
      curr$logit_p = logit_p
      curr$log_likelihood = log_likelihood
      curr$accept = curr$accept + 1
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

  # proposal
  logit_p_proposal_sigma = 1

  # accept/reject
  accept_reject = numeric(2)

  n_sources = ncol(phi)

  # posterior
  posterior_logit_p = matrix(NA, iterations/thinning, n_sources - 1)

  #' sample logit_p from the prior
  logit_p = rnorm(n_sources-1, 0, logit_p_sigma)

  #' convert from logit_p to P
  p = inverse_logit(logit_p)

  log_likelihood = log_lik(humans, phi, p)

  # storage for the current iteration
  curr = list(logit_p = logit_p, log_likelihood = log_likelihood, accept=0)

  # main MCMC loop
  post_i = 0;
  for (i in seq_len(iterations+burnin)) {

    # update the p's
    curr = update_p(curr, humans, phi)

    # sample
    if (i %% 1000 == 0) {
      cat("Up to iteration", i, "of", burnin + iterations, "\n")
    }
    if (i > burnin && i %% thinning == 0) {
      post_i = post_i + 1;
      posterior_logit_p[post_i,] = curr$logit_p
    }
  }

  list(post = posterior_logit_p, ar = curr$accept / (iterations + burnin))
}
