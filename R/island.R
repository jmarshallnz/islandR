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

update_priors = function(curr) {
  # hyper-priors
  alpha_mu   = 0
  alpha_prec = 0.1
  tau_shape  = 0.1
  tau_rate   = 0.1

  # first, our means
  for (i in 1:ncol(curr$logit_p)) {
    # posterior will be alpha_mu*alpha_prec + n * tau * mean(p) / alpha_prec + n * tau
    mu_hat = mean(curr$logit_p[,i])
    n      = nrow(curr$logit_p)
    prec = alpha_prec + n * curr$tau
    mu   = (alpha_mu*alpha_prec + n*mu_hat*curr$tau) / prec
    curr$mu[i] = rnorm(1, mu, 1/sqrt(prec))
  }

  # update TAU
  # posterior will be (tau_shape + n/2, tau_rate + sum(p - mu)^2/2)
  rss = 0;
  for (i in 1:ncol(curr$logit_p))
    rss = rss + sum((curr$logit_p[,i] - curr$mu[i])^2)

  n     = nrow(curr$logit_p) * ncol(curr$logit_p)
  shape = tau_shape + n / 2
  rate  = tau_rate  + rss / 2
  curr$tau   = rgamma(1, shape, rate=rate)

  curr
}

update_p = function(curr, humans, phi) {

  # proposal distribution
  logit_p_proposal_sigma_2 = 1

  # TODO: Could do this at random instead
  # for each covariate patterns
  for (x in 1:nrow(curr$logit_p)) {
    # we can optimise this quite a bit, as we need only update
    # the log likelihood corresponding to this covariate pattern
    # as every other param is staying the same

    h = humans[[x]]

    # update each of the logit_p's
    for (id in 1:ncol(curr$logit_p)) {
      logit_p = curr$logit_p[x,];
      logit_p[id] = rnorm(1, curr$logit_p[x,id], logit_p_proposal_sigma);

      # inverse logit to new P
      p = inverse_logit(logit_p)

      # compute prior-hastings ratio
      # Prior-Hastings ratio = Proposal(f,f')/Proposal(f',f) * Prior(f')/Prior(f)
      # Proposal is normal distribution so is symmetric, so this drops down to the prior.

      # Prior in our hierarchical model is determined by mu + tau
      # exp(((p-mu)^2-(p'-mu)^2)/2*tau)
      log_hastings_ratio = ((curr$logit_p[x,id] - curr$mu[id])^2 -
                            (logit_p[id] - curr$mu[id])^2)*0.5*curr$tau;

      # compute likelihood ratio
      log_likelihood = log_lik(h, phi, p)
      log_likelihood_ratio = log_likelihood - curr$log_likelihood[x]

      # accept/reject
      log_alpha = log_likelihood_ratio + log_hastings_ratio
      if (log_alpha > 0 || runif(1) < exp(log_alpha)) {
        # accept new values, copy to old
        curr$logit_p[x,] = logit_p
        curr$log_likelihood[x] = log_likelihood
        curr$accept = curr$accept + 1
      } else {
        curr$reject = curr$reject + 1
      }
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
  n_times   = length(humans)

  # posterior
  posterior_logit_p = list()

  #' sample logit_p from the prior
  logit_p = matrix(rnorm(n_times * (n_sources-1), 0, logit_p_sigma), n_times)
  mu      = rep(0, n_sources - 1)
  tau     = 1

  #' compute log-likelihood for each covariate pattern
  log_likelihood = numeric(n_times)
  for (x in 1:n_times) {
    p = inverse_logit(logit_p[x])
    log_likelihood[x] = log_lik(h[[x]], phi, p)
  }

  # storage for the current iteration
  curr = list(logit_p = logit_p, mu = mu, tau = tau, log_likelihood = log_likelihood, accept=0, reject=0)

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
      posterior_logit_p[[post_i]] = list(logit_p = curr$logit_p,
                                         mu      = curr$mu,
                                         tau     = curr$tau)
    }
  }

  list(post = posterior_logit_p, ar = curr$accept / (curr$accept + curr$reject))
}
