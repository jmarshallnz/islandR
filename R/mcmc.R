#' @importFrom MASS mvrnorm
NULL

# Log-likelihood function
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

update_hyper_theta = function(curr) {
  # hyper-priors
  theta_0    = matrix(0, nrow(curr$theta), ncol(curr$theta))
  theta_prec = diag(0.1, nrow(curr$theta)) # same prior across all sources
  tau_shape  = 1 #0.1
  tau_rate   = 10 #0.1
  rho_0      = 0
  rho_prec   = 4

  ntimes     = max(curr$t, 1)

  if (ntimes > 1) {

    ## AR1 model

  # The sampling is based on
  #
  # Bayes regression with autoregressive errors: A Gibbs sampling approach
  # by Siddhartha Gibb (1993)
  #
  # However, it is generalised so that it works for multivariate responses and
  # parameter vectors.
  #
  # 1. Transform to eliminate the auto-correlation
  #
  # p_star = rho(L)p
  # X_star = rho(L)X
  #
  # NOTE: We assume here that the data are ordered in time, and that
  #       subsequent rows have the same non-time covariates. Thus
  #       any location for example needs to be blocked so that within
  #       block the time is sequential.
  #
  t1 = curr$t != 1
  t2 = curr$t != ntimes
  p = curr$p[t1,] - curr$rho*curr$p[t2,]
  X = curr$X[t1,] - curr$rho*curr$X[t2,]

  #
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

  var_hat   = solve(theta_prec + curr$tau * t(X) %*% X)
  theta_hat = var_hat %*% (theta_prec %*% theta_0 + curr$tau * t(X) %*% p)

  # now, to do the sampling we need to iterate over the columns of theta_hat
  # as it is a matrix
  for (j in 1:ncol(curr$theta))
    curr$theta[,j] = mvrnorm(1, theta_hat[,j], var_hat)

  # 2. Sample tau
  #
  # prior tau ~ Gamma(tau_shape, tau_rate)
  #
  # posterior will be (tau_shape + n/2, tau_rate + (p - X theta)'(p - X theta)/2)

  e = p - X %*% curr$theta
  rss = sum(e^2)

  n     = nrow(p)*ncol(p)
  shape = tau_shape + n / 2
  rate  = tau_rate  + rss / 2
  curr$tau   = rgamma(1, shape, rate=rate)

  # 3. Sample rho
  #
  # Prior is TN(rho_0, rho_prec)
  #
  # Given y, tau, theta, the error e_t = p_t - x_t theta becomes degenerate. Thus
  #
  # e_t = rho e_{t-1} + eta_t
  #
  # so the posterior is a truncated normal

  e   = curr$p - curr$X %*% curr$theta
  E   = e[t2,]
  eps = e[t1,]

  var = 1 / (rho_prec + curr$tau * sum(E^2))
  mu  = var * (rho_prec * rho_0 + curr$tau * sum(E*eps))

  rho = rnorm(1, mu, sqrt(var))
  if (rho < 1 && rho > -1)
    curr$rho = rho


  } else {

    # no time, so just solve directly

#    curr$theta = matrix(0, nrow(curr$theta), ncol(curr$theta))

    # deterministic update as there is no hierachical model without time involved.
    # if there was a hierarchical model, we'd need the bit commented out as well,
    # but I dunno what sort of model that would be. I guess we'd need a random effect
    # for this, but there's no obvious random effects to include in our model, right?

    # maybe a different p for the same covariate pattern across multiple people?
    # which I guess may be soaking up additional variation?

#    curr$theta = solve(t(X) %*% X) %*% t(X) %*% p
#    curr$tau = 1
#
     X = curr$X
     p = curr$p

     var_hat   = solve(theta_prec + curr$tau * t(X) %*% X)
     theta_hat = var_hat %*% (theta_prec %*% theta_0 + curr$tau * t(X) %*% p)

     # now, to do the sampling we need to iterate over the columns of theta_hat
     # as it is a matrix
     for (j in 1:ncol(curr$theta))
       curr$theta[,j] = mvrnorm(1, theta_hat[,j], var_hat)

     # 2. Sample tau
     #
     # prior tau ~ Gamma(tau_shape, tau_rate)
     #
     # posterior will be (tau_shape + n/2, tau_rate + (p - X theta)'(p - X theta)/2)
     e = p - X %*% curr$theta
     rss = sum(e^2)

     n     = nrow(p)*ncol(p)
     shape = tau_shape + n / 2
     rate  = tau_rate  + rss / 2
     curr$tau   = rgamma(1, shape, rate=rate)
  }
  curr
}

log_lik_full = function(humans, phi, p) {
  ll = numeric(length(humans))
  for (i in 1:length(humans)) {
    ll[i] = log_lik(humans[[i]], phi, p[i,])
  }
  ll
}

# non-hierarchical update for theta
update_theta = function(curr, humans, phi) {

  # proposal distribution
  theta_proposal_sigma = 2

  # prior distribution
  theta_0    = matrix(0, nrow(curr$theta), ncol(curr$theta))
  theta_prec = diag(0.1, nrow(curr$theta)) # same prior across all sources

  # in this updater, p is just X %*% theta

  # for each coefficient, sampled at random
  rows = sample(nrow(curr$theta))
  for (i in rows) {

    # for each source, sampled at random
    cols = sample(ncol(curr$theta))
    for (j in cols) {
      # update theta[i,j] with RWMH
      theta = curr$theta
      theta[i,j] = rnorm(1, curr$theta[i,j], theta_proposal_sigma)
      # compute new p
      p = curr$X %*% theta

      # proposal ratio is symmetric, so need only the prior ratio
      t0_prop = theta[,j]      - theta_0[,j]
      t0_curr = curr$theta[,j] - theta_0[,j]
      # exp(-0.5*(t(t0_prop) %*% theta_prec %*% t0_prop)) / exp(-0.5*(t(t0_curr) %*% theta_prec %*% t0_curr))

      log_hastings_ratio = -0.5*(t(t0_prop) %*% theta_prec %*% t0_prop -
                                 t(t0_curr) %*% theta_prec %*% t0_curr)

      # compute likelihood ratio.
      # (partial log-likelihood could be used here, as we need only update
      #  the pattern affected by the theta that was changed. However, this is
      #  probably more cumbersome to compute, so we instead just update
      #  the entire p)

      log_likelihood = log_lik_full(humans, phi, p)
      log_likelihood_ratio = sum(log_likelihood - curr$log_likelihood)

      # accept/reject
      log_alpha = log_likelihood_ratio + log_hastings_ratio
      if (log_alpha > 0 || runif(1) < exp(log_alpha)) {
        # accept new values, copy to old
        curr$p[,j] = p[,j]
        curr$theta[i,j] = theta[i,j]
        curr$log_likelihood = log_likelihood
        curr$accept = curr$accept + 1
      } else {
        curr$reject = curr$reject + 1
      }
    }
  }
  return(curr)
}

update_ranef = function(curr, humans, phi) {

  # proposal distribution
  p_proposal_sigma = 2

  # first up, compute the fitted values and residuals
  mu = curr$X %*% curr$theta

  n_times = max(curr$t, 1)

  # for each covariate pattern (sampled at random)
  # TODO: if we decide to block update, we'd need to change this
  rows = sample(nrow(curr$p))
  for (i in rows) {

    cols = sample(ncol(curr$p))
    for (j in cols) {

    # note that the likelihood only changes for the particular
    # covariate pattern, so we needn't compute the entire thing

    # update the corresponding p

    # The below assumes we're proposing from the conditional prior
    #
    # Thus the prior-hastings ratio is 1.
    #
    # The conditional prior is
    #
    # e[t] ~ rho*e[t-1] + Normal(0, tau)
    #
    # which can be rewritten in various forms. We could potentially
    # block sample using this technique.
    #
    # we assume throughout that the arrangement of the data are such
    # that the previous and/or next time point for this covariate
    # pattern is the next or previous observation respectively

    p = curr$p[i,]

    if (n_times > 1) {

      ## AR1 model...

    updater = sample(2, 1)

    t = curr$t[i]
    if (updater == 1) {
      e = curr$p - mu
      en = e[i,]
      en[j] = rnorm(1, e[i,j], 2*p_proposal_sigma)
      p[j] = mu[i,j] + en[j]
      if (t == 1) {
        # prior ratio
        # p[1] - mu[1] - rho*(p[2] - mu[2]) ~ Normal(0, tau)
        log_hastings_ratio = -0.5*curr$tau*((en[j] - curr$rho*e[i+1,j])^2
                                           -(e[i,j] - curr$rho*e[i+1,j])^2)
      } else if (t == n_times) {
        log_hastings_ratio = -0.5*curr$tau*((en[j] - curr$rho*e[i-1,j])^2
                                           -(e[i,j] - curr$rho*e[i-1,j])^2)
      } else {
        log_hastings_ratio = -0.5*0.5*curr$tau*((en[j] - 0.5*curr$rho*(e[i-1,j] + e[i+1,j]))^2
                                               -(e[i,j] - 0.5*curr$rho*(e[i-1,j] + e[i+1,j]))^2)
      }
    } else {

      if (t == 1) {
        # rearranging to find the conditional prior on e[1] gives
        #
        #  e[1] ~ Normal(1/rho*e[2], tau/rho^2)
        #
        # is what you'd think. However, it turns out (mysteriously)
        # that it is actually
        #
        #  e[1] ~ Normal(rho*e[2], tau)
        #
        # TODO: not exactly sure why though??? Google backcasting?
        #
        # as we're proposing from the prior, the hastings ratio is 1
        #
        # however, if rho is close to 0, we want to sample independently
        p[j] = rnorm(1, mu[i,j] + (curr$p[i+1,j] - mu[i+1,j])*curr$rho, 1/sqrt(curr$tau))
      } else if (t == n_times) {
        #  e[t] ~ Normal(rho*e[t-1], tau)
        p[j] = rnorm(1, mu[i,j] + (curr$p[i-1,j] - mu[i-1,j])*curr$rho, 1/sqrt(curr$tau))
      } else {
        # Here we're proposing from the prior. This means proposal is prior
        # so hastings ratio is 1.

        # Alternate is to propose from some other symmetric distribution, so
        # hastings ratio is prior ratio.

        # We can use the above and below to get a better proposal conditioned
        # on all the other values.
        #   p[t] - mu[t] - rho*(p[t-1] - mu[t-1]) ~ Normal(0, tau)
        #   1/rho(p[t+1] - mu[t+1]) - (p[t] - mu[t]) ~ Normal(0, tau/rho^2)

        #   p[t] - mu[t] - rho*(p[t-1] - mu[t-1]) ~ Normal(0, tau)
        #   -1/rho(p[t+1] - mu[t+1]) + (p[t] - mu[t]) ~ Normal(0, tau/rho^2)
        #
        #   p[t] - mu[t] - rho/2*(p[t-1] - mu[t-1]) - 1/rho/2*(p[t+1] - mu[t+1]) ~ Normal(0, tau*(1 + 1/rho^2) / 4)

        # BUT: Mysteriously, the backcasting suggests this is not the case???
        #     Again, google backcasting AR(1) processes
        # f[t] ~ Normal(mu + rho*(f[t-1] + f[t+1] - 2*mu)/2, tau*2)
        p[j] = rnorm(1, mu[i,j] + 0.5*(curr$p[i-1,j] - mu[i-1,j] + curr$p[i+1,j] - mu[i+1,j])*curr$rho, 1/sqrt(2*curr$tau))
      }
      log_hastings_ratio = 0;
    }

    } else {
      ## No time aspect -> use standard M-H
      e = curr$p - mu
      en = e[i,]
      en[j] = rnorm(1, e[i,j], p_proposal_sigma)
      p[j] = mu[i,j] + en[j]

      # prior ratio
      # p[1] - mu[1] ~ Normal(0, tau)
      log_hastings_ratio = -0.5*curr$tau*(en[j]^2 - e[i,j]^2)
    }

    # compute likelihood ratio
    log_likelihood = log_lik(humans[[i]], phi, p)
    log_likelihood_ratio = log_likelihood - curr$log_likelihood[i]

    # accept/reject
    log_alpha = log_likelihood_ratio + log_hastings_ratio
    if (log_alpha > 0 || runif(1) < exp(log_alpha)) {
      # accept new values, copy to old
      curr$p[i,] = p
      curr$log_likelihood[i] = log_likelihood
      curr$accept = curr$accept + 1
    } else {
      curr$reject = curr$reject + 1
    }
    }
  }
  return(curr)
}

mcmc_no_ar1 = function(humans, X, phi, iterations = 10000, burnin = 1000, thinning = 100) {

  # priors
  logit_p_sigma = 1

  # accept/reject
  accept_reject = numeric(2)

  # posterior
  posterior = list()

  # parameter vector
  n_sources = ncol(phi)
  theta   = matrix(0, ncol(X), n_sources-1)
  rownames(theta) <- colnames(X)
  colnames(theta) <- colnames(phi)[-n_sources]

  # initialise p
  p = matrix(0, nrow(X), n_sources-1)

  # compute log-likelihood for each covariate pattern
  log_likelihood = numeric(length(humans))
  for (x in seq_along(humans)) {
    log_likelihood[x] = log_lik(humans[[x]], phi, p[x,])
  }

  # storage for the current iteration
  curr = list(p = p, X = X, t = NULL, theta = theta, tau = NULL, rho = NULL, log_likelihood = log_likelihood, accept=0, reject=0)

  # main MCMC loop
  post_i = 0;
  for (i in seq_len(iterations+burnin)) {

    # update the theta's + p's
    curr = update_theta(curr, humans, phi)

    # sample
    if (i %% 1000 == 0) {
      cat("Up to iteration", i, "of", burnin + iterations, "\n")
    }
    if (i > burnin && i %% thinning == 0) {
      post_i = post_i + 1;
      posterior[[post_i]] = list(p     = curr$p,
                                 theta = curr$theta)
    }
  }

  list(post = posterior, ar = c(curr$accept, curr$reject))
}

mcmc = function(humans, x, formula, phi, iterations = 10000, burnin = 1000) {

  # MCMC parameters
  thinning   = 100

  # priors
  logit_p_sigma = 1

  # accept/reject
  accept_reject = numeric(2)

  # posterior
  posterior = list()

  # time column (in case time is important)
  t = x$Time

  # design matrix
  if (is.null(x))
    x = data.frame(dummy=1)
  X = model.matrix(formula, data=x)

  # parameter vector
  n_sources = ncol(phi)
  theta   = matrix(0, ncol(X), n_sources-1)
  rownames(theta) <- colnames(X)

  ## HACK
#  theta   = matrix(0, ncol(X), n_sources)
  ## ENDHACK

  # precision, auto-correlation
  tau     = 1
  rho     = 0

  hierarchical = FALSE
  # initialise p
  p = matrix(0, nrow(X), n_sources-1)

  ## HACK
#  p = matrix(0, nrow(X), n_sources)
  ## ENDHACK

  # compute log-likelihood for each covariate pattern
  log_likelihood = numeric(length(humans))
  for (x in seq_along(humans)) {
    log_likelihood[x] = log_lik(humans[[x]], phi, p[x,])
  }

  # storage for the current iteration
  curr = list(p = p, X = X, t = t, theta = theta, tau = tau, rho = rho, log_likelihood = log_likelihood, accept=0, reject=0)

  # main MCMC loop
  post_i = 0;
  for (i in seq_len(iterations+burnin)) {

    if (hierarchical) {
      # update priors
      curr = update_hyper_theta(curr)

      # update the p's
      curr = update_ranef(curr, humans, phi)
    } else {
      # update the theta's + p's
      curr = update_theta(curr, humans, phi)
    }

    # sample
    if (i %% 1000 == 0) {
      cat("Up to iteration", i, "of", burnin + iterations, "\n")
    }
    if (i > burnin && i %% thinning == 0) {
      post_i = post_i + 1;
      posterior[[post_i]] = list(p     = curr$p,
                                 theta = curr$theta,
                                 tau   = curr$tau,
                                 rho   = curr$rho)
    }
  }

  list(post = posterior, ar = c(curr$accept, curr$reject))
}
