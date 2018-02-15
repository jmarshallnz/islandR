water_attribution = function(curr, water, phi) {

  # proposal distribution
  theta_proposal_sigma = 2

  # prior distribution for water
  theta_w0    = matrix(0, nrow(curr$theta[[1]]), ncol(curr$theta_w))
  theta_w_prec = diag(0.1, nrow(curr$theta_w)) # same prior across all sources

  # prior distribution for humans
  theta_h0    = matrix(0, nrow(curr$theta_h), ncol(curr$theta_h))
  theta_h_prec = diag(0.1, nrow(curr$theta_h)) # same prior across all sources

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

}

# basic idea for tree based attribution.
# each source and sink is a node in a directed tree (not sure if it can be a DAG?)
# each has a current realisation from a probability distribution p=p(ST|node)
# each has a set of data (STs) - not needed for sources
# each has a list of nodes it is linked to.
# each has a design matrix of other data - not needed for sources, only sinks. This has dimension
#     equal to the number of STs by number of parameters.
# each has a current realisation of parameters - dimension number of sources minus 1 by #params per obs.
#
# We then iterate over nodes, updating them.
#  For sinks, we update the parameter and then compute the new probability distribution.
#  For sources, we leave everything fixed - nothing changes as we assume we know everything.
# We know we have a source as it doesn't link to anything.

# THEN, later on,
