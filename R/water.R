log_likelihood <- function(node, F) {
  # run over data in this and any downstream nodes, multiplying shit up
  sinks <- nodes[[node$sinks]]
  p <- numeric(length(node$y))
  for (i in 1:length(node$sinks))
    p <- p + nodes[[i]]$p * F[i]
  node$y * log(p)
}

delogit <- function(p) {
  P <- exp(cbind(p, 0))
  P/rowSums(P)
}

# The model object is a list of sources/sinks. Each node
# has data, current probability dist, links to sources and
# links to sinks. The links to sources are associated with
# parameter values and a design matrix.
update_nodes <- function(nodes) {
  # run through the model nodes, updating each one
  for (i in 1:length(nodes)) {
    # if we're a sink, run over the parameters
    node <- nodes[[i]]
    rows <- sample(nrow(node$theta))
    for (r in rows) {
      # for each source, sampled at random
      cols = sample(ncol(node$theta))
      for (c in cols) {
        # update theta[r,c] with RWMH
        theta = node$theta
        theta[r,c] = rnorm(1, theta[r,c], theta_proposal_sigma)
        # compute new f
        f = node$X %*% theta

        # proposal ratio is symmetric, so need only the prior ratio
        t0_prop = theta[,j]      - theta_0[,j]
        t0_curr = node$theta[,j] - theta_0[,j]
        # exp(-0.5*(t(t0_prop) %*% theta_prec %*% t0_prop)) / exp(-0.5*(t(t0_curr) %*% theta_prec %*% t0_curr))

        log_hastings_ratio = -0.5*(t(t0_prop) %*% theta_prec %*% t0_prop -
                                     t(t0_curr) %*% theta_prec %*% t0_curr)

        # transform f from logit scale
        F = delogit(f)

        # use source mixing to compute ST distribution for this sink.

        # Hmm, how do covariates at a sink/source level work? The next sink
        # down doesn't care about those covariates, right?? So if we get an
        # improved attribution for subclasses of water how does that flow
        # through to the end result? What we normally end up with is the
        # ST dist given the covariate pattern (i.e. covariate pattern alters
        # the weights of the source ST dists). What we want is the ST dist
        # marginal. P(ST) = sum_i P(ST | covariate_i) P(covariate_i)

        # p will be a matrix of size #ST times #covariate patterns (X)
        # which is also the nrow of F. This is also the same size as y
        # I guess? (ATM we compress it down to table up only what we have
        # but this expands it out with zeros). This is really inefficient
        # in the case where we have lots of unique covariate patterns maybe?
        # in that case, X will basically be # of isolates where we work out
        # p per one (but we also need p's for the ones we don't observe?)
        p = matrix(0, nrow(node$y), ncol(node$y))
        # of vectors with length equal to number
        # of covariate patterns,
        p = rep(0, length(node$y))
        for (i in 1:length(node$sinks))
          p <- p + nodes[[i]]$p * F[i]

        # compute likelihood. This will involve p from above plus
        ll = log_likelihood(node, F)

      }
    }
  }
}

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
