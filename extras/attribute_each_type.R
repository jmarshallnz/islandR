# here is the code for computing per-human attribution

# TODO: Tidy this up so it's much faster so we can include it directly in the package. ATM we're basically
#       running a loop on what should be sort-of (due to log/exp) a matrix multiply.

attribute_types <- function(mod, types = NULL, prior_p = NULL) {
  genotype_dist <- mod$genotype_data
  genotype_iters <- iterations(genotype_dist)

  post <- predict(mod, FUN=identity) %>%
    mutate(iteration = (identity-1)%/%genotype_iters+1)

  # Human posterior attribution per type: Combine human attribution with type
  # distribution for each source via Bayes.
  # i.e. P(Source | ST) = P(ST | Source) * P(Source) / sum(P(ST|Source) * P(Source))
  #
  # Note that Animal attribution will use the same, but with P(Source) being constant
  # for each source. In this case, can just use P(Source) = 1
  #
  # run through for each TYPE.
  posterior_prob <- function(log_p_given_source, p) {
    # This should do:
    #  num <- exp(log_p_given_source) * p
    #  num / sum(num)
    # but, log_p_given_source is often very small, leading to underflows (i.e. exp = 0)
    # whereby sum(num) = 0. Thus, we first factor out the maximum log to avoid this
    max_log_p_given_source <- max(log_p_given_source)
    num <- exp(log_p_given_source - max_log_p_given_source) * p
    #    num <- p_given_source * p
    num / sum(num)
  }

  p_type <- as.data.frame(genotype_dist, types)
  grouped <- p_type %>% left_join(post, by=c("iteration", "source"="Source")) %>%
    group_by(iteration, identity, X, type)
  if (is.null(prior_p)) {
    done <- grouped %>%
      mutate(p_source_given_st = posterior_prob(log_p, p))
  } else {
    done <- grouped %>%
      mutate(p_source_given_st = posterior_prob(log_p, prior_p))
  }
  done %>% ungroup() %>% rename(Type = type, Source=source)
}

