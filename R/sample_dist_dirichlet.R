rdirichlet<-function(n,a)
{
  l  <- length(a);
  x  <- matrix(rgamma(l*n,a), ncol=l, byrow=TRUE);
  sm <- x %*% rep(1,l);
  x / as.vector(sm);
}

#' Fit the sampling distribution of genotypes to sources using a Dirichlet/multinomial distribution.
#'
#' This routine estimates the sampling distribution of the genotypes on each
#' source by assuming the observed counts arise from a multinomial distribution,
#' with a Dirichlet prior on the parameter vector, giving a Dirichlet posterior.
#' Once estimated, the sampling distribution may be used with
#' \code{\link{attribution}} to attribute cases to sources based on genotype data in addition to
#' other covariates on the cases.
#'
#' @export
#' @param formula A formula of the form Source ~ Genotype
#' @param prior   The Dirichlet prior. Currently only a constant. Defaults to 1.
#' @param non_primary one or more sources that should be considered 'output' sources. No genotype distribution is computed for these,
#'        but P(ST | source) will be computed for these STs in addition to those observed on other sources.
#' @param samples the number of iterations to sample. Defaults to 100.
#' @param data optional data frame from which to take variables in \code{formula} and \code{sequence}.
#' @return an object of class dirichlet which derives from sampling_dist.
#' @seealso \code{\link{st_fit}}, \code{\link{print.sample_dist}}, \code{\link{plot.sample_dist}}, \code{\link{summary.sample_dist}}
st_fit_dirichlet <- function(formula, prior = 1, non_primary = "Human", samples = 100, data) {
  mod.terms = terms(formula, data=data)
  mod.frame = model.frame(formula, data=data)

  if (attr(mod.terms, "response") == 0)
    stop("formula needs a left hand side")
  if (length(labels(mod.terms)) != 1)
    stop("formula needs a single entry on the right hand side")

  mod.vars = attr(mod.terms, "variables")
  response = all.vars(mod.terms[[2]]) # response variable

  # find unique types, as they're what we want to give attribution for
  # (technically only need non-primary types, but it's easiest to return
  # all, and they're interesting anyway)
  type = labels(mod.terms)

  # We basically just need a matrix of unique types by their counts
  counts = table(data[,type], data[,response])

  # then we just need to pick out the source columns
  source_names = colnames(counts)[!colnames(counts) %in% non_primary]
  source_counts = counts[,source_names]

  # and now we add on prior
  alpha = source_counts + prior

  # and generate samples from a Dirichlet per column
  dirichlet_matrix <- function(x) {
    a <- apply(x, 2, function(y) { rdirichlet(1, y) })
    rownames(a) <- rownames(x)
    colnames(a) <- colnames(x)
    log(a)
  }

  out <- replicate(samples, dirichlet_matrix(alpha))
  # righto, now construct a useful object...
  x = list(types = rownames(counts), sources = source_names,
           sampling_distribution = out, model = "dirichlet")
  class(x) = c("dirichlet", "sample_dist")
  x
}
