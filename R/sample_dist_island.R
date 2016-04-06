#' @useDynLib islandR
#' @importFrom Rcpp sourceCpp
NULL

#' Fit the sampling distribution of genotypes to sources using the assymmetric island model.
#'
#' This routine uses the asymmetric island model (See \href{http://dx.doi.org/10.1371/journal.pgen.1000203}{D. Wilson et al 2008}) to estimate the sampling distribution
#' of genotypes on each source. Once estimated, the sampling distribution may be used with
#' \code{\link{attribution}} to attribute cases to sources based on genotype data in addition to
#' other covariates on the cases.
#'
#' @export
#' @param formula A formula of the form Genotype ~ Source
#' @param sequences A formula of the form ~ASP + GLN + GLT specifying columns for the allelic profile of each type
#' @param non_primary one or more sources that should be considered 'output' sources. No genotype distribution is computed for these,
#'        but P(ST | source) will be computed for these STs in addition to those observed on other sources.
#' @param iters the number of iterations to sample. An additional burnin period of 10% of this value will be used, and
#'        100 samples will be taken (or every sample, whichever is smaller).
#' @param data optional data frame from which to take variables in \code{formula} and \code{sequence}.
#' @return an object of class island which derives from sampling_dist.
#' @seealso \code{\link{st_fit}}, \code{\link{print.sample_dist}}, \code{\link{plot.sample_dist}}, \code{\link{summary.sample_dist}}
st_fit_island <- function(formula, sequences, non_primary = "Human", iters = 10000, data) {
  mod.terms = terms(formula, data=data)
  mod.frame = model.frame(formula, data=data)
  allele.frame = model.frame(sequences, data=data)

  if (attr(mod.terms, "response") == 0)
    stop("formula needs a left hand side")
  if (length(labels(mod.terms)) != 1)
    stop("formula needs a single entry on the right hand side")
  if (nrow(mod.frame) != nrow(allele.frame))
    stop("mismatch in number of rows between source/type and allelic profile")

  mod.vars = attr(mod.terms, "variables")
  response = all.vars(mod.terms[[2]]) # response variable

  # find unique types, as they're what we want to give attribution for
  # (techinically only need non-primary types, but it's easiest to return
  # all, and they're interesting anyway)
  type = labels(mod.terms)
  unique = !duplicated(data[,type])

  # filter out the non-primaries, and unique types
  sources = data[!(data[,response] %in% non_primary),]
  types = data[unique,]

  # convert reponse variable to a factor, and then to numeric
  sources[,response] = droplevels(as.factor(sources[,response]))
  source_names = levels(sources[,response])
  sources[,response] = as.numeric(sources[,response])

  sources.frame = model.frame(formula, data = sources)
  allele.frame  = model.frame(sequences, data = sources)

  source.frame = cbind(Type=sources.frame[,type], allele.frame, Source=sources.frame[,response])

  # add on the unique types to estimate the sample distribution on
  types.frame = model.frame(formula, data = types)
  allele.frame  = model.frame(sequences, data = types)

  type.frame = cbind(Type=types.frame[,type], allele.frame, Source=0)

  # convert to a matrix
  island.mat = as.matrix(rbind(source.frame, type.frame))

  # run the island model
  out = island(island.mat, niter = iters)

  # now regenerate useful summary information... (in future assign names in C++ land?) Also
  # need to change output in C++ land so it gives all types rather than duplicates and only
  # human STs (as in future, other non-primary sts will be needed, like water)

  set_names <- function(x) {
    rownames(x) = type.frame$Type
    colnames(x) = source_names
    x
  }
  out = lapply(out, set_names)

  # righto, now construct a useful object...
  x = list(types = type.frame$Type, sequences = type.frame[,-c(1,ncol(type.frame))], sources = source_names, sampling_distribution = simplify2array(out), model = "island")
  class(x) = c("island", "sample_dist")
  x
}
