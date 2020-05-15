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
#' @param formula A formula of the form Source ~ Genotype
#' @param sequences A formula of the form ~ASP + GLN + GLT specifying columns for the allelic profile of each type
#' @param non_primary one or more sources that should be considered 'output' sources. No genotype distribution is computed for these,
#'        but P(ST | source) will be computed for these STs in addition to those observed on other sources.
#' @param samples the number of samples required from the posterior after burnin. Defaults to 100.
#' @param burnin  the number of samples to discard for burnin. Defaults to 10.
#' @param thin    the number of iterations per sample taken. Defaults to 100.
#' @param priors - a list of priors for the fit. Defaults to 1,c(1,1),c(1,1).
#' @param control - a list of model control information. Defaults to mutation="source", recombination="source". Specify "constant" for
#'        a single probability across all sources.
#' @param data optional data frame from which to take variables in \code{formula} and \code{sequence}.
#' @return an object of class island which derives from sampling_dist.
#' @seealso \code{\link{st_fit}}, \code{\link{print.sample_dist}}, \code{\link{plot.sample_dist}}, \code{\link{summary.sample_dist}}
st_fit_island <- function(formula,
                          sequences,
                          non_primary = "Human",
                          samples = 100, burnin = 10, thin = 100,
                          priors = list(migration=1,mutation=c(1,1),recombination=c(1,1)),
                          control = list(mutation='source', recombination='source'),
                          data) {
  mod.terms = terms(formula, data=data)
  mod.frame = model.frame(formula, data=data)
  allele.frame = model.frame(sequences, data=data)

  if (length(priors$recombination) != 2) {
    stop("Beta prior on priors$recombination must be length 2")
  }
  if (length(priors$mutation) != 2) {
    stop("Beta prior on priors$mutation must be length 2")
  }

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
  unique = !duplicated(data[[type]])

  # the island model assumes all alleles (and all types) are non-negative integers, as it uses them to index stuff.
  # let's assume that they're just character strings, convert them to factors and then back again. We'll then have
  # a map for them.
  sequence_cols <- all.vars(terms(sequences))
  allele.factor <- lapply(data[sequence_cols], as.factor)
  allele.map    <- lapply(allele.factor, levels)
  allele.numeric <- lapply(allele.factor, as.numeric)

  # replace our sequence data with numeric. We can use allele.map to map back again
  data[sequence_cols] <- allele.numeric

  # filter out the non-primaries, and unique types
  sources = data[!(data[[response]] %in% non_primary),]
  types = data[unique,]

  # convert reponse variable to a factor, and then to numeric
  sources[[response]] = droplevels(as.factor(sources[[response]]))
  source_names = levels(sources[[response]])
  sources[[response]] = as.numeric(sources[[response]])

  sources.frame = model.frame(formula, data = sources)
  allele.frame  = model.frame(sequences, data = sources)

  source.frame = cbind(Type=sources.frame[[type]], allele.frame, Source=sources.frame[[response]])

  # add on the unique types to estimate the sample distribution on.

  # NOTE: This only works for unique *HUMAN* types. It doesn't make sense for unique *SOURCE* types
  #       due to the way the island model leave-one out works. i.e. Human types need to use a different
  #       computation to source types for uniqueness (see human_unique and beast_unique in island model code)
  #
  #       I think the only way to fix this easily without massively hacking the island model code (e.g.
  #       refactoring it so that it uses distance between types) is to have a flag for actual human type
  #       and beast type. i.e. get the island model to spit out the attribution for the beast types as
  #       well and then patch it up afterwards?
  types.frame = model.frame(formula, data = types)
  allele.frame  = model.frame(sequences, data = types)

  type.frame = cbind(Type=types.frame[[type]], allele.frame, Source=0)

  # convert to a matrix
  island.mat = as.matrix(rbind(source.frame, type.frame))

  # extract the mutation/recombination parameter matrices
  num_sources = length(source_names)
  X_mutation = diag(1, num_sources)
  if (!is.null(control$mutation) && control$mutation == 'constant')
    X_mutation = matrix(1, num_sources, 1)
  X_recombination = diag(1, num_sources)
  if (!is.null(control$recombination) && control$recombination == 'constant')
    X_recombination = matrix(1, num_sources, 1)

  # run the island model
  out = island(island.mat,
               beta_migration = priors$migration,
               gamma_mutation = priors$mutation,
               gamma_recombination = priors$recombination,
               X_mutation = X_mutation,
               X_recombination = X_recombination,
               samples = samples,
               burnin = burnin,
               thin = thin)

  # now regenerate useful summary information... (in future assign names in C++ land?) Also
  # need to change output in C++ land so it gives all types rather than duplicates and only
  # human STs (as in future, other non-primary sts will be needed, like water)

  set_names <- function(x) {
    rownames(x) = type.frame$Type
    colnames(x) = source_names
    x
  }
  hum_lik = lapply(out$hum_lik, set_names)

  # assign names to the evolution traces
  colnames(out$evolution) <- c("Iteration",
                               apply(as.matrix(expand.grid("A",1:length(source_names),1:length(source_names)))[,c(1,3,2)], 1, paste0, collapse=""),
                               paste0("M", 1:length(source_names)),
                               paste0("R", 1:length(source_names)),
                               "Likelihood")

  # assign names to the acceptance rates
  colnames(out$accept) <- c("Accept", "Reject")
  move_types <- c("A swap", "A rw", "M swap", "M rw", "R swap", "R rw")

  # righto, now construct a useful object...
  sequences = type.frame[,-c(1,ncol(type.frame))]
  # convert our sequences back
  x <- lapply(sequence_cols, function(col) { allele.map[[col]][ sequences[[col]] ] } )
  sequences[sequence_cols] <- x

  x = list(types = type.frame$Type,
           sequences = sequences,
           sources = source_names,
           sampling_distribution = hum_lik,
           evolution_params = data.frame(out$evolution[-1,]),
           acceptance = data.frame(Type = move_types, out$accept),
           model = "island")
  class(x) = c("island", "sample_dist")
  x
}
