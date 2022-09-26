library(islandR)
library(tidyverse)

# fit AI model to manawatu data
st = st_fit(formula = Source ~ ST,
            non_primary = "Human",
            data = manawatu,
            method="island",
            sequences = ~ ASP + GLN + GLT + GLY + PGM + TKT + UNC)

# prior for eaah sources
prior <- tibble(source = setdiff(unique(manawatu$Source), "Human"),
                prior = 1/4)

# the st_fit function returns log_p, and for precision, when computing
# p/sum(p) it is best to first factor out the maximal exponent:
exp_over_exp_sum_log <- function(x) { y = x - max(x); exp(y)/sum(exp(y)) }

# compute posterior
st_df |> as.data.frame() |> # convert to data.frame
  left_join(prior) |> # join our prior info in
  mutate(log_pp = log(prior) + log_p) |> # multiply prior*p (numerator/denominator of Bayes' Theorem)
  group_by(type, iteration) |>
  mutate(posterior = exp_over_exp_sum_log(log_pp)) |> # compute exp(log_pp)/sum(exp(log_pp))
  select(type, iteration, posterior)

