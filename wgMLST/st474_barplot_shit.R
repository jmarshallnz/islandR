library(readxl)
library(islandR)
library(dplyr)
library(tidyr)
library(forcats)
library(ggplot2)
library(purrr)

load("wgMLST/474_attribution_update.RData")

# now do prediction on each isolate. The way to do this is to combine the predictions from
# above with the predictions from st_fit. Idea is that st_fit gives P(ST | source)
# and the attribution gives P(source), so P(source | ST) = P(ST | source)*P(source) / sum_j P(ST | j) P(j)
# using Bayes Thm.
posterior <- predict(mod, FUN=identity)
post <- posterior %>% mutate(iteration = (identity-1)%/%100+1)

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
  num / sum(num)
}

human_post_types <- list()
for (i in seq_along(st$types)) {
  cat("Processing type", i, "\n")
  p_type <- st$sampling_distribution[i,,]
  colnames(p_type) <- 1:ncol(p_type)
  # for each iteration from st_fit, we have an MCMC run for the attribution, so
  # pull these out and compute via Bayes Thm
  p_type <- p_type %>% as.data.frame %>% tibble::rownames_to_column("Source") %>%
    gather(iteration, p_given_source, -Source, convert=TRUE)

  human_post_types[[i]] <- p_type %>% left_join(post, by=c("iteration", "Source")) %>%
    group_by(iteration, identity) %>%
    mutate(p_source_given_st = posterior_prob(p_given_source, p)) %>%
    group_by(Source) %>% summarize(v_source_given_st = var(p_source_given_st),
                                   p_source_given_st = mean(p_source_given_st),
                                   Type = st$types[i])
}
human_type_attr <- do.call(rbind, human_post_types)
human_type_attr <- human_type_attr %>% left_join(final %>% select(FILE, Type=ST, source))

# right, now generate our barplots
human_types <- human_type_attr %>% filter(source == "Human") %>%
  group_by(Type) %>% mutate(v_source_given_st = sum(v_source_given_st)) %>% ungroup %>%
  spread(Source, p_source_given_st) %>% select(FILE, Uncertainty=v_source_given_st, EnvWater:`Supplier A`)

write.csv(human_types, "wgMLST/human_barplots_with_uncertainty.csv", row.names=FALSE)

gg_types <- human_types %>% arrange(`Supplier A`) %>% group_by(FILE) %>% gather(Source, P, EnvWater:`Supplier A`) %>%
  ungroup %>%
  mutate(FILE = fct_reorder2(FILE, Source, P, function(x, y) { y[x == "Supplier A"] }))

ggplot(gg_types, aes(x=FILE, y=P, fill=Source)) + geom_col()
